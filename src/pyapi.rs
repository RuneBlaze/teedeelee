use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple, PyType};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use bincode;
use zstd;

use crate::porcelain::{
    CompleteTreeView, DistanceMatrixView, MSCTreeView, TopologyTreeView, TreeViewBase, TreeViewType,
};
use crate::tree::{parse_newick, MSCTreeCollection, SortCriterion, SortCriterionKind, SortOrder, TaxonSet, Tree, TreeCollection};

#[pyclass(name = "TreeSet", module = "teedeelee")]
#[derive(Clone)]
pub struct PyCompleteTreeView {
    inner: CompleteTreeView,
}

/// A Python wrapper for TopologyTreeView
#[pyclass(name = "TopologySet", module = "teedeelee")]
#[derive(Clone)]
pub struct PyTopologyTreeView {
    inner: TopologyTreeView,
}

/// A Python wrapper for MSCTreeView
#[pyclass(name = "MSCTreeSet", module = "teedeelee")]
#[derive(Clone)]
pub struct PyMSCTreeView {
    inner: MSCTreeView,
}

/// A Python wrapper for DistanceMatrixView
#[pyclass(name = "DistanceMatrix", module = "teedeelee")]
#[derive(Clone)]
pub struct PyDistanceMatrixView {
    inner: DistanceMatrixView,
}

/// A Python wrapper for a single tree with its taxon set
#[pyclass(name = "Tree", module = "teedeelee")]
#[derive(Clone)]
pub struct PySingleTree {
    taxon_set: Arc<TaxonSet>,
    tree: Arc<Tree>,
}

/// A Python wrapper for TaxonSet providing read-only access
#[pyclass(name = "TaxonSet", module = "teedeelee")]
#[derive(Clone)]
pub struct PyTaxonSet {
    inner: Arc<TaxonSet>,
}

/// A Python wrapper for MSCTreeView in topology-only mode
#[pyclass(name = "MSCTopologySet", module = "teedeelee")]
#[derive(Clone)]
pub struct PyMSCTopologyTreeView {
    inner: MSCTreeView,
}

/// A Python wrapper for a collection of MSC topology tree views
#[pyclass(name = "FamilyOfMSC", module = "teedeelee")]
#[derive(Clone)]
pub struct PyFamilyOfMSC {
    views: Vec<Arc<PyMSCTopologyTreeView>>,
}

/// Add these enum definitions for Python
#[pyclass(name = "SortCriterionKind", eq, eq_int, module = "teedeelee")]
#[derive(Clone, Copy, Eq, PartialEq)]
pub enum PySortCriterionKind {
    LexicographicalOrder,
    ChildCount,
    DescendantCount,
}

#[pyclass(name = "SortOrder", eq, eq_int, module = "teedeelee")]
#[derive(Clone, Copy, Eq, PartialEq)]
pub enum PySortOrder {
    Ascending,
    Descending,
}

/// Add this struct for Python
#[pyclass(name = "SortCriterion", eq, module = "teedeelee")]
#[derive(Clone, Eq, PartialEq)]
pub struct PySortCriterion {
    kind: PySortCriterionKind,
    order: PySortOrder,
}

/// Add these as proper Python enums
#[pyclass(name = "SortBy", eq, eq_int, module = "teedeelee")]
#[derive(Clone, Copy, Eq, PartialEq)]
pub enum PySortBy {
    LexicographicalOrder,
    ChildCount,
    DescendantCount,
}

#[pymethods]
impl PyCompleteTreeView {
    #[new]
    fn new(arg: Bound<'_, PyAny>) -> PyResult<Self> {
        // If input is a string, use existing newick parser
        if let Ok(newick_str) = arg.extract::<String>() {
            return CompleteTreeView::from_newick_string(&newick_str)
                .map(|inner| PyCompleteTreeView { inner })
                .map_err(|e| PyValueError::new_err(e));
        }

        // If input is a list, combine the trees
        if let Ok(tree_list) = arg.downcast::<PyList>() {
            let mut collection = TreeCollection::new(false);

            // First pass: merge all taxon sets
            for item in tree_list.iter() {
                let tree: PySingleTree = item.extract()?;
                for (name, _) in tree.taxon_set.to_id.iter() {
                    collection.taxon_set.request(name.clone());
                }
            }

            // Second pass: add trees with remapped taxa
            for item in tree_list.iter() {
                let tree: PySingleTree = item.extract()?;
                let mut id_map = HashMap::new();
                for (name, &old_id) in tree.taxon_set.to_id.iter() {
                    let new_id = collection.taxon_set.to_id[name];
                    id_map.insert(old_id, new_id);
                }
                let mut new_tree = (*tree.tree).clone();
                new_tree.remap_taxa(&id_map);
                collection.trees.push(new_tree);
            }

            return Ok(PyCompleteTreeView {
                inner: CompleteTreeView(TreeViewBase {
                    taxon_set: Arc::new(collection.taxon_set),
                    trees: collection.trees.into(),
                    topology_only: false,
                }),
            });
        }

        Err(PyValueError::new_err(
            "Argument must be either a Newick string or a list of Tree objects",
        ))
    }

    #[getter]
    fn ngenes(&self) -> usize {
        self.inner.ngenes()
    }

    #[getter]
    fn ntaxa(&self) -> usize {
        self.inner.ntaxa()
    }

    fn slice(&self, start: usize, end: usize) -> PyResult<Self> {
        Ok(PyCompleteTreeView {
            inner: self.inner.slice(start, end),
        })
    }

    fn get_distance_matrix(&self, index: usize) -> Option<PyDistanceMatrixView> {
        self.inner
            .get_distance_matrix(index)
            .map(|inner| PyDistanceMatrixView { inner })
    }

    fn get_taxon_name(&self, id: usize) -> Option<String> {
        self.inner.get_taxon_name(id).map(String::from)
    }

    fn to_newick_vec(&self) -> Vec<String> {
        self.inner.to_newick_vec()
    }

    fn to_newick_string(&self) -> String {
        self.inner.to_newick_string()
    }

    fn as_topology(&self) -> PyTopologyTreeView {
        PyTopologyTreeView {
            inner: self.inner.clone().as_topology(),
        }
    }

    fn __getitem__(&self, index: usize) -> PyResult<PySingleTree> {
        if index >= self.inner.ngenes() {
            return Err(PyIndexError::new_err(format!(
                "Tree index {} out of bounds. Valid range: 0..{}",
                index,
                self.inner.ngenes()
            )));
        }

        Ok(PySingleTree {
            taxon_set: Arc::clone(&self.inner.0.taxon_set),
            tree: Arc::new(self.inner.0.trees[index].clone()),
        })
    }

    #[getter]
    fn taxon_set(&self) -> PyTaxonSet {
        PyTaxonSet {
            inner: Arc::clone(&self.inner.0.taxon_set),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.ngenes()
    }

    fn __str__(&self) -> String {
        format!(
            "TreeSet with {} trees, {} taxa",
            self.inner.ngenes(),
            self.inner.ntaxa()
        )
    }

    fn __repr__(&self) -> String {
        format!("TreeSet('{}')", self.inner.to_newick_string())
    }

    fn newick(&self) -> String {
        self.inner.to_newick_string()
    }

    #[classmethod]
    fn from_file(_cls: Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let contents = std::fs::read_to_string(path)
            .map_err(|e| PyValueError::new_err(format!("Failed to read file: {}", e)))?;
        CompleteTreeView::from_newick_string(&contents)
            .map(|inner| PyCompleteTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    /// Relabels taxa according to the provided mapping
    fn remap(&self, mapping: HashMap<String, String>) -> PyResult<Self> {
        self.inner
            .remap(&mapping)
            .map(|inner| PyCompleteTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        bincode::serialize(&self.inner)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        self.inner = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        Ok(())
    }
}

#[pymethods]
impl PyTopologyTreeView {
    #[new]
    fn new(newick_str: &str) -> PyResult<Self> {
        TopologyTreeView::from_newick_string(newick_str)
            .map(|inner| PyTopologyTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    #[getter]
    fn ngenes(&self) -> usize {
        self.inner.ngenes()
    }

    #[getter]
    fn ntaxa(&self) -> usize {
        self.inner.ntaxa()
    }

    fn slice(&self, start: usize, end: usize) -> PyResult<Self> {
        Ok(PyTopologyTreeView {
            inner: self.inner.slice(start, end),
        })
    }

    fn get_distance_matrix(&self, index: usize) -> Option<PyDistanceMatrixView> {
        self.inner
            .get_distance_matrix(index)
            .map(|inner| PyDistanceMatrixView { inner })
    }

    fn get_taxon_name(&self, id: usize) -> Option<String> {
        self.inner.get_taxon_name(id).map(String::from)
    }

    fn to_newick_vec(&self) -> Vec<String> {
        self.inner.to_newick_vec()
    }

    fn to_newick_string(&self) -> String {
        self.inner.to_newick_string()
    }

    fn restriction(&self, taxa: Vec<String>) -> PyResult<Self> {
        let taxa_refs: Vec<&str> = taxa.iter().map(|s| s.as_str()).collect();
        self.inner
            .restriction(&taxa_refs)
            .map(|inner| PyTopologyTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    fn __getitem__(&self, index: usize) -> PyResult<PySingleTree> {
        if index >= self.inner.ngenes() {
            return Err(PyIndexError::new_err(format!(
                "Tree index {} out of bounds. Valid range: 0..{}",
                index,
                self.inner.ngenes()
            )));
        }

        Ok(PySingleTree {
            taxon_set: Arc::clone(&self.inner.0.taxon_set),
            tree: Arc::new(self.inner.0.trees[index].clone()),
        })
    }

    #[getter]
    fn taxon_set(&self) -> PyTaxonSet {
        PyTaxonSet {
            inner: Arc::clone(&self.inner.0.taxon_set),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.ngenes()
    }

    fn __str__(&self) -> String {
        format!(
            "TopologySet with {} trees, {} taxa",
            self.inner.ngenes(),
            self.inner.ntaxa()
        )
    }

    fn __repr__(&self) -> String {
        format!("TopologySet('{}')", self.inner.to_newick_string())
    }

    fn newick(&self) -> String {
        self.inner.to_newick_string()
    }

    #[classmethod]
    fn from_file(_cls: Bound<'_, PyType>, path: &str) -> PyResult<Self> {
        let contents = std::fs::read_to_string(path)
            .map_err(|e| PyValueError::new_err(format!("Failed to read file: {}", e)))?;
        TopologyTreeView::from_newick_string(&contents)
            .map(|inner| PyTopologyTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    /// Relabels taxa according to the provided mapping
    fn remap(&self, mapping: HashMap<String, String>) -> PyResult<Self> {
        self.inner
            .remap(&mapping)
            .map(|inner| PyTopologyTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        bincode::serialize(&self.inner)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        self.inner = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        Ok(())
    }
}

#[pymethods]
impl PyMSCTreeView {
    #[new]
    fn new(gene_trees: &str, species_tree: &str) -> PyResult<Self> {
        let mut collection = MSCTreeCollection::new(false);

        // Parse gene trees
        let gene_collection = TreeCollection::from_newick_string(gene_trees, false)
            .map_err(|e| PyValueError::new_err(e))?;
        collection.taxon_set = gene_collection.taxon_set;
        collection.gene_trees = gene_collection.trees;

        // Parse species tree directly using parse_newick
        let species_tree = parse_newick(&mut collection.taxon_set, species_tree);
        
        // Verify the species tree is valid
        if species_tree.ntaxa == 0 {
            return Err(PyValueError::new_err("No species tree found in input string"));
        }

        collection.species_tree = species_tree;

        Ok(PyMSCTreeView {
            inner: MSCTreeView::new(collection),
        })
    }

    #[getter]
    fn ngenes(&self) -> usize {
        self.inner.ngenes()
    }

    #[getter]
    fn ntaxa(&self) -> usize {
        self.inner.ntaxa()
    }

    fn slice(&self, start: usize, end: usize) -> PyResult<Self> {
        Ok(PyMSCTreeView {
            inner: self.inner.slice(start, end),
        })
    }

    fn get_distance_matrix(&self, index: usize) -> Option<PyDistanceMatrixView> {
        self.inner
            .get_distance_matrix(index)
            .map(|inner| PyDistanceMatrixView { inner })
    }

    fn get_species_distance_matrix(&self) -> PyDistanceMatrixView {
        PyDistanceMatrixView {
            inner: self.inner.get_species_distance_matrix(),
        }
    }

    fn get_taxon_name(&self, id: usize) -> Option<String> {
        self.inner.get_taxon_name(id).map(String::from)
    }

    fn to_newick_vec(&self) -> Vec<String> {
        self.inner.to_newick_vec()
    }

    fn to_newick_string(&self) -> String {
        self.inner.to_newick_string()
    }

    fn species_tree_to_newick(&self) -> String {
        self.inner.species_tree_to_newick()
    }

    fn restriction(&self, taxa: Vec<String>) -> PyResult<Self> {
        let taxa_refs: Vec<&str> = taxa.iter().map(|s| s.as_str()).collect();
        self.inner
            .restriction(&taxa_refs)
            .map(|inner| PyMSCTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    fn __getitem__(&self, index: usize) -> PyResult<PySingleTree> {
        if index >= self.inner.ngenes() {
            return Err(PyIndexError::new_err(format!(
                "Tree index {} out of bounds. Valid range: 0..{}",
                index,
                self.inner.ngenes()
            )));
        }

        let tree = self
            .inner
            .get_tree(index)
            .expect("Index was checked to be in bounds");

        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) => &base.taxon_set,
            TreeViewType::Topology(base) => &base.taxon_set,
        };

        Ok(PySingleTree {
            taxon_set: Arc::clone(taxon_set),
            tree: Arc::new(tree.clone()),
        })
    }

    fn get_species_tree(&self) -> PySingleTree {
        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) | TreeViewType::Topology(base) => &base.taxon_set,
        };

        PySingleTree {
            taxon_set: Arc::clone(taxon_set),
            tree: Arc::new(self.inner.species_tree.clone()),
        }
    }

    #[getter]
    fn taxon_set(&self) -> PyTaxonSet {
        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) => &base.taxon_set,
            TreeViewType::Topology(base) => &base.taxon_set,
        };
        PyTaxonSet {
            inner: Arc::clone(taxon_set),
        }
    }

    fn as_topology(&self) -> PyMSCTopologyTreeView {
        PyMSCTopologyTreeView {
            inner: self.inner.clone().as_topology(),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.ngenes()
    }

    fn __str__(&self) -> String {
        format!(
            "MSCTreeSet with {} gene trees, {} taxa",
            self.inner.ngenes(),
            self.inner.ntaxa()
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "MSCTreeSet('{}', '{}')",
            self.inner.to_newick_string(),
            self.inner.species_tree_to_newick()
        )
    }

    fn newick(&self) -> String {
        self.inner.to_newick_string()
    }

    #[classmethod]
    fn from_files(
        _cls: Bound<'_, PyType>,
        gene_trees_path: &str,
        species_tree_path: &str,
    ) -> PyResult<Self> {
        let gene_trees = std::fs::read_to_string(gene_trees_path)
            .map_err(|e| PyValueError::new_err(format!("Failed to read gene trees file: {}", e)))?;
        let species_tree = std::fs::read_to_string(species_tree_path).map_err(|e| {
            PyValueError::new_err(format!("Failed to read species tree file: {}", e))
        })?;

        PyMSCTreeView::new(&gene_trees, &species_tree)
    }

    /// Relabels taxa according to the provided mapping
    fn remap(&self, mapping: HashMap<String, String>) -> PyResult<Self> {
        self.inner
            .remap(&mapping)
            .map(|inner| PyMSCTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        bincode::serialize(&self.inner)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        self.inner = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        Ok(())
    }
}

#[pymethods]
impl PyDistanceMatrixView {
    fn get(&self, i: usize, j: usize) -> f64 {
        self.inner.get(i, j)
    }

    fn get_by_name(&self, taxon1: &str, taxon2: &str) -> PyResult<f64> {
        self.inner
            .get_by_name(taxon1, taxon2)
            .map_err(|e| PyValueError::new_err(e))
    }

    #[getter]
    fn ntaxa(&self) -> usize {
        self.inner.ntaxa()
    }

    fn get_taxon_name(&self, id: usize) -> Option<String> {
        self.inner.get_taxon_name(id).map(String::from)
    }

    fn __getitem__(&self, key: Bound<'_, PyTuple>) -> PyResult<f64> {
        if key.len() != 2 {
            return Err(PyValueError::new_err("Expected a tuple of two values"));
        }

        let item0 = key.get_item(0)?;
        let item1 = key.get_item(1)?;

        // If both items are strings, use get_by_name
        if item0.is_instance_of::<pyo3::types::PyString>()
            && item1.is_instance_of::<pyo3::types::PyString>()
        {
            let taxon1: String = item0.extract()?;
            let taxon2: String = item1.extract()?;
            return self
                .inner
                .get_by_name(&taxon1, &taxon2)
                .map_err(|e| PyValueError::new_err(e));
        }

        // If both items are integers, use get
        if item0.is_instance_of::<pyo3::types::PyInt>()
            && item1.is_instance_of::<pyo3::types::PyInt>()
        {
            let i: usize = item0.extract()?;
            let j: usize = item1.extract()?;

            if i >= self.inner.ntaxa() || j >= self.inner.ntaxa() {
                return Err(PyIndexError::new_err(format!(
                    "Index ({}, {}) out of bounds. Valid range: 0..{}",
                    i,
                    j,
                    self.inner.ntaxa()
                )));
            }

            return Ok(self.inner.get(i, j));
        }

        Err(PyValueError::new_err(
            "Index must be either (int, int) or (str, str)",
        ))
    }

    #[getter]
    fn taxon_set(&self) -> PyTaxonSet {
        PyTaxonSet {
            inner: Arc::clone(&self.inner.taxon_set),
        }
    }

    fn __str__(&self) -> String {
        let n = self.inner.ntaxa();
        let mut result = String::new();
        
        // Add header row with taxon names
        result.push_str("      "); // Padding for row labels
        for j in 0..n {
            if let Some(name) = self.inner.get_taxon_name(j) {
                result.push_str(&format!("{:>8}", name));
            }
        }
        result.push('\n');
        
        // Add matrix rows
        for i in 0..n {
            if let Some(name) = self.inner.get_taxon_name(i) {
                result.push_str(&format!("{:<6}", name));
            }
            for j in 0..n {
                let val = self.inner.get(i, j);
                // Check if the value is effectively an integer
                if (val - val.round()).abs() < f64::EPSILON {
                    result.push_str(&format!("{:8.0}", val));
                } else {
                    result.push_str(&format!("{:8.3}", val));
                }
            }
            result.push('\n');
        }
        
        result
    }

    fn __repr__(&self) -> String {
        format!("DistanceMatrix(size={})", self.inner.ntaxa())
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        bincode::serialize(&self.inner)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        self.inner = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        Ok(())
    }
}

#[pymethods]
impl PySingleTree {
    #[new]
    fn new(newick_str: &str) -> PyResult<Self> {
        let collection = TreeCollection::from_newick_string(newick_str, false)
            .map_err(|e| PyValueError::new_err(e))?;

        if collection.trees.is_empty() {
            return Err(PyValueError::new_err("No tree found in input string"));
        }

        Ok(PySingleTree {
            taxon_set: Arc::new(collection.taxon_set),
            tree: Arc::new(collection.trees[0].clone()),
        })
    }

    #[getter]
    fn ntaxa(&self) -> usize {
        self.taxon_set.len()
    }

    fn get_taxon_name(&self, id: usize) -> Option<String> {
        self.taxon_set.names.get(id).map(String::from)
    }

    fn get_distance_matrix(&self) -> PyDistanceMatrixView {
        PyDistanceMatrixView {
            inner: DistanceMatrixView {
                matrix: Arc::new(self.tree.distance_matrix()),
                taxon_set: Arc::clone(&self.taxon_set),
            },
        }
    }

    fn to_newick(&self) -> String {
        // Helper function to write tree recursively
        fn write_subtree(tree: &Tree, node: usize, taxon_set: &TaxonSet, result: &mut String) {
            if tree.is_leaf(node) {
                result.push_str(&taxon_set.names[tree.taxa[node] as usize]);
            } else {
                result.push('(');
                let mut first = true;
                for child in tree.children(node) {
                    if !first {
                        result.push(',');
                    }
                    first = false;
                    write_subtree(tree, child, taxon_set, result);
                }
                result.push(')');
            }
        }

        let mut result = String::new();
        write_subtree(&self.tree, 0, &self.taxon_set, &mut result);
        result.push(';');
        result
    }

    fn restriction(&self, taxa: Vec<String>) -> PyResult<Self> {
        // Convert taxa strings to a HashSet of taxon IDs
        let mut keep_taxa = HashSet::new();
        for taxon in &taxa {
            if let Some(id) = self.taxon_set.to_id.get(taxon.as_str()) {
                keep_taxa.insert(*id);
            } else {
                return Err(PyValueError::new_err(format!(
                    "Taxon '{}' not found in tree",
                    taxon
                )));
            }
        }

        // Create new taxon set and mapping
        let mut new_taxon_set = TaxonSet::new();
        let mut id_map = HashMap::new();
        for taxon in taxa {
            let old_id = self.taxon_set.to_id[&taxon];
            let new_id = new_taxon_set.request(taxon);
            id_map.insert(old_id, new_id);
        }

        // Restrict the tree
        let new_tree = self.tree.restrict(&keep_taxa, &id_map);

        Ok(PySingleTree {
            taxon_set: Arc::new(new_taxon_set),
            tree: Arc::new(new_tree),
        })
    }

    fn get_taxa(&self) -> Vec<String> {
        self.taxon_set.names.clone()
    }

    fn __str__(&self) -> String {
        self.to_newick()
    }

    fn __repr__(&self) -> String {
        format!("SingleTree('{}')", self.to_newick())
    }

    #[getter]
    fn taxon_set(&self) -> PyTaxonSet {
        PyTaxonSet {
            inner: Arc::clone(&self.taxon_set),
        }
    }

    fn __len__(&self) -> usize {
        self.tree.num_nodes()
    }

    /// Relabels taxa according to the provided mapping
    fn remap(&self, mapping: HashMap<String, String>) -> PyResult<Self> {
        // Create a TreeViewBase to use its remap functionality
        let view_base = TreeViewBase {
            taxon_set: Arc::clone(&self.taxon_set),
            trees: Arc::new([(*self.tree).clone()]),
            topology_only: false,
        };

        view_base.remap(&mapping).map(|new_base| PySingleTree {
            taxon_set: Arc::clone(&new_base.taxon_set),
            tree: Arc::new(new_base.trees[0].clone()),
        }).map_err(|e| PyValueError::new_err(e))
    }

    /// Sort the tree by a single criterion
    #[pyo3(signature = (criterion, ascending=None))]
    fn sort_by(&self, criterion: PySortBy, ascending: Option<bool>) -> PyResult<Self> {
        let order = ascending.unwrap_or(true);
        let criterion = SortCriterion {
            kind: match criterion {
                PySortBy::LexicographicalOrder => SortCriterionKind::LexicographicalOrder,
                PySortBy::ChildCount => SortCriterionKind::ChildCount,
                PySortBy::DescendantCount => SortCriterionKind::DescendantCount,
            },
            order: if order { SortOrder::Ascending } else { SortOrder::Descending },
        };

        let mut new_tree = (*self.tree).clone();
        new_tree.sort(&self.taxon_set, &[criterion]);

        Ok(PySingleTree {
            taxon_set: Arc::clone(&self.taxon_set),
            tree: Arc::new(new_tree),
        })
    }

    /// Sort the tree by multiple criteria
    fn sort_by_multiple(&self, criteria: Vec<(PySortBy, Option<bool>)>) -> PyResult<Self> {
        let rust_criteria: Vec<SortCriterion> = criteria
            .into_iter()
            .map(|(criterion, ascending)| {
                let order = ascending.unwrap_or(true);
                SortCriterion {
                    kind: match criterion {
                        PySortBy::LexicographicalOrder => SortCriterionKind::LexicographicalOrder,
                        PySortBy::ChildCount => SortCriterionKind::ChildCount,
                        PySortBy::DescendantCount => SortCriterionKind::DescendantCount,
                    },
                    order: if order { SortOrder::Ascending } else { SortOrder::Descending },
                }
            })
            .collect();

        let mut new_tree = (*self.tree).clone();
        new_tree.sort(&self.taxon_set, &rust_criteria);

        Ok(PySingleTree {
            taxon_set: Arc::clone(&self.taxon_set),
            tree: Arc::new(new_tree),
        })
    }

    /// Sort the tree according to the given criteria
    fn sort(&self, criteria: Vec<PySortCriterion>) -> PyResult<Self> {
        // Convert Python criteria to Rust criteria
        let rust_criteria: Vec<SortCriterion> = criteria
            .into_iter()
            .map(|c| SortCriterion {
                kind: match c.kind {
                    PySortCriterionKind::LexicographicalOrder => SortCriterionKind::LexicographicalOrder,
                    PySortCriterionKind::ChildCount => SortCriterionKind::ChildCount,
                    PySortCriterionKind::DescendantCount => SortCriterionKind::DescendantCount,
                },
                order: match c.order {
                    PySortOrder::Ascending => SortOrder::Ascending,
                    PySortOrder::Descending => SortOrder::Descending,
                },
            })
            .collect();

        // Create a new tree and sort it
        let mut new_tree = (*self.tree).clone();
        new_tree.sort(&self.taxon_set, &rust_criteria);

        // Return new PySingleTree with sorted tree
        Ok(PySingleTree {
            taxon_set: Arc::clone(&self.taxon_set),
            tree: Arc::new(new_tree),
        })
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        let state = (Arc::as_ref(&self.taxon_set), Arc::as_ref(&self.tree));
        bincode::serialize(&state)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        let (taxon_set, tree): (TaxonSet, Tree) = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        self.taxon_set = Arc::new(taxon_set);
        self.tree = Arc::new(tree);
        Ok(())
    }

    fn newick(&self) -> String {
        self.to_newick()
    }
}

#[pymethods]
impl PyTaxonSet {
    #[getter]
    fn len(&self) -> usize {
        self.inner.len()
    }

    fn get_id(&self, taxon_name: &str) -> Option<usize> {
        self.inner.to_id.get(taxon_name).copied()
    }

    fn get_name(&self, id: usize) -> Option<String> {
        if id < self.inner.names.len() {
            Some(self.inner.names[id].clone())
        } else {
            None
        }
    }

    fn names(&self) -> Vec<String> {
        self.inner.names.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __str__(&self) -> String {
        format!("TaxonSet with {} taxa", self.inner.len())
    }

    fn __repr__(&self) -> String {
        format!("TaxonSet({})", self.inner.names.join(", "))
    }

    fn __contains__(&self, item: &str) -> bool {
        self.inner.to_id.contains_key(item)
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        bincode::serialize(Arc::as_ref(&self.inner))
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        let taxon_set: TaxonSet = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        self.inner = Arc::new(taxon_set);
        Ok(())
    }
}

#[pymethods]
impl PyMSCTopologyTreeView {
    #[new]
    fn new(gene_trees: &str, species_tree: &str) -> PyResult<Self> {
        let mut collection = MSCTreeCollection::new(true); // Note: true for topology-only

        // Parse gene trees
        let gene_collection = TreeCollection::from_newick_string(gene_trees, true)
            .map_err(|e| PyValueError::new_err(e))?;
        collection.taxon_set = gene_collection.taxon_set;
        collection.gene_trees = gene_collection.trees;

        // Parse species tree directly
        let species_tree = parse_newick(&mut collection.taxon_set, species_tree);
        
        // Verify the species tree is valid
        if species_tree.ntaxa == 0 {
            return Err(PyValueError::new_err("No species tree found in input string"));
        }

        collection.species_tree = species_tree;

        Ok(PyMSCTopologyTreeView {
            inner: MSCTreeView::new(collection).as_topology(),
        })
    }

    #[getter]
    fn ngenes(&self) -> usize {
        self.inner.ngenes()
    }

    #[getter]
    fn ntaxa(&self) -> usize {
        self.inner.ntaxa()
    }

    fn slice(&self, start: usize, end: usize) -> PyResult<Self> {
        Ok(PyMSCTopologyTreeView {
            inner: self.inner.slice(start, end),
        })
    }

    fn get_distance_matrix(&self, index: usize) -> Option<PyDistanceMatrixView> {
        self.inner
            .get_distance_matrix(index)
            .map(|inner| PyDistanceMatrixView { inner })
    }

    fn get_species_distance_matrix(&self) -> PyDistanceMatrixView {
        PyDistanceMatrixView {
            inner: self.inner.get_species_distance_matrix(),
        }
    }

    fn get_taxon_name(&self, id: usize) -> Option<String> {
        self.inner.get_taxon_name(id).map(String::from)
    }

    fn to_newick_vec(&self) -> Vec<String> {
        self.inner.to_newick_vec()
    }

    fn to_newick_string(&self) -> String {
        self.inner.to_newick_string()
    }

    fn species_tree_to_newick(&self) -> String {
        self.inner.species_tree_to_newick()
    }

    fn restriction(&self, taxa: Vec<String>) -> PyResult<Self> {
        let taxa_refs: Vec<&str> = taxa.iter().map(|s| s.as_str()).collect();
        self.inner
            .restriction(&taxa_refs)
            .map(|inner| PyMSCTopologyTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    fn __getitem__(&self, index: usize) -> PyResult<PySingleTree> {
        self.inner.get_tree(index)
            .ok_or_else(|| PyIndexError::new_err("Index out of bounds"))
            .map(|tree| PySingleTree {
                taxon_set: Arc::new(self.inner.taxon_set().clone()),
                tree: Arc::new(tree.clone()),
            })
    }

    fn get_species_tree(&self) -> PySingleTree {
        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) | TreeViewType::Topology(base) => &base.taxon_set,
        };

        PySingleTree {
            taxon_set: Arc::clone(taxon_set),
            tree: Arc::new(self.inner.species_tree.clone()),
        }
    }

    #[getter]
    fn taxon_set(&self) -> PyTaxonSet {
        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) => &base.taxon_set,
            TreeViewType::Topology(base) => &base.taxon_set,
        };
        PyTaxonSet {
            inner: Arc::clone(taxon_set),
        }
    }

    fn __len__(&self) -> usize {
        self.inner.ngenes()
    }

    fn __str__(&self) -> String {
        format!(
            "MSCTopologySet with {} gene trees, {} taxa",
            self.inner.ngenes(),
            self.inner.ntaxa()
        )
    }

    fn __repr__(&self) -> String {
        format!(
            "MSCTopologySet('{}', '{}')",
            self.inner.to_newick_string(),
            self.inner.species_tree_to_newick()
        )
    }

    fn newick(&self) -> String {
        self.inner.to_newick_string()
    }

    #[classmethod]
    fn from_files(
        _cls: Bound<'_, PyType>,
        gene_trees_path: &str,
        species_tree_path: &str,
    ) -> PyResult<Self> {
        let gene_trees = std::fs::read_to_string(gene_trees_path)
            .map_err(|e| PyValueError::new_err(format!("Failed to read gene trees file: {}", e)))?;
        let species_tree = std::fs::read_to_string(species_tree_path).map_err(|e| {
            PyValueError::new_err(format!("Failed to read species tree file: {}", e))
        })?;

        let view = PyMSCTopologyTreeView::new(&gene_trees, &species_tree)?;
        Ok(view)
    }

    /// Relabels taxa according to the provided mapping
    fn remap(&self, mapping: HashMap<String, String>) -> PyResult<Self> {
        self.inner
            .remap(&mapping)
            .map(|inner| PyMSCTopologyTreeView { inner })
            .map_err(|e| PyValueError::new_err(e))
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        bincode::serialize(&self.inner)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        self.inner = bincode::deserialize(state)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        Ok(())
    }
}

#[pymethods]
impl PyFamilyOfMSC {
    #[new]
    fn new() -> Self {
        PyFamilyOfMSC { views: Vec::new() }
    }

    fn add(&mut self, view: PyMSCTopologyTreeView) {
        self.views.push(Arc::new(view));
    }

    fn __len__(&self) -> usize {
        self.views.len()
    }

    fn __str__(&self) -> String {
        format!("FamilyOfMSC with {} MSC topology sets", self.views.len())
    }

    fn __repr__(&self) -> String {
        format!("FamilyOfMSC with {} members", self.views.len())
    }

    #[classmethod]
    fn from_files(_cls: Bound<'_, PyType>, gene_trees_paths: Vec<String>, species_trees_paths: Vec<String>) -> PyResult<Self> {
        if gene_trees_paths.len() != species_trees_paths.len() {
            return Err(PyValueError::new_err(
                "Number of gene trees files must match number of species trees files",
            ));
        }

        let mut family = PyFamilyOfMSC::new();
        for (gene_path, species_path) in gene_trees_paths.iter().zip(species_trees_paths.iter()) {
            let view = PyMSCTopologyTreeView::from_files(_cls.clone(), gene_path, species_path)?;
            family.add(view);
        }
        Ok(family)
    }

    fn get_sizes(&self) -> Vec<(usize, usize)> {
        self.views
            .iter()
            .map(|view| (view.ngenes(), view.ntaxa()))
            .collect()
    }

    fn __getitem__(&self, index: usize) -> PyResult<PyMSCTopologyTreeView> {
        self.views
            .get(index)
            .map(|view| (**view).clone())
            .ok_or_else(|| PyIndexError::new_err("Index out of bounds"))
    }

    pub fn __getstate__(&self) -> PyResult<Vec<u8>> {
        // Convert Arc<PyMSCTopologyTreeView> to Vec<MSCTreeView>
        let views: Vec<MSCTreeView> = self.views.iter()
            .map(|view| view.inner.clone())
            .collect();
        
        // First serialize with bincode
        let serialized = bincode::serialize(&views)
            .map_err(|e| PyValueError::new_err(format!("Serialization error: {}", e)))?;
        
        // Then compress with zstd
        zstd::encode_all(serialized.as_slice(), 3) // compression level 3 for good balance
            .map_err(|e| PyValueError::new_err(format!("Compression error: {}", e)))
    }

    pub fn __setstate__(&mut self, state: &[u8]) -> PyResult<()> {
        // First decompress with zstd
        let decompressed = zstd::decode_all(state)
            .map_err(|e| PyValueError::new_err(format!("Decompression error: {}", e)))?;
        
        // Then deserialize with bincode
        let views: Vec<MSCTreeView> = bincode::deserialize(&decompressed)
            .map_err(|e| PyValueError::new_err(format!("Deserialization error: {}", e)))?;
        
        // Convert back to Arc<PyMSCTopologyTreeView>
        self.views = views.into_iter()
            .map(|view| Arc::new(PyMSCTopologyTreeView { inner: view }))
            .collect();
        Ok(())
    }
}

#[pymethods]
impl PySortCriterionKind {
    #[classmethod]
    fn lexicographical(_cls: Bound<'_, PyType>) -> Self {
        Self::LexicographicalOrder
    }

    #[classmethod]
    fn child_count(_cls: Bound<'_, PyType>) -> Self {
        Self::ChildCount
    }

    #[classmethod]
    fn descendant_count(_cls: Bound<'_, PyType>) -> Self {
        Self::DescendantCount
    }
}

#[pymethods]
impl PySortOrder {
    #[classmethod]
    fn ascending(_cls: Bound<'_, PyType>) -> Self {
        Self::Ascending
    }

    #[classmethod]
    fn descending(_cls: Bound<'_, PyType>) -> Self {
        Self::Descending
    }
}

#[pymethods]
impl PySortCriterion {
    #[new]
    fn new(kind: PySortCriterionKind, order: PySortOrder) -> Self {
        Self { kind, order }
    }
}

#[pymethods]
impl PySortBy {
    fn __str__(&self) -> String {
        match self {
            Self::LexicographicalOrder => "lexicographical".to_string(),
            Self::ChildCount => "child_count".to_string(),
            Self::DescendantCount => "descendant_count".to_string(),
        }
    }
}

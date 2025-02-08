use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::sync::Arc;

use crate::porcelain::{
    CompleteTreeView, DistanceMatrixView, MSCTreeView, TopologyTreeView, TreeViewBase, TreeViewType,
};
use crate::tree::{MSCTreeCollection, TaxonSet, Tree, TreeCollection};

/// A Python wrapper for CompleteTreeView
#[pyclass(name = "TreeSet")]
#[derive(Clone)]
pub struct PyCompleteTreeView {
    inner: CompleteTreeView,
}

/// A Python wrapper for TopologyTreeView
#[pyclass(name = "TopologySet")]
#[derive(Clone)]
pub struct PyTopologyTreeView {
    inner: TopologyTreeView,
}

/// A Python wrapper for MSCTreeView
#[pyclass(name = "MSCTreeSet")]
#[derive(Clone)]
pub struct PyMSCTreeView {
    inner: MSCTreeView,
}

/// A Python wrapper for DistanceMatrixView
#[pyclass(name = "DistanceMatrix")]
#[derive(Clone)]
pub struct PyDistanceMatrixView {
    inner: DistanceMatrixView,
}

/// A Python wrapper for a single tree with its taxon set
#[pyclass(name = "Tree")]
#[derive(Clone)]
pub struct PySingleTree {
    taxon_set: Arc<TaxonSet>,
    tree: Arc<Tree>,
}

/// A Python wrapper for TaxonSet providing read-only access
#[pyclass(name = "TaxonSet")]
#[derive(Clone)]
pub struct PyTaxonSet {
    inner: Arc<TaxonSet>,
}

/// A Python wrapper for MSCTreeView in topology-only mode
#[pyclass(name = "MSCTopologySet")]
#[derive(Clone)]
pub struct PyMSCTopologyTreeView {
    inner: MSCTreeView,
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
            return Err(PyValueError::new_err("Tree index out of bounds"));
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
            return Err(PyValueError::new_err("Tree index out of bounds"));
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

        // Parse species tree
        let species_collection = TreeCollection::from_newick_string(species_tree, false)
            .map_err(|e| PyValueError::new_err(e))?;
        if !species_collection.trees.is_empty() {
            collection.species_tree = species_collection.trees[0].clone();
        } else {
            return Err(PyValueError::new_err(
                "No species tree found in input string",
            ));
        }

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
            return Err(PyValueError::new_err("Tree index out of bounds"));
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
        let tree = self
            .inner
            .species_tree
            .get_tree(0)
            .expect("Species tree should always have exactly one tree");

        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) => &base.taxon_set,
            TreeViewType::Topology(base) => &base.taxon_set,
        };

        PySingleTree {
            taxon_set: Arc::clone(taxon_set),
            tree: Arc::new(tree.clone()),
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
                return Err(PyValueError::new_err("Index out of bounds"));
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
                    if tree.lengths[child] >= 0.0 {
                        result.push(':');
                        result.push_str(&tree.lengths[child].to_string());
                    }
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

    fn get_names(&self) -> Vec<String> {
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
}

#[pymethods]
impl PyMSCTopologyTreeView {
    #[new]
    fn new(gene_trees: &str, species_tree: &str) -> PyResult<Self> {
        let msc_view = PyMSCTreeView::new(gene_trees, species_tree)?;
        Ok(PyMSCTopologyTreeView {
            inner: msc_view.inner.as_topology(),
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
        if index >= self.inner.ngenes() {
            return Err(PyValueError::new_err("Tree index out of bounds"));
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
        let tree = self
            .inner
            .species_tree
            .get_tree(0)
            .expect("Species tree should always have exactly one tree");

        let taxon_set = match &self.inner.base {
            TreeViewType::Complete(base) => &base.taxon_set,
            TreeViewType::Topology(base) => &base.taxon_set,
        };

        PySingleTree {
            taxon_set: Arc::clone(taxon_set),
            tree: Arc::new(tree.clone()),
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
}

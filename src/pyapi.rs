use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyList, PyTuple};
use std::sync::Arc;

use crate::porcelain::{CompleteTreeView, DistanceMatrixView, MSCTreeView, TopologyTreeView};
use crate::tree::{MSCTreeCollection, TreeCollection};

/// A Python wrapper for CompleteTreeView
#[pyclass(name = "CompleteTreeView")]
#[derive(Clone)]
pub struct PyCompleteTreeView {
    inner: CompleteTreeView,
}

/// A Python wrapper for TopologyTreeView
#[pyclass(name = "TopologyTreeView")]
#[derive(Clone)]
pub struct PyTopologyTreeView {
    inner: TopologyTreeView,
}

/// A Python wrapper for MSCTreeView
#[pyclass(name = "MSCTreeView")]
#[derive(Clone)]
pub struct PyMSCTreeView {
    inner: MSCTreeView,
}

/// A Python wrapper for DistanceMatrixView
#[pyclass(name = "DistanceMatrixView")]
#[derive(Clone)]
pub struct PyDistanceMatrixView {
    inner: DistanceMatrixView,
}

#[pymethods]
impl PyCompleteTreeView {
    #[new]
    fn new(newick_str: &str) -> PyResult<Self> {
        CompleteTreeView::from_newick_string(newick_str)
            .map(|inner| PyCompleteTreeView { inner })
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
}

#[pymethods]
impl PyMSCTreeView {
    #[new]
    fn new(gene_trees: &str, species_tree: &str) -> PyResult<Self> {
        let mut collection = MSCTreeCollection::new();

        // Parse gene trees
        let gene_collection =
            TreeCollection::from_newick_string(gene_trees).map_err(|e| PyValueError::new_err(e))?;
        collection.taxon_set = gene_collection.taxon_set;
        collection.gene_trees = gene_collection.trees;

        // Parse species tree
        let species_collection = TreeCollection::from_newick_string(species_tree)
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
}

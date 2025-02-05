use pyo3::prelude::*;

mod porcelain;
mod pyapi;
mod tree;

use pyapi::{PyCompleteTreeView, PyDistanceMatrixView, PyMSCTreeView, PyTopologyTreeView};

/// A Python module implemented in Rust.
#[pymodule]
fn teedeelee(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Add the Python classes to the module
    m.add_class::<PyCompleteTreeView>()?;
    m.add_class::<PyTopologyTreeView>()?;
    m.add_class::<PyMSCTreeView>()?;
    m.add_class::<PyDistanceMatrixView>()?;

    Ok(())
}

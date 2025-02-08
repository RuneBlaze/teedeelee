use crate::tree::{
    DistanceMatrix, MSCTreeCollection, TaxonSet, Tree, TreeCollection, TreeCollectionTrait,
};
use either::Either;
use std::cmp::{max, min};
use std::collections::{HashMap, HashSet};
use std::ops::Index;
use std::sync::Arc;

// Common functionality shared between Complete and Topology views
#[derive(Clone, Debug)]
pub(crate) struct TreeViewBase {
    pub(crate) taxon_set: Arc<TaxonSet>,
    pub(crate) trees: Arc<[Tree]>,
}

/// A thread-safe, shareable view of a tree collection with complete information
/// (including branch lengths and support values).
#[derive(Clone)]
pub struct CompleteTreeView(pub(crate) TreeViewBase);

/// A thread-safe, shareable view of a tree collection with topology-only information
#[derive(Clone)]
pub struct TopologyTreeView(pub(crate) TreeViewBase);

/// A thread-safe, shareable view of a multi-species coalescent tree collection
#[derive(Clone)]
pub struct MSCTreeView {
    pub(crate) base: TreeViewBase,
    pub(crate) species_tree: Arc<Tree>,
}

/// A thread-safe, shareable view of a distance matrix
#[derive(Clone)]
pub struct DistanceMatrixView {
    pub(crate) matrix: Arc<DistanceMatrix>,
    pub(crate) taxon_set: Arc<TaxonSet>,
}

impl DistanceMatrixView {
    /// Get distance between taxa by their indices
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.matrix.get(i, j)
    }

    /// Get distance between taxa by their names
    pub fn get_by_name(&self, taxon1: &str, taxon2: &str) -> Result<f64, String> {
        let i = self
            .taxon_set
            .to_id
            .get(taxon1)
            .ok_or_else(|| format!("Taxon '{}' not found", taxon1))?;
        let j = self
            .taxon_set
            .to_id
            .get(taxon2)
            .ok_or_else(|| format!("Taxon '{}' not found", taxon2))?;
        Ok(self.matrix.get(*i, *j))
    }

    /// Get the number of taxa
    pub fn ntaxa(&self) -> usize {
        self.matrix.size
    }

    /// Get the name of a taxon by its index
    pub fn get_taxon_name(&self, id: usize) -> Option<&str> {
        self.taxon_set.names.get(id).map(|s| s.as_str())
    }
}

// Add implementation for indexing with tuples
impl Index<(usize, usize)> for DistanceMatrixView {
    type Output = f64;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        static ZERO: f64 = 0.0;
        if i == j {
            &ZERO
        } else {
            &self.matrix.distances[self.matrix.get_index(min(i, j), max(i, j))]
        }
    }
}

// Modify the macro to use TreeCollectionTrait
macro_rules! impl_common_methods {
    ($type:ty) => {
        impl $type {
            /// Returns the number of trees in the collection
            pub fn ngenes(&self) -> usize {
                self.0.trees.len()
            }

            /// Returns the number of taxa in the collection
            pub fn ntaxa(&self) -> usize {
                self.0.taxon_set.len()
            }

            /// Returns a slice of trees as a new TreeView
            pub fn slice(&self, start: usize, end: usize) -> Self {
                Self(TreeViewBase {
                    taxon_set: Arc::clone(&self.0.taxon_set),
                    trees: Arc::from(&self.0.trees[start..end]),
                })
            }

            /// Returns an iterator over the trees
            pub fn iter(&self) -> impl Iterator<Item = &Tree> {
                self.0.trees.iter()
            }

            /// Gets a reference to a specific tree
            pub fn get_tree(&self, index: usize) -> Option<&Tree> {
                self.0.trees.get(index)
            }

            /// Gets the taxon name for a given ID
            pub fn get_taxon_name(&self, id: usize) -> Option<&str> {
                self.0.taxon_set.names.get(id).map(|s| s.as_str())
            }

            /// Convert the tree collection to a vector of Newick strings
            pub fn to_newick_vec(&self) -> Vec<String> {
                let collection = TreeCollectionView {
                    taxon_set: &self.0.taxon_set,
                    trees: &self.0.trees,
                };
                collection.to_newick_vec()
            }

            /// Convert the tree collection to a single string with semicolon+newline-separated Newick trees
            pub fn to_newick_string(&self) -> String {
                let collection = TreeCollectionView {
                    taxon_set: &self.0.taxon_set,
                    trees: &self.0.trees,
                };
                collection.to_newick_string()
            }

            /// Returns the distance matrix for the i-th tree
            pub fn get_distance_matrix(&self, index: usize) -> Option<DistanceMatrixView> {
                self.0.trees.get(index).map(|tree| DistanceMatrixView {
                    matrix: Arc::new(tree.distance_matrix()),
                    taxon_set: Arc::clone(&self.0.taxon_set),
                })
            }
        }
    };
}

// Add a helper struct to implement TreeCollectionTrait for views
struct TreeCollectionView<'a> {
    taxon_set: &'a TaxonSet,
    trees: &'a [Tree],
}

impl<'a> TreeCollectionTrait for TreeCollectionView<'a> {
    fn taxon_set(&self) -> &TaxonSet {
        self.taxon_set
    }

    fn trees(&self) -> &[Tree] {
        self.trees
    }
}

impl_common_methods!(CompleteTreeView);
impl_common_methods!(TopologyTreeView);

impl CompleteTreeView {
    /// Creates a new CompleteTreeView from a TreeCollection
    pub fn new(collection: TreeCollection) -> Self {
        CompleteTreeView(TreeViewBase {
            taxon_set: Arc::new(collection.taxon_set),
            trees: Arc::from(collection.trees),
        })
    }

    /// Converts this view to a topology-only view
    pub fn as_topology(self) -> TopologyTreeView {
        TopologyTreeView(TreeViewBase {
            taxon_set: self.0.taxon_set,
            trees: self.0.trees,
        })
    }

    /// Creates a new CompleteTreeView from a string containing semicolon-separated Newick trees
    pub fn from_newick_string(newick_str: &str) -> Result<Self, String> {
        TreeCollection::from_newick_string(newick_str)
            .map(Self::new)
            .map_err(|e| e.to_string())
    }
}

impl TopologyTreeView {
    /// Creates a new TopologyTreeView from a TreeCollection
    pub fn new(collection: TreeCollection) -> Self {
        TopologyTreeView(TreeViewBase {
            taxon_set: Arc::new(collection.taxon_set),
            trees: Arc::from(collection.trees),
        })
    }

    /// Creates a new topology-only TreeView containing only the specified taxa
    pub fn restriction(&self, taxa: &[&str]) -> Result<Self, String> {
        // Convert taxa strings to a HashSet of taxon IDs
        let mut keep_taxa = HashSet::new();
        for taxon in taxa {
            if let Some(id) = self.0.taxon_set.to_id.get(*taxon) {
                keep_taxa.insert(*id);
            } else {
                return Err(format!("Taxon '{}' not found in tree", taxon));
            }
        }

        // Create new taxon set and mapping
        let mut new_taxon_set = TaxonSet::new();
        let mut id_map = HashMap::new();
        for taxon in taxa {
            let old_id = self.0.taxon_set.to_id[*taxon];
            let new_id = new_taxon_set.request(taxon.to_string());
            id_map.insert(old_id, new_id);
        }

        // Restrict all trees
        let new_trees: Vec<Tree> = self
            .0
            .trees
            .iter()
            .map(|tree| tree.restrict(&keep_taxa, &id_map))
            .collect();

        Ok(TopologyTreeView(TreeViewBase {
            taxon_set: Arc::new(new_taxon_set),
            trees: Arc::from(new_trees),
        }))
    }

    /// Creates a new TopologyTreeView from a string containing semicolon-separated Newick trees
    pub fn from_newick_string(newick_str: &str) -> Result<Self, String> {
        TreeCollection::from_newick_string(newick_str)
            .map(Self::new)
            .map_err(|e| e.to_string())
    }
}

// Implement common methods for MSCTreeView
impl MSCTreeView {
    /// Creates a new MSCTreeView from an MSCTreeCollection
    pub fn new(collection: MSCTreeCollection) -> Self {
        MSCTreeView {
            base: TreeViewBase {
                taxon_set: Arc::new(collection.taxon_set),
                trees: Arc::from(collection.gene_trees),
            },
            species_tree: Arc::new(collection.species_tree),
        }
    }

    /// Returns a reference to the species tree
    pub fn species_tree(&self) -> &Tree {
        &self.species_tree
    }

    // Implement the common methods directly since we can't use the macro
    /// Returns the number of gene trees in the collection
    pub fn ngenes(&self) -> usize {
        self.base.trees.len()
    }

    /// Returns the number of taxa in the collection
    pub fn ntaxa(&self) -> usize {
        self.base.taxon_set.len()
    }

    /// Returns a slice of gene trees as a new MSCTreeView
    pub fn slice(&self, start: usize, end: usize) -> Self {
        MSCTreeView {
            base: TreeViewBase {
                taxon_set: Arc::clone(&self.base.taxon_set),
                trees: Arc::from(&self.base.trees[start..end]),
            },
            species_tree: Arc::clone(&self.species_tree),
        }
    }

    /// Returns an iterator over the gene trees
    pub fn iter(&self) -> impl Iterator<Item = &Tree> {
        self.base.trees.iter()
    }

    /// Gets a reference to a specific gene tree
    pub fn get_tree(&self, index: usize) -> Option<&Tree> {
        self.base.trees.get(index)
    }

    /// Gets the taxon name for a given ID
    pub fn get_taxon_name(&self, id: usize) -> Option<&str> {
        self.base.taxon_set.names.get(id).map(|s| s.as_str())
    }

    /// Convert the gene trees to a vector of Newick strings
    pub fn to_newick_vec(&self) -> Vec<String> {
        let collection = TreeCollectionView {
            taxon_set: &self.base.taxon_set,
            trees: &self.base.trees,
        };
        collection.to_newick_vec()
    }

    /// Convert the gene trees to a single string with semicolon+newline-separated Newick trees
    pub fn to_newick_string(&self) -> String {
        let collection = TreeCollectionView {
            taxon_set: &self.base.taxon_set,
            trees: &self.base.trees,
        };
        collection.to_newick_string()
    }

    /// Convert the species tree to Newick format
    pub fn species_tree_to_newick(&self) -> String {
        let collection = TreeCollectionView {
            taxon_set: &self.base.taxon_set,
            trees: std::slice::from_ref(&self.species_tree),
        };
        collection.tree_to_newick(&self.species_tree)
    }

    /// Creates a new MSCTreeView containing only the specified taxa
    pub fn restriction(&self, taxa: &[&str]) -> Result<Self, String> {
        // Convert taxa strings to a HashSet of taxon IDs
        let mut keep_taxa = HashSet::new();
        for taxon in taxa {
            if let Some(id) = self.base.taxon_set.to_id.get(*taxon) {
                keep_taxa.insert(*id);
            } else {
                return Err(format!("Taxon '{}' not found in tree", taxon));
            }
        }

        // Create new taxon set and mapping
        let mut new_taxon_set = TaxonSet::new();
        let mut id_map = HashMap::new();
        for taxon in taxa {
            let old_id = self.base.taxon_set.to_id[*taxon];
            let new_id = new_taxon_set.request(taxon.to_string());
            id_map.insert(old_id, new_id);
        }

        // Restrict all trees
        let new_gene_trees: Vec<Tree> = self
            .base
            .trees
            .iter()
            .map(|tree| tree.restrict(&keep_taxa, &id_map))
            .collect();

        // Restrict species tree
        let new_species_tree = self.species_tree.restrict(&keep_taxa, &id_map);

        Ok(MSCTreeView {
            base: TreeViewBase {
                taxon_set: Arc::new(new_taxon_set),
                trees: Arc::from(new_gene_trees),
            },
            species_tree: Arc::new(new_species_tree),
        })
    }

    /// Returns the distance matrix for the i-th gene tree
    pub fn get_distance_matrix(&self, index: usize) -> Option<DistanceMatrixView> {
        self.base.trees.get(index).map(|tree| DistanceMatrixView {
            matrix: Arc::new(tree.distance_matrix()),
            taxon_set: Arc::clone(&self.base.taxon_set),
        })
    }

    /// Returns the distance matrix for the species tree
    pub fn get_species_distance_matrix(&self) -> DistanceMatrixView {
        DistanceMatrixView {
            matrix: Arc::new(self.species_tree.distance_matrix()),
            taxon_set: Arc::clone(&self.base.taxon_set),
        }
    }
}

/// A builder for creating TreeViews with specific configurations
pub struct TreeViewBuilder {
    topology_only: bool,
    msc: bool, // New field to indicate MSC tree collection
}

impl TreeViewBuilder {
    /// Creates a new TreeViewBuilder with default settings
    pub fn new() -> Self {
        TreeViewBuilder {
            topology_only: false,
            msc: false,
        }
    }

    /// Sets whether the view should be topology-only
    pub fn topology_only(mut self, topology_only: bool) -> Self {
        self.topology_only = topology_only;
        self
    }

    /// Sets whether to build an MSC tree view
    pub fn msc(mut self, msc: bool) -> Self {
        self.msc = msc;
        self
    }

    /// Builds a TreeView from a TreeCollection or MSCTreeCollection
    pub fn build<T>(
        self,
        collection: T,
    ) -> Result<Either<Either<CompleteTreeView, TopologyTreeView>, MSCTreeView>, String>
    where
        T: Into<Either<TreeCollection, MSCTreeCollection>>,
    {
        match collection.into() {
            Either::Left(tree_collection) => {
                if self.msc {
                    Err("Cannot create MSCTreeView from regular TreeCollection".to_string())
                } else if self.topology_only {
                    Ok(Either::Left(Either::Right(TopologyTreeView::new(
                        tree_collection,
                    ))))
                } else {
                    Ok(Either::Left(Either::Left(CompleteTreeView::new(
                        tree_collection,
                    ))))
                }
            }
            Either::Right(msc_collection) => {
                if !self.msc {
                    Err("Expected regular TreeCollection but got MSCTreeCollection".to_string())
                } else {
                    Ok(Either::Right(MSCTreeView::new(msc_collection)))
                }
            }
        }
    }
}

impl Default for TreeViewBuilder {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::{TreeCollection, MSCTreeCollection};

    #[test]
    fn test_complete_tree_view_basics() {
        let newick = "(A:1.0,B:2.0,(C:1.5,D:0.5):2.0);";
        let view = CompleteTreeView::from_newick_string(newick).unwrap();
        
        assert_eq!(view.ngenes(), 1);
        assert_eq!(view.ntaxa(), 4);
        assert!(view.get_tree(0).is_some());
        assert!(view.get_tree(1).is_none());
        
        // Test taxon names
        assert_eq!(view.get_taxon_name(0), Some("A"));
        assert_eq!(view.get_taxon_name(1), Some("B"));
        assert_eq!(view.get_taxon_name(4), None);
    }

    #[test]
    fn test_topology_tree_view_restriction() {
        let newick = "(A,B,(C,D));(A,C,(B,D));";
        let view = TopologyTreeView::from_newick_string(newick).unwrap();
        
        assert_eq!(view.ngenes(), 2);
        
        // Restrict to taxa A, B
        let restricted = view.restriction(&["A", "B"]).unwrap();
        assert_eq!(restricted.ntaxa(), 2);
        assert_eq!(restricted.ngenes(), 2);
        
        // Test invalid taxon
        assert!(view.restriction(&["A", "X"]).is_err());
    }

    #[test]
    fn test_distance_matrix_view() {
        let newick = "(A:1.0,B:2.0,(C:1.5,D:0.5):2.0);";
        let view = CompleteTreeView::from_newick_string(newick).unwrap();
        
        let matrix = view.get_distance_matrix(0).unwrap();
        
        // Test matrix access
        assert_eq!(matrix.get(0, 0), 0.0); // Distance to self is 0
        assert!(matrix.get_by_name("A", "B").is_ok());
        assert!(matrix.get_by_name("A", "X").is_err());
        
        // Test indexing
        assert_eq!(matrix[(0, 1)], matrix.get(0, 1));
        assert_eq!(matrix[(1, 0)], matrix.get(0, 1)); // Symmetric
    }

    #[test]
    fn test_msc_tree_view() {
        let mut collection = MSCTreeCollection::new();
        let gene_trees = TreeCollection::from_newick_string("(A,B,(C,D));(A,C,(B,D));").unwrap();
        let species_tree = gene_trees.trees[0].clone();
        collection.gene_trees = gene_trees.trees;
        collection.taxon_set = gene_trees.taxon_set;
        collection.species_tree = species_tree;
        let view = MSCTreeView::new(collection);
        
        assert_eq!(view.ngenes(), 2);
        assert_eq!(view.ntaxa(), 4);
        
        // Test restriction
        let restricted = view.restriction(&["A", "B"]).unwrap();
        assert_eq!(restricted.ntaxa(), 2);
        assert_eq!(restricted.ngenes(), 2);
    }

    #[test]
    fn test_tree_view_slicing() {
        let newick = "(A,B);(C,D);(E,F);";
        let view = CompleteTreeView::from_newick_string(newick).unwrap();
        let view_clone = view.clone();  // Clone before moving
        
        let sliced = view.slice(1, 3);
        assert_eq!(sliced.ngenes(), 2);
        assert_eq!(sliced.ntaxa(), view_clone.ntaxa());
        
        let topology_view = view_clone.as_topology();
        assert_eq!(topology_view.ngenes(), sliced.ngenes());
    }

    #[test]
    fn test_newick_conversion() {
        let newick = "(A:1.0,B:2.0);(C:1.5,D:0.5);";
        let view = CompleteTreeView::from_newick_string(newick).unwrap();
        
        let newick_vec = view.to_newick_vec();
        assert_eq!(newick_vec.len(), 2);
        
        let newick_string = view.to_newick_string();
        assert!(newick_string.contains("A:1.0"));
        assert!(newick_string.contains("D:0.5"));
    }
}

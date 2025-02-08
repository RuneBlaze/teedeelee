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
    pub(crate) topology_only: bool,
}

/// A thread-safe, shareable view of a tree collection with complete information
/// (including branch lengths and support values).
#[derive(Clone)]
pub struct CompleteTreeView(pub(crate) TreeViewBase);

/// A thread-safe, shareable view of a tree collection with topology-only information
#[derive(Clone)]
pub struct TopologyTreeView(pub(crate) TreeViewBase);

// Add this enum to represent the tree view type
#[derive(Clone)]
pub enum TreeViewType {
    Complete(TreeViewBase),
    Topology(TreeViewBase),
}

/// A thread-safe, shareable view of a multi-species coalescent tree collection
#[derive(Clone)]
pub struct MSCTreeView {
    pub(crate) base: TreeViewType,
    pub(crate) species_tree: TreeViewType,
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
                    topology_only: self.0.topology_only,
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

// Add a helper struct for topology views
struct TopologyCollectionView<'a> {
    taxon_set: &'a TaxonSet,
    trees: &'a [Tree],
}

impl<'a> TreeCollectionTrait for TopologyCollectionView<'a> {
    fn taxon_set(&self) -> &TaxonSet {
        self.taxon_set
    }

    fn trees(&self) -> &[Tree] {
        self.trees
    }

    fn include_branch_lengths(&self) -> bool {
        false // Topology views don't include branch lengths
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
            topology_only: collection.topology_only,
        })
    }

    /// Converts this view to a topology-only view
    pub fn as_topology(self) -> TopologyTreeView {
        TopologyTreeView(TreeViewBase {
            taxon_set: self.0.taxon_set,
            trees: self.0.trees,
            topology_only: self.0.topology_only,
        })
    }

    /// Creates a new CompleteTreeView from a string containing semicolon-separated Newick trees
    pub fn from_newick_string(newick_str: &str) -> Result<Self, String> {
        TreeCollection::from_newick_string(newick_str, false)
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
            topology_only: true,
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
            topology_only: true,
        }))
    }

    /// Creates a new TopologyTreeView from a string containing semicolon-separated Newick trees
    pub fn from_newick_string(newick_str: &str) -> Result<Self, String> {
        TreeCollection::from_newick_string(newick_str, true)
            .map(Self::new)
            .map_err(|e| e.to_string())
    }
}

// Add methods to convert between Complete and Topology views
impl TreeViewType {
    fn as_topology(self) -> Self {
        match self {
            Self::Complete(base) => Self::Topology(base),
            Self::Topology(base) => Self::Topology(base),
        }
    }

    fn ngenes(&self) -> usize {
        match self {
            Self::Complete(base) | Self::Topology(base) => base.trees.len(),
        }
    }

    fn ntaxa(&self) -> usize {
        match self {
            Self::Complete(base) | Self::Topology(base) => base.taxon_set.len(),
        }
    }

    // Additional common methods
    pub fn get_tree(&self, index: usize) -> Option<&Tree> {
        match self {
            Self::Complete(base) | Self::Topology(base) => base.trees.get(index),
        }
    }

    fn get_taxon_name(&self, id: usize) -> Option<&str> {
        match self {
            Self::Complete(base) | Self::Topology(base) => {
                base.taxon_set.names.get(id).map(|s| s.as_str())
            }
        }
    }

    fn iter(&self) -> impl Iterator<Item = &Tree> {
        match self {
            Self::Complete(base) | Self::Topology(base) => base.trees.iter(),
        }
    }

    fn slice(&self, start: usize, end: usize) -> Self {
        match self {
            Self::Complete(base) => Self::Complete(TreeViewBase {
                taxon_set: Arc::clone(&base.taxon_set),
                trees: Arc::from(&base.trees[start..end]),
                topology_only: base.topology_only,
            }),
            Self::Topology(base) => Self::Topology(TreeViewBase {
                taxon_set: Arc::clone(&base.taxon_set),
                trees: Arc::from(&base.trees[start..end]),
                topology_only: true,
            }),
        }
    }

    fn get_distance_matrix(&self, index: usize) -> Option<DistanceMatrixView> {
        self.get_tree(index).map(|tree| DistanceMatrixView {
            matrix: Arc::new(tree.distance_matrix()),
            taxon_set: Arc::clone(match self {
                Self::Complete(base) | Self::Topology(base) => &base.taxon_set,
            }),
        })
    }

    fn to_newick_vec(&self) -> Vec<String> {
        let collection = TreeCollectionView {
            taxon_set: match self {
                Self::Complete(base) | Self::Topology(base) => &base.taxon_set,
            },
            trees: match self {
                Self::Complete(base) | Self::Topology(base) => &base.trees,
            },
        };
        collection.to_newick_vec()
    }

    fn to_newick_string(&self) -> String {
        let collection = TreeCollectionView {
            taxon_set: match self {
                Self::Complete(base) | Self::Topology(base) => &base.taxon_set,
            },
            trees: match self {
                Self::Complete(base) | Self::Topology(base) => &base.trees,
            },
        };
        collection.to_newick_string()
    }

    fn restriction(&self, taxa: &[&str]) -> Result<Self, String> {
        // Convert taxa strings to a HashSet of taxon IDs
        let base = match self {
            Self::Complete(base) | Self::Topology(base) => base,
        };

        let mut keep_taxa = HashSet::new();
        for taxon in taxa {
            if let Some(id) = base.taxon_set.to_id.get(*taxon) {
                keep_taxa.insert(*id);
            } else {
                return Err(format!("Taxon '{}' not found in tree", taxon));
            }
        }

        // Create new taxon set and mapping
        let mut new_taxon_set = TaxonSet::new();
        let mut id_map = HashMap::new();
        for taxon in taxa {
            let old_id = base.taxon_set.to_id[*taxon];
            let new_id = new_taxon_set.request(taxon.to_string());
            id_map.insert(old_id, new_id);
        }

        // Restrict all trees
        let new_trees: Vec<Tree> = base
            .trees
            .iter()
            .map(|tree| tree.restrict(&keep_taxa, &id_map))
            .collect();

        Ok(match self {
            Self::Complete(_) => Self::Complete(TreeViewBase {
                taxon_set: Arc::new(new_taxon_set),
                trees: Arc::from(new_trees),
                topology_only: base.topology_only,
            }),
            Self::Topology(_) => Self::Topology(TreeViewBase {
                taxon_set: Arc::new(new_taxon_set),
                trees: Arc::from(new_trees),
                topology_only: true,
            }),
        })
    }
}

impl MSCTreeView {
    pub fn new(collection: MSCTreeCollection) -> Self {
        let base = TreeViewType::Complete(TreeViewBase {
            taxon_set: Arc::new(collection.taxon_set),
            trees: Arc::from(collection.gene_trees),
            topology_only: collection.topology_only,
        });

        let species_tree = TreeViewType::Complete(TreeViewBase {
            taxon_set: Arc::clone(match &base {
                TreeViewType::Complete(b) => &b.taxon_set,
                TreeViewType::Topology(b) => &b.taxon_set,
            }),
            trees: Arc::from(vec![collection.species_tree]),
            topology_only: collection.topology_only,
        });

        MSCTreeView { base, species_tree }
    }

    pub fn as_topology(self) -> Self {
        MSCTreeView {
            base: self.base.as_topology(),
            species_tree: self.species_tree.as_topology(),
        }
    }

    // Update existing methods to handle both variants
    pub fn ngenes(&self) -> usize {
        self.base.ngenes()
    }

    pub fn ntaxa(&self) -> usize {
        self.base.ntaxa()
    }

    pub fn slice(&self, start: usize, end: usize) -> Self {
        MSCTreeView {
            base: self.base.slice(start, end),
            species_tree: self.species_tree.clone(), // Species tree doesn't get sliced
        }
    }

    pub fn get_distance_matrix(&self, index: usize) -> Option<DistanceMatrixView> {
        self.base.get_distance_matrix(index)
    }

    pub fn get_species_distance_matrix(&self) -> DistanceMatrixView {
        // Assuming species_tree always has exactly one tree
        self.species_tree.get_distance_matrix(0).unwrap()
    }

    pub fn get_taxon_name(&self, id: usize) -> Option<&str> {
        self.base.get_taxon_name(id)
    }

    pub fn to_newick_vec(&self) -> Vec<String> {
        self.base.to_newick_vec()
    }

    pub fn to_newick_string(&self) -> String {
        self.base.to_newick_string()
    }

    pub fn species_tree_to_newick(&self) -> String {
        self.species_tree.to_newick_string()
    }

    pub fn restriction(&self, taxa: &[&str]) -> Result<Self, String> {
        Ok(MSCTreeView {
            base: self.base.restriction(taxa)?,
            species_tree: self.species_tree.restriction(taxa)?,
        })
    }

    pub fn get_tree(&self, index: usize) -> Option<&Tree> {
        self.base.get_tree(index)
    }

    pub fn iter(&self) -> impl Iterator<Item = &Tree> {
        self.base.iter()
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

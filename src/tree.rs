use ahash::AHashSet;
use rand::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::{cmp::min, collections::BTreeMap, collections::HashMap};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonSet {
    pub to_id: BTreeMap<String, usize>,
    pub names: Vec<String>, // TODO: no need to keep two copies of the same string
    last: usize,
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

impl TaxonSet {
    pub fn request(&mut self, taxon_name: String) -> usize {
        *self.to_id.entry(taxon_name.clone()).or_insert_with(|| {
            self.names.push(taxon_name);
            self.last += 1;
            self.last - 1
        })
    }

    pub fn retrieve(&self, taxon_name: &str) -> usize {
        *self.to_id.get(taxon_name).unwrap_or_else(|| {
            panic!(
                "Taxon '{}' not found. Available taxa: {:?}",
                taxon_name, self.names
            )
        })
    }

    pub fn new() -> Self {
        TaxonSet {
            to_id: BTreeMap::new(),
            names: Vec::new(),
            last: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.last
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tree {
    pub taxa: Vec<i32>,
    pub parents: Vec<i32>,
    pub support: Vec<f64>, // branch support
    pub lengths: Vec<f64>, // branch lengths
    pub firstchild: Vec<i32>,
    pub nextsib: Vec<i32>,
    pub childcount: Vec<u32>,
    pub fake_root: bool,
    pub ntaxa: usize,
}

impl Tree {
    pub fn children(&self, node: usize) -> ChildrenIterator {
        ChildrenIterator::new(self, node)
    }

    pub fn postorder(&self) -> PostorderIterator {
        PostorderIterator::new(self)
    }

    pub fn ancestors(&self, node: usize) -> AncestorsIterator {
        AncestorsIterator::new(self, node)
    }

    pub fn postorder_from(&self, node: usize) -> PostorderIterator {
        PostorderIterator::from_node(self, node)
    }

    pub fn is_leaf(&self, node: usize) -> bool {
        if self.childcount[node] == 0 {
            return true;
        }
        false
    }

    pub fn is_root(&self, node: usize) -> bool {
        node == 0
    }

    pub fn num_nodes(&self) -> usize {
        self.taxa.len() - if self.fake_root { 1 } else { 0 }
    }

    pub fn xor_clades(&self, transl: &Vec<u64>, universe: u64) -> HashSet<u64> {
        let mut bips = HashSet::new();
        let mut bit_reprs = vec![0u64; self.taxa.len()];
        let mut calculated_fake_root = false;
        for node in self.postorder() {
            if self.is_leaf(node) {
                bit_reprs[node] = transl[self.taxa[node as usize] as usize];
            } else {
                if self.is_root(node) {
                    continue;
                }
                if self.is_root(self.parents[node as usize] as usize) & self.fake_root {
                    if calculated_fake_root
                        || self
                            .children(self.parents[node as usize] as usize)
                            .any(|it| self.is_leaf(it))
                    {
                        continue;
                    } else {
                        calculated_fake_root = true;
                    }
                }
                let clade = self
                    .children(node)
                    .map(|c| bit_reprs[c])
                    .reduce(|acc, i| acc ^ i)
                    .unwrap();
                bit_reprs[node] = clade;
                bips.insert(min(clade, universe ^ clade));
            }
        }
        bips
    }

    /// Suppresses unifurcations in the tree by removing nodes with exactly one child
    /// Returns None if there are no unifurcations to suppress
    fn suppress_unifurcations(&self) -> Option<Tree> {
        // Fast path: check if there are any unifurcations
        let has_unifurcations = (1..self.taxa.len()).any(|node| self.childcount[node] == 1);
        if !has_unifurcations {
            return None;
        }

        let mut new_taxa = vec![-1];
        let mut new_parents = vec![0];
        let mut new_firstchild = vec![-1];
        let mut new_nextsib = vec![-1];
        let mut new_childcount = vec![0];
        let mut new_support = vec![0.0];
        let mut new_lengths = vec![0.0];

        // Map from old node indices to new node indices
        let mut node_map = HashMap::new();
        node_map.insert(0, 0);
        let mut next_idx = 1;

        // Process all nodes except root
        for node in 1..self.taxa.len() {
            // Skip nodes with exactly one child
            if self.childcount[node] == 1 {
                continue;
            }

            // Find effective parent (skip unifurcation nodes)
            let mut parent = self.parents[node] as usize;
            while parent != 0 && self.childcount[parent] == 1 {
                parent = self.parents[parent] as usize;
            }

            let parent_idx = node_map[&parent];
            let new_idx = next_idx;
            next_idx += 1;
            node_map.insert(node, new_idx);

            // Add node to new tree
            new_taxa.push(self.taxa[node]);
            new_parents.push(parent_idx as i32);
            new_firstchild.push(-1);
            new_nextsib.push(-1);
            new_childcount.push(0);

            // Accumulate branch lengths through unifurcations
            let mut length = self.lengths[node];
            let mut curr = node;
            while self.parents[curr] as usize != parent {
                curr = self.parents[curr] as usize;
                length += self.lengths[curr];
            }
            new_lengths.push(length);

            // Use support from lowest unifurcation
            new_support.push(self.support[node]);

            // Update parent's child pointers
            new_childcount[parent_idx] += 1;
            if new_firstchild[parent_idx] == -1 {
                new_firstchild[parent_idx] = new_idx as i32;
            } else {
                let mut sibling = new_firstchild[parent_idx] as usize;
                while new_nextsib[sibling] != -1 {
                    sibling = new_nextsib[sibling] as usize;
                }
                new_nextsib[sibling] = new_idx as i32;
            }
        }

        let suppressed_tree = Tree {
            taxa: new_taxa,
            parents: new_parents,
            support: new_support,
            lengths: new_lengths,
            firstchild: new_firstchild,
            nextsib: new_nextsib,
            childcount: new_childcount,
            fake_root: self.fake_root,
            ntaxa: self.ntaxa,
        };

        // If the suppressed tree's root (node 0) has a single child then unwrap it by re-rooting.
        if suppressed_tree.childcount[0] == 1 {
            let new_root = suppressed_tree.firstchild[0] as usize;
            return Some(suppressed_tree.re_root(new_root));
        }

        Some(suppressed_tree)
    }

    /// Creates a restricted version of the tree containing only the specified taxa,
    /// suppressing unifurcations (i.e. collapsing internal nodes with only one kept descendant).
    pub fn restrict(&self, keep_taxa: &HashSet<usize>, id_map: &HashMap<usize, usize>) -> Tree {
        // First pass: compute, for each node, the number of kept taxa in its subtree.
        let mut kept_counts = vec![0usize; self.taxa.len()];
        let post_order: Vec<usize> = self.postorder().collect();
        for &node in &post_order {
            if self.is_leaf(node) {
                kept_counts[node] = if keep_taxa.contains(&(self.taxa[node] as usize)) {
                    1
                } else {
                    0
                };
            } else {
                let mut sum = 0;
                for child in self.children(node) {
                    sum += kept_counts[child];
                }
                kept_counts[node] = sum;
            }
        }

        // Determine essential nodes:
        // Keep the node if it is the root,
        // OR if it is a kept leaf,
        // OR if it has at least 2 kept taxa in its subtree.
        let mut essential_nodes = HashSet::new();
        for node in 0..self.taxa.len() {
            if node == 0
                || (self.is_leaf(node) && keep_taxa.contains(&(self.taxa[node] as usize)))
                || kept_counts[node] > 1
            {
                essential_nodes.insert(node);
            }
        }

        // Build new tree structure using only essential nodes.
        let mut new_taxa = vec![-1];
        let mut new_parents = vec![0];
        let mut new_firstchild = vec![-1];
        let mut new_nextsib = vec![-1];
        let mut new_childcount = vec![0];

        // Map the root (node 0) to the new tree's root.
        let mut node_map = HashMap::new();
        node_map.insert(0, 0);
        let mut next_idx = 1;

        // Helper: find the nearest essential ancestor (it might be several levels above)
        let nearest_essential = |node: usize| -> usize {
            let mut p = self.parents[node] as usize;
            while !essential_nodes.contains(&p) {
                p = self.parents[p] as usize;
            }
            p
        };

        // Second pass: add only those nodes that are essential.
        // (Start from 1 because root is already in the map.)
        for node in 1..self.taxa.len() {
            if !essential_nodes.contains(&node) {
                continue;
            }
            // Get the effective parent: if the node's immediate parent is not essential,
            // then collapse (skip) that unifurcation by finding the nearest essential ancestor.
            let effective_parent = if essential_nodes.contains(&(self.parents[node] as usize)) {
                self.parents[node] as usize
            } else {
                nearest_essential(node)
            };

            // (Because we traverse in increasing order, effective_parent should have already been processed.)
            if !node_map.contains_key(&effective_parent) {
                continue; // Should not happen with proper traversal
            }
            let parent_idx = node_map[&effective_parent];
            let new_idx = next_idx;
            next_idx += 1;
            node_map.insert(node, new_idx);

            // For leaves, use the taxon id from id_map. For internal nodes, keep it as -1.
            let taxon_val = if self.is_leaf(node) {
                id_map
                    .get(&(self.taxa[node] as usize))
                    .map(|&id| id as i32)
                    .unwrap_or(-1)
            } else {
                -1
            };

            new_taxa.push(taxon_val);
            new_parents.push(parent_idx as i32);
            new_firstchild.push(-1);
            new_nextsib.push(-1);
            new_childcount.push(0);

            // Update parent's child pointers.
            new_childcount[parent_idx] += 1;
            if new_firstchild[parent_idx] == -1 {
                new_firstchild[parent_idx] = new_idx as i32;
            } else {
                let mut sibling = new_firstchild[parent_idx] as usize;
                while new_nextsib[sibling] != -1 {
                    sibling = new_nextsib[sibling] as usize;
                }
                new_nextsib[sibling] = new_idx as i32;
            }
        }

        // Create initial restricted tree
        let restricted = Tree {
            taxa: new_taxa,
            parents: new_parents,
            support: vec![0.0; next_idx],
            lengths: vec![0.0; next_idx],
            firstchild: new_firstchild,
            nextsib: new_nextsib,
            childcount: new_childcount,
            fake_root: false,
            ntaxa: keep_taxa.len(),
        };

        // Suppress any remaining unifurcations
        restricted.suppress_unifurcations().unwrap_or(restricted)
    }

    pub fn distance_matrix(&self) -> DistanceMatrix {
        DistanceMatrix::compute_topological(self)
    }

    pub fn remap_taxa(&mut self, id_map: &HashMap<usize, usize>) {
        for taxon_id in self.taxa.iter_mut() {
            if let Some(&new_id) = id_map.get(&(*taxon_id as usize)) {
                *taxon_id = new_id as i32;
            }
        }
    }

    /// Sorts the children of each node according to the given criteria.
    /// Earlier criteria in the list take precedence over later ones.
    pub fn sort(&mut self, taxon_set: &TaxonSet, criteria: &[SortCriterion]) {
        // First compute all the data we need for comparisons
        let descendant_counts = if criteria.iter().any(|c| c.kind == SortCriterionKind::DescendantCount) {
            let mut counts = vec![0usize; self.taxa.len()];
            for node in self.postorder() {
                counts[node] = if self.is_leaf(node) {
                    1
                } else {
                    self.children(node).map(|c| counts[c]).sum()
                }
            }
            Some(counts)
        } else {
            None
        };

        // Pre-collect sorted children vectors for each node
        let mut sorted_children: Vec<Vec<usize>> = vec![Vec::new(); self.taxa.len()];
        for node in self.postorder() {
            sorted_children[node] = self.children(node).collect();
        }

        // Define a helper function for recursive comparison
        fn compare_nodes_recursive(
            a: usize,
            b: usize,
            taxon_set: &TaxonSet,
            tree: &Tree,
            criteria: &[SortCriterion],
            sorted_children: &Vec<Vec<usize>>,
            descendant_counts: Option<&Vec<usize>>,
        ) -> std::cmp::Ordering {
            for criterion in criteria {
                let mut cmp = match criterion.kind {
                    SortCriterionKind::LexicographicalOrder => {
                        if tree.is_leaf(a) && tree.is_leaf(b) {
                            let name_a = &taxon_set.names[tree.taxa[a] as usize];
                            let name_b = &taxon_set.names[tree.taxa[b] as usize];
                            name_a.cmp(name_b)
                        } else if tree.is_leaf(a) || tree.is_leaf(b) {
                            let val_a = if tree.is_leaf(a) { &taxon_set.names[tree.taxa[a] as usize] } else { "" };
                            let val_b = if tree.is_leaf(b) { &taxon_set.names[tree.taxa[b] as usize] } else { "" };
                            val_a.cmp(val_b)
                        } else {
                            let children_a = &sorted_children[a];
                            let children_b = &sorted_children[b];
                            for (child_a, child_b) in children_a.iter().zip(children_b.iter()) {
                                let child_cmp = compare_nodes_recursive(*child_a, *child_b, taxon_set, tree, criteria, sorted_children, descendant_counts);
                                if child_cmp != std::cmp::Ordering::Equal {
                                    return child_cmp;
                                }
                            }
                            children_a.len().cmp(&children_b.len())
                        }
                    },
                    SortCriterionKind::ChildCount => tree.childcount[a].cmp(&tree.childcount[b]),
                    SortCriterionKind::DescendantCount => {
                        if let Some(ref counts) = descendant_counts {
                            counts[a].cmp(&counts[b])
                        } else {
                            std::cmp::Ordering::Equal
                        }
                    }
                };
                if criterion.order == SortOrder::Descending {
                    cmp = cmp.reverse();
                }
                if cmp != std::cmp::Ordering::Equal {
                    return cmp;
                }
            }
            std::cmp::Ordering::Equal
        }

        // Sort children of each node
        for node in 0..self.taxa.len() {
            if self.childcount[node] <= 1 {
                continue;
            }

            // Create an immutable snapshot of sorted_children to avoid borrow conflicts
            let snapshot = sorted_children.clone();

            // Sort children according to criteria using the snapshot
            sorted_children[node].sort_by(|&a, &b| {
                compare_nodes_recursive(a, b, taxon_set, self, criteria, &snapshot, descendant_counts.as_ref())
            });
            
            // Update the tree structure with the new ordering
            self.firstchild[node] = sorted_children[node][0] as i32;
            for i in 0..sorted_children[node].len() - 1 {
                self.nextsib[sorted_children[node][i]] = sorted_children[node][i + 1] as i32;
            }
            self.nextsib[sorted_children[node][sorted_children[node].len() - 1]] = -1;
        }
    }

    /// Re-roots the tree at `new_root`, doing a DFS re-indexing so that the node formerly at `new_root`
    /// becomes index 0 and all pointers (parents, children) are updated accordingly.
    fn re_root(&self, new_root: usize) -> Tree {
        use std::collections::HashMap;
        // Map from old indices to new indices (in DFS order)
        let mut new_index_map: HashMap<usize, usize> = HashMap::new();
        let mut stack = vec![new_root];
        let mut order = Vec::new(); // Order in which nodes are visited

        while let Some(n) = stack.pop() {
            if new_index_map.contains_key(&n) {
                continue;
            }
            let new_idx = new_index_map.len();
            new_index_map.insert(n, new_idx);
            order.push(n);

            // Collect children of node `n`
            let mut children = Vec::new();
            let mut child = self.firstchild[n];
            while child != -1 {
                children.push(child as usize);
                child = self.nextsib[child as usize];
            }
            // Push children in reverse order so that the left-most child is processed first.
            for &c in children.iter().rev() {
                stack.push(c);
            }
        }

        let new_n = order.len();
        let mut taxa = vec![0; new_n];
        let mut parents = vec![0; new_n];
        let mut support = vec![0.0; new_n];
        let mut lengths = vec![0.0; new_n];
        let mut firstchild = vec![-1; new_n];
        let mut nextsib = vec![-1; new_n];
        let mut childcount = vec![0; new_n];

        // Copy over the field values and update parent pointers using the mapping.
        for (new_idx, &old_idx) in order.iter().enumerate() {
            taxa[new_idx] = self.taxa[old_idx];
            support[new_idx] = self.support[old_idx];
            lengths[new_idx] = self.lengths[old_idx];
            if new_idx == 0 {
                // New root: no parent.
                parents[new_idx] = 0;
            } else {
                let old_parent = self.parents[old_idx] as usize;
                parents[new_idx] = *new_index_map.get(&old_parent).unwrap() as i32;
            }
        }

        // Reconstruct child pointers by gathering children using the new indices.
        let mut children_map: Vec<Vec<usize>> = vec![Vec::new(); new_n];
        for i in 1..new_n {
            let p = parents[i] as usize;
            children_map[p].push(i);
        }
        for i in 0..new_n {
            if !children_map[i].is_empty() {
                firstchild[i] = children_map[i][0] as i32;
                childcount[i] = children_map[i].len() as u32;
                for (j, &child) in children_map[i].iter().enumerate() {
                    if j < children_map[i].len() - 1 {
                        nextsib[child] = children_map[i][j + 1] as i32;
                    } else {
                        nextsib[child] = -1;
                    }
                }
            } else {
                childcount[i] = 0;
            }
        }

        // Recompute number of leaves (taxa) if needed.
        let ntaxa = taxa.iter().filter(|&&t| t != -1).count();
        Tree {
            taxa,
            parents,
            support,
            lengths,
            firstchild,
            nextsib,
            childcount,
            fake_root: false,
            ntaxa,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TreeCollection {
    pub taxon_set: TaxonSet,
    pub trees: Vec<Tree>,
    pub topology_only: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MSCTreeCollection {
    pub taxon_set: TaxonSet,
    pub gene_trees: Vec<Tree>,
    pub species_tree: Tree,
    pub topology_only: bool,
}

/// Common functionality for collections of trees
pub trait TreeCollectionTrait {
    fn taxon_set(&self) -> &TaxonSet;
    fn trees(&self) -> &[Tree];

    // Add this method to determine if we should include branch lengths
    fn include_branch_lengths(&self) -> bool {
        true // Default implementation includes branch lengths
    }

    fn ngenes(&self) -> usize {
        self.trees().len()
    }

    fn ntaxa(&self) -> usize {
        self.taxon_set().len()
    }

    /// Convert a single tree to Newick format
    fn tree_to_newick(&self, tree: &Tree) -> String {
        fn write_subtree(
            tree: &Tree,
            node: usize,
            taxon_set: &TaxonSet,
            include_lengths: bool,
            result: &mut String,
        ) {
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
                    write_subtree(tree, child, taxon_set, include_lengths, result);
                    if include_lengths && tree.lengths[child] >= 0.0 {
                        result.push(':');
                        result.push_str(&tree.lengths[child].to_string());
                    }
                }
                result.push(')');
            }
        }

        let mut result = String::new();
        write_subtree(
            tree,
            0,
            self.taxon_set(),
            self.include_branch_lengths(),
            &mut result,
        );
        result.push(';');
        result
    }

    /// Convert the tree collection to a vector of Newick strings
    fn to_newick_vec(&self) -> Vec<String> {
        self.trees()
            .iter()
            .map(|tree| self.tree_to_newick(tree))
            .collect()
    }

    /// Convert the tree collection to a single string with semicolon+newline-separated Newick trees
    fn to_newick_string(&self) -> String {
        self.trees()
            .iter()
            .map(|tree| self.tree_to_newick(tree))
            .collect::<Vec<_>>()
            .join("\n")
    }
}

impl TreeCollectionTrait for TreeCollection {
    fn taxon_set(&self) -> &TaxonSet {
        &self.taxon_set
    }

    fn trees(&self) -> &[Tree] {
        &self.trees
    }

    fn include_branch_lengths(&self) -> bool {
        !self.topology_only
    }
}

impl MSCTreeCollection {
    pub fn new(topology_only: bool) -> Self {
        MSCTreeCollection {
            taxon_set: TaxonSet::new(),
            gene_trees: Vec::new(),
            species_tree: Tree {
                taxa: vec![-42],
                parents: vec![0],
                support: vec![-1.0],
                lengths: vec![-1.0],
                firstchild: vec![-1],
                nextsib: vec![-1],
                childcount: vec![0],
                fake_root: false,
                ntaxa: 0,
            },
            topology_only,
        }
    }

    pub fn species_tree(&self) -> &Tree {
        &self.species_tree
    }
}

impl TreeCollectionTrait for MSCTreeCollection {
    fn taxon_set(&self) -> &TaxonSet {
        &self.taxon_set
    }

    fn trees(&self) -> &[Tree] {
        &self.gene_trees
    }

    fn include_branch_lengths(&self) -> bool {
        !self.topology_only
    }
}

pub struct PostorderIterator {
    // s1 : Vec<usize>,
    s2: Vec<usize>,
}

pub struct AncestorsIterator<'a> {
    tree: &'a Tree,
    current: i32,
}

impl<'a> AncestorsIterator<'a> {
    pub fn new(tree: &'a Tree, taxon: usize) -> Self {
        let n = taxon;
        if n == 0usize {
            AncestorsIterator { tree, current: 0 }
        } else {
            AncestorsIterator {
                tree,
                current: tree.parents[n] as i32,
            }
        }
    }
}

impl<'a> Iterator for AncestorsIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current == 0 {
            None
        } else {
            let res = self.current as usize;
            self.current = self.tree.parents[res] as i32;
            Some(res)
        }
    }
}

pub struct ChildrenIterator<'a> {
    tree: &'a Tree,
    current: i32,
}

impl<'a> ChildrenIterator<'a> {
    pub fn new(tree: &'a Tree, node: usize) -> Self {
        ChildrenIterator {
            tree,
            current: tree.firstchild[node],
        }
    }
}

impl<'a> Iterator for ChildrenIterator<'a> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current == -1 {
            None
        } else {
            let res = self.current as usize;
            self.current = self.tree.nextsib[self.current as usize];
            Some(res)
        }
    }
}

impl PostorderIterator {
    pub fn new(tree: &Tree) -> Self {
        Self::from_node(tree, 0)
    }

    pub fn from_node(tree: &Tree, node: usize) -> Self {
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();
        s1.push(node);
        while let Some(n) = s1.pop() {
            s2.push(n);
            tree.children(n).for_each(|c| s1.push(c));
        }
        PostorderIterator {
            // s1,
            s2,
        }
    }

    pub fn from_node_excluding(tree: &Tree, node: usize, exclusion: &AHashSet<usize>) -> Self {
        let mut s1 = Vec::new();
        let mut s2 = Vec::new();
        s1.push(node);
        while let Some(n) = s1.pop() {
            s2.push(n);
            tree.children(n).for_each(|c| {
                if !exclusion.contains(&c) {
                    s1.push(c);
                }
            });
        }
        PostorderIterator { s2 }
    }
}

impl Iterator for PostorderIterator {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.s2.pop()
    }
}

pub fn parse_newick(taxon_set: &mut TaxonSet, newick: &str) -> Tree {
    let mut taxa: Vec<i32> = vec![-42];
    let mut parents: Vec<i32> = vec![0];
    let mut support: Vec<f64> = vec![-1.0];
    let mut lengths: Vec<f64> = vec![-1.0];
    let mut childcount: Vec<u32> = vec![0];
    let mut firstchild: Vec<i32> = vec![-1];
    let mut nextsib: Vec<i32> = vec![-1];
    let mut ntaxa: usize = 0;
    // we just reuse TreeSwift's logic
    let mut n: usize = 0; // the current node
    let mut chars = newick.chars().fuse().peekable();
    while let Some(c) = chars.next() {
        // println!("{}", c);
        if c == ';' {
            break;
        } else if c == ' ' {
            continue;
        } else if c == '(' {
            taxa.push(-1);
            childcount[n as usize] += 1;
            parents.push(n as i32);
            support.push(0.0);
            lengths.push(0.0);
            childcount.push(0);
            firstchild.push(-1);
            nextsib.push(-1);
            firstchild[n] = (taxa.len() - 1) as i32;
            n = taxa.len() - 1;
        } else if c == ')' {
            n = parents[n] as usize;
        } else if c == ',' {
            nextsib[n] = (taxa.len()) as i32;
            n = parents[n] as usize;
            taxa.push(-1);
            childcount[n as usize] += 1;
            parents.push(n as i32);
            support.push(0.0);
            lengths.push(0.0);
            childcount.push(0);
            firstchild.push(-1);
            nextsib.push(-1);
            n = taxa.len() - 1;
        } else if c == '[' {
            let mut cnt = 0usize;
            loop {
                match chars.next() {
                    Some(']') => {
                        if cnt == 0 {
                            break;
                        } else {
                            cnt -= 1;
                        }
                    }
                    Some('[') => {
                        cnt += 1;
                    }
                    Some(_) | None => {}
                }
            }
        } else if c == ':' {
            let mut ls = "".to_string();
            loop {
                match chars.peek() {
                    Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ls.push(chars.next().unwrap());
                    }
                }
            }
            if ls.len() > 0 {
                lengths[n as usize] = ls.parse::<f64>().unwrap();
            }
        } else {
            let mut ts = c.to_string();
            loop {
                match chars.peek() {
                    Some(':') | Some(',') | Some(')') | Some(';') | Some('[') | None => {
                        break;
                    }
                    Some(_) => {
                        ts.push(chars.next().unwrap());
                    }
                }
            }
            if childcount[n] == 0 {
                // println!("{}", ts);
                taxa[n] = taxon_set.request(ts) as i32;
                ntaxa += 1;
            }
        }
    }

    Tree {
        taxa,
        parents,
        support,
        lengths,
        firstchild,
        nextsib,
        childcount,
        fake_root: false,
        ntaxa,
    }
}

impl TreeCollection {
    pub fn new(topology_only: bool) -> Self {
        TreeCollection {
            taxon_set: TaxonSet::new(),
            trees: Vec::new(),
            topology_only,
        }
    }

    pub fn from_newick_string(newick_str: &str, topology_only: bool) -> Result<Self, String> {
        let mut taxon_set = TaxonSet::new();
        let mut trees = Vec::new();

        for newick in newick_str.split(';').filter(|s| !s.trim().is_empty()) {
            trees.push(parse_newick(&mut taxon_set, &format!("{};", newick.trim())));
        }

        if trees.is_empty() {
            return Err("No valid trees found in input string".to_string());
        }

        Ok(TreeCollection {
            taxon_set,
            trees,
            topology_only,
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DistanceMatrix {
    pub(crate) distances: Vec<f64>,
    pub(crate) size: usize,
}

impl DistanceMatrix {
    pub fn new(size: usize) -> Self {
        // For n taxa, we need n*(n-1)/2 entries for upper triangular matrix
        let n_entries = (size * (size - 1)) / 2;
        DistanceMatrix {
            distances: vec![0.0; n_entries],
            size,
        }
    }

    /// Convert (i,j) coordinates to vector index
    pub fn get_index(&self, i: usize, j: usize) -> usize {
        debug_assert!(i < j);
        debug_assert!(j < self.size);
        // Formula for accessing upper triangular matrix stored as vector
        (self.size * i) - ((i * (i + 1)) / 2) + j - i - 1
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            return 0.0;
        }
        if i < j {
            self.distances[self.get_index(i, j)]
        } else {
            self.distances[self.get_index(j, i)]
        }
    }

    pub fn set(&mut self, i: usize, j: usize, value: f64) {
        if i == j {
            return;
        }
        let idx = if i < j {
            self.get_index(i, j)
        } else {
            self.get_index(j, i)
        };
        self.distances[idx] = value;
    }

    pub fn add(&mut self, i: usize, j: usize, value: f64) {
        if i == j {
            return;
        }
        let idx = if i < j {
            self.get_index(i, j)
        } else {
            self.get_index(j, i)
        };
        self.distances[idx] += value;
    }

    /// Compute topological distance matrix for a tree
    pub fn compute_topological(tree: &Tree) -> Self {
        let n = tree.ntaxa;
        let mut matrix = DistanceMatrix::new(n);

        // Store leaf-to-root distances for each leaf
        let mut leaf_dists = Vec::<Vec<(usize, f64)>>::new();
        leaf_dists.resize(tree.taxa.len(), Vec::new());

        for node in tree.postorder() {
            if tree.is_leaf(node) {
                leaf_dists[node].push((node, 0.0));
            } else {
                // Process children
                for c in tree.children(node) {
                    // Add 1.0 to all distances (topological distance)
                    for e in leaf_dists.get_mut(c).unwrap() {
                        e.1 += 1.0;
                    }
                }

                // Calculate pairwise distances between leaves in different subtrees
                let node_children: Vec<usize> = tree.children(node).collect();
                for c1 in 0..(node_children.len() - 1) {
                    let leaves_c1 = &leaf_dists[node_children[c1]];
                    for c2 in (c1 + 1)..node_children.len() {
                        let leaves_c2 = &leaf_dists[node_children[c2]];
                        for (u, ud) in leaves_c1.iter() {
                            for (v, vd) in leaves_c2.iter() {
                                let dist = ud + vd;
                                let u_leaf = tree.taxa[*u] as usize;
                                let v_leaf = tree.taxa[*v] as usize;
                                let (l, r) = if u_leaf < v_leaf {
                                    (u_leaf, v_leaf)
                                } else {
                                    (v_leaf, u_leaf)
                                };
                                matrix.set(l, r, dist);
                            }
                        }
                    }
                }

                // Combine leaf lists for the next iteration
                for (i, e) in tree.children(node).enumerate() {
                    if i == 0 {
                        leaf_dists.swap(node, tree.firstchild[node] as usize);
                    } else {
                        leaf_dists.push(vec![]);
                        let mut v = leaf_dists.swap_remove(e);
                        leaf_dists[node].append(&mut v);
                    }
                }
            }
        }

        matrix
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortOrder {
    Ascending,
    Descending,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SortCriterion {
    pub kind: SortCriterionKind,
    pub order: SortOrder,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortCriterionKind {
    LexicographicalOrder,
    ChildCount,
    DescendantCount,
}

impl SortCriterion {
    pub fn ascending(kind: SortCriterionKind) -> Self {
        SortCriterion { kind, order: SortOrder::Ascending }
    }

    pub fn descending(kind: SortCriterionKind) -> Self {
        SortCriterion { kind, order: SortOrder::Descending }
    }
}

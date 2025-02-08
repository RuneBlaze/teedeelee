import pytest
from hypothesis import given, settings, strategies as st
import treeswift as ts
from teedeelee import MSCTreeSet, FamilyOfMSC, TreeSet, Tree
import glob
import os
import logging
from typing import List, Tuple
import numpy as np
from itertools import combinations, islice, chain

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def load_test_data(
    base_dir: str = "/Users/lbq/Downloads/S101", num_datasets: int = 5
) -> Tuple[List[str], List[str]]:
    """
    Load test data from S101 directory structure.

    Args:
        base_dir: Base directory containing the numbered subdirectories
        num_datasets: Number of datasets to load

    Returns:
        Tuple of (gene trees files list, species tree files list)
    """
    subdirs = sorted(glob.glob(os.path.join(base_dir, "*")))[:num_datasets]
    gene_trees_files = []
    species_tree_files = []

    for subdir in subdirs:
        if not os.path.isdir(subdir):
            continue

        gene_trees_file = os.path.join(subdir, "truegenetrees")
        species_tree_file = os.path.join(subdir, "s_tree.trees")

        if not (os.path.exists(gene_trees_file) and os.path.exists(species_tree_file)):
            logger.warning(f"Missing required files in {subdir}")
            continue

        gene_trees_files.append(gene_trees_file)
        species_tree_files.append(species_tree_file)

    return gene_trees_files, species_tree_files


def check_tree_faithfulness(
    species_tree_file: str, gene_trees_file: str, subset: List[str]
) -> bool:
    """
    Check if the restriction operation preserves faithfulness properties

    Returns:
        bool: True if faithfulness is preserved
    """
    # Create MSC tree set directly from files using from_files
    msctree = MSCTreeSet.from_files(gene_trees_file, species_tree_file)

    # Get original and restricted trees
    original_trees = msctree.as_topology()
    restricted_trees = original_trees.restriction(subset)

    # Check basic properties
    if len(restricted_trees) != msctree.ngenes:
        return False

    for tree in restricted_trees:
        tree_ts = ts.read_tree_newick(str(tree))
        leaf_labels = {node.label for node in tree_ts.traverse_leaves()}
        if leaf_labels != set(subset):
            return False

    return True


def test_restriction_faithfulness():
    """Test that restriction operations preserve faithfulness properties"""
    gene_trees_files, species_tree_files = load_test_data()

    for gene_trees_file, species_tree_file in zip(gene_trees_files, species_tree_files):
        logger.info(f"\nTesting dataset:")
        logger.info(f"Gene trees: {gene_trees_file}")
        logger.info(f"Species tree: {species_tree_file}")

        # Load species tree using TreeSet.from_file
        species_tree_set = TreeSet.from_file(species_tree_file)
        species_tree = species_tree_set[0]  # Get first tree

        # Get all taxa
        all_taxa = species_tree.get_taxa()

        # Test with different subset sizes
        for subset_size in [3, len(all_taxa) // 2, len(all_taxa)]:
            if subset_size > len(all_taxa):
                continue

            subset = all_taxa[:subset_size]  # Take first n taxa
            logger.info(f"Testing with subset size {subset_size}: {subset[:5]}...")

            assert check_tree_faithfulness(
                species_tree_file, gene_trees_file, subset
            ), f"Faithfulness not preserved for subset size {subset_size}"


def test_family_msc_faithfulness():
    """Test faithfulness properties of FamilyOfMSC operations"""
    gene_trees_files, species_tree_files = load_test_data()

    # Create family of MSC directly using from_files
    family = FamilyOfMSC.from_files(gene_trees_files, species_tree_files)

    # Basic property checks
    assert len(family) == len(
        gene_trees_files
    ), f"Family should contain {len(gene_trees_files)} MSC tree sets"

    # Test each MSC set in the family individually
    for i, msc_set in enumerate(family):
        logger.info(f"\nTesting MSC set {i} in family")

        species_tree = msc_set.get_species_tree()
        all_taxa = species_tree.get_taxa()

        # Test restriction with different subset sizes
        for subset_size in [3, len(all_taxa) // 2, len(all_taxa)]:
            if subset_size > len(all_taxa):
                continue

            subset = all_taxa[:subset_size]
            logger.info(f"Testing restriction with subset size {subset_size}")

            # Get topology set and restrict it
            topology_set = msc_set
            restricted_set = topology_set.restriction(subset)

            # Verify taxa in restricted trees
            for tree in restricted_set:
                tree_ts = ts.read_tree_newick(str(tree))
                restricted_taxa = {node.label for node in tree_ts.traverse_leaves()}
                assert (
                    restricted_taxa == set(subset)
                ), f"Restricted tree contains incorrect taxa. Expected {subset}, got {restricted_taxa}"


def test_edge_cases():
    """Test edge cases for faithfulness properties"""
    gene_trees_files, species_tree_files = load_test_data(num_datasets=1)

    # Load trees using from_file
    species_tree_set = TreeSet.from_file(species_tree_files[0])
    species_tree = species_tree_set[0]
    all_taxa = species_tree.get_taxa()

    # Test minimum subset size
    min_subset = all_taxa[:2]
    assert check_tree_faithfulness(
        species_tree_files[0], gene_trees_files[0], min_subset
    ), "Faithfulness not preserved for minimum subset"

    # Test subset preserving internal node
    internal_subset = all_taxa[:3]
    assert check_tree_faithfulness(
        species_tree_files[0], gene_trees_files[0], internal_subset
    ), "Faithfulness not preserved when keeping internal node"

    # Test with all taxa
    assert check_tree_faithfulness(
        species_tree_files[0], gene_trees_files[0], all_taxa
    ), "Faithfulness not preserved with all taxa"


def test_topology_preservation():
    """Test that topological relationships are preserved after restriction"""
    gene_trees_files, species_tree_files = load_test_data(num_datasets=1)

    # Load trees using from_file
    gene_trees_set = TreeSet.from_file(gene_trees_files[0])

    # Get original topology set
    topology_set = gene_trees_set.as_topology()

    # Get taxa from first tree
    all_taxa = gene_trees_set[0].get_taxa()

    # Test with different subset sizes
    for subset_size in [3, len(all_taxa) // 2, len(all_taxa)]:
        if subset_size > len(all_taxa):
            continue

        subset = all_taxa[:subset_size]
        logger.info(f"\nTesting topology preservation with subset size {subset_size}")

        # Get restricted topology
        restricted_set = topology_set.restriction(subset)

        for tree in restricted_set:
            tree_ts = ts.read_tree_newick(str(tree))

            # Check basic tree properties
            assert (
                tree_ts.num_nodes() >= 2 * len(subset) - 1
            ), f"Invalid tree structure after restriction: {tree}"

            # Verify all leaves are from subset
            leaf_labels = {node.label for node in tree_ts.traverse_leaves()}
            assert leaf_labels == set(
                subset
            ), f"Restricted tree contains incorrect taxa: {leaf_labels} vs {subset}"


def normalize_tree_str(tree_str: str) -> str:
    """Normalize a tree string by removing rooted markers and standardizing format"""
    # Remove rooted marker and any whitespace
    tree_str = tree_str.replace("[&R]", "").strip()
    # Remove edge lengths
    tree_str = tree_str.replace(":1", "")
    # Remove semicolon
    tree_str = tree_str.strip(";")

    # Use TreeSwift to standardize the tree format
    ts_tree = ts.read_tree_newick(tree_str + ";")
    ts_tree.order(
        "num_descendants_then_label"
    )  # This helps standardize the internal structure
    # ts_tree.suppress_unifurcations()
    normalized = ts_tree.newick()

    # Clean up the normalized stringe
    normalized = normalized.strip(";")
    normalized = normalized.replace(":1.0", "")

    while normalized.startswith("((") and normalized.endswith("))"):
        normalized = normalized[1:-1]

    return normalized.strip()


def test_distance_matrix_faithfulness():
    """Test that distance matrices are faithful between TreeSwift and teedeelee"""
    gene_files, species_files = load_test_data(num_datasets=1)

    # Load trees using both libraries
    td_trees = TreeSet.from_file(gene_files[0]).as_topology()
    with open(gene_files[0], "r") as f:
        ts_trees = [ts.read_tree_newick(line.strip()) for line in f]

    # Set all TreeSwift edge lengths to 1 for topological distance
    for tree in ts_trees:
        for node in tree.traverse_postorder():
            if not node.is_root():
                node.edge_length = 1

    MAX_COMPARISONS = 10  # choose(5,2)
    import itertools

    oracle = iter(itertools.chain([4] * 5, itertools.count(5)))

    # Test each corresponding pair of trees
    for idx, (td_tree, ts_tree) in enumerate(zip(td_trees, ts_trees)):
        # Get distance matrices
        td_dists = td_trees.get_distance_matrix(idx)
        ts_dists = ts_tree.distance_matrix(leaf_labels=True)

        # Get sorted taxa list for consistent ordering
        taxa = sorted(ts_dists.keys())
        n_taxa = len(taxa)

        # Compare distances between pairs, taking at most MAX_COMPARISONS
        pairs = list(islice(combinations(taxa, 2), MAX_COMPARISONS))

        for tax1, tax2 in pairs:
            ts_dist = ts_dists[tax1][tax2]
            td_dist = td_dists.get_by_name(tax1, tax2)
            assert abs(ts_dist - td_dist) < 1e-10, (
                f"Distance mismatch for {tax1}->{tax2}: ts={ts_dist}, td={td_dist}; "
                f"tree (normalized): {normalize_tree_str(str(td_tree))}"
            )

        # Test after restriction with gradually increasing subset sizes
        # size =
        size = min(next(oracle), n_taxa // 4)
        subset = taxa[:size]
        logger.info(f"Testing restriction with {size} taxa")

        # Get restricted trees and their distance matrices
        td_sub = td_trees.restriction(subset)[idx]
        td_sub_dists = Tree(td_sub.newick()).get_distance_matrix()

        ts_sub = ts_tree.extract_tree_with(set(subset))
        for node in ts_sub.traverse_postorder():
            if not node.is_root():
                node.edge_length = 1
        ts_sub_dists = ts_sub.distance_matrix(leaf_labels=True)

        # Compare distances in restricted trees, taking at most MAX_COMPARISONS
        pairs = list(islice(combinations(subset, 2), MAX_COMPARISONS))

        for tax1, tax2 in pairs:
            ts_dist = ts_sub_dists[tax1][tax2]
            td_dist = td_sub_dists.get_by_name(tax1, tax2)
            assert abs(ts_dist - td_dist) < 1e-10, (
                f"Distance mismatch after restriction for {tax1}->{tax2}:\n"
                f"ts={ts_dist}, td={td_dist}\n"
                f"Original tree: {normalize_tree_str(str(td_tree))}\n"
                f"Restricted (ts): {normalize_tree_str(ts_sub.newick())}\n"
                f"Restricted (td): {normalize_tree_str(str(td_sub))}"
            )


def create_ladder_tree(n: int) -> str:
    """Create a ladder tree with n leaves in Newick format."""
    if n < 2:
        raise ValueError("Ladder tree must have at least 2 leaves")

    # Start with two leaves
    tree = "(1,2)"

    # Add remaining leaves one by one
    for i in range(3, n + 1):
        tree = f"({tree},{i})"

    return tree + ";"


def test_ladder_tree_restriction():
    """Test restriction faithfulness using ladder trees and all possible 4-leaf subsets"""
    # Create a ladder tree with 6 leaves
    ladder_newick = create_ladder_tree(6)

    # Create TreeSwift tree
    ts_tree = ts.read_tree_newick(ladder_newick)

    # Create teedeelee tree
    td_tree = Tree(ladder_newick)

    # Get all taxa
    all_taxa = sorted([str(i) for i in range(1, 7)])

    # Test all possible 4-leaf subsets
    for subset in combinations(all_taxa, 4):
        subset = list(subset)
        logger.info(f"Testing subset: {subset}")

        # Get restricted trees
        ts_restricted = ts_tree.extract_tree_with(set(subset))
        td_restricted = td_tree.restriction(subset)

        # Compare tree structures using clades
        assert trees_have_same_clades(ts_restricted.newick(), str(td_restricted)), (
            f"Tree mismatch for subset {subset}:\n"
            f"Original tree: {ladder_newick}\n"
            f"TreeSwift: {ts_restricted.newick()}\n"
            f"teedeelee: {td_restricted}"
        )


def create_random_tree(n: int, seed: int = None) -> str:
    """Create a random binary tree with n leaves in Newick format."""
    if seed is not None:
        np.random.seed(seed)

    # Start with all leaves as separate nodes
    nodes = [[str(i)] for i in range(1, n + 1)]

    # Randomly join nodes until we have a single tree
    while len(nodes) > 1:
        # Pick two random nodes to join
        idx1, idx2 = np.random.choice(len(nodes), size=2, replace=False)
        node1 = nodes.pop(idx1)
        node2 = nodes.pop(idx2 if idx2 < idx1 else idx2 - 1)

        # Join the nodes
        new_node = [f"({','.join(node1)},{','.join(node2)})"]
        nodes.append(new_node)

    return nodes[0][0] + ";"


def test_random_tree_restriction():
    """Test restriction faithfulness using random trees and all possible 4-leaf subsets"""
    # Create multiple random trees with 6 leaves
    for seed in range(10000):  # Test 100 different random trees
        random_newick = create_random_tree(6, seed=seed)

        # Create TreeSwift tree
        ts_tree = ts.read_tree_newick(random_newick)

        # Create teedeelee tree
        td_tree = Tree(random_newick)

        # Get all taxa
        all_taxa = sorted([str(i) for i in range(1, 7)])

        # Test all possible 5-leaf subsets
        for subset in combinations(all_taxa, 5):
            subset = list(subset)
            logger.info(f"Testing random tree {seed} with subset: {subset}")

            # Get restricted trees
            ts_restricted = ts_tree.extract_tree_with(set(subset))
            td_restricted = td_tree.restriction(subset)

            # Compare tree structures using clades
            assert trees_have_same_clades(ts_restricted.newick(), str(td_restricted)), (
                f"Tree mismatch for random tree {seed}, subset {subset}:\n"
                f"Original tree: {random_newick}\n"
                f"TreeSwift: {ts_restricted.newick()}\n"
                f"teedeelee: {td_restricted}"
            )


def trees_have_same_clades(tree1_str: str, tree2_str: str) -> bool:
    """
    Compare two trees by converting them to TreeSwift trees and checking if they have identical clades.

    Args:
        tree1_str: First tree in Newick format
        tree2_str: Second tree in Newick format

    Returns:
        bool: True if trees have identical clades, False otherwise
    """
    try:
        # Parse trees with TreeSwift
        ts_tree1 = ts.read_tree_newick(tree1_str)
        ts_tree2 = ts.read_tree_newick(tree2_str)

        # Get clades for each tree
        clades1 = set()
        clades2 = set()

        # Helper to get clades from a tree
        def get_clades(tree, clades_set):
            for node in tree.traverse_postorder():
                if not node.is_leaf():
                    # Get leaf labels under this node
                    leaves = frozenset(leaf.label for leaf in node.traverse_leaves())
                    clades_set.add(leaves)

        # Get clades for both trees
        get_clades(ts_tree1, clades1)
        get_clades(ts_tree2, clades2)

        # Compare clade sets
        return clades1 == clades2

    except Exception as e:
        logger.error(f"Error comparing trees: {e}")
        logger.error(f"Tree1: {tree1_str}")
        logger.error(f"Tree2: {tree2_str}")
        return False

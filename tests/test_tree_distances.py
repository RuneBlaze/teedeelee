import pytest
import treeswift
from teedeelee import TreeSet, TopologySet
from hypothesis import given, settings, strategies as st
import random


def parse_treeswift(newick):
    """Helper to parse a newick string with treeswift"""
    return treeswift.read_tree_newick(newick)


def test_basic_topology_distances():
    """Test that basic topological distances match between libraries"""
    # Simple balanced tree - using unit branch lengths for topology
    newick = "((A:1,B:1):1,(C:1,D:1):1);"

    # Parse with both libraries
    ts_tree = parse_treeswift(newick)
    td_tree = TopologySet(newick)

    # Get distance matrices
    td_dist = td_tree.get_distance_matrix(0)

    # Get node mapping
    ts_nodes = ts_tree.label_to_node(selection="leaves")

    # Compare distances between all pairs
    taxa = ["A", "B", "C", "D"]
    for i, taxon1 in enumerate(taxa):
        for j, taxon2 in enumerate(taxa[i + 1 :], i + 1):
            node1 = ts_nodes[taxon1]
            node2 = ts_nodes[taxon2]
            ts_dist = ts_tree.distance_between(node1, node2)
            td_dist_val = td_dist.get_by_name(taxon1, taxon2)
            assert (
                abs(ts_dist - td_dist_val) < 1e-10
            ), f"Distance mismatch for {taxon1}->{taxon2}: treeswift={ts_dist}, teedeelee={td_dist_val}"


def test_branch_length_distances():
    """Test that branch length distances match between libraries"""
    # Tree with unit branch lengths for topological comparison
    newick = "(A:1,B:1,(C:1,D:1):1);"

    # Parse with both libraries
    ts_tree = parse_treeswift(newick)
    td_tree = TreeSet(newick)

    # Get distance matrices and node mapping
    td_dist = td_tree.get_distance_matrix(0)
    ts_nodes = ts_tree.label_to_node(selection="leaves")

    # Compare distances between all pairs
    taxa = ["A", "B", "C", "D"]
    for i, taxon1 in enumerate(taxa):
        for j, taxon2 in enumerate(taxa[i + 1 :], i + 1):
            node1 = ts_nodes[taxon1]
            node2 = ts_nodes[taxon2]
            ts_dist = ts_tree.distance_between(node1, node2)
            td_dist_val = td_dist.get_by_name(taxon1, taxon2)
            assert (
                abs(ts_dist - td_dist_val) < 1e-10
            ), f"Distance mismatch for {taxon1}->{taxon2}: treeswift={ts_dist}, teedeelee={td_dist_val}"


def retopologize_ts_tree(ts_tree):
    """Helper to set all non-root edge lengths to 1 in a TreeSwift tree"""
    for node in ts_tree.traverse_postorder():
        if not node.is_root():
            node.edge_length = 1.0
    return ts_tree


def test_topology_restriction():
    """Test that tree restriction operations match between libraries"""
    # Original tree with unit branch lengths
    newick = "((A:1,B:1):1,(C:1,D:1):1);"

    # Parse with both libraries
    ts_tree = parse_treeswift(newick)
    td_tree = TopologySet(newick)

    # Restrict to subset of taxa
    subset = ["A", "C", "D"]

    # TreeSwift restriction with retopologizing
    ts_restricted = retopologize_ts_tree(ts_tree.extract_tree_with(set(subset)))
    ts_nodes = ts_restricted.label_to_node(selection="leaves")

    # Teedeelee restriction
    td_restricted = td_tree.restriction(subset)
    td_dist = td_restricted.get_distance_matrix(0)

    # Compare distances in restricted trees
    for i, taxon1 in enumerate(subset):
        for j, taxon2 in enumerate(subset[i + 1 :], i + 1):
            node1 = ts_nodes[taxon1]
            node2 = ts_nodes[taxon2]
            ts_dist = ts_restricted.distance_between(node1, node2)
            td_dist_val = td_dist.get_by_name(taxon1, taxon2)
            assert (
                abs(ts_dist - td_dist_val) < 1e-10
            ), f"Distance mismatch after restriction for {taxon1}->{taxon2}: treeswift={ts_dist}, teedeelee={td_dist_val}"


def test_complex_topology():
    """Test distances on a more complex tree structure"""
    # Using unit branch lengths for topology comparison
    newick = "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);"

    # Parse with both libraries
    ts_tree = parse_treeswift(newick)
    td_tree = TopologySet(newick)

    # Get distance matrices and node mapping
    td_dist = td_tree.get_distance_matrix(0)
    ts_nodes = ts_tree.label_to_node(selection="leaves")

    # Compare distances between some interesting pairs
    test_pairs = [
        ("A", "B"),  # Close pair
        ("A", "H"),  # Distant pair
        ("C", "D"),  # Close pair in different subtree
        ("E", "G"),  # Moderately distant pair
    ]

    for taxon1, taxon2 in test_pairs:
        node1 = ts_nodes[taxon1]
        node2 = ts_nodes[taxon2]
        ts_dist = ts_tree.distance_between(node1, node2)
        td_dist_val = td_dist.get_by_name(taxon1, taxon2)
        assert (
            abs(ts_dist - td_dist_val) < 1e-10
        ), f"Distance mismatch for {taxon1}->{taxon2}: treeswift={ts_dist}, teedeelee={td_dist_val}"


def test_error_handling():
    """Test error handling for invalid inputs"""
    newick = "((A:1,B:1):1,(C:1,D:1):1);"
    td_tree = TopologySet(newick)

    # Test restriction with invalid taxon
    with pytest.raises(Exception):
        td_tree.restriction(["A", "X"])  # X doesn't exist

    # Test distance matrix with invalid taxon
    td_dist = td_tree.get_distance_matrix(0)
    with pytest.raises(Exception):
        td_dist.get_by_name("A", "X")  # X doesn't exist


def test_unrooted_tree():
    """Test distances on an unrooted tree"""
    # Unrooted tree with unit branch lengths
    newick = "(A:1,B:1,C:1,D:1);"

    ts_tree = parse_treeswift(newick)
    td_tree = TreeSet(newick)

    td_dist = td_tree.get_distance_matrix(0)
    ts_nodes = ts_tree.label_to_node(selection="leaves")

    taxa = ["A", "B", "C", "D"]
    for i, taxon1 in enumerate(taxa):
        for j, taxon2 in enumerate(taxa[i + 1 :], i + 1):
            node1 = ts_nodes[taxon1]
            node2 = ts_nodes[taxon2]
            ts_dist = ts_tree.distance_between(node1, node2)
            td_dist_val = td_dist.get_by_name(taxon1, taxon2)
            assert (
                abs(ts_dist - td_dist_val) < 1e-10
            ), f"Distance mismatch for {taxon1}->{taxon2}: treeswift={ts_dist}, teedeelee={td_dist_val}"


def generate_random_newick(taxa_list, current_depth=0, max_depth=3):
    """Helper to generate a random newick string with unit branch lengths"""
    if len(taxa_list) <= 2:
        # Base case: create a simple subtree
        random.shuffle(taxa_list)
        return f"({','.join(f'{taxon}:1' for taxon in taxa_list)})"

    if current_depth >= max_depth:
        # At max depth, create a multifurcating node with all remaining taxa
        random.shuffle(taxa_list)
        return f"({','.join(f'{taxon}:1' for taxon in taxa_list)})"

    # Randomly decide between bifurcation and multifurcation
    n_children = random.randint(2, min(4, len(taxa_list)))

    # Split taxa into n_children groups
    random.shuffle(taxa_list)
    split_points = sorted(
        [0]
        + [random.randint(1, len(taxa_list) - 1) for _ in range(n_children - 1)]
        + [len(taxa_list)]
    )
    groups = [
        taxa_list[split_points[i] : split_points[i + 1]]
        for i in range(len(split_points) - 1)
    ]

    # Filter out empty groups
    groups = [g for g in groups if g]

    # Ensure at least 2 groups
    while len(groups) < 2:
        # Find largest group to split
        largest_idx = max(range(len(groups)), key=lambda i: len(groups[i]))
        if len(groups[largest_idx]) > 1:
            split = random.randint(1, len(groups[largest_idx]) - 1)
            new_group = groups[largest_idx][split:]
            groups[largest_idx] = groups[largest_idx][:split]
            groups.append(new_group)

    # Recursively build subtrees
    subtrees = []
    for group in groups:
        if len(group) == 1:
            # Single taxon - add directly without recursion
            subtrees.append(f"{group[0]}:1")
        else:
            # Multiple taxa - recurse
            subtrees.append(generate_random_newick(group, current_depth + 1, max_depth))

    return f"({','.join(subtrees)})"


@st.composite
def tree_taxa_strategy(draw):
    """Strategy to generate a list of unique taxa names"""
    n_taxa = draw(st.integers(min_value=3, max_value=15))
    return [f"t{i}" for i in range(n_taxa)]


@settings(max_examples=100)  # Adjust number of test cases as needed
@given(taxa=tree_taxa_strategy())
def test_random_tree_distances(taxa):
    """Test that distances match between libraries for randomly generated trees"""
    # Generate a random tree structure
    newick = generate_random_newick(taxa) + ";"

    # Parse with both libraries
    ts_tree = retopologize_ts_tree(parse_treeswift(newick))
    td_tree = TreeSet(newick)

    # Get distance matrices and node mapping
    td_dist = td_tree.get_distance_matrix(0)
    ts_nodes = ts_tree.label_to_node(selection="leaves")

    # Compare distances between all pairs
    for i, taxon1 in enumerate(taxa):
        for j, taxon2 in enumerate(taxa[i + 1 :], i + 1):
            node1 = ts_nodes[taxon1]
            node2 = ts_nodes[taxon2]
            ts_dist = ts_tree.distance_between(node1, node2)
            td_dist_val = td_dist.get_by_name(taxon1, taxon2)
            assert abs(ts_dist - td_dist_val) < 1e-10, (
                f"Distance mismatch for {taxon1}->{taxon2} in tree {newick}: "
                f"treeswift={ts_dist}, teedeelee={td_dist_val}; treeswift_newick={ts_tree.newick()}, teedeelee_newick={td_tree.to_newick_string()}"
            )


@settings(max_examples=50)  # Reduce number of examples since each test is more complex
@given(taxa=tree_taxa_strategy())
def test_random_tree_restriction_distances(taxa):
    """Test that distances are preserved after restriction operations"""
    # Generate a random tree structure
    newick = generate_random_newick(taxa) + ";"

    # Parse with both libraries
    ts_tree = parse_treeswift(newick)
    td_tree = TopologySet(newick)

    # Choose a random subset of taxa (at least 3 taxa)
    subset_size = random.randint(3, len(taxa))
    subset = random.sample(taxa, subset_size)

    # Get restricted trees
    ts_restricted = retopologize_ts_tree(ts_tree.extract_tree_with(set(subset)))
    td_restricted = td_tree.restriction(subset)

    # Get distance matrices for restricted trees
    td_dist = td_restricted.get_distance_matrix(0)
    ts_nodes = ts_restricted.label_to_node(selection="leaves")

    # Compare distances between all pairs in the restricted trees
    for i, taxon1 in enumerate(subset):
        for j, taxon2 in enumerate(subset[i + 1 :], i + 1):
            node1 = ts_nodes[taxon1]
            node2 = ts_nodes[taxon2]
            ts_dist = ts_restricted.distance_between(node1, node2)
            td_dist_val = td_dist.get_by_name(taxon1, taxon2)
            assert abs(ts_dist - td_dist_val) < 1e-10, (
                f"Distance mismatch after restriction for {taxon1}->{taxon2} in tree {newick}: "
                f"treeswift={ts_dist}, teedeelee={td_dist_val}; treeswift_newick={ts_restricted.newick()}, teedeelee_newick={td_restricted.to_newick_string()}"
            )


def test_restriction_edge_cases():
    """Test edge cases for tree restriction operations"""
    newick = "((A:1,B:1):1,(C:1,D:1):1);"
    td_tree = TopologySet(newick)

    # Test restriction to minimum size (2 taxa)
    min_subset = ["A", "B"]
    restricted = td_tree.restriction(min_subset)
    dist_matrix = restricted.get_distance_matrix(0)
    assert (
        dist_matrix.get_by_name("A", "B") == 2.0
    ), "Incorrect distance after restriction to minimum subset"

    # Test restriction preserving internal node
    subset_with_internal = ["A", "B", "C"]
    restricted = td_tree.restriction(subset_with_internal)
    dist_matrix = restricted.get_distance_matrix(0)
    assert (
        dist_matrix.get_by_name("A", "B") == 2.0
    ), "Distance between close pair incorrect"
    assert dist_matrix.get_by_name("A", "C") == 3.0, "Distance across tree incorrect"
    assert dist_matrix.get_by_name("B", "C") == 3.0, "Distance across tree incorrect"

from teedeelee import (
    MSCTreeSet,
    SortCriterion,
    SortCriterionKind,
    SortOrder,
    SortBy,
    FamilyOfMSC,
)
from smallperm import PseudoRandomPermutation as PRP
import treeswift as ts
import logging

import cramjam

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

MAX_TREES_TO_LOG = 2  # Constant to control number of trees logged


def only_topology(tree: ts.Tree) -> ts.Tree:
    """Remove branch lengths and internal node labels to show only topology."""
    for node in tree.traverse_postorder():
        node.edge_length = None
        if node.num_children() > 0:
            node.label = None
    return tree


def normalize_tree(tree: ts.Tree) -> ts.Tree:
    """
    Normalize a tree by ordering nodes using TreeSwift's ordering method.

    Args:
        tree: TreeSwift Tree object to normalize

    Returns:
        Normalized TreeSwift Tree
    """
    tree.order("num_descendants_then_label")
    tree.suppress_unifurcations()
    return tree


def rust_tree_to_newick(tree_str: str) -> str:
    """
    Convert Rust tree string to proper Newick format.

    Args:
        tree_str: Tree string from Rust implementation

    Returns:
        Properly formatted Newick string
    """
    # Ensure the string ends with a semicolon
    if not tree_str.endswith(";"):
        tree_str += ";"
    return tree_str


def random_projection_rust(
    gene_trees_file: str, species_tree_file: str, target_dimension: int, seed: int = 0
):
    """
    Performs random projection using the Rust implementation (teedeelee)

    Args:
        gene_trees_file: Path to file containing gene trees
        species_tree_file: Path to file containing species tree
        target_dimension: Number of taxa to keep in projection
        seed: Random seed for reproducibility

    Returns:
        Tuple of (projected gene trees, projected species tree)
    """
    logger.info("=== Rust Implementation Steps ===")

    # Load trees using Rust implementation
    msctree = MSCTreeSet.from_files(gene_trees_file, species_tree_file)
    family = FamilyOfMSC.from_files([gene_trees_file], [species_tree_file])
    logger.info(f"Loaded {msctree.ngenes} gene trees with {msctree.ntaxa} taxa")

    import pickle as pkl

    # Number of bytes in the pickled object
    dumped = pkl.dumps(msctree)
    logger.info(
        f"Pickled object size: {len(dumped)} bytes. This is {len(dumped) / 1024 / 1024} MB."
    )
    # Get all taxon names
    names = msctree.taxon_set.names()
    N = len(names)
    logger.info(f"Total taxa: {N}")
    logger.info(f"Original taxa names: {names[:5]}...")

    logger.info(f"Family of MSC: {family}")
    logger.info(f"Family of MSC size: {len(family)}")
    dumped = pkl.dumps(family)
    logger.info(
        f"Pickled family of MSC size: {len(dumped)} bytes. This is {len(dumped) / 1024 / 1024} MB."
    )

    sorted_names = sorted(names)
    logger.info(f"Sorted taxa names: {sorted_names[:5]}...")

    # Generate random subset of taxa
    prp = PRP(N, seed)
    subset = [sorted_names[i] for _, i in zip(range(N // 10), prp)]
    logger.info(f"Selected taxa subset: {subset}")

    # Create index mapping for taxa
    mapping = {str(j): str(i) for i, j in enumerate(subset)}
    logger.info(f"Taxa mapping: {mapping}")

    # Get projected trees
    projected_trees = msctree.as_topology().restriction(subset)
    logger.info(f"After restriction - first {MAX_TREES_TO_LOG} trees:")
    for i in range(min(MAX_TREES_TO_LOG, len(projected_trees))):
        logger.info(f"  Tree {i}: {projected_trees[i]}")

    # logger.info("Projected species tree:")
    logger.info(
        "Projected species tree: "
        + normalize_tree(
            ts.read_tree_newick(projected_trees.species_tree_to_newick())
        ).newick()
    )

    remapped_trees = projected_trees.remap(mapping)
    logger.info(f"After remapping - first {MAX_TREES_TO_LOG} trees:")
    for i in range(min(MAX_TREES_TO_LOG, len(remapped_trees))):
        logger.info(f"  Tree {i}: {remapped_trees[i]}")

    # Convert to TreeSwift trees and normalize them
    swift_trees = [
        normalize_tree(ts.read_tree_newick(rust_tree_to_newick(str(tree))))
        for _, tree in zip(range(MAX_TREES_TO_LOG), remapped_trees)
    ]
    logger.info(f"After normalization - first {MAX_TREES_TO_LOG} trees:")
    for i in range(min(MAX_TREES_TO_LOG, len(swift_trees))):
        logger.info(f"  Tree {i}: {swift_trees[i].newick()}")

    # Get and sort species tree
    species_tree = remapped_trees.get_species_tree()
    logger.info(f"Species tree before sorting: {species_tree}")

    sorted_stree = species_tree.sort_by_multiple(
        [(SortBy.DescendantCount, True), (SortBy.LexicographicalOrder, True)]
    )
    logger.info(f"Final sorted species tree: {sorted_stree}")

    return swift_trees, sorted_stree


def random_projection_treeswift(
    gene_trees_file: str, species_tree_file: str, target_dimension: int, seed: int = 0
):
    """
    Performs random projection using TreeSwift implementation

    Args:
        gene_trees_file: Path to file containing gene trees
        species_tree_file: Path to file containing species tree
        target_dimension: Number of taxa to keep in projection
        seed: Random seed for reproducibility

    Returns:
        Tuple of (projected gene trees, projected species tree)
    """
    logger.info("\n=== TreeSwift Implementation Steps ===")

    # Read input files
    with open(gene_trees_file) as f:
        gene_trees = [line.strip() for line in f if line.strip()]

    with open(species_tree_file) as f:
        species_tree = f.read().strip()

    logger.info(f"Loaded {len(gene_trees)} gene trees")
    logger.info(f"Original species tree: {species_tree}")

    # Parse Newick strings into TreeSwift objects and clean them
    gene_tree_objects = [
        only_topology(ts.read_tree_newick(tree)) for tree in gene_trees
    ]
    species_tree_object = only_topology(ts.read_tree_newick(species_tree))

    # Collect all taxa from species tree
    taxa = {leaf.label for leaf in species_tree_object.traverse_leaves()}
    taxa_sorted = sorted(list(taxa))
    logger.info(f"Total taxa: {len(taxa)}")
    logger.info(f"Original taxa names: {list(taxa)[:5]}...")

    # Select random subset of taxa
    prp = PRP(len(taxa), seed)
    selected_taxa = [taxa_sorted[i] for i in prp[:target_dimension]]
    logger.info(f"Selected taxa subset: {selected_taxa}")

    taxa_to_index = {taxon: str(i) for i, taxon in enumerate(selected_taxa)}
    logger.info(f"Taxa mapping: {taxa_to_index}")

    # Project gene trees
    projected_gene_trees = []
    for i, tree in enumerate(gene_tree_objects):
        subtree = tree.extract_tree_with(selected_taxa)
        if i < MAX_TREES_TO_LOG:
            logger.info(f"Tree {i} after restriction: {subtree.newick()}")
        subtree.rename_nodes(taxa_to_index)
        if i < MAX_TREES_TO_LOG:
            logger.info(f"Tree {i} after remapping: {subtree.newick()}")
        # Normalize the tree before adding it
        normalize_tree(subtree)
        projected_gene_trees.append(subtree)

    # Project species tree
    projected_species_tree = species_tree_object.extract_tree_with(selected_taxa)
    logger.info(f"Species tree after restriction: {projected_species_tree.newick()}")

    projected_species_tree.rename_nodes(taxa_to_index)
    logger.info(f"Species tree after remapping: {projected_species_tree.newick()}")

    projected_species_tree.order("num_descendants_then_label")
    logger.info(f"Final sorted species tree: {projected_species_tree.newick()}")

    return projected_gene_trees, projected_species_tree


if __name__ == "__main__":
    # Example usage
    gtrees_file = "examples/gtrees.tre"
    stree_file = "examples/stree.tre"
    target_dim = 10

    # Compare both implementations
    rust_trees, rust_stree = random_projection_rust(gtrees_file, stree_file, target_dim)
    swift_trees, swift_stree = random_projection_treeswift(
        gtrees_file, stree_file, target_dim
    )

    print("\nFinal Comparison:")
    print("Rust implementation results:")
    print("Species tree:")
    print(rust_stree)
    print("\nFirst few gene trees:")
    for i in range(min(MAX_TREES_TO_LOG, len(rust_trees))):
        print(f"Tree {i}: {rust_trees[i].newick()}")

    print("\nTreeSwift implementation results:")
    print("Species tree:")
    print(swift_stree.newick())
    print("\nFirst few gene trees:")
    for i in range(min(MAX_TREES_TO_LOG, len(swift_trees))):
        print(f"Tree {i}: {swift_trees[i].newick()}")

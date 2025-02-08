import treeswift as ts
from teedeelee import Tree, TreeSet, TopologySet
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


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
    ts_tree.order("num_descendants_then_label")
    ts_tree.suppress_unifurcations()
    normalized = ts_tree.newick()

    # Clean up the normalized string
    normalized = normalized.strip(";")
    normalized = normalized.replace(":1.0", "")

    while normalized.startswith("((") and normalized.endswith("))"):
        normalized = normalized[1:-1]

    return normalized.strip()


def main():
    # A simple test tree
    # tree_str = "(0,((100,10),1));"
    tree_str = "(((1,10),100),0);"

    # Create both TreeSwift and teedeelee trees
    ts_tree = ts.read_tree_newick(tree_str)
    td_tree = Tree(tree_str)

    print(ts_tree.newick())
    print(td_tree.newick())

    print("Original tree (normalized):", normalize_tree_str(tree_str))

    # Calculate distances
    # TreeSwift distances
    for node in ts_tree.traverse_postorder():
        if not node.is_root():
            node.edge_length = 1
    ts_dist_dict = ts_tree.distance_matrix(leaf_labels=True)

    # teedeelee distances
    td_dist = td_tree.get_distance_matrix()

    print(td_dist)

    print("\nDistances between 0 and 1:")
    print(f"TreeSwift distance: {ts_dist_dict['0']['1']}")
    print(f"teedeelee distance: {td_dist.get_by_name('0', '1')}")


if __name__ == "__main__":
    main()

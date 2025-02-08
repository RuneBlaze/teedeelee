import treeswift as ts
from teedeelee import Tree, SortBy
from demo import normalize_tree

def test_tree_normalization():
    # Define the two test trees
    tree1_newick = "(7,(9,((1,(5,6)),(8,((2,4),(0,3))))));";
    tree2_newick = "(7,(9,((1,(5,6)),(8,((0,3),(2,4))))));";

    # Parse and normalize both trees
    tree1 = normalize_tree(ts.read_tree_newick(tree1_newick))
    tree2 = normalize_tree(ts.read_tree_newick(tree2_newick))

    # Get Newick strings for comparison
    newick1 = tree1.newick()
    newick2 = tree2.newick()

    print("Tree 1 normalized:", newick1)
    print("Tree 2 normalized:", newick2)
    print("\nTrees are identical after normalization:", newick1 == newick2)

def test_tree_normalization_teedeelee():
    # Define the two test trees
    tree1_newick = "(7,(9,((1,(5,6)),(8,((2,4),(0,3))))));";
    tree2_newick = "(7,(9,((1,(5,6)),(8,((0,3),(2,4))))));";

    # Parse trees using teedeelee and sort them
    tree1 = Tree(tree1_newick).sort_by_multiple([
        (SortBy.DescendantCount, True),
        (SortBy.LexicographicalOrder, True)
    ])
    tree2 = Tree(tree2_newick).sort_by_multiple([
        (SortBy.DescendantCount, True),
        (SortBy.LexicographicalOrder, True)
    ])

    # Get string representations for comparison
    str1 = str(tree1)
    str2 = str(tree2)

    print("\nTeedeelee implementation:")
    print("Tree 1 normalized:", str1)
    print("Tree 2 normalized:", str2)
    print("Trees are identical after normalization:", str1 == str2)

if __name__ == "__main__":
    test_tree_normalization()
    test_tree_normalization_teedeelee() 
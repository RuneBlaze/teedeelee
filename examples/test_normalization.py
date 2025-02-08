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

def test_tree_normalization_random_jiggles():
    import random

    def tuple_to_newick(tree):
        """Convert a tuple-based tree to Newick format."""
        if isinstance(tree, tuple):
            return f"({','.join(tuple_to_newick(child) for child in tree)})"
        return str(tree)

    def random_tuple_tree(leaves):
        """Generate a random binary tree structure with given leaf labels."""
        if len(leaves) == 1:
            return leaves[0]
        if len(leaves) == 2:
            return (leaves[0], leaves[1]) if random.random() < 0.5 else (leaves[1], leaves[0])
        
        # Randomly split leaves into two groups
        split_point = random.randint(1, len(leaves)-1)
        left_leaves = leaves[:split_point]
        right_leaves = leaves[split_point:]
        
        return (random_tuple_tree(left_leaves), random_tuple_tree(right_leaves))

    def random_jiggle(tree):
        """Randomly swap children in a tuple-based tree."""
        if not isinstance(tree, tuple):
            return tree
        
        # Convert to list for modification
        children = list(tree)
        # Recursively jiggle children
        children = [random_jiggle(child) for child in children]
        # Randomly swap if binary node
        if len(children) == 2 and random.random() < 0.5:
            children[0], children[1] = children[1], children[0]
        
        return tuple(children)

    # Generate base tree structure
    leaves = list(range(10))  # 0-9 as leaf labels
    base_tree = random_tuple_tree(leaves)
    
    num_variants = 5
    variants = []
    
    print("\nTuple-based random jiggle test:")
    print(f"Base tree: {tuple_to_newick(base_tree)}")
    
    for i in range(num_variants):
        # Create a jiggled variant
        jiggled = random_jiggle(base_tree)
        variants.append(jiggled)
        print(f"Variant {i+1}: {tuple_to_newick(jiggled)}")
        
        # Convert to Tree object and normalize
        tree = Tree(tuple_to_newick(jiggled))
        normalized = tree.sort_by_multiple([
            (SortBy.DescendantCount, True),
            (SortBy.LexicographicalOrder, True)
        ])
        print(f"Normalized {i+1}: {normalized}")

    # Convert all variants to Tree objects and normalize
    normalized_variants = [
        str(Tree(tuple_to_newick(v)).sort_by_multiple([
            (SortBy.DescendantCount, True),
            (SortBy.LexicographicalOrder, True)
        ]))
        for v in variants
    ]
    
    all_same = all(variant == normalized_variants[0] for variant in normalized_variants)
    print("\nAll variants normalized to same tree:", all_same)

if __name__ == "__main__":
    test_tree_normalization()
    test_tree_normalization_teedeelee()
    test_tree_normalization_random_jiggles() 
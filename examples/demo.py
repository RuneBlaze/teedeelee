from teedeelee import TreeSet, TopologySet, MSCTreeSet

def demo_complete_tree():
    print("\n=== Tree Set Demo ===")
    # Example Newick string with branch lengths
    newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    
    # Create a TreeSet
    trees = TreeSet(newick)
    
    print(f"Number of taxa: {trees.ntaxa}")
    print(f"Original tree: {trees.to_newick_string()}")
    
    # Get distance matrix
    dist_matrix = trees.get_distance_matrix(0)
    if dist_matrix:
        print("\nPairwise distances:")
        for i in range(dist_matrix.ntaxa):
            for j in range(i + 1, dist_matrix.ntaxa):
                name1 = dist_matrix.get_taxon_name(i)
                name2 = dist_matrix.get_taxon_name(j)
                distance = dist_matrix[i, j]
                print(f"{name1} -> {name2}: {distance}")

def demo_topology_tree():
    print("\n=== Topology Set Demo ===")
    # Example Newick string
    newick = "((A,B),(C,D));"
    
    # Create a TopologySet
    trees = TopologySet(newick)
    
    print(f"Original tree: {trees.to_newick_string()}")
    
    # Demonstrate restriction to subset of taxa
    restricted = trees.restriction(["A", "C", "D"])
    print(f"Restricted to [A,C,D]: {restricted.to_newick_string()}")

def demo_msc_tree():
    print("\n=== MSC Tree Set Demo ===")
    # Example gene trees and species tree
    gene_trees = "((A,B),(C,D));((A,C),(B,D));"
    species_tree = "((A,B),(C,D));"
    
    # Create MSC tree set
    msc = MSCTreeSet(gene_trees, species_tree)
    
    print(f"Number of gene trees: {msc.ngenes}")
    print(f"Species tree: {msc.species_tree_to_newick()}")
    print(f"Gene trees: {msc.to_newick_vec()}")
    
    # Get species distance matrix
    sp_dist = msc.get_species_distance_matrix()
    print("\nSpecies pairwise distances:")
    for i in range(sp_dist.ntaxa):
        for j in range(i + 1, sp_dist.ntaxa):
            name1 = sp_dist.get_taxon_name(i)
            name2 = sp_dist.get_taxon_name(j)
            distance = sp_dist[i, j]
            print(f"{name1} -> {name2}: {distance}")

def main():
    print("Teedeelee Demo\n")
    
    try:
        demo_complete_tree()
        demo_topology_tree()
        demo_msc_tree()
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main() 
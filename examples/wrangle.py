import glob
import os
import pickle
import time
from contextlib import contextmanager
from teedeelee import FamilyOfMSC
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@contextmanager
def timer(description: str):
    """Context manager for timing code blocks and logging the duration.

    Args:
        description: Description of the operation being timed
    """
    start = time.perf_counter()
    yield
    elapsed = time.perf_counter() - start
    logger.info(f"{description}: {elapsed:.2f} seconds")


def get_pickle_size(obj) -> tuple[int, float]:
    """Get the size of a pickled object in bytes and MB.

    Args:
        obj: Object to pickle

    Returns:
        Tuple of (size in bytes, size in MB)
    """
    pickled = pickle.dumps(obj)
    size_bytes = len(pickled)
    size_mb = size_bytes / (1024 * 1024)
    return size_bytes, size_mb


def process_s101_data(base_dir: str = "/Users/lbq/Downloads/S101"):
    """
    Process the S101 data directory structure to create a single FamilyOfMSC object
    containing all datasets.

    Args:
        base_dir: Base directory containing the numbered subdirectories

    Returns:
        FamilyOfMSC object containing all datasets
    """
    # Get all subdirectories (01, 02, etc.)
    subdirs = sorted(glob.glob(os.path.join(base_dir, "*")))
    gene_trees_files = []
    species_tree_files = []

    for subdir in subdirs:
        if not os.path.isdir(subdir):
            continue

        # Find the required files
        gene_trees_file = os.path.join(subdir, "truegenetrees")
        species_tree_file = os.path.join(subdir, "s_tree.trees")

        if not (os.path.exists(gene_trees_file) and os.path.exists(species_tree_file)):
            logger.warning(f"Missing required files in {subdir}")
            continue

        logger.info(f"\nFound dataset in: {subdir}")
        logger.info(f"Gene trees file: {gene_trees_file}")
        logger.info(f"Species tree file: {species_tree_file}")

        gene_trees_files.append(gene_trees_file)
        species_tree_files.append(species_tree_file)

    logger.info(f"\nCreating FamilyOfMSC from {len(gene_trees_files)} datasets...")

    # Create single FamilyOfMSC object with timing
    with timer("Creating FamilyOfMSC"):
        family = FamilyOfMSC.from_files(gene_trees_files, species_tree_files)

    # Get pickle size
    bytes_size, mb_size = get_pickle_size(family)
    logger.info(f"Pickled family size: {bytes_size:,} bytes ({mb_size:.2f} MB)")
    logger.info(f"Created FamilyOfMSC with {len(family)} trees")

    return family


def save_family(family: FamilyOfMSC, output_path: str = "processed_family.pkl"):
    """Save processed family to a pickle file."""
    with timer(f"Saving family to {output_path}"):
        with open(output_path, "wb") as f:
            pickle.dump(family, f)


if __name__ == "__main__":
    with timer("Total processing time"):
        family = process_s101_data()

        # Get and log size information
        bytes_size, mb_size = get_pickle_size(family)
        logger.info(f"\nFinal FamilyOfMSC:")
        logger.info(f"Number of trees: {len(family)}")
        logger.info(f"Pickled size: {bytes_size:,} bytes ({mb_size:.2f} MB)")

        # Save the processed family
        save_family(family)

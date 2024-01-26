import os
import tempfile
from markerrepo.annotation import write_annotation

def test_write_annotation():
    # Creating a sample ct_dict with various scenarios:
    # - Cluster1 has a clear highest score.
    # - Cluster2 has two cell types with the same score but different match counts.
    # - Cluster3 is empty to test handling of clusters with no cell types.
    ct_dict = {
        "Cluster1": {
            "CellTypeA": [100, 5, 10, 8],
            "CellTypeB": [95, 4, 10, 7]
        },
        "Cluster2": {
            "CellTypeC": [90, 5, 10, 6],
            "CellTypeD": [90, 6, 10, 6]  # Same score as CellTypeC but higher match_count
        },
        "Cluster3": {
            # No cell types for this cluster
        }
    }

    # Using a temporary directory to store output files
    with tempfile.TemporaryDirectory() as tmpdir:
        output = tmpdir
        write_annotation(ct_dict, output)

        # Checking the main annotation file
        with open(os.path.join(output, "annotation.txt"), "r") as f:
            lines = f.readlines()
            assert lines[0].strip() == "Cluster1\tCellTypeA"
            assert lines[1].strip() == "Cluster2\tCellTypeD"
            assert len(lines) == 2  # Cluster3 should not be present since it has no cell types

        # Checking the detail files for each cluster
        for cluster in ct_dict:
            if ct_dict[cluster]:
                with open(os.path.join(output, f"ranks/cluster_{cluster}.txt"), "r") as f:
                    detail_lines = f.readlines()
                    assert len(detail_lines) == len(ct_dict[cluster])
                    # Ensuring each line in the detail file starts with the correct cell type
                    # in the order they are expected to be in
                    for i, cell_type in enumerate(sorted(ct_dict[cluster], key=lambda x: (-ct_dict[cluster][x][0], -ct_dict[cluster][x][1]))):
                        assert detail_lines[i].startswith(cell_type)

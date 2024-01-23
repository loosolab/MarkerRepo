import pytest
from markerrepo import wrappers as wrap

# Check if an exception is raised when neither marker_repo nor SCSA is enabled
def test_no_method_selected():
    with pytest.raises(ValueError, match="At least one of 'marker_repo' or 'SCSA' must be True."):
        wrap.run_annotation(adata=None, marker_repo=False, SCSA=False)

# Check if an exception is raised when no marker lists are provided
def test_no_marker_lists():
    with pytest.raises(ValueError, match="No marker lists provided. Please provide a list of marker list paths."):
        wrap.run_annotation(adata=None, marker_repo=True, SCSA=True, marker_lists=[])

#  Check if an exception is raised when a non-existent marker list file is provided
def test_nonexistent_marker_list():
    with pytest.raises(FileNotFoundError, match=r"Marker list file not found: .+"):
        wrap.run_annotation(adata=None, marker_repo=True, SCSA=True, marker_lists=["nonexistent_file.csv"])

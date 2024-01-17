import pytest
from markerrepo.wrappers import create_marker_lists

# Test for successful list creation with valid inputs
def test_create_marker_lists_success():
    repo_path = '/mnt/workspace/mkessle/projects/annotate_by_marker_and_features'
    output_path = f'{repo_path}/tests/data'
    filename = 'test_create_marker_lists_success'   
    column_specific_terms = {"Source": "panglao.se"}
    result = create_marker_lists(organism='mouse', repo_path=repo_path, style='score', column_specific_terms=column_specific_terms, 
                                 path=output_path, file_name=filename)

    assert isinstance(result, list), "Expected a list of paths"

# Test for handling invalid style
def test_invalid_style():
    with pytest.raises(ValueError) as excinfo:
        create_marker_lists(style='invalid_style')
    assert "The parameter 'style' must be one of 'two_column', 'score', 'ui', or 'panglao'" in str(excinfo.value)

# Test for file not found
def test_repo_path_not_found():
    with pytest.raises(FileNotFoundError) as excinfo:
        create_marker_lists(repo_path='/invalid/path')
    assert "The specified repository path '/invalid/path' does not exist." in str(excinfo.value)
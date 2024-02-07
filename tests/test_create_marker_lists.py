import pytest
from markerrepo import wrappers as wrap

repo_path = "."

# Test for successful list creation with valid inputs
def test_create_marker_lists_success():
    output_path = f'{repo_path}/tests/data'
    filename = 'test_create_marker_lists_success'   
    column_specific_terms = {"Source": "panglao.se"}
    result = wrap.create_marker_lists(organism='mouse', repo_path=repo_path, style='score', column_specific_terms=column_specific_terms, 
                                 path=output_path, file_name=filename)

    assert isinstance(result, list), "Expected a list of paths"

# Test for handling invalid style
def test_invalid_style():
    with pytest.raises(ValueError):
        wrap.create_marker_lists(style='invalid_style')

# Test for file not found
def test_repo_path_not_found():
    with pytest.raises(FileNotFoundError):
        wrap.create_marker_lists(repo_path='/invalid/path')

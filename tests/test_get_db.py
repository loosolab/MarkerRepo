import pytest
from unittest import mock
import pandas as pd
from markerrepo import marker_repo as mr

repo_path = "."

def test_get_db_success():
    # Assuming you have a valid path with YAML files
    valid_path = repo_path
    df = mr.get_db(repo_path=valid_path)
    assert isinstance(df, pd.DataFrame), "Expected a DataFrame"


def test_get_db_path_not_found():
    with pytest.raises(FileNotFoundError):
        mr.get_db(repo_path="/invalid/path")


def test_get_db_no_yaml_found():
    valid_path = repo_path
    
    # Mock os.path.exists to return True
    with mock.patch('os.path.exists', return_value=True):
        # Mock os.walk to return an empty list, simulating no YAML files found
        with mock.patch('os.walk', return_value=[]):
            with pytest.raises(FileNotFoundError):
                mr.get_db(repo_path=valid_path)
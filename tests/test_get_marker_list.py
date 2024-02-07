import pytest
from unittest import mock
import pandas as pd
from  markerrepo import marker_repo as mr

def test_get_marker_list():
    # Simulated YAML data object
    simulated_yaml_data = {
        "marker_list": [
            {"name": "Cell type A", "markers": ["Gene1", "Gene2"]},
            {"name": "Cell type B", "markers": ["Gene3"]}
        ]
    }

    # Mock 'open' and 'yaml.safe_load'
    with mock.patch("builtins.open", mock.mock_open()):
        with mock.patch("yaml.safe_load", return_value=simulated_yaml_data):
            df = mr.get_marker_list("irrelevant/file/path.yaml")

    # Check the structure of the expected DataFrame
    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == ["Marker", "Info"]
    assert len(df) == 3  # Three marker entries
    assert df.iloc[0]["Marker"] == "Gene1"
    assert df.iloc[0]["Info"] == "Cell type A"
    assert df.iloc[1]["Marker"] == "Gene2"
    assert df.iloc[1]["Info"] == "Cell type A"
    assert df.iloc[2]["Marker"] == "Gene3"
    assert df.iloc[2]["Info"] == "Cell type B"

import pytest
from unittest.mock import MagicMock
from pandas import DataFrame
from markerrepo.annotation import update_adata_with_markers

# Create a Mock AnnData object
class MockAdata:
    def __init__(self):
        self.uns = {}

# Test functions
def test_add_new_entry():
    # Test adding a new entry to an empty AnnData object
    mock_adata = MockAdata()
    df = DataFrame(index=["marker1", "marker2"])
    update_adata_with_markers(mock_adata, "new_list", df)

    assert "new_list" in mock_adata.uns["MarkerRepo"]["marker_lists"]
    assert mock_adata.uns["MarkerRepo"]["marker_lists"]["new_list"] == ["marker1", "marker2"]


def test_ignore_overwrite():
    # Test the ignore_overwrite functionality where existing entry is overwritten without prompt
    mock_adata = MockAdata()
    mock_adata.uns["MarkerRepo"] = {"marker_lists": {"existing_list": ["old_marker"]}}
    df = DataFrame(index=["new_marker"])
    update_adata_with_markers(mock_adata, "existing_list", df, ignore_overwrite=True)

    assert mock_adata.uns["MarkerRepo"]["marker_lists"]["existing_list"] == ["new_marker"]


def test_user_interaction_overwrite(monkeypatch):
    # Test for overwrite interaction
    mock_adata = MockAdata()
    mock_adata.uns["MarkerRepo"] = {"marker_lists": {"existing_list": ["old_marker"]}}
    df = DataFrame(index=["new_marker"])
    monkeypatch.setattr('builtins.input', lambda _: "O")
    update_adata_with_markers(mock_adata, "existing_list", df)
    assert mock_adata.uns["MarkerRepo"]["marker_lists"]["existing_list"] == ["new_marker"]


def test_user_interaction_update(monkeypatch):
    # Test for update interaction
    mock_adata = MockAdata()
    mock_adata.uns["MarkerRepo"] = {"marker_lists": {"existing_list": ["old_marker"]}}
    df = DataFrame(index=["new_marker"])
    monkeypatch.setattr('builtins.input', lambda _: "U")
    update_adata_with_markers(mock_adata, "existing_list", df)
    assert set(mock_adata.uns["MarkerRepo"]["marker_lists"]["existing_list"]) == {"old_marker", "new_marker"}


def test_user_interaction_leave_unchanged(monkeypatch):
    # Test for leaving unchanged interaction
    mock_adata = MockAdata()
    mock_adata.uns["MarkerRepo"] = {"marker_lists": {"existing_list": ["old_marker"]}}
    df = DataFrame(index=["new_marker"])
    monkeypatch.setattr('builtins.input', lambda _: "L")
    original_list = mock_adata.uns["MarkerRepo"]["marker_lists"]["existing_list"].copy()
    update_adata_with_markers(mock_adata, "existing_list", df)

    assert mock_adata.uns["MarkerRepo"]["marker_lists"]["existing_list"] == original_list


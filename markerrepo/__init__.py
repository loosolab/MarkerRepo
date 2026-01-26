
# import with prefix _ to hide them
import importlib as _importlib

submodules = [
    "annotation",
    #"check_uids",  # command line tool
    "generate_metafile",
    "homology",
    "marker_repo",
    "parsing",
    "plotting",
    "scoring",
    "update_uids",
    "utils",
    "validate_yaml",
    "wrappers"
]

# define what is exported in this module
__all__ = submodules

def __dir__():
    """Return the defined submodules."""
    return __all__


def __getattr__(name):
    """Lazyload modules (inspired by scipy)."""
    if name in submodules:
        # return import to make it directly available
        return _importlib.import_module(f"markerrepo.{name}")
    else:
        raise AttributeError(f"Module 'markerrepo' does not contain a module named '{name}'.")

import importlib
import pytest

SUBMODULES = [
    "npl", "npl.core", "npl.descriptors", "npl.calculators",
    "npl.optimization", "npl.visualize", "npl.utils",
]


@pytest.mark.parametrize("module_name", SUBMODULES)
def test_submodule_imports(module_name):
    importlib.import_module(module_name)


def test_no_monte_carlo_module():
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("npl.monte_carlo")


def test_declared_exports_resolve():
    import npl.descriptors as d
    import npl.optimization as o
    for name in d.__all__:
        assert getattr(d, name) is not None
    for name in o.__all__:
        assert getattr(o, name) is not None

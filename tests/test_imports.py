import importlib
import pytest

SUBMODULES = [
    "npl", "npl.core", "npl.descriptors", "npl.calculators",
    "npl.optimization", "npl.monte_carlo", "npl.visualize", "npl.utils",
]


@pytest.mark.parametrize("module_name", SUBMODULES)
def test_submodule_imports(module_name):
    importlib.import_module(module_name)


def test_monte_carlo_entry_point_importable():
    from npl.monte_carlo import run_monte_carlo
    assert callable(run_monte_carlo)


def test_declared_exports_resolve():
    import npl.descriptors as d
    import npl.optimization as o
    import npl.monte_carlo as mc
    for name in d.__all__:
        assert getattr(d, name) is not None
    for name in o.__all__:
        assert getattr(o, name) is not None
    for name in mc.__all__:
        assert getattr(mc, name) is not None

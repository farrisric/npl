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


def test_import_and_use_leave_logging_and_files_alone(tmp_path):
    # Fresh interpreter: the test session has already imported npl.
    import subprocess
    import sys
    code = (
        'import logging\n'
        'import npl.monte_carlo, npl.optimization\n'
        'from npl.calculators import TOPCalculator\n'
        'from npl.descriptors import ExtendedTopologicalFeaturesClassifier as E\n'
        "TOPCalculator('ETOP', stoichiometry='Pt151Cu50', feature_classifier=E)\n"
        'assert logging.getLogger().handlers == [], logging.getLogger().handlers\n'
    )
    result = subprocess.run([sys.executable, '-W', 'error', '-c', code], cwd=tmp_path,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert list(tmp_path.iterdir()) == []


def test_setup_logging_configures_only_the_npl_logger(tmp_path):
    import logging
    import npl
    root_handlers = list(logging.getLogger().handlers)
    log_file = tmp_path / 'run.log'
    try:
        npl.setup_logging(file=str(log_file))
        logger = npl.setup_logging(file=str(log_file))  # a second call replaces, not stacks
        assert len([h for h in logger.handlers if getattr(h, '_npl_managed', False)]) == 2
        logging.getLogger('npl.monte_carlo').info('hello')
        assert 'INFO - hello' in log_file.read_text()
        assert logging.getLogger().handlers == root_handlers
    finally:
        npl.setup_logging(stream=False)  # removes and closes the handlers
        logging.getLogger('npl').setLevel(logging.NOTSET)
        logging.getLogger('npl').propagate = True

# Contributing to npl

Thanks for your interest in improving `npl`. Contributions of all kinds are
welcome: bug reports, feature requests, documentation, and code.

## Reporting bugs and requesting features

Please open an issue on the
[GitHub issue tracker](https://github.com/farrisric/npl/issues). For bugs,
include what you expected, what happened, a minimal reproducing example, and
your OS, Python version, and `npl` version.

## Development setup

```sh
git clone https://github.com/farrisric/npl.git
cd npl
pip install -e ".[test]"
```

## Running the tests

```sh
pytest
```

CI runs the test suite on Python 3.10-3.13.

## Linting

```sh
flake8 npl/
```

## Pull requests

Keep pull requests focused, add tests for new behavior, and make sure
`pytest` and `flake8` pass before requesting review.

from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="npl",
    version="1.0.3",
    description="Nanoparticle Library for computational analysis of nanoparticles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Riccardo Farris",
    author_email="rfarris@ub.edu",
    url="",
    packages=find_packages(),
    install_requires=[
        "acat>=1.7.1",
        "scikit_learn>=1.5.0",
        "scipy>=1.10.0",
        "sortedcontainers>=2.4.0"
    ],
)

from setuptools import setup, find_packages

setup(
    name="npl",
    version="1.0.0",
    description="Nanoparticle Library for computational analysis of nanoparticles",
    long_description="This is a long description for the Nanoparticle Library.",
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
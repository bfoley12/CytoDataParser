# setup.py
from setuptools import setup, find_packages

setup(
    name="cytodataparser",
    version="0.3.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)

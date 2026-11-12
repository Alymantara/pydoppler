"""
PyDoppler packaging.

Build and upload a release:
    python -m build
    python -m twine check dist/*
    python -m twine upload dist/*
"""
import os
import re

from setuptools import setup

HERE = os.path.abspath(os.path.dirname(__file__))
GITHUB = "https://github.com/alymantara/pydoppler"
GITHUB_BRANCH = "main"   # branch the README images are served from on PyPI


def read_version():
    with open(os.path.join(HERE, "pydoppler", "__init__.py")) as f:
        return re.search(r'^__version__\s*=\s*"([^"]+)"', f.read(), re.M).group(1)


def read_long_description():
    with open(os.path.join(HERE, "README.md"), encoding="utf-8") as f:
        text = f.read()
    # Relative image/link paths work on GitHub but not on PyPI: make them absolute.
    raw = f"https://raw.githubusercontent.com/alymantara/pydoppler/{GITHUB_BRANCH}/"
    blob = f"{GITHUB}/blob/{GITHUB_BRANCH}/"
    text = text.replace('src="docs/images/', f'src="{raw}docs/images/')
    text = text.replace("](notebooks/", f"]({blob}notebooks/")
    return text


setup(
    name="pydoppler",
    version=read_version(),
    author="Juan V. Hernandez Santisteban",
    author_email="jvhs1@st-andrews.ac.uk",
    description=("Doppler tomography in pure Python, based on Henk Spruit's "
                 "maximum-entropy Doppler mapping code"),
    long_description=read_long_description(),
    long_description_content_type="text/markdown",
    url=GITHUB,
    project_urls={
        "Source": GITHUB,
        "Bug Tracker": f"{GITHUB}/issues",
        "Tutorial": f"{GITHUB}/blob/{GITHUB_BRANCH}/notebooks/PyDoppler_tutorial.ipynb",
        "ASCL": "https://ascl.net/2106.003",
    },
    license="MIT",
    packages=["pydoppler"],
    include_package_data=True,
    package_data={
        "pydoppler": ["test_data/sample_script.py", "test_data/ugem99/*"],
    },
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20",
        "scipy>=1.7",
        "matplotlib>=3.5",
    ],
    keywords=["astronomy", "doppler tomography", "maximum entropy",
              "cataclysmic variables", "accretion discs", "spectroscopy"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Operating System :: OS Independent",
    ],
)

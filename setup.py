from os import path
from setuptools import setup, find_packages
import sys
import versioneer


# NOTE: This file must remain Python 2 compatible for the foreseeable future,
# to ensure that we error out properly for people with outdated setuptools
# and/or pip.
min_version = (3, 7)
if sys.version_info < min_version:
    error = """
Lightshow does not support Python {0}.{1}.
Python {2}.{3} and above is required. Check your Python version like so:

python3 --version

This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
Upgrade pip like so:

pip install --upgrade pip
""".format(
        *(sys.version_info[:2] + min_version)
    )
    sys.exit(error)

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.rst"), encoding="utf-8") as readme_file:
    readme = readme_file.read()

with open(path.join(here, "requirements.txt")) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [
        line
        for line in requirements_file.read().splitlines()
        if not line.startswith("#")
    ]

with open(path.join(here, "requirements-dev.txt")) as requirements_file:
    requirements_dev = [
        line
        for line in requirements_file.read().splitlines()
        if not line.startswith("#")
    ]

# pip install -e ".[dev]"
setup(
    name="lightshow",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A one-stop-shop for building input files for computational "
    " spectroscopy.",
    long_description=readme,
    author="Brookhaven National Laboratory",
    author_email="mcarbone@bnl.gov, john.vinson@nist.gov, dlu@bnl.gov, "
    "vorwerk@physik.hu-berlin.de, xiaqu@bnl.gov",
    url="https://github.com/AI-multimodal/Lightshow",
    python_requires=">={}".format(".".join(str(n) for n in min_version)),
    packages=find_packages(exclude=["docs", "lightshow/_tests"]),
    entry_points={},
    include_package_data=True,
    package_data={
        "xanes_bench": [
            # When adding files here, remember to update MANIFEST.in as well,
            # or else they will not be included in the distribution on PyPI!
            # 'path/to/data_file',
        ]
    },
    install_requires=requirements,
    extras_require={"dev": requirements_dev},
    license="Modified BSD (3-clause)",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)

#!/bin/bash

bash scripts/install.sh build
echo "__version__ = '$(dunamai from any --style=pep440 --no-metadata)'" >lightshow/_version.py
flit build
git checkout lightshow/_version.py

#!/bin/bash

build_docs() {

	if [[ "${GITHUB_ACTION_IS_RUNNING}" = 1 ]]; then
		bash scripts/install.sh doc
	fi

	make -C docs/ html

	# Helper when running on local. If not running in a GitHub Actions
	# environment, this will attempt to open index.html with the users'
	# default program
	if [[ -z "${GITHUB_ACTION_IS_RUNNING}" ]]; then
		open docs/build/html/index.html
	fi

}

pip install toml
bash scripts/install.sh
bash scripts/install.sh doc
echo "__version__ = '$(dunamai from any --style=pep440 --no-metadata)'" >lightshow/_version.py
build_docs
git checkout lightshow/_version.py

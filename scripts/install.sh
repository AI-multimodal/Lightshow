#!/bin/bash

install() {
	if [ ! "$1" ]; then
		echo "No argument supplied - installing core dependencies"
		TARGET='["project"]["dependencies"]'
	else
		echo "Target supplied: $1"
		TARGET="['project']['optional-dependencies']['$1']"
	fi

	echo "installing $TARGET"

	python -c "import toml; c = toml.load('pyproject.toml'); print('\n'.join(c$TARGET))" |
		pip install -r /dev/stdin
}

pip install toml
install "$1"

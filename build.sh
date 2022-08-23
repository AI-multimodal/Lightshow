#!/bin/bash

replace_version_in_init () {
    version="$(dunamai from git --no-metadata --style semver)"
    dunamai check "$version" --style semver
    sed_command="s/...  # semantic-version-placeholder/'$version'/g"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "$sed_command" lightshow/__init__.py
    else
        sed -i "$sed_command" lightshow/__init__.py
    fi
    echo "__init__ version set to" "$version"
    export _TMP_VERSION="$version"
}

reverse_replace_version_in_init () {
    sed_command="s/'$_TMP_VERSION'/...  # semantic-version-placeholder/g"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "$sed_command" lightshow/__init__.py
    else
        sed -i "$sed_command" lightshow/__init__.py
    fi
    echo "__init__ version" "$version" "reset to placeholder"
    unset _TMP_VERSION
}

# Good stuff. The poor man's toml parser
# https://github.com/pypa/pip/issues/8049

install_doc_requirements_only () {
    python3 -c 'import toml; c = toml.load("pyproject.toml"); print("\n".join(c["project"]["optional-dependencies"]["doc"]))' | pip install -r /dev/stdin
}

install_test_requirements_only () {
    python3 -c 'import toml; c = toml.load("pyproject.toml"); print("\n".join(c["project"]["optional-dependencies"]["test"]))' | pip install -r /dev/stdin
}

install_requirements() {
    python3 -c 'import toml; c = toml.load("pyproject.toml"); print("\n".join(c["project"]["dependencies"]))' | pip install -r /dev/stdin
}

install_development_requirements () {
    install_requirements
    install_test_requirements_only
    install_doc_requirements_only
}

build_docs () {

    if [[ "${GITHUB_ACTION_IS_RUNNING}" = 1 ]]; then
        flit install --deps=production --extras=doc
    fi

    make -C docs/ html

    if [[ -z "${GITHUB_ACTION_IS_RUNNING}" ]]; then
        open docs/build/html/index.html
    fi
    
}

install_flit_dunamai () {
    pip install flit~=3.7
    pip install dunamai~=1.12
}

echo "Running build script with CLI args:" "$@"

# Check if this is a GitHub Action. This changes the behavior of the scripts
if [[ "${GITHUB_ACTION_IS_RUNNING}" = 1 ]]; then
    install_flit_dunamai  # Assume we have this installed for local tests
    echo "GitHub Action is running"
else
    echo "GitHub Action is not running - assuming local behavior"
fi


# Iterate over all arguments in order
for var in "$@"
do
    echo "Found" "$var"

    if [ "$var" = "docs" ]; then
        replace_version_in_init
        build_docs
        reverse_replace_version_in_init

    elif [ "$var" = "test" ]; then
        flit install --deps=production --extras=test

    elif [ "$var" = "install-dev-requirements" ]; then
        pip install toml
        install_flit_dunamai
        install_development_requirements

    elif [ "$var" = "_CI-test-requirements" ]; then
        pip install toml
        install_requirements
        install_test_requirements_only

    elif [ "$var" = "_CI-docs-requirements" ]; then
        pip install toml
        install_requirements
        install_doc_requirements_only

    # elif [ "$var" = "publish" ]; then
    #     flit_publish
    fi
done



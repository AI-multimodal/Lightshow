#!/bin/bash

replace_version_in_init () {
    version="$(dunamai from any)"
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

# Replace the version string placeholder in the __init__ with the real thing
replace_version_in_init

# Iterate over all arguments in order
for var in "$@"
do
    echo "Found" "$var"
    if [ "$var" = "docs" ]; then
        build_docs
    elif [ "$var" = "test" ]; then
        flit install --deps=production --extras=test
    # elif [ "$var" = "publish" ]; then
    #     flit_publish
    fi
done

# Put the version string placeholder back
reverse_replace_version_in_init

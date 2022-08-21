#!/bin/bash

replace_version_in_init () {
    version="$(dunamai from any)"
    sed_command="s/...  # semantic-version-placeholder/'$version'/g"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "$sed_command" lightshow/__init__.py
    else
        sed -i "$sed_command" lightshow/__init__.py
    fi
    export _TMP_VERSION="$version"
}

reverse_replace_version_in_init () {
    sed_command="s/'$_TMP_VERSION'/...  # semantic-version-placeholder/g"
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sed -i '' "$sed_command" lightshow/__init__.py
    else
        sed -i "$sed_command" lightshow/__init__.py
    fi
    unset _TMP_VERSION
}

build_docs () {

    replace_version_in_init
    make -C docs/ html
    reverse_replace_version_in_init

    # If the GitHub token is not in the environment variables, we just open
    # the built docs in the browser. The GitHub token is present during CI
    # runs via GitHub Actions
    if [[ -z "${GITHUB_TOKEN}" ]]; then
        open docs/build/html/index.html
    fi
    
}

install_flit_dunamai () {
    pip install flit~=3.7
    pip install dunamai~=1.12
}

# Builds the docs from scratch
build_CI_docs () {
    install_flit_dunamai
    replace_version_in_init
    flit install --deps=production --extras=doc
    reverse_replace_version_in_init
    build_docs
}

build_CI_tests () {
    install_flit_dunamai
    replace_version_in_init
    flit install --deps=production --extras=test
    reverse_replace_version_in_init
    build_docs
}

flit_publish () {
    install_flit_dunamai
    replace_version_in_init
    flit publish --deps=production
    reverse_replace_version_in_init
}

# Iterate over all arguments in order
for var in "$@"
do
    echo "Found" "$var"
    if [ "$var" = "docs" ]; then
        build_CI_docs
    elif [ "$var" = "test" ]; then
        build_CI_tests
    # elif [ "$var" = "publish" ]; then
    #     flit_publish
    fi
done

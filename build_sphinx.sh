#!/bin/bash

sphinx-apidoc -f -o docs/source lightshow
make -C docs/ html
open docs/build/html/index.html

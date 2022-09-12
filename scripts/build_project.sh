#!/bin/bash

pip install flit~=3.7
bash scripts/install.sh
flit build

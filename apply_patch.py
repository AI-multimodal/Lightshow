#! /usr/bin/env python

import argparse
import ase
import os
import time
import shutil
import xanes_bench


def main():
    warn_text = "Warning: This script will overwrite the ase/io/espresso.py."
    parser = argparse.ArgumentParser(description=warn_text)
    parser.add_argument("-f", "--force", action='store_true', help="Enforce overwrite")
    options = parser.parse_args()

    if options.force:
        dest_fn = os.path.join(os.path.dirname(ase.__file__), "io", "espresso.py")
        src_fn = os.path.join(os.path.dirname(xanes_bench.__file__), "Xspectra", "espresso.py")
        ts = time.ctime().replace(' ', '_').replace(":", "_")
        shutil.move(dest_fn, f"{dest_fn}.{ts}")
        shutil.copy(src_fn, dest_fn)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()

#!/bin/sh
# Builds the two shared libraries used by pyslitdec_ox.py
set -e
cd "$(dirname "$0")"
cc -O2 -fPIC -Wall -shared -o libslitdec_ox.so slitdec_ox.c -lm
cc -O2 -fPIC -shared -o libslitdec_orig.so slitdec_orig.c -lm

#!/usr/bin/env bash

for i in ./GoGraph/classes/*.pyx; do mv -v "${i::-1}"{x,} ; done
rm -vf ./GoGraph/classes/*.cpython*\.so
rm -vf ./GoGraph/classes/*\.c

for i in ./utilities/*.pyx; do mv -v "${i::-1}"{x,} ; done
rm -vf ./utilities/*.cpython*\.so
rm -vf ./utilities/*\.c

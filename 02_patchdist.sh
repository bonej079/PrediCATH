#!/usr/bin/env bash

if [ "$1" = "classes" ]; then
    mv  *.cpython-37m-x86_64-linux-gnu.so ./GoGraph/classes/
fi

if [ "$1" = "utilities" ]; then
    mv  *.cpython-37m-x86_64-linux-gnu.so ./utilities
fi
#!/bin/bash

# # === Create build directory and build ===
mkdir -p build
cd build
cmake ..
make -j4



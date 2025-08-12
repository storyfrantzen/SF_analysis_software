#!/bin/bash

# Create build directory and build
mkdir -p build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/install ..
make -j4

# Install the built files (library + pcm) to the install prefix
make install

# Remove intermediate ROOT dictionary and shared library files in build root to avoid duplicates
# rm -f libBranchVarsDict.so libBranchVarsDict_rdict.pcm


export LD_LIBRARY_PATH=/work/clas12/storyf/SF_analysis_software/build/install/lib:$LD_LIBRARY_PATH
#!/bin/bash
# Build Sundials.

# Paths
SCRIPT_PATH=$(dirname "${BASH_SOURCE[@]}")
ROOT_PATH=$(cd "${SCRIPT_PATH}/.." && pwd)
SUNDIALS_BUILD_PATH="${ROOT_PATH}/thirdparty/sundials/build/"
SUNDIALS_INSTALL_PATH="${ROOT_PATH}/thirdparty/sundials/install/"

# Options
BUILD_TYPE="Release"

cmake=${CMAKE:-cmake}
make=${MAKE:-make}

# Make build directory and cd to it
mkdir -p "${SUNDIALS_BUILD_PATH}"
cd "${SUNDIALS_BUILD_PATH}" || exit

#   -DSUNDIALS_INDEX_SIZE=32 \
# Install sundials with 32 bit index integers as used in Eigen.
${cmake} \
  -DCMAKE_INSTALL_PREFIX="${SUNDIALS_INSTALL_PATH}" \
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_ARKODE=OFF \
  -DBUILD_CVODE=OFF \
  -DBUILD_IDA=OFF \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_STATIC_LIBS=ON \
  -DBUILD_NVECTOR_MANYVECTOR=OFF \
  -DBUILD_SUNNONLINSOL_PETSCSNES=OFF \
  -DEXAMPLES_INSTALL=OFF \
  -DEXAMPLES_ENABLE_C=OFF \
  ..

${make}
${make} install

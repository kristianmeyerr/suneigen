# Powershell script to install SUNDIALS on Windows Machines

# Paths
$ROOT_PATH="$PSScriptRoot" + "\.."
$SUNDIALS_BUILD_PATH="$ROOT_PATH" + "\thirdparty\sundials\build"
$SUNDIALS_INSTALL_PATH="$ROOT_PATH" + "\thirdparty\sundials\install"

# Options
$BUILD_TYPE="Release"

# Create build and install directories
New-Item -Path "$SUNDIALS_BUILD_PATH" -ItemType Directory
New-Item -Path "$SUNDIALS_INSTALL_PATH" -ItemType Directory

# cd to the build directory ..
cd $SUNDIALS_BUILD_PATH

# Build and install sundials
cmake `
    -DCMAKE_INSTALL_PREFIX="${SUNDIALS_INSTALL_PATH}" `
    -DCMAKE_BUILD_TYPE="${BUILD_TYPE}" `
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON `
    -DBUILD_ARKODE=OFF `
    -DBUILD_CVODE=OFF `
    -DBUILD_IDA=OFF `
    -DBUILD_SHARED_LIBS=OFF `
    -DBUILD_STATIC_LIBS=ON `
    -DBUILD_NVECTOR_MANYVECTOR=OFF `
    -DBUILD_SUNNONLINSOL_PETSCSNES=OFF `
    -DEXAMPLES_INSTALL=OFF `
    -DEXAMPLES_ENABLE_C=OFF `
    ..

cmake --build . --config $BUILD_TYPE --target all_build

cmake --build . --config $BUILD_TYPE --target install

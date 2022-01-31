# Powershell script to install SUNDIALS on Windows Machines

# Paths
$ROOT_PATH="$PSScriptRoot" + "\.."
$SUNDIALS_BUILD_PATH="$ROOT_PATH" + "\thirdparty\sundials\build"
$SUNDIALS_INSTALL_PATH="$ROOT_PATH" + "\thirdparty\sundials\install"

# Options
$BUILD_TYPE="Release"

# Create build and install directories
ni $SUNDIALS_BUILD_PATH
ni $SUNDIALS_INSTALL_PATH

# cd to the build directory ..
cd $SUNDIALS_BUILD_PATH

# Build and install sundials
cmake -DCMAKE_INSTALL_PREFIX="$SUNDIALS_INSTALL_PATH" -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DEXAMPLES_ENABLE_C=OFF ..

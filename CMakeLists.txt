cmake_minimum_required(VERSION 2.6)
PROJECT(linalg)

#
# Register own CMake extensions:
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -frtti -fexceptions -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -Wall -frtti -fexceptions")

#
# Search for libraries
#
FIND_PACKAGE(LAPACK REQUIRED)

FIND_PACKAGE(PythonLibs REQUIRED)

#
# Make sure libfluctuator can be found
#
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})


#
# Compiler flags:
#


#
# traverse into source tree:
#
add_subdirectory(src)
add_subdirectory(python)
add_subdirectory(test)



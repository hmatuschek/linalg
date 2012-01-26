cmake_minimum_required(VERSION 2.6)
PROJECT(linalg)

#
# Register own CMake extensions:
#
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)


#
# MacOS X specific stuff:
#
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
SET(CMAKE_OSX_DEPLOYMENT_TARGET 10.5)
#SET(CMAKE_OSX_ARCHITECTURES "i386")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

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



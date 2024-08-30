# Main CMakeLists file for OpenCAEPoroX
#
# Sample usages:
#   cmake .                             // build in default configuration
#   cmake -DCMAKE_BUILD_TYPE=Release .  // build in Release configuration
#   cmake -DCMAKE_BUILD_TYPE=Debug .    // build in Debug configuration
#   cmake -DCMAKE_CXX_COMPILER=g++ .    // build with GNU C++ compiler
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON . // build with verbose on
#   cmake -DUSE_FASPCPR=ON .            // build with FASPCPR support
#   cmake -DUSE_FASP4BLKOIL=ON .        // build with FASP4BLKOIL support
#   cmake -DUSE_FASP4CUDA=ON .          // build with FASP4CUDA support
#   cmake -DUSE_UMFPACK=ON .            // build with UMFPACK support
#   cmake -DUSE_OPENMP=ON .             // build with OpenMP support

###############################################################################
## General environment setting
###############################################################################

# Minimum cmake version needed
cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0076 NEW)

# Load user and defaults configure file
set(USER_CONFIG "${CMAKE_CURRENT_SOURCE_DIR}/config/user.cmake" CACHE PATH "Path to optional user configuration file.")
include("${USER_CONFIG}" OPTIONAL)
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

# Name of the project: should appear after find_program
project(OpenCAEPoroX C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Set module lookup path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config/cmake/modules)

# Add math and stdc++ libraries if not included by default
if (UNIX OR LINUX)
    set(ADD_STDLIBS m stdc++)
endif()

# Do not skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
 
# When building, don't use the install RPATH already
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
 
# Add the automatically determined parts of the RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(OCP_ENABLE_TESTING)
    enable_testing()
endif()


# LIBNAME is used in the external dependencies, to make it easier if we 
# decide to change the name in the future
set(LIBNAME OpenCAEPoroX CACHE INTERNAL "The library name")
add_library(OpenCAEPoroX STATIC)

###############################################################################
## External dependencies
###############################################################################

add_subdirectory(external)

#find_package(MPI)
#target_include_directories(${PROJECT_NAME} PUBLIC "${MPI_CXX_INCLUDE_PATH}")
#target_link_libraries(${LIBNAME} PUBLIC MPI::MPI_CXX)

###############################################################################
## Set configuration types and default compiler flags
###############################################################################

set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Configs" FORCE)
if (MSVC)
    # Compiler flags for Windows MSVC
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /DDEBUG /W4")
else()
    # Compiler flags for Linux or Mac OS X
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 ")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -Wall -g")
endif()


if (OCP_USE_MPI)
    find_package(MPI REQUIRED)
    set(MPI_CXX_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH})
endif()

if (OCP_USE_LAPACK)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

if (OCP_USE_METIS)
    find_package(METIS REQUIRED)
endif()

if (OCP_USE_PARMETIS)
    find_package(PARMETIS REQUIRED)
endif()

if (OCP_USE_FASP)
    find_package(FASP REQUIRED)
endif()

if (OCP_USE_FASP4BLKOIL)
    find_package(FASP4BLKOIL REQUIRED)
endif()


# Find third-party libraries
set(OCP_TPLS OPENMP LAPACK BLAS METIS PARMETIS PETSC PETSCSOLVER FASP FASP4BLKOIL MPI_CXX)
#
set(TPL_LIBRARIES "")
set(TPL_INCLUDE_DIRS "")
foreach(TPL IN LISTS OCP_TPLS)
    if (${TPL}_FOUND)
        message(STATUS "OCP: using package ${TPL}")
        list(APPEND TPL_LIBRARIES ${${TPL}_LIBRARIES})
        list(APPEND TPL_INCLUDE_DIRS ${${TPL}_INCLUDE_DIRS})
    endif()
endforeach(TPL)
#
list(REMOVE_DUPLICATES TPL_LIBRARIES)
list(REMOVE_DUPLICATES TPL_INCLUDE_DIRS)
#
include_directories(${TPL_INCLUDE_DIRS})


###############################################################################
## Project specific parameters
###############################################################################

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(main)

# Link third-party libraries
target_link_libraries(${LIBNAME} PUBLIC ${TPL_LIBRARIES})

add_subdirectory(config)

#
configure_file(
        "${PROJECT_SOURCE_DIR}/config/cmake/tpl_config.hpp.in"
        "${PROJECT_BINARY_DIR}/config/tpl_config.hpp")
#
file(WRITE "${PROJECT_BINARY_DIR}/config.hpp"
"// Auto-generated file.
#define OCP_CONFIG_FILE \"${PROJECT_BINARY_DIR}/config/tpl_config.hpp\"
")

target_compile_definitions(OpenCAEPoroX PUBLIC "OCP_CONFIG_FILE=\"${PROJECT_BINARY_DIR}/config/tpl_config.hpp\"")


# set default build type: Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "the build type" FORCE)
endif()
message(STATUS "Set build type to ${CMAKE_BUILD_TYPE} configuration")


if(WIN32)
	MESSAGE(STATUS "Windows Environments Now")
	find_package(MPI REQUIRED)
	target_link_libraries(${LIBNAME} PUBLIC MPI::MPI_CXX)
endif()


###############################################################################
## Installtion targets for lib and executable files
###############################################################################

# Set install location for OpenCAEPoroX library
set(CMAKE_INSTALL_LIBDIR "${PROJECT_SOURCE_DIR}/lib/")
install(TARGETS OpenCAEPoroX
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# CMake uninstall target
if(NOT TARGET uninstall)
    configure_file("${CMAKE_MODULE_PATH}/cmake_uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                   IMMEDIATE @ONLY)
    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

# Once done this will define
#  FASPXX_FOUND        - System has FASPXX
#  FASPXX_DIR          - The FASPXX directory
#  FASPXX_INCLUDE_DIRS - The FASPXX include directories
#  FASPXX_LIBRARIES    - The libraries needed to use FASPXX
#
#  Li Zhao
#  04/27/2024

# message(STATUS "Looking for FASPXX")

set(FASPXX_DIR "${FASPXX_DIR}")

# Check for header file
find_path(FASPXX_INCLUDE_DIRS faspxx.h
   HINTS ${FASPXX_DIR}/include $ENV{FASPXX_DIR}/include 
   DOC "Directory where the FASPXX header is located")
mark_as_advanced(FASPXX_INCLUDE_DIRS)

# Check for FASPXX library
find_library(FASPXX_LIBRARIES FASPXX
    HINTS ${FASPXX_DIR}/lib $ENV{FASPXX_DIR}/lib 
    DOC "The FASPXX library")
mark_as_advanced(FASPXX_LIBRARIES)

set(FASPXX_LIBRARIES ${FASPXX_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASPXX
    "FASPXX could not be found. Check FASPXX_DIR."
    FASPXX_INCLUDE_DIRS FASPXX_LIBRARIES)

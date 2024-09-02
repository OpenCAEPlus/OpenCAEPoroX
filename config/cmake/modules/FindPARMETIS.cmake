# - Try to find PARMETIS
# 
#  OUTPUT:
#  PARMETIS_FOUND        - system has PARMETIS
#  PARMETIS_INCLUDE_DIRS - include directories for PARMETIS
#  PARMETIS_LIBRARIES    - libraries for PARMETIS
#

message(STATUS "Checking for package 'PARMETIS'")

# Check for header file
find_path(PARMETIS_INCLUDE_DIRS parmetis.h
  HINTS ${PARMETIS_DIR} ${PARMETIS_DIR}/include ${PARMETIS_DIR}/PARMETIS/include $ENV{PARMETIS_DIR}/include $ENV{PARMETIS_DIR}/PARMETIS/include ${PARMETIS_DIR}/lib ${PARMETIS_DIR}/PARMETIS/lib $ENV{PARMETIS_DIR}/lib $ENV{PARMETIS_DIR}/PARMETIS/lib
  PATH_SUFFIXES suitesparse ufsparse
  DOC "Directory where the PARMETIS header is located"
  )
mark_as_advanced(PARMETIS_INCLUDE_DIRS)

# Check for GKlib library, should be installed at the same location as PARMETIS
# find_library(GKLIB_LIB GKlib
#   HINTS ${PARMETIS_DIR} ${PARMETIS_DIR}/PARMETIS ${PARMETIS_BUILD_DIR}/lib ${PARMETIS_DIR}/PARMETIS/lib $ENV{PARMETIS_DIR}/lib $ENV{PARMETIS_DIR}/PARMETIS/lib
#   DOC "The GKlib library"
#   )

# Check for PARMETIS library
find_library(PARMETIS_LIB parmetis
  HINTS ${PARMETIS_DIR} ${PARMETIS_DIR}/PARMETIS ${PARMETIS_DIR}/lib ${PARMETIS_DIR}/PARMETIS/lib $ENV{PARMETIS_BUILD_DIR}/libparmetis $ENV{PARMETIS_DIR}/PARMETIS/lib
  DOC "The PARMETIS library"
  )

set(PARMETIS_LIBRARIES "${GKLIB_LIB}" "${PARMETIS_LIB}") 
mark_as_advanced(PARMETIS_LIBRARIES)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS
  "PARMETIS could not be found. Be sure to set PARMETIS_DIR."
  PARMETIS_LIBRARIES PARMETIS_INCLUDE_DIRS)

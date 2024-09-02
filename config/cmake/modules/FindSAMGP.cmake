# ##############################################################################
# For SAMGSOLVER
# ##############################################################################


set(SAMG_DIR "${SAMG_DIR}")

# Check for header file
find_path(SAMGP_INCLUDE_DIRS samg.h
   HINTS ${SAMG_DIR}/include $ENV{SAMG_DIR}/include ${PROJECT_SOURCE_DIR}/SAMG/include
   DOC "Directory where the SAMG header is located")
mark_as_advanced(SAMGP_INCLUDE_DIRS)

# Check for SAMG library
find_library(SAMGP_LIBRARIES samgp
    HINTS ${SAMG_DIR}/lib $ENV{SAMG_DIR}/lib ${PROJECT_SOURCE_DIR}/SAMG/lib
    DOC "The SAMG library")
mark_as_advanced(SAMGP_LIBRARIES)

# Collect libraries
set(SAMGP_LIBRARIES ${SAMGP_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SAMGP
    "SAMGP could not be found. Check SAMG_DIR."
    SAMGP_LIBRARIES SAMGP_INCLUDE_DIRS)

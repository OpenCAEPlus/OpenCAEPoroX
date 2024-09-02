# ##############################################################################
# For AS (adaptive_solver)
# ##############################################################################


set(AS_DIR "${AS_DIR}")

# Check for header file
find_path(AS_INCLUDE_DIRS as.h
   HINTS ${AS_DIR}/include $ENV{AS_DIR}/include ${PROJECT_SOURCE_DIR}/AS/include
   DOC "Directory where the AS header is located")
mark_as_advanced(AS_INCLUDE_DIRS)

# Check for AS library
find_library(AS_LIBRARIES adaptive_solver
    HINTS ${AS_DIR}/lib $ENV{AS_DIR}/lib ${PROJECT_SOURCE_DIR}/AS/lib
    DOC "The AS library")
mark_as_advanced(AS_LIBRARIES)

# Collect libraries
set(AS_LIBRARIES ${AS_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AS
    "AS could not be found. Check AS_DIR."
    AS_LIBRARIES AS_INCLUDE_DIRS)

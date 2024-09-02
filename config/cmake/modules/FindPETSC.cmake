# ##############################################################################
# For PETSCSOLVER
# ##############################################################################


set(PETSC_DIR "${PETSC_DIR}")

# Check for header file
find_path(PETSC_INCLUDE_DIRS petsc.h
   HINTS ${PETSC_DIR}/include $ENV{PETSC_DIR}/include ${PROJECT_SOURCE_DIR}/PETSC/include
   DOC "Directory where the PETSC header is located")
mark_as_advanced(PETSC_INCLUDE_DIRS)

# Check for PETSC library
find_library(PETSC_LIBRARIES petsc
    HINTS ${PETSC_DIR}/lib $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib ${PROJECT_SOURCE_DIR}/PETSC/lib
    DOC "The PETSC library")
mark_as_advanced(PETSC_LIBRARIES)

# Collect libraries
set(PETSC_LIBRARIES ${PETSC_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSC
    "PETSC could not be found. Check PETSC_DIR."
    PETSC_LIBRARIES PETSC_INCLUDE_DIRS)

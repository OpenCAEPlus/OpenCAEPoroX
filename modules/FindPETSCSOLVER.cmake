# ##############################################################################
# For PETSCSOLVERSOLVER
# ##############################################################################


set(PETSCSOLVER_DIR "${PETSCSOLVER_DIR}")

# Check for header file
find_path(PETSCSOLVER_INCLUDE_DIRS PETScBSolverPS.h
   HINTS ${PETSCSOLVER_DIR}/include $ENV{PETSCSOLVER_DIR}/include ${PROJECT_SOURCE_DIR}/PETSCSOLVER/include
   DOC "Directory where the PETSCSOLVER header is located")
mark_as_advanced(PETSCSOLVER_INCLUDE_DIRS)

# Check for PETSCSOLVER library
find_library(PETSCSOLVER_LIBRARIES petsc_solver
    HINTS ${PETSCSOLVER_DIR}/lib $ENV{PETSCSOLVER_DIR}/lib ${PROJECT_SOURCE_DIR}/PETSCSOLVER/lib
    DOC "The PETSCSOLVER library")
mark_as_advanced(PETSCSOLVER_LIBRARIES)

# Collect libraries
set(PETSCSOLVER_LIBRARIES ${PETSCSOLVER_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSCSOLVER
    "PETSCSOLVER could not be found. Check PETSCSOLVER_DIR."
    PETSCSOLVER_LIBRARIES PETSCSOLVER_INCLUDE_DIRS)

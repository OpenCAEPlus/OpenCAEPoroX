# ##############################################################################
# For PETSCSOLVER
# ##############################################################################
option(USE_PETSCSOLVER "Use PETSCSOLVER" OFF)

if(USE_PETSCSOLVER)

  find_package(PETSCSOLVER)

  if (PETSCSOLVER_FOUND)
    message(STATUS "INFO: PETSCSOLVER found")
    add_library(petsc_solver STATIC IMPORTED GLOBAL)
    set_property(
    TARGET petsc_solver 
    APPEND
    PROPERTY IMPORTED_LOCATION ${PETSCSOLVER_LIBRARIES})
    set_property(
    TARGET petsc_solver 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSCSOLVER_INCLUDE_DIRS})
    set_property(
      TARGET petsc_solver 
      APPEND 
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_PETSCSOLVER=1")
    target_link_libraries(${LIBNAME} PUBLIC petsc_solver)
  else (PETSCSOLVER_FOUND)
    message(STATUS "INFO: PETSCSOLVER was requested but not found!")
  endif(PETSCSOLVER_FOUND)

endif(USE_PETSCSOLVER)



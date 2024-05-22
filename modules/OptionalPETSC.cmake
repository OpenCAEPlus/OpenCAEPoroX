# ##############################################################################
# For PETSC
# ##############################################################################
option(USE_PETSC "Use PETSC" OFF)

if(USE_PETSC)

  find_package(PETSC)

  if (PETSC_FOUND)
    message(STATUS "INFO: PETSC found")
    add_library(petsc SHARED IMPORTED GLOBAL)
    set_property(
    TARGET petsc 
    APPEND
    PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
    set_property(
    TARGET petsc 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
    # set_property(
    #   TARGET petsc 
    #   APPEND 
    #   PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_PETSC=1")
    target_link_libraries(${LIBNAME} PUBLIC petsc)
  else (PETSC_FOUND)
    message(STATUS "INFO: PETSCSOLVER was requested but not found!")
  endif(PETSC_FOUND)

endif(USE_PETSC)



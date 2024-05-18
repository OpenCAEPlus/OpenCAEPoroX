# ##############################################################################
# For AS (adaptive_solver)
# ##############################################################################
option(USE_AS "Use AS" OFF)

if(USE_AS)
  # try to find AS
  find_package(AS)
  find_package(PETSC)

  if (AS_FOUND)
    message(STATUS "INFO: AS found")
    add_library(adaptive_solver STATIC IMPORTED GLOBAL)
    # add_library(petsc SHARED IMPORTED GLOBAL)
    set_property(
    TARGET adaptive_solver 
    APPEND
    PROPERTY IMPORTED_LOCATION ${AS_LIBRARIES})
    # set_property(
    # TARGET petsc 
    # APPEND
    # PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
    set_property(
    TARGET adaptive_solver 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${AS_INCLUDE_DIRS})
    # set_property(
    # TARGET petsc 
    # APPEND
    # PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
    set_property(
      TARGET adaptive_solver 
      APPEND 
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_PETSCSOLVER=1" "WITH_AS=1")
      # PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_AS=1")
    target_link_libraries(${LIBNAME} PUBLIC adaptive_solver)
    # target_link_libraries(${LIBNAME} PUBLIC petsc)
  else (AS_FOUND)
    message(STATUS "INFO: AS was requested but not found!")
  endif(AS_FOUND)

endif(USE_AS)



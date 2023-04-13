# ##############################################################################
# For PARMETISSOLVER
# ##############################################################################
option(USE_PARMETIS "Use PARMETIS" OFF)

if(USE_PARMETIS)

  # set the path to find specific modules
  set(PARMETIS_DIR "${PARMETIS_DIR}")

  find_package(PARMETIS)

  if(PARMETIS_FOUND)
    # include_directories(${PARMETIS_INCLUDE_DIRS}) add_definitions(-D__SOLVER_PARMETIS__)
    add_library(parmetis STATIC IMPORTED GLOBAL)
    set_property(
      TARGET parmetis 
      PROPERTY IMPORTED_LOCATION ${PARMETIS_LIBRARIES})
    set_property(
      TARGET parmetis 
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PARMETIS_INCLUDE_DIRS})
      target_link_libraries(${LIBNAME} PUBLIC parmetis)
  else(PARMETIS_FOUND)
    message(STATUS "INFO: PARMETIS was requested but not found!")
    message(STATUS "INFO: Going to try download and install from git repo")
    # FetchContent_Declare(
    #   PARMETIS GIT_REPOSITORY https://github.com/PARMETISDevTeam/PARMETISsolver.git)
    # FetchContent_MakeAvailable(PARMETIS)
  endif(PARMETIS_FOUND)

endif(USE_PARMETIS)



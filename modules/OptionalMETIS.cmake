# ##############################################################################
# For METISSOLVER
# ##############################################################################
option(USE_METIS "Use METIS" OFF)

if(USE_METIS)

  # set the path to find specific modules
  set(METIS_DIR "${METIS_DIR}")

  find_package(METIS)

  if(METIS_FOUND)
    # include_directories(${METIS_INCLUDE_DIRS}) add_definitions(-D__SOLVER_METIS__)
    add_library(metis STATIC IMPORTED GLOBAL)
    set_property(
      TARGET metis 
      PROPERTY IMPORTED_LOCATION ${METIS_LIBRARIES})
    set_property(
      TARGET metis 
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIRS})
      target_link_libraries(${LIBNAME} PUBLIC metis)
  else(METIS_FOUND)
    message(STATUS "INFO: METIS was requested but not found!")
    message(STATUS "INFO: Going to try download and install from git repo")
    # FetchContent_Declare(
    #   METIS GIT_REPOSITORY https://github.com/METISDevTeam/METISsolver.git)
    # FetchContent_MakeAvailable(METIS)
  endif(METIS_FOUND)

endif(USE_METIS)



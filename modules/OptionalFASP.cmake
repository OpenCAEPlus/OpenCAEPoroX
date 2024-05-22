# ##############################################################################
# For FASPSOLVER
# ##############################################################################
option(USE_FASP "Use FASP" OFF)

if(USE_FASP)

  # set the path to find specific modules
  set(FASP_DIR "${FASP_DIR}")

  find_package(FASP)

  if(FASP_FOUND)
    # include_directories(${FASP_INCLUDE_DIRS}) add_definitions(-D__SOLVER_FASP__)
    add_library(fasp STATIC IMPORTED GLOBAL)
    set_property(
      TARGET fasp 
      PROPERTY IMPORTED_LOCATION ${FASP_LIBRARIES})
    set_property(
        TARGET fasp
        APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_FASP=1")      
    set_property(
      TARGET fasp 
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FASP_INCLUDE_DIRS})
      target_link_libraries(${LIBNAME} PUBLIC fasp)
  else(FASP_FOUND)
    message(STATUS "INFO: FASP was requested but not found!")
    message(STATUS "INFO: Going to try download and install from git repo")
    FetchContent_Declare(
      fasp GIT_REPOSITORY https://github.com/FaspDevTeam/faspsolver.git)
    FetchContent_MakeAvailable(fasp)
  endif(FASP_FOUND)

endif(USE_FASP)




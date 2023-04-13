# ##############################################################################
# For SAMGSOLVER
# ##############################################################################
option(USE_SAMGP "Use SAMGP" OFF)

if(USE_SAMGP)

  # try to find SAMG
  find_package(SAMGP)

  if (SAMGP_FOUND)
    message(STATUS "INFO: SAMGP found")
    add_library(samgp SHARED IMPORTED GLOBAL)
    set_property(
    TARGET samgp 
    APPEND
    PROPERTY IMPORTED_LOCATION ${SAMGP_LIBRARIES})
    set_property(
    TARGET samgp 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SAMGP_INCLUDE_DIRS})
    set_property(
      TARGET samgp 
      APPEND 
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_SAMG=1" "SAMG_LCASE_USCORE=1")
    target_link_libraries(${LIBNAME} PUBLIC samgp)
  else (SAMGP_FOUND)
    message(STATUS "INFO: SAMGP was requested but not found!")
  endif(SAMGP_FOUND)

endif(USE_SAMGP)



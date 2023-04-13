# ##############################################################################
# For SAMGSOLVER
# ##############################################################################
option(USE_SAMGP "Use SAMGP" OFF)

if(USE_SAMGP)

  # try to find SAMG
  find_package(SAMGPINTERFACE)

  if (SAMGPINTERFACE_FOUND)
    message(STATUS "INFO: SAMGPINTERFACE found")
    add_library(samgp_interface SHARED IMPORTED GLOBAL)
    set_property(
    TARGET samgp_interface 
    APPEND
    PROPERTY IMPORTED_LOCATION ${SAMGPINTERFACE_LIBRARIES})
    set_property(
    TARGET samgp_interface 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SAMGPINTERFACE_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC samgp_interface)
  else (SAMGPINTERFACE_FOUND)
    message(STATUS "INFO: SAMGP was requested but not found!")
  endif(SAMGPINTERFACE_FOUND)

endif(USE_SAMGP)



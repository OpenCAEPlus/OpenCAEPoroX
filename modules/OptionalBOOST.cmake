# ##############################################################################
# For Graph in BOOST library
# ##############################################################################
option(USE_BOOST "Use BOOST" OFF)

if(USE_BOOST)

  # set the path to find specific modules
  set(BOOST_DIR "${BOOST_DIR}")

  find_package(BOOST)

  if(BOOST_FOUND)
    add_library(boost_graph STATIC IMPORTED)
    add_library(boost_serialization STATIC IMPORTED)
    add_library(boost_filesystem STATIC IMPORTED)

    set_property(TARGET boost_graph PROPERTY IMPORTED_LOCATION ${BOOST_LIBRARIES})
    set_property(TARGET boost_serialization PROPERTY IMPORTED_LOCATION ${BOOST_SERIALIZATION_LIBRARY})
    set_property(TARGET boost_filesystem PROPERTY IMPORTED_LOCATION ${BOOST_FILESYSTEM_LIBRARY})
    # set_property(TARGET boost_graph APPEND PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_BGL=1" "WITH_PBGL=1")      
   
    set_property(TARGET boost_graph PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BOOST_INCLUDE_DIRS})
    set_property(TARGET boost_serialization PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BOOST_INCLUDE_DIRS})
    set_property(TARGET boost_filesystem PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${BOOST_INCLUDE_DIRS})

    set_property(TARGET ${LIBNAME} APPEND PROPERTY COMPILE_DEFINITIONS "WITH_BGL=1") 
    target_link_libraries(${LIBNAME} PUBLIC boost_graph boost_serialization boost_filesystem)

  else(BOOST_FOUND)
    message(
      WARNING
        "WARNING: BOOST was requested but not found! Continue without it."
    )
  endif(BOOST_FOUND)

endif(USE_BOOST)






# Once done this will define
#  BOOST_FOUND        - System has BOOST
#  BOOST_DIR          - The BOOST directory
#  BOOST_INCLUDE_DIRS - The BOOST include directories
#  BOOST_LIBRARIES    - The libraries needed to use BOOST
#
#  Li Zhao
#  05/07/2024

set(BOOST_DIR "${BOOST_DIR}")

# Check for header file
find_path(BOOST_INCLUDE_DIRS boost/graph/adjacency_list.hpp #boost/graph/connected_components.hpp
   HINTS ${BOOST_DIR}/include $ENV{BOOST_DIR} $ENV{BOOST_DIR}/include
   DOC "Directory where the BOOST header is located")
mark_as_advanced(BOOST_INCLUDE_DIRS)

# Check for BOOST library
find_library(BOOST_LIBRARIES boost_graph
    HINTS ${BOOST_DIR}/lib $ENV{BOOST_DIR}/lib $ENV{BOOST_DIR}/stage/lib
    DOC "The BOOST library")


set(_BOOST_LIBRARY_SEARCH_DIRS
  ${BOOST_DIR}
  $ENV{BOOST_DIR}
)
set(_BOOST_LIBRARY_DIR_SUFFIXES "lib")

## boost_serialization
set(_SER_LIBRARY boost_serialization)
find_library(SER_LIBRARY
    NAMES ${_SER_LIBRARY}
    PATHS ${_BOOST_LIBRARY_SEARCH_DIRS}
    PATH_SUFFIXES ${_BOOST_LIBRARY_DIR_SUFFIXES}
    DOC "Path to boost_filesystem library"
)
if (NOT SER_LIBRARY)
    message(STATUS "boost_serialization could not be found. Check BOOST_DIR/lib")
else()
    # list(APPEND BOOST_LIBRARIES ${SER_LIBRARY})
    set(BOOST_SERIALIZATION_LIBRARY ${SER_LIBRARY}) 
endif()


## boost_filesystem
set(_FS_LIBRARY boost_filesystem)
find_library(FS_LIBRARY
    NAMES ${_FS_LIBRARY}
    PATHS ${_BOOST_LIBRARY_SEARCH_DIRS}
    PATH_SUFFIXES ${_BOOST_LIBRARY_DIR_SUFFIXES}
    DOC "Path to boost_filesystem library"
)
if (NOT FS_LIBRARY)
    message(STATUS "boost_filesystem could not be found. Check BOOST_DIR/lib")
else()
    # list(APPEND BOOST_LIBRARIES ${FS_LIBRARY})
    set(BOOST_FILESYSTEM_LIBRARY ${FS_LIBRARY})  
endif()



mark_as_advanced(BOOST_LIBRARIES)
mark_as_advanced(BOOST_SERIALIZATION_LIBRARY)
mark_as_advanced(BOOST_FILESYSTEM_LIBRARY)


# Collect libraries
set(BOOST_LIBRARIES ${BOOST_LIBRARIES})
# set(BOOST_SERIALIZATION_LIBRARY ${BOOST_SERIALIZATION_LIBRARY})
# set(BOOST_FILESYSTEM_LIBRARY ${BOOST_FILESYSTEM_LIBRARY})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BOOST
    "BOOST could not be found. Check BOOST_DIR."
    BOOST_LIBRARIES BOOST_INCLUDE_DIRS)

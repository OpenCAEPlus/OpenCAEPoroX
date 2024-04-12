# ##############################################################################
# For PETSCSOLVERSOLVER
# ##############################################################################
option(USE_PETSCSOLVER "Use PETSCSOLVER" ON)

if(USE_PETSCSOLVER)

  # # try to find PETSCSOLVER_origin

  # include_directories($ENV{PETSCSOLVER_origin_DIR}/include/)

  # # Check for header file
  # find_path(PETSCSOLVER_origin_INCLUDE_DIRS petsc.h
  #   HINTS ${PETSCSOLVER_DIR}/include $ENV{PETSCSOLVER_origin_DIR}/include ${PROJECT_SOURCE_DIR}/petsc/include
  #   DOC "Directory where the PETSCSOLVER header is located")

  # # Check for petsc library
  # find_library(PETSCSOLVER_origin_LIBRARIES petsc
  #     HINTS $ENV{PETSCSOLVER_origin_DIR}/arch-linux2-c-debug/lib $ENV{PETSCSOLVER_DIR}/lib ${PROJECT_SOURCE_DIR}/PETSCSOLVER/lib
  #     DOC "The PETSCSOLVER library")

  # set(PETSCSOLVER_origin_LIBRARIES ${PETSCSOLVER_origin_LIBRARIES})
  # set(PETSCSOLVER_origin_INCLUDE_DIRS ${PETSCSOLVER_origin_INCLUDE_DIRS})

  # add_library(petsc SHARED)

  # target_link_libraries(${LIBNAME} PUBLIC ${PETSCSOLVER_origin_LIBRARIES})

  # try to find PETSCSOLVER
  find_package(PETSCSOLVER)
  find_package(PETSC)
  target_include_directories(${PROJECT_NAME} PUBLIC "${PETSC_INCLUDE_DIRS}")
  target_include_directories(${PROJECT_NAME} PUBLIC "${PETSC_DIR}/${PETSC_ARCH}/include")
  message("------> PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}, HAHAHAHAH ${PETSC_DIR}/${PETSC_ARCH}/include")

  if (PETSCSOLVER_FOUND)
    message(STATUS "INFO: PETSCSOLVER found")
    add_library(petsc_solver STATIC IMPORTED GLOBAL)
    add_library(petsc SHARED IMPORTED GLOBAL)
    set_property(
    TARGET petsc_solver 
    APPEND
    PROPERTY IMPORTED_LOCATION ${PETSCSOLVER_LIBRARIES})
    set_property(
    TARGET petsc 
    APPEND
    PROPERTY IMPORTED_LOCATION ${PETSC_LIBRARIES})
    set_property(
    TARGET petsc_solver 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSCSOLVER_INCLUDE_DIRS})
    set_property(
    TARGET petsc 
    APPEND
    PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PETSC_INCLUDE_DIRS})
    set_property(
      TARGET petsc_solver 
      APPEND 
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_PETSCSOLVER=1")
    target_link_libraries(${LIBNAME} PUBLIC petsc_solver)
    target_link_libraries(${LIBNAME} PUBLIC petsc)
  else (PETSCSOLVER_FOUND)
    message(STATUS "INFO: PETSCSOLVER was requested but not found!")
  endif(PETSCSOLVER_FOUND)

endif(USE_PETSCSOLVER)



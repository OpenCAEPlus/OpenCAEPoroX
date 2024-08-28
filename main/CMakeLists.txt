# Main executable target: testOpenCAEPoro
add_executable(testOpenCAEPoro)
target_sources(testOpenCAEPoro PRIVATE Main.cpp)
target_link_libraries(testOpenCAEPoro PUBLIC OpenCAEPoroX ${ADD_STDLIBS})
install(TARGETS testOpenCAEPoro DESTINATION ${PROJECT_SOURCE_DIR})

if(BUILD_TEST)
  include(CTest)

  add_test(NAME spe1a
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe1a/spe1a.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

  add_test(NAME spe3
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe3/spe3.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

  add_test(NAME spe4
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe4/spe4/spe4.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

  add_test(NAME spe4m
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe4/spe4m/spe4m.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

  add_test(NAME spe5
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe5/spe5.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

#  add_test(NAME spe6
#          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe6/spe6.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
#          )

  add_test(NAME spe9_FIM
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe9/spe9_FIM.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

  add_test(NAME spe9_IMPEC
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe9/spe9_IMPEC.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

  add_test(NAME spe10
          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe10/spe10.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
          )

#  add_test(NAME spe11
#          COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/data/spe11/spe11a/spe11a.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1
#          )

endif()

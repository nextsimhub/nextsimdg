function(add_nextdimdg_test)
  # Define the supported set of keywords
  set(prefix ARG)
  set(noValues ENABLE_OPENMP)
  set(singleValues TARGET)
  set(multiValues SOURCES)

  # Process the arguments passed in
  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix}
    "${noValues}"
    "${singleValues}"
    "${multiValues}"
    ${ARGN})

  add_executable(${ARG_TARGET} ${ARG_SOURCES})
  #target_link_libraries(${ARG_TARGET} PRIVATE ${PROJECT_NAME}::${NextsimDGLib})
  target_link_libraries(${ARG_TARGET} PRIVATE ${NextsimDGLib})  
  foreach(numProcs RANGE 1 2)
    add_test(NAME ${ARG_TARGET}_n${numProcs} COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${numProcs} ${MPIEXEC_PREFLAGS} ./${ARG_TARGET})
    if(ARG_ENABLE_OPENMP)
      set_tests_properties(${ARG_TARGET}_n${numProcs} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2")
    endif()
  endforeach()
endfunction()

# Define tests
if(ENABLE_OPENMP)
  add_nextdimdg_test(TARGET multithreading_square SOURCES multithreading_square.cpp ENABLE_OPENMP)
  add_nextdimdg_test(TARGET multithreading_rectangle1 SOURCES multithreading_rectangle1.cpp ENABLE_OPENMP)
  add_nextdimdg_test(TARGET multithreading_rectangle2 SOURCES multithreading_rectangle2.cpp ENABLE_OPENMP)
else()
  add_nextdimdg_test(TARGET multithreading_square SOURCES multithreading_square.cpp)
  add_nextdimdg_test(TARGET multithreading_rectangle1 SOURCES multithreading_rectangle1.cpp)
  add_nextdimdg_test(TARGET multithreading_rectangle2 SOURCES multithreading_rectangle2.cpp)
endif()
add_nextdimdg_test(TARGET example1 SOURCES example1.cpp ../../src/Timer.cpp)
add_nextdimdg_test(TARGET example2a SOURCES example2a.cpp)
add_nextdimdg_test(TARGET example2b SOURCES example2b.cpp)
add_nextdimdg_test(TARGET example2c SOURCES example2c.cpp)
add_nextdimdg_test(TARGET symlaplace SOURCES symlaplace.cpp)



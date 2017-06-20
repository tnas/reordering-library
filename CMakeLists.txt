cmake_minimum_required(VERSION 2.8)

project(reordering-library)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(PROFILER_ENABLED ON)

set(SCOREP_PROFILER "scorep")
set(COMPILER_FLAGS "-Wall -fopenmp -o3")
set(CXX_COMPILER_FLAGS "g++ ${COMPILER_FLAGS}")
set(C_COMPILER_FLAGS "gcc ${COMPILER_FLAGS}")

if(PROFILER_ENABLED)
	foreach(LANG C CXX)
		message(STATUS "Enabling scorep for ${LANG}")
		set(CMAKE_${LANG}_COMPILER "${SCOREP_PROFILER}")
		set(CMAKE_${LANG}_FLAGS "${${LANG}_COMPILER_FLAGS}")
	endforeach()
else()
	message(STATUS "Default compilers")
	foreach(LANG C CXX)
		set(CMAKE_${LANG}_FLAGS "${COMPILER_FLAGS}")
	endforeach()
endif()


add_executable(reordering-library program.c 
	./CommonFiles/util.c 
	./CommonFiles/matrix.c 
	./CommonFiles/linked_list.c 
	./CommonFiles/graph.c 
	./CommonFiles/util_parallel.c
	./CommonFiles/matrix_parallel.c
	./CommonFiles/graph_parallel.c 
	./CommonFiles/graph_hsl.c
	./Reorderings/spectral.c 
	./Reorderings/rcm.c 
	./Reorderings/rcm_hsl.c
	./Reorderings/rcm_parallel.c
	./Reorderings/rcm_boost.cpp
	./Reorderings/sloan.c
	./Reorderings/sloan_hsl.c
	./Reorderings/sloan_parallel.c
	./Reorderings/sloan_boost.cpp
	./UnitTests/test_suite_reordering.c
	./UnitTests/test_suite_prefixsum.c
	./UnitTests/test_suite_matrix.c
	./UnitTests/test_graph_parallel.c
	./UnitTests/test_util.c
	./UnitTests/test_linked_list.c)
	
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(reordering-library ${Boost_LIBRARIES})
target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/hsl_mc73d.o)
target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/libhsl_mc73.a)
target_link_libraries(reordering-library m gfortran blas)

install(TARGETS reordering-library RUNTIME DESTINATION bin)

add_subdirectory(CommonFiles)
add_subdirectory(Matrices)
add_subdirectory(Reorderings)
add_subdirectory(UnitTests)
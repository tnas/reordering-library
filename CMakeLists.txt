cmake_minimum_required(VERSION 2.8)

# project(reordering-library C CXX Fortran)
project(reordering-library)

# set(CMAKE_VERBOSE_MAKEFILE TRUE)

add_executable(reordering-library  program.c 
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
	./Reorderings/rcm_octave.cpp
	./Reorderings/sloan.c
	./Reorderings/sloan_hsl.c
	./Reorderings/sloan_parallel.c
	./Reorderings/sloan_boost.cpp
	./UnitTests/test_suite_reordering.c
	./UnitTests/test_suite_prefixsum.c
	./UnitTests/test_suite_matrix.c
	./UnitTests/test_graph_parallel.c
	./UnitTests/test_util.c)
	
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
# file(GLOB_RECURSE OCTAVE_INCLUDE LIST_DIRECTORIES true RELATIVE /usr/include/octave-4.0.2/octave octave.h)
# message("OCTAVE INCLUDE: " ${OCTAVE_INCLUDE})
# include_directories(/usr/include)
include_directories(/usr/include/octave-4.0.2/octave)
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# set(OCTAVE_INCLUDE /usr/include/octave-4.0.0/octave)
# set(MKOCTFILE "mkoctfile")
# set(OCTAVE_RCM_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/rcm_octave.cpp")

# add_custom_command(
# 	TARGET reordering-library
# 	PRE_BUILD
# 	COMMAND ${MKOCTFILE} -v ${OCTAVE_RCM_SOURCE}
# 	COMMENT "-- Compiling Octave RCM algorithm.")

# add_custom_command(
# 	OUTPUT rcm_octave.o
# 	COMMAND ${MKOCTFILE} -v ${OCTAVE_RCM_SOURCE}
# 	COMMENT "-- Compiling Octave RCM algorithm.")
# add_library(rcm_octave STATIC rcm_octave.o)
# set_source_files_properties(rcm_octave.o 
# 	PROPERTIES 
# 	EXTERNAL_OBJECT true
# 	GENERATED true)
# set_target_properties(rcm_octave
# 	PROPERTIES
# 	LINKER_LANGUAGE CXX)
# target_link_libraries(reordering-library rcm_octave)

target_link_libraries(reordering-library ${Boost_LIBRARIES})
target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/hsl_mc73d.o)
target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/libhsl_mc73.a)
target_link_libraries(reordering-library m gfortran blas octave)

install(TARGETS reordering-library RUNTIME DESTINATION bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -o3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -fopenmp -o3")

# set_property(TARGET reordering-library PROPERTY LINKER_LANGUAGE CXX)

add_subdirectory(CommonFiles)
add_subdirectory(Matrices)
add_subdirectory(Reorderings)
add_subdirectory(UnitTests)
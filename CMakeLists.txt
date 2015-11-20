cmake_minimum_required(VERSION 2.8)

project(reordering-library)

add_executable(reordering-library program.c ./CommonFiles/blas.c ./CommonFiles/matrix.c ./CommonFiles/linked_list.c ./CommonFiles/graph.c ./CommonFiles/graph_parallel.c ./Reorderings/rcm.c ./Reorderings/rcm_unordered.c ./Reorderings/sloan.c ./Reorderings/spectral.c)

target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/hsl_mc73d.o)
target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/libhsl_mc73.a)

install(TARGETS reordering-library RUNTIME DESTINATION bin)

set(CMAKE_C_FLAGS "-lm -Wall -lgfortran -lblas")
set(CMAKE_C_COMPILER "/usr/bin/gcc")
#set(CMAKE_VERBOSE_MAKEFILE ON)

add_subdirectory(CommonFiles)
add_subdirectory(Reorderings)
add_subdirectory(Matrices)
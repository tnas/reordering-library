cmake_minimum_required(VERSION 2.8)

project(reordering-library)

add_executable(reordering-library program.c ./CommonFiles/blas.c ./CommonFiles/matrix.c ./CommonFiles/linked_list.c ./CommonFiles/graph.c ./CommonFiles/graph_parallel.c ./Reorderings/rcm.c ./Reorderings/rcm_unordered.c ./Reorderings/sloan.c)

#target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/hsl_mc73d.o)
#target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/libhsl_mc73.a)

install(TARGETS reordering-library RUNTIME DESTINATION bin)

set(CMAKE_C_FLAGS "-lm -Wall -lgfortran -lblas -fopenmp")
set(CMAKE_C_COMPILER "/usr/bin/gcc")

add_subdirectory(CommonFiles)
add_subdirectory(Matrices)
add_subdirectory(Reorderings)
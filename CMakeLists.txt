cmake_minimum_required(VERSION 2.8)

project(reordering-library C Fortran)

add_executable(reordering-library program.c ./Reorderings/hsl_mc60/mc60d.f ./CommonFiles/util.c ./CommonFiles/blas.c ./CommonFiles/matrix.c ./CommonFiles/linked_list.c ./CommonFiles/graph.c ./CommonFiles/graph_parallel.c ./Reorderings/rcm.c Reorderings/rcm_parallel.c ./Reorderings/sloan.c ./UnitTests/rcm_parallel_test.c)

target_link_libraries(reordering-library m)
target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc60/libmc60.a)

#target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/hsl_mc73d.o)
#target_link_libraries(reordering-library ${CMAKE_CURRENT_SOURCE_DIR}/Reorderings/hsl_mc73/libhsl_mc73.a)

install(TARGETS reordering-library RUNTIME DESTINATION bin)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -lblas -fopenmp")
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set_property(TARGET reordering-library PROPERTY LINKER_LANGUAGE C)

add_subdirectory(CommonFiles)
add_subdirectory(Matrices)
add_subdirectory(Reorderings)
add_subdirectory(UnitTests)
add_subdirectory(Big-Matrices)
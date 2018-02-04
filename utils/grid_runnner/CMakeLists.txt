cmake_minimum_required(VERSION 3.5)
include(../../src/CMakeLists.txt)

project(run_param_grid)

set(BUILD_SHARED_LIBS FALSE)

set(SRC
    run_param_grid.c
    )

add_executable(run_param_grid ${SRC})
target_link_libraries(run_param_grid meraxes_lib)
install(TARGETS run_param_grid DESTINATION bin COMPONENT bin)

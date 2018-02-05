cmake_minimum_required(VERSION 3.5)
project(run_param_grid)
set(BUILD_SHARED_LIBS OFF)

find_package(Meraxes PATHS ${CMAKE_CURRENT_SOURCE_DIR}/../../share/meraxes)

set(SRC
    run_param_grid.c
    )

add_executable(run_param_grid ${SRC})
target_link_libraries(run_param_grid meraxes_lib)

# install(TARGETS run_param_grid DESTINATION bin COMPONENT bin)

set_property(TARGET run_param_grid PROPERTY C_STANDARD 99)

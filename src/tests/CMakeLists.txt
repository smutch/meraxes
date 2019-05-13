find_path(CRITERION_INCLUDE_DIR
    NAMES criterion/criterion.h
    PATHS ${CRITERION_ROOT_DIR}/include
    )

find_library(CRITERION_LIBRARY
    NAMES criterion
    PATHS ${CRITERION_ROOT_DIR}/lib
    )

# add_executable(test_init test_init.c)
# set_property(TARGET test_init PROPERTY C_STANDARD 99)
# target_include_directories(test_init PRIVATE ${CRITERION_INCLUDE_DIR} ${CMAKE_SOURCE_DIR})
# target_link_libraries(test_init PRIVATE ${CRITERION_LIBRARY} meraxes_lib)
# add_test(NAME test_init COMMAND test_init)

add_executable(test_parse_snaplist test_parse_snaplist.c)
set_property(TARGET test_parse_snaplist PROPERTY C_STANDARD 99)
target_include_directories(test_parse_snaplist PRIVATE ${CRITERION_INCLUDE_DIR} ${CMAKE_SOURCE_DIR})
target_link_libraries(test_parse_snaplist PRIVATE ${CRITERION_LIBRARY} meraxes_lib)
add_test(NAME test_parse_snaplist COMMAND test_parse_snaplist)

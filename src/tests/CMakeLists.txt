add_executable(test_init test_init.c)

target_link_libraries(test_init meraxes_lib)
target_link_libraries(test_init criterion)

add_test(NAME test_init COMMAND test_init)

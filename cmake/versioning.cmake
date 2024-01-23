find_package(Git REQUIRED)

execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    OUTPUT_VARIABLE GITREF
    OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${GIT_EXECUTABLE} --no-pager diff --no-color
    OUTPUT_VARIABLE GITDIFF
    OUTPUT_STRIP_TRAILING_WHITESPACE)

if (NOT GITDIFF STREQUAL "")
    string(REPLACE "\\\"" "\\\\\"" GITDIFF ${GITDIFF})
    string(REPLACE "\n" "\\n" GITDIFF ${GITDIFF})
endif()

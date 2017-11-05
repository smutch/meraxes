# project setup
cmake_minimum_required(VERSION 3.5)
project(meraxes)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/local.cmake OPTIONAL)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set_property(DIRECTORY PROPERTY EP_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/3rd_party")
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE CACHE BOOL "Export compile commands (for clang-tidy)" FORCE)
find_package(Git REQUIRED)

# Build type
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -Wall -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -Wall -march=native")
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type.")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/.. CACHE PATH "..." FORCE)
endif()

# ensure RPATH is used for builds
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Compiler options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")

# need to set extra options for clang
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()

# source files
file(GLOB SOURCES core/*.[ch] *.h physics/*.[ch])
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core/meraxes.c)

# output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# build options
set(N_HISTORY_SNAPS 5 CACHE VALUE "Number of snapshots for delayed SN feedback")
add_definitions(-DN_HISTORY_SNAPS=${N_HISTORY_SNAPS})
set(USE_JeMalloc TRUE CACHE BOOL "Build with jemalloc (may improve performance on some systems)")
set(USE_BUNDLED_JeMalloc FALSE CACHE BOOL "Download and build against bundled jemalloc")
set(USE_BUNDLED_FFTW FALSE CACHE BOOL "Download and build against bundled fftw")

# library
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared library version (required for Mhysa)")
add_library(meraxes_lib ${SOURCES})
set_target_properties(meraxes_lib PROPERTIES OUTPUT_NAME meraxes)
install(TARGETS meraxes_lib DESTINATION lib COMPONENT lib)
install(FILES meraxes.h ${CMAKE_BINARY_DIR}/git.h DESTINATION include COMPONENT lib)
add_custom_target(install.lib
    ${CMAKE_COMMAND} 
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE} 
    -DCOMPONENT=lib
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake) 
add_dependencies(install.lib meraxes_lib)

# executable
if(NOT BUILD_SHARED_LIBS)
    add_executable(meraxes ${CMAKE_CURRENT_SOURCE_DIR}/core/meraxes.c)
    target_link_libraries(meraxes meraxes_lib)
    install(TARGETS meraxes DESTINATION bin COMPONENT bin)
    add_custom_target(install.meraxes
        ${CMAKE_COMMAND} 
        -DBUILD_TYPE=${CMAKE_BUILD_TYPE} 
        -DCOMPONENT=bin
        -P ${CMAKE_BINARY_DIR}/cmake_install.cmake) 
    add_dependencies(install.meraxes meraxes)

    add_custom_target(setuprun
        COMMENT "Setting up run with minimal required files..."
        COMMAND test -n \"\${RUNDIR}\" || (echo "*** Must set RUNDIR to execute target 'setuprun'. ***" && exit 1)
        COMMAND echo "RUNDIR = \${RUNDIR}"
        COMMAND mkdir -p \${RUNDIR}
        COMMAND cp -n ${CMAKE_CURRENT_SOURCE_DIR}/../input/input.par \${RUNDIR}
        COMMAND cp -n ${CMAKE_CURRENT_SOURCE_DIR}/../input/snaplist.txt \${RUNDIR}
        COMMAND cp ${CMAKE_BINARY_DIR}/bin/meraxes \${RUNDIR}
        )
    add_dependencies(setuprun meraxes)
endif()

# Meraxes headers
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
    OUTPUT_VARIABLE GITREF
    OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${GIT_EXECUTABLE} --no-pager diff --no-color
    OUTPUT_VARIABLE GITDIFF
    OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE "\\n" "\\\\n" GITDIFF GITDIFF)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/git.h.in ${CMAKE_BINARY_DIR}/git.h ESCAPE_QUOTES @ONLY)

add_custom_target(cleandist
    COMMENT "Clean everything, including 3rd party libs."
    COMMAND ${CMAKE_MAKE_PROGRAM} clean
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning 3rd party libs..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party
    )

if(EXISTS ${CLANG_FORMAT_BINARY})
    add_custom_target(format
        COMMENT "Use clang-format to format the source code consistently."
        COMMAND ${CLANG_FORMAT_BINARY} -style=WebKit -i ${SOURCES}
        )
endif()

##############
# DEPENDENCIES
##############

# MATH
target_link_libraries(meraxes_lib m)

# JEMALLOC
if(USE_JeMalloc OR USE_BUNDLED_JeMalloc)
    if (NOT USE_BUNDLED_JeMalloc)
        find_package(JeMalloc REQUIRED)
        target_link_libraries(meraxes_lib ${JEMALLOC_LIBRARIES})
    else()
        include(cmake/BundleJeMalloc.cmake)
        bundle_jemalloc()
    endif()
endif()

# MPI
find_package(MPI REQUIRED)
target_include_directories(meraxes_lib PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(meraxes_lib ${MPI_C_LIBRARIES})

# MLOG
find_library(MLOG_LIBRARIES libmlog.a PATHS ${MLOG_ROOT})
if((NOT MLOG_LIBRARIES) OR (NOT EXISTS ${MLOG_LIBRARIES}))
    # message(STATUS "Unable to find mlog. Building it as part of Meraxes")
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mlog/mlog.c")
        message(STATUS "Checking out mlog submodule")
        execute_process(COMMAND git submodule update --init
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    set(MLOG_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/mlog")
    add_library(mlog STATIC "${MLOG_ROOT}/mlog.c" "${MLOG_ROOT}/mlog.h")
    set_property(TARGET mlog PROPERTY POSITION_INDEPENDENT_CODE ON)
    target_include_directories(meraxes_lib PUBLIC ${MLOG_ROOT})
    target_link_libraries(meraxes_lib mlog)
else()
    find_path(MLOG_INCLUDE_DIRS mlog.h PATHS ${MLOG_ROOT})
    target_include_directories(meraxes_lib PUBLIC ${MLOG_INCLUDE_DIRS})
    target_link_libraries(meraxes_lib ${MLOG_LIBRARIES})
endif()

# GSL
find_package(GSL REQUIRED)
target_include_directories(meraxes_lib PUBLIC ${GSL_INCLUDE_DIRS})
target_link_libraries(meraxes_lib ${GSL_LIBRARIES})

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)
if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Meraxes requires HDF5 with parallel support.")
endif()
target_include_directories(meraxes_lib PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(meraxes_lib ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
add_definitions(${HDF5_DEFINITIONS})

# THIS IS A TEMPORARY HACK TO OVERCOME A BUG IN SOME VERSIONS OF CMAKE
if(NOT (${HDF5_C_HL_LIBRARIES} MATCHES "libhdf5_hl"))
    list(GET HDF5_C_LIBRARIES 0 HDF5_LIB)
    get_filename_component(HDF5_LIBDIR ${HDF5_LIB} DIRECTORY)
    target_link_libraries(meraxes_lib "${HDF5_LIBDIR}/libhdf5_hl.so")
endif()

# FFTW
if(NOT USE_BUNDLED_FFTW)
    find_package(FFTW REQUIRED)
    target_include_directories(meraxes_lib PUBLIC ${FFTW_INCLUDE_DIRS})
    target_link_libraries(meraxes_lib ${FFTW_LIBRARIES})
else()
    include(cmake/BundleFFTW.cmake)
    bundle_fftw()
endif()

# OpenMP
option(USE_OpenMP "Use OpenMP" OFF)
if(USE_OpenMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Profiling
option(ENABLE_PROFILING "Enable profiling of executable with gperftools." OFF)
if(ENABLE_PROFILING)
    find_library(GPERF_LIBRARIES profiler PATHS "${GPERF_ROOT}/*")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

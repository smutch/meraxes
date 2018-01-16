# project setup
cmake_minimum_required(VERSION 3.5)
# Set default compilers (must come before project line)
set(CMAKE_C_COMPILER   "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(meraxes)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/local.cmake OPTIONAL)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set_property(DIRECTORY PROPERTY EP_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/3rd_party")
option(CMAKE_EXPORT_COMPILE_COMMANDS "Export compile commands (for clang-tidy)" ON)
find_package(Git REQUIRED)

# build options
set(N_HISTORY_SNAPS 5 CACHE VALUE "Number of snapshots for delayed SN feedback")
add_definitions(-DN_HISTORY_SNAPS=${N_HISTORY_SNAPS})
option(USE_BUNDLED_FFTW "Download and build against bundled fftw" OFF)
option(BUILD_TESTS "Build criterion tests" OFF)

# Build type
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -Wall -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -Wall -march=native")
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type.")

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/.. CACHE PATH "..." FORCE)
endif()

# ensure RPATH is used for builds
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Compiler options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
set(CMAKE_CXX_FLAGS "-std=c++11")

# The following line is needed to deal with this annoying problem:
#    https://developers.redhat.com/blog/2015/02/05/gcc5-and-the-c11-abi/
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

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

# library
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared library version (required for Mhysa)")

# CUDA & CUFFT (& add executable: meraxes)
#  Note that the CUDA version should be >= 7.5
#  because c++11 support is needed!
option(USE_CUDA "Use CUDA" OFF)
if(USE_CUDA)
    find_package(CUDA REQUIRED)
    add_definitions(-DUSE_CUDA)
    option(USE_CUFFT "Use cuFFT" OFF)
    if(USE_CUFFT)
        add_definitions(-DUSE_CUFFT)
    endif()

    # Collect CUDA source code
    file(GLOB SOURCES_CU core/*.cu *.hh)
    list(APPEND SOURCES ${SOURCES_CU})

    # Pass options to NVCC
    set(CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(CUDA_PROPAGATE_HOST_FLAGS ON)
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_30")

    ################################################################################
    # This ugly hack is needed to work around a bug in the fftw3.h header file
    #    see the following for a discussion: https://github.com/FFTW/fftw3/issues/18
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -D__INTEL_COMPILER")
    ################################################################################

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CUDA_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CUDA_EXE_LINKER_FLAGS}")

    # Add the Meraxes lib to the project
    cuda_add_library(meraxes_lib ${SOURCES})

    # CUFFT
    if(USE_CUFFT)
        cuda_add_cufft_to_target(meraxes_lib)
    endif()
else()
    add_library(meraxes_lib ${SOURCES})
endif()

set_target_properties(meraxes_lib PROPERTIES OUTPUT_NAME meraxes)
install(TARGETS meraxes_lib DESTINATION lib COMPONENT lib)
install(FILES ${CMAKE_BINARY_DIR}/meraxes.h ${CMAKE_BINARY_DIR}/git.h DESTINATION include COMPONENT lib)
add_custom_target(install.lib
    ${CMAKE_COMMAND} 
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE} 
    -DCOMPONENT=lib
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake) 
add_dependencies(install.lib meraxes_lib)

# executable
if(NOT BUILD_SHARED_LIBS)
    if(USE_CUDA)
        cuda_add_executable(meraxes ${CMAKE_CURRENT_SOURCE_DIR}/core/meraxes.c)
    else()
        add_executable(meraxes ${CMAKE_CURRENT_SOURCE_DIR}/core/meraxes.c)
    endif()
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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/meraxes.h.in ${CMAKE_BINARY_DIR}/meraxes.h ESCAPE_QUOTES @ONLY)
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

#######
# TESTS
#######

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

##############
# DEPENDENCIES
##############

# MATH
target_link_libraries(meraxes_lib m)

# MPI
find_package(MPI REQUIRED)
target_include_directories(meraxes_lib PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(meraxes_lib ${MPI_C_LIBRARIES})

# MLOG
find_library(MLOG_LIBRARIES libmlog.a PATHS ${MLOG_ROOT})
if((NOT MLOG_LIBRARIES) OR (NOT EXISTS ${MLOG_LIBRARIES}))
    message(STATUS "Unable to find mlog. Building it as part of Meraxes")
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

# Profiling
option(ENABLE_PROFILING "Enable profiling of executable with gperftools." OFF)
if(ENABLE_PROFILING)
    find_library(GPERF_LIBRARIES profiler PATHS "${GPERF_ROOT}/*")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

# project setup
cmake_minimum_required(VERSION 3.5)

# Set default compilers (must come before project line)
set(CMAKE_C_COMPILER   "gcc")
set(CMAKE_CXX_COMPILER "g++")

# Set project name
project(meraxes)
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# Git
find_package(Git REQUIRED)

# Uncomment for a verbose run
set(CMAKE_VERBOSE_MAKEFILE ON)

# Build type
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -DDEBUG -Wall -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -Wall -march=native")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type.")

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

# Code options
if(N_HISTORY_SNAPS)
    add_definitions(-DN_HISTORY_SNAPS=${N_HISTORY_SNAPS})
endif()

# Compiler options
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -std=gnu99")
set(CMAKE_CXX_FLAGS "-std=c++11")

# The following line is needed to deal with this annoying problem:
#    https://developers.redhat.com/blog/2015/02/05/gcc5-and-the-c11-abi/
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

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

# Source files
file(GLOB SOURCES_C core/*.[ch] *.h physics/*.[ch])

# MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
link_libraries(${MPI_C_LIBRARIES})

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
    add_custom_command(TARGET mlog PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${MLOG_ROOT}/mlog.h $<TARGET_FILE_DIR:mlog>)
    link_libraries(mlog)
else()
    find_path(MLOG_INCLUDE_DIRS mlog.h PATHS ${MLOG_ROOT})
    include_directories(${MLOG_INCLUDE_DIRS})
    link_libraries(${MLOG_LIBRARIES})
endif()

# GSL
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
link_libraries(${GSL_LIBRARIES})

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)
if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Meraxes requires HDF5 with parallel support.")
endif()
include_directories(${HDF5_INCLUDE_DIRS})
link_libraries(${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
add_definitions(${HDF5_DEFINITIONS})

# THIS IS A TEMPORARY HACK TO OVERCOME A BUG IN SOME VERSIONS OF CMAKE
if(NOT (${HDF5_C_HL_LIBRARIES} MATCHES "libhdf5_hl"))
    list(GET HDF5_C_LIBRARIES 0 HDF5_LIB)
    get_filename_component(HDF5_LIBDIR ${HDF5_LIB} DIRECTORY)
    link_libraries("${HDF5_LIBDIR}/libhdf5_hl.so")
endif()

# FFTW
find_path(FFTW_INCLUDE_DIRS fftw3.h PATHS "${FFTW_ROOT}/*" NO_DEFAULT_PATH)
find_path(FFTW_INCLUDE_DIRS fftw3.h)
find_library(FFTW_LIBRARIES fftw3f PATHS "${FFTW_ROOT}/*" NO_DEFAULT_PATH)
find_library(FFTW_LIBRARIES fftw3f)
find_library(FFTW_MPI_LIBRARIES fftw3f_mpi PATHS "${FFTW_ROOT}/*" NO_DEFAULT_PATH)
find_library(FFTW_MPI_LIBRARIES fftw3f_mpi)
include_directories(${FFTW_INCLUDE_DIRS})
link_libraries(${FFTW_LIBRARIES} ${FFTW_MPI_LIBRARIES})

# OpenMP
option(USE_OpenMP "Use OpenMP" OFF)
if(USE_OpenMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# CUDA & CUFFT (& add executable: meraxes)
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
    set(SOURCES ${SOURCES_C} ${SOURCES_CU})

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
    cuda_add_executable(meraxes ${SOURCES})

    # CUFFT
    if(USE_CUFFT)
        CUDA_ADD_CUFFT_TO_TARGET(meraxes)
    endif()

else()
    set(SOURCES ${SOURCES_C})
    add_executable(meraxes ${SOURCES})
endif()

# Profiling
option(ENABLE_PROFILING "Enable profiling with gperftools." OFF)
if(ENABLE_PROFILING)
    find_library(GPERF_LIBRARIES profiler PATHS "${GPERF_ROOT}/*")
    link_libraries(${GPERF_LIBRARIES})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

# math
link_libraries(m)

# executable
install(TARGETS meraxes DESTINATION bin COMPONENT bin)
add_custom_target(install.meraxes
    ${CMAKE_COMMAND} 
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE} 
    -DCOMPONENT=bin
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake) 
add_dependencies(install.meraxes meraxes)

# library
add_library(lib SHARED ${SOURCES})
set_target_properties(lib PROPERTIES OUTPUT_NAME meraxes)
set_property(TARGET lib PROPERTY POSITION_INDEPENDENT_CODE ON)
install(TARGETS lib DESTINATION lib COMPONENT lib)
install(FILES meraxes.h ${CMAKE_BINARY_DIR}/git.h DESTINATION include COMPONENT lib)
add_custom_target(install.lib
    ${CMAKE_COMMAND} 
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE} 
    -DCOMPONENT=lib
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake) 
add_dependencies(install.lib lib)

add_custom_target(setuprun
    COMMENT "Setting up run with minimal required files..."
    COMMAND echo "RUNDIR = \${RUNDIR}"
    COMMAND mkdir -p \${RUNDIR}
    COMMAND cp -n ${CMAKE_CURRENT_SOURCE_DIR}/../input/input.par \${RUNDIR}
    COMMAND cp -n ${CMAKE_CURRENT_SOURCE_DIR}/../input/snaplist.txt \${RUNDIR}
    COMMAND cp ${CMAKE_BINARY_DIR}/meraxes \${RUNDIR}
    )
add_dependencies(setuprun meraxes)

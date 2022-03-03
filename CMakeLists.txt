###################
#  project setup  #
###################

cmake_minimum_required(VERSION 3.18)

# include local custom settings (must come before project line to allow
# compilers to be set)
include(${CMAKE_CURRENT_SOURCE_DIR}/local.cmake OPTIONAL)

project(meraxes C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build type
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
if(CMAKE_BUILD_TYPE EQUAL "Debug" AND NOT BUILD_SHARED_LIBS)
    set(SANITIZE_FLAGS "-fsanitize=address,leak,undefined -fno-common -fno-omit-frame-pointer")
else()
    set(SANITIZE_FLAGS "")
endif()


# We don't want to install to /usr/... by default
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/target/${CMAKE_BUILD_TYPE}" CACHE PATH "..." FORCE)
endif()

# output paths for builds and installs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ensure RPATH is used for builds so we can move the binary around without
# installing first
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# The following line is needed to deal with this annoying problem:
#    https://developers.redhat.com/blog/2015/02/05/gcc5-and-the-c11-abi/
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0")

# need to set extra options for clang
if("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()


#################
# build options #
#################

set(N_HISTORY_SNAPS 5 CACHE STRING "Number of snapshots for delayed SN feedback")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared library version (required for Mhysa)")
option(CALC_MAGS "Calculate magnitude" OFF)
set(MAGS_N_SNAPS 3 CACHE STRING "Number of snapshots to compute magnitudes")
set(MAGS_N_BANDS 6 CACHE STRING "Number of bands to compute")
set(SECTOR_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/sector" CACHE PATH "Base directory of sector library")
option(BUILD_TESTS "Build criterion tests" OFF)
option(GDB "Drop into GDB with mpi_debug_here() calls" OFF)
option(ENABLE_PROFILING "Enable profiling of executable with gperftools." OFF)
option(USE_CUDA "Build with CUDA support for reionization calculations" OFF)

add_subdirectory(src)


#####################
#  Custom commands  #
#####################

# Provide the setuprun build command
if(NOT BUILD_SHARED_LIBS)
    set(INPUT_FILE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/input")
    configure_file(${INPUT_FILE_DIR}/input.par ${CMAKE_BINARY_DIR}/input.par ESCAPE_QUOTES @ONLY)

    add_custom_target(setuprun
        COMMENT "Setting up run with minimal required files..."
        COMMAND test -n \"\${RUNDIR}\" || (echo "*** Must set RUNDIR to execute target 'setuprun'. ***" && exit 1)
        COMMAND echo "RUNDIR = \${RUNDIR}"
        COMMAND mkdir -p \${RUNDIR}
        COMMAND cp -n ${CMAKE_BINARY_DIR}/input.par \${RUNDIR} || (exit 0)
        COMMAND cp ${CMAKE_BINARY_DIR}/bin/meraxes \${RUNDIR}
        )
    add_dependencies(setuprun meraxes)
endif()

if(EXISTS ${CLANG_FORMAT_BINARY})
    add_custom_target(format
        COMMENT "Use clang-format to format the source code consistently."
        COMMAND ${CLANG_FORMAT_BINARY} -i ${SOURCES}
        )
endif()


#########
# TESTS #
#########

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()



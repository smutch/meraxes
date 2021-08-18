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

# We don't want to install to /usr/... by default
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/.. CACHE PATH "..." FORCE)
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


########################################
#  Meraxes lib and executable targets  #
########################################

# TODO: This is bad practice, we should really list the sources explicitly
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core/*.[ch] *.h ${CMAKE_CURRENT_SOURCE_DIR}/physics/*.[ch])
list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/core/meraxes.c)
    
# THE MERAXES TARGETS
if(BUILD_SHARED_LIBS)
    add_library(meraxes_lib SHARED ${SOURCES})
else()
    add_library(meraxes_lib STATIC ${SOURCES})
endif()
set_target_properties(meraxes_lib PROPERTIES
    OUTPUT_NAME meraxes
    PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/meraxes.h;${CMAKE_CURRENT_BINARY_DIR}/meraxes_conf.h"
    C_STANDARD 99
    )
add_executable(meraxes ${CMAKE_CURRENT_SOURCE_DIR}/core/meraxes.c)
set_target_properties(meraxes PROPERTIES C_STANDARD 99)
target_link_libraries(meraxes PRIVATE meraxes_lib)

target_compile_definitions(meraxes_lib PUBLIC
    N_HISTORY_SNAPS=${N_HISTORY_SNAPS}
    MAGS_N_SNAPS=${MAGS_N_SNAPS}
    MAGS_N_BANDS=${MAGS_N_BANDS}
    $<$<BOOL:GDB>:GDB>
    $<$<CONFIG:Debug>:DEBUG>
    )

target_compile_options(meraxes_lib PRIVATE
    $<$<CONFIG:Debug>:-Wall -Wextra -Werror=implicit ${SANITIZE_FLAGS}>
    $<$<CONFIG:Release>:-march=native -ffast-math>
    )

# Profiling
if(ENABLE_PROFILING)
    find_library(GPERF_LIBRARIES profiler PATHS "${GPERF_ROOT}/*")
    target_compile_options(meraxes_lib PRIVATE -pg)
endif()

# Make sure we get all the files we need for an install of the library
install(TARGETS meraxes_lib EXPORT MeraxesConfig
    DESTINATION lib
    PUBLIC_HEADER DESTINATION include
    COMPONENT lib)
add_custom_target(install.lib
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCOMPONENT=lib
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
add_dependencies(install.lib meraxes_lib)

install(TARGETS meraxes DESTINATION bin COMPONENT bin)
add_custom_target(install.meraxes
    ${CMAKE_COMMAND}
    -DBUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCOMPONENT=bin
    -P ${CMAKE_BINARY_DIR}/cmake_install.cmake)
add_dependencies(install.meraxes meraxes)

# Templated headers
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/versioning.cmake)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/meraxes_conf.h.in ${CMAKE_BINARY_DIR}/meraxes_conf.h ESCAPE_QUOTES @ONLY)
target_include_directories(meraxes_lib BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(meraxes BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})


##########
#  CUDA  #
##########

if(USE_CUDA)
    enable_language(CUDA)
    set_target_properties(meraxes_lib meraxes PROPERTIES
        CUDA_ARCHITECTURES 60
        CUDA_STANDARD 11
        CXX_STANDARD 11
        )

    find_package(CUDAToolkit REQUIRED)

    # Add CUDA source code and link libraries
    # TODO: This is bad practice, we should really list the sources explicitly
    file(GLOB SOURCES_CU core/*.cu *.hh meraxes_gpu.h)
    target_sources(meraxes_lib PRIVATE ${SOURCES_CU})

    # link libraries
    target_link_libraries(meraxes_lib PUBLIC CUDA::cudart_static)
    target_compile_definitions(meraxes_lib PUBLIC USE_CUDA)

    # This hack is needed to work around a bug in the fftw3.h header file
    # see the following for a discussion: https://github.com/FFTW/fftw3/issues/18
    target_compile_definitions(meraxes_lib PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:__INTEL_COMPILER>)
endif()


#########
# TESTS #
#########

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


################
# DEPENDENCIES # 
################

# MATH
target_link_libraries(meraxes_lib PUBLIC m)

# MPI
find_package(MPI REQUIRED)
target_link_libraries(meraxes_lib PUBLIC MPI::MPI_C)

# MAGNITUDES
if(CALC_MAGS)
    add_definitions(-DCALC_MAGS)
    file(GLOB SECTOR_SOURCES ${SECTOR_ROOT}/*.c ${SECTOR}/sector.h)
    add_library(sector STATIC ${SECTOR_SOURCES})
    target_link_libraries(sector PRIVATE MPI::MPI_C)
    target_link_libraries(meraxes_lib PUBLIC sector)
    set_target_properties(sector PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        PUBLIC_HEADER "${SECTOR_ROOT}/sector.h")
    target_include_directories(meraxes_lib PUBLIC
        $<BUILD_INTERFACE:${SECTOR_ROOT}>
        $<INSTALL_INTERFACE:include>)
    add_custom_command(
        TARGET sector
        POST_BUILD
        COMMAND objcopy
        ARGS --redefine-sym interp=sector_interp --redefine-sym trapz_table=sector_trapz_table $<TARGET_FILE:sector>
        )
    install(TARGETS sector EXPORT MeraxesConfig
        DESTINATION lib
        PUBLIC_HEADER DESTINATION include
        COMPONENT lib)
endif()

# MLOG
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/mlog/mlog.c")
    message(STATUS "Checking out mlog submodule")
    execute_process(COMMAND git submodule update --init
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
set(MLOG_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/mlog")
add_library(mlog STATIC "${MLOG_ROOT}/mlog.c" "${MLOG_ROOT}/mlog.h")
target_link_libraries(mlog PRIVATE MPI::MPI_C)
target_include_directories(mlog PRIVATE ${MPI_INCLUDE_PATH})
set_target_properties(mlog PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    PUBLIC_HEADER "${MLOG_ROOT}/mlog.h")
target_include_directories(meraxes_lib PUBLIC
    $<BUILD_INTERFACE:${MLOG_ROOT}>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(meraxes_lib PUBLIC mlog)

install(TARGETS mlog EXPORT MeraxesConfig
    DESTINATION lib
    PUBLIC_HEADER DESTINATION include
    COMPONENT lib)

# GSL
find_package(GSL REQUIRED)
target_link_libraries(meraxes_lib PRIVATE GSL::gsl GSL::gslcblas)

# HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)
if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Meraxes requires HDF5 with parallel support.")
endif()
target_include_directories(meraxes_lib PUBLIC ${HDF5_INCLUDE_DIRS})
target_link_libraries(meraxes_lib PUBLIC ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES})
add_definitions(${HDF5_DEFINITIONS})

# FFTW
find_package(FFTW REQUIRED)
target_include_directories(meraxes_lib PUBLIC ${FFTW_INCLUDE_DIRS})
target_link_libraries(meraxes_lib PUBLIC ${FFTW_LIBRARIES})

# Exports for other codes to build against this version of Meraxes
if(CALC_MAGS)
    export(TARGETS meraxes_lib mlog sector FILE MeraxesConfig.cmake)
else()
    export(TARGETS meraxes_lib mlog FILE MeraxesConfig.cmake)
endif()
install(EXPORT MeraxesConfig DESTINATION share/meraxes/cmake)


#####################
#  Custom commands  #
#####################

# Provide the setuprun build command
if(NOT BUILD_SHARED_LIBS)
    set(INPUT_FILE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../input")
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

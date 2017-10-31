# - Try to find fftw
# Once done this will define
#  FFTW_FOUND - System has fftw
#  FFTW_INCLUDE_DIRS - The fftw include directories
#  FFTW_LIBRARIES - The libraries needed to use fftw

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PC_FFTW QUIET fftw)
endif()

set(FFTW_DEFINITIONS ${PC_FFTW_CFLAGS_OTHER})

find_path(FFTW_INCLUDE_DIR fftw3.h
          PATHS ${PC_FFTW_INCLUDEDIR} ${PC_FFTW_INCLUDE_DIRS} "${FFTW_ROOT}/include")

function(find_fftw_lib name dst_var)
    # If we're asked to use static linkage, add it as a preferred library name.
    if(FFTW_USE_STATIC)
        list(APPEND FFTW_NAMES
            "${CMAKE_STATIC_LIBRARY_PREFIX}lib${name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
        list(INSERT FFTW_NAMES 0
            "${CMAKE_STATIC_LIBRARY_PREFIX}lib${name}${CMAKE_STATIC_LIBRARY_SUFFIX}")
    endif()

    list(APPEND FFTW_NAMES ${name})

    find_library(${dst_var} NAMES ${FFTW_NAMES}
        PATHS "${FFTW_ROOT}/lib"
        HINTS ${PC_FFTW_LIBDIR} ${PC_FFTW_LIBRARY_DIRS})
endfunction()

find_fftw_lib(fftw3f FFTW_LIBRARY)
find_fftw_lib(fftw3f_mpi FFTW_MPI_LIBRARY)

set(FFTW_LIBRARIES ${FFTW_LIBRARY} ${FFTW_MPI_LIBRARY})
set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(FFTW DEFAULT_MSG
    FFTW_LIBRARY FFTW_MPI_LIBRARY FFTW_INCLUDE_DIR)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY FFTW_MPI_LIBRARY)

# - Try to find fftw
# Once done this will define
#  FFTW_FOUND - System has fftw
#  FFTW_INCLUDE_DIRS - The fftw include directories
#  FFTW_LIBRARIES - The libraries needed to use fftw

find_package(PkgConfig)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PC_FFTW QUIET fftw3f)
endif()

set(FFTW_DEFINITIONS ${PC_FFTW_CFLAGS_OTHER})

find_path(FFTW_INCLUDE_DIR fftw3.h
          PATHS ${PC_FFTW_INCLUDEDIR} ${PC_FFTW_INCLUDE_DIRS} "${FFTW_ROOT}/include")

find_library(FFTW_MPI_LIBRARY NAME fftw3f_mpi
    PATHS "${FFTW_ROOT}/lib"
    HINTS ${PC_FFTW_LIBDIR} ${PC_FFTW_LIBRARY_DIRS})

find_library(FFTW_LIBRARY NAME fftw3f
    PATHS "${FFTW_ROOT}/lib"
    HINTS ${PC_FFTW_LIBDIR} ${PC_FFTW_LIBRARY_DIRS})

set(FFTW_LIBRARIES ${FFTW_LIBRARY} ${FFTW_MPI_LIBRARY})
set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_MPI_LIBRARY FFTW_INCLUDE_DIR)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY FFTW_MPI_LIBRARY)

function(bundle_fftw)
    ExternalProject_Add(fftw-bundle
        URL http://fftw.org/fftw-3.3.7.tar.gz
        URL_MD5 0d5915d7d39b3253c1cc05030d79ac47
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --enable-type-prefix --enable-mpi --enable-float --enable-shared
        BUILD_COMMAND make -j4
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS 0
        INSTALL_COMMAND make install
        EXCLUDE_FROM_ALL 1
        )
    ExternalProject_Get_Property(fftw-bundle INSTALL_DIR)

    set(FFTW_INCLUDE_DIRS "${INSTALL_DIR}/include")
    set(FFTW_LIBRARY "${INSTALL_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX}")
    set(FFTW_MPI_LIBRARY "${INSTALL_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}fftw3f_mpi${CMAKE_SHARED_LIBRARY_SUFFIX}")

    add_library(fftw SHARED IMPORTED)
    set_target_properties(fftw PROPERTIES IMPORTED_LOCATION ${FFTW_LIBRARY})
    add_library(fftw-mpi SHARED IMPORTED)
    set_target_properties(fftw-mpi PROPERTIES IMPORTED_LOCATION ${FFTW_MPI_LIBRARY})

    if(NOT EXISTS ${FFTW_LIBRARY})
        add_dependencies(fftw fftw-bundle)
        add_dependencies(fftw-mpi fftw-bundle)
    endif()

    set(FFTW_LIBRARIES ${FFTW_MPI_LIBRARY} ${FFTW_LIBRARY})
    target_include_directories(meraxes_lib PRIVATE ${FFTW_INCLUDE_DIRS})
    target_include_directories(meraxes PRIVATE ${FFTW_INCLUDE_DIRS})
    target_link_libraries(meraxes_lib PRIVATE fftw fftw-mpi)
    target_link_libraries(meraxes PRIVATE fftw fftw-mpi)
endfunction()

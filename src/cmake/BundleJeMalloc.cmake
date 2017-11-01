function(bundle_jemalloc)
    ExternalProject_Add(jemalloc-bundle
        URL https://github.com/jemalloc/jemalloc/releases/download/5.0.1/jemalloc-5.0.1.tar.bz2
        URL_MD5 507f7b6b882d868730d644510491d18f
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
        BUILD_COMMAND make -j4
        BUILD_IN_SOURCE 1
        BUILD_ALWAYS 0
        INSTALL_COMMAND make install
        EXCLUDE_FROM_ALL 1
        )
    ExternalProject_Get_Property(jemalloc-bundle INSTALL_DIR)
    add_library(jemalloc SHARED IMPORTED)
    set(JEMALLOC_LIBRARIES ${INSTALL_DIR}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}jemalloc.so)
    set_target_properties(jemalloc PROPERTIES IMPORTED_LOCATION ${JEMALLOC_LIBRARIES})
    if(NOT EXISTS ${JEMALLOC_LIBRARIES})
        add_dependencies(jemalloc jemalloc-bundle)
    endif()
    target_link_libraries(meraxes_lib PRIVATE jemalloc)
    target_link_libraries(meraxes PRIVATE jemalloc)
endfunction()

if(EIGEN3_INCLUDE_DIR)
	set(EIGEN3_FOUND TRUE)
else()
    find_path(EIGEN3_INCLUDE_DIR NAMES 
		Eigen/Core
		PATHS
		${CMAKE_INSTALL_PREFIX}/include
		${CMAKE_CURRENT_SOURCE_DIR}/../eigen
        PATH_SUFFIXES eigen3 eigen        
    )
    include( FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Eigen3 DEFAULT_MSG EIGEN3_INCLUDE_DIR)   
    mark_as_advanced(EIGEN3_INCLUDE_DIR)
endif()
cmake_minimum_required(VERSION 3.17)

set(BGAL_PQP_SRC        
        Build.cpp
        BV.cpp
        PQP.cpp
		TriDist.cpp
		)

# Get static lib
add_library(PQP STATIC ${BGAL_PQP_SRC})
set_target_properties(PQP PROPERTIES VERSION ${VERSION})
set_target_properties(PQP PROPERTIES CLEAN_DIRECT_OUTPUT 1)
target_include_directories(PQP PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)


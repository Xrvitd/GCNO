cmake_minimum_required(VERSION 3.17)

set(BGAL_Tessellation3D_SRC        
        Side3D.cpp
		Tessellation3D.cpp
		)
# Get static lib
add_library(Tessellation3D STATIC ${BGAL_Tessellation3D_SRC})
target_link_libraries(Tessellation3D Algorithm BaseShape Model ${Boost_LIBRARIES})
set_target_properties(Tessellation3D PROPERTIES VERSION ${VERSION})
set_target_properties(Tessellation3D PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(Tessellation3D PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

cmake_minimum_required(VERSION 3.17)

set(BGAL_Tessellation2D_SRC        
        Side2D.cpp
		Tessellation2D.cpp
		)

# Get static lib
add_library(Tessellation2D STATIC ${BGAL_Tessellation2D_SRC})
target_link_libraries(Tessellation2D Algorithm BaseShape ${Boost_LIBRARIES})
set_target_properties(Tessellation2D PROPERTIES VERSION ${VERSION})
set_target_properties(Tessellation2D PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(Tessellation2D PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

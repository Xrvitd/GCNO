cmake_minimum_required(VERSION 3.17)

set(BGAL_Draw_SRC        
        DrawPS.cpp
        )


# Get static lib
add_library(Draw STATIC ${BGAL_Draw_SRC})
target_link_libraries(Draw Algorithm BaseShape ${Boost_LIBRARIES})
set_target_properties(Draw PROPERTIES VERSION ${VERSION})
set_target_properties(Draw PROPERTIES CLEAN_DIRECT_OUTPUT 1)

target_include_directories(Draw PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
	
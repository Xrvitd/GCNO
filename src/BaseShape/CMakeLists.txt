cmake_minimum_required(VERSION 3.17)

set(BGAL_BaseShape_SRC        
        KDTree.cpp
        Line.cpp
        Point.cpp
        Polygon.cpp
        Triangle.cpp
        )
# Get static lib
add_library(BaseShape STATIC ${BGAL_BaseShape_SRC})
target_link_libraries(BaseShape Algorithm ${Boost_LIBRARIES})
set_target_properties(BaseShape PROPERTIES VERSION ${VERSION})
set_target_properties(BaseShape PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(BaseShape PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
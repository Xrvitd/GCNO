cmake_minimum_required(VERSION 3.17)

set(BGAL_Geodesic_SRC        
        Dijkstra/Dijkstra.cpp
        AbstractMethod.cpp
        )

# Get static lib
add_library(Geodesic STATIC ${BGAL_Geodesic_SRC})
target_link_libraries(Geodesic BaseShape Algorithm Model PQP ${Boost_LIBRARIES})
set_target_properties(Geodesic PROPERTIES VERSION ${VERSION})
set_target_properties(Geodesic PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(Geodesic PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
	
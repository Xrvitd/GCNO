cmake_minimum_required(VERSION 3.17)

set(BGAL_Integral_SRC        
        Integral.cpp
        Tetrahedron_arbq_rule.cpp
        )

# Get static lib
add_library(Integral STATIC ${BGAL_Integral_SRC})
target_link_libraries(Integral Algorithm BaseShape ${Boost_LIBRARIES})
set_target_properties(Integral PROPERTIES VERSION ${VERSION})
set_target_properties(Integral PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(Integral PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)
	
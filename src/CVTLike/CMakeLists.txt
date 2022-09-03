cmake_minimum_required(VERSION 3.17)

set(BGAL_CVTLike_SRC        
		CPD.cpp
		CVT.cpp
		)
# Get static lib
add_library(CVTLike STATIC ${BGAL_CVTLike_SRC})
target_link_libraries(CVTLike Algorithm BaseShape Model Tessellation2D Tessellation3D Optimization ${Boost_LIBRARIES})
set_target_properties(CVTLike PROPERTIES VERSION ${VERSION})
set_target_properties(CVTLike PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(CVTLike PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

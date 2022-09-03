cmake_minimum_required(VERSION 3.17)

set(BGAL_Reconstruction_SRC        
        MarchingTerahedra/MarchingTetrahedra.cpp
		)

# Get static lib
add_library(Reconstruction STATIC ${BGAL_Reconstruction_SRC})
target_link_libraries(Reconstruction Algorithm BaseShape Model PQP ${Boost_LIBRARIES})
set_target_properties(Reconstruction PROPERTIES VERSION ${VERSION})
set_target_properties(Reconstruction PROPERTIES CLEAN_DIRECT_OUTPUT 1)


target_include_directories(Reconstruction PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

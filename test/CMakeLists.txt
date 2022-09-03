cmake_minimum_required(VERSION 3.17)
project(mytest)
SET(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH})

find_package(Eigen3 REQUIRED)
if (Eigen3_FOUND)
    message(STATUS "${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif ()

# Get Boost
find_package(Boost REQUIRED)
if (Boost_FOUND)
    message(STATUS "BOOST FOUNDED")
    include_directories(${Boost_INCLUDE_DIRS})
endif ()

# Get CGAL
find_package(CGAL REQUIRED)
if (CGAL_FOUND)
    include(${CGAL_USE_FILE})
else ()
    message("ERROR: this program requires CGAL and will not be compiled.")
endif ()

find_package(BGAL REQUIRED)
if(BGAL_FOUND)
include_directories(${BGAL_DIR})
endif()
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC BGAL::Algorithm BGAL::BaseShape BGAL::Draw BGAL::Geodesic BGAL::Integral BGAL::Model BGAL::Optimization BGAL::PointCloudProcessing BGAL::PQP BGAL::Reconstruction BGAL::Tessellation2D BGAL::Tessellation3D BGAL::CVTLike ${Boost_LIBRARIES})

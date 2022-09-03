if (BGAL_FOUND)
    return()
endif ()

find_path(BGAL_INCLUDE_DIR BGAL/LBFGS.h
        HINTS
        ${BGAL}
        ENV BGAL
        PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/BGAL
        ${CMAKE_SOURCE_DIR}/../BGAL
        ${CMAKE_SOURCE_DIR}/../../BGAL
        /usr/local
        /usr/local/BGAL
        PATH_SUFFIXES include
        )
find_path(BGAL_LIB_DIR libBGAL.dylib
        HINTS
        ${BGAL}
        ENV BGAL
        PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/BGAL
        ${CMAKE_SOURCE_DIR}/../BGAL
        ${CMAKE_SOURCE_DIR}/../../BGAL
        /usr/local
        /usr/local/BGAL
        PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BGAL
        "\nBGAL not found ${CMAKE_SOURCE_DIR}/../BGAL"
        BKHao_INCLUDE_DIR)
mark_as_advanced(BGAL_INCLUDE_DIR)

#include(BGAL)
# Get Eigen3
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
cmake_minimum_required(VERSION 3.17)

set(BGAL_Optimization_SRC        
        ALGLIB/alglibinternal.cpp
        ALGLIB/alglibmisc.cpp
        ALGLIB/ap.cpp
        ALGLIB/dataanalysis.cpp
        ALGLIB/diffequations.cpp
        ALGLIB/fasttransforms.cpp
        ALGLIB/integration.cpp
        ALGLIB/interpolation.cpp
        ALGLIB/linalg.cpp
        ALGLIB/optimization.cpp
        ALGLIB/solvers.cpp
        ALGLIB/specialfunctions.cpp
        ALGLIB/statistics.cpp
        LBFGS/LBFGS.cpp
        GradientDescent/GradientDescent.cpp
        LinearSystem/LinearSystem.cpp
        LBFGS/LBFGS.cpp
        )

# Get static lib
add_library(Optimization STATIC ${BGAL_Optimization_SRC})
set_target_properties(Optimization PROPERTIES VERSION ${VERSION})
set_target_properties(Optimization PROPERTIES CLEAN_DIRECT_OUTPUT 1)
target_include_directories(Optimization PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

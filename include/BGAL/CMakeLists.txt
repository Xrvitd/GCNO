cmake_minimum_required(VERSION 3.17)

set(BGAL_HEADERS
        # Algorithm
        Algorithm/BOC/BOC.h
        # BK_BaseShape
        BaseShape/KDTree.h
        BaseShape/Line.h
        BaseShape/Point.h
        BaseShape/Polygon.h
        BaseShape/Triangle.h
        # Draw
        Draw/DrawPS.h
        # Geodesic
        Geodesic/Dijkstra/Dijkstra.h
        Geodesic/AbstractMethod.h
        # Intergral
        Integral/Integral.h
        Integral/Tetrahedron_arbq_rule.h
        # BK_Model
        Model/ManifoldModel.h
        Model/Model.h
        Model/Model_Iterator.h
        # Optimization
        Optimization/ALGLIB/alglibinternal.h
        Optimization/ALGLIB/alglibmisc.h
        Optimization/ALGLIB/ap.h
        Optimization/ALGLIB/dataanalysis.h
        Optimization/ALGLIB/diffequations.h
        Optimization/ALGLIB/fasttransforms.h
        Optimization/ALGLIB/integration.h
        Optimization/ALGLIB/interpolation.h
        Optimization/ALGLIB/linalg.h
        Optimization/ALGLIB/optimization.h
        Optimization/ALGLIB/solvers.h
        Optimization/ALGLIB/specialfunctions.h
        Optimization/ALGLIB/statistics.h
        Optimization/LBFGS/LBFGS.h
        Optimization/GradientDescent/GradientDescent.h
        Optimization/LinearSystem/LinearSystem.h
        # PointCloudProcessing
        PointCloudProcessing/Registration/ICP/ICP.h
        # PQP
        PQP/Build.h
        PQP/BV.h
        PQP/BVTQ.h
        PQP/MatVec.h
        PQP/OBB_Disjoint.h
        PQP/PQP.h
        PQP/PQP_Compile.h
        PQP/PQP_Internal.h
        PQP/RectDist.h
        PQP/Tri.h
        PQP/TriDist.h
        PQP/GetTime.h
        # Reconstruction
        Reconstruction/MarchingTetrahedra/MarchingTetrahedra.h
        # Tessellation2D
        Tessellation2D/Side2D.h
        Tessellation2D/Tessellation2D.h
        # Tessellation3D
        Tessellation3D/Side3D.h
        Tessellation3D/Tessellation3D.h
        # CVTLike
        CVTLike/CPD.h
        )

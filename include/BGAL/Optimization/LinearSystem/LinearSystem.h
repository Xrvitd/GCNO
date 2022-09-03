#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace BGAL
{
    class _LinearSystem
    {
    public:
        static Eigen::VectorXd solve_ldlt(const Eigen::SparseMatrix<double> &CoeffMat,
                                          const Eigen::VectorXd &right);
        static Eigen::VectorXd solve_ldlt(const Eigen::SparseMatrix<double> &CoeffMat,
                                          const Eigen::VectorXd &right,
                                          const double &pinvtoler);
        static Eigen::VectorXd solve_llt(const Eigen::SparseMatrix<double> &CoeffMat,
                                         const Eigen::VectorXd &right);
    };
} // namespace BGAL
#include "BGAL/Optimization/LinearSystem/LinearSystem.h"
namespace BGAL
{
  Eigen::VectorXd _LinearSystem::solve_ldlt(const Eigen::SparseMatrix<double> &CoeffMat, const Eigen::VectorXd &right)
  {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(CoeffMat);
    Eigen::VectorXd res = ldlt.solve(right);
    return res;
  }
  Eigen::VectorXd _LinearSystem::solve_ldlt(const Eigen::SparseMatrix<double> &CoeffMat,
                                            const Eigen::VectorXd &right,
                                            const double &pinvtoler)
  {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(CoeffMat);
    Eigen::VectorXd X_0 = ldlt.permutationP() * right;
    Eigen::VectorXd X_1 = ldlt.matrixL().solve(X_0);
    Eigen::VectorXd X_2(ldlt.vectorD().size());
    X_2.setZero();
    for (int i = 0; i < ldlt.vectorD().size(); ++i)
      if (abs(ldlt.vectorD()(i)) > pinvtoler)
        X_2[i] = X_1[i] / ldlt.vectorD()(i);
    //else
    //	X_2[i] = X_1[i] / pinvtoler;
    Eigen::VectorXd X_3 = ldlt.matrixU().solve(X_2);
    return ldlt.permutationPinv() * X_3;
  }
  Eigen::VectorXd _LinearSystem::solve_llt(const Eigen::SparseMatrix<double> &CoeffMat, const Eigen::VectorXd &right)
  {
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> cho(CoeffMat);
    Eigen::VectorXd res = cho.solve(right);
    return res;
  }
} // namespace BGAL
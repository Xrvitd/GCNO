#include "BGAL/PointCloudProcessing/Registration/ICP/ICP.h"
#include "BGAL/BaseShape/KDTree.h"
#include <Eigen/SVD>
namespace BGAL
{
  _ICP::_ICP()
  {
    _static_points.clear();
    _max_iteration = 400;
    _num_random_sample = 400;
    _epsilon = 1e-8;
  }
  _ICP::_ICP(const std::vector<_Point3> &static_points)
      : _static_points(static_points)
  {
    _max_iteration = 400;
    _num_random_sample = 400 > (int)(_static_points.size()) / 5 ? 400 : (int)(_static_points.size());
    _epsilon = 1e-8;
  }
  Eigen::Matrix4d _ICP::registration_(std::vector<_Point3> dynamic_points) const
  {
    int __num_random_sample = _num_random_sample;
    if (_num_random_sample > dynamic_points.size())
      __num_random_sample = dynamic_points.size();
    Eigen::Matrix3d rotationM;
    _Point3 translationV(0, 0, 0);
    rotationM.setIdentity();
    Eigen::Matrix4d RTM;
    RTM.setIdentity();
    _Point3 dynamic_mid(0.0, 0.0, 0.0);
    _Point3 static_mid = dynamic_mid;
    for (int i = 0; i < dynamic_points.size(); ++i)
    {
      dynamic_mid += dynamic_points[i];
    }
    dynamic_mid /= dynamic_points.size();
    for (int i = 0; i < _static_points.size(); ++i)
    {
      static_mid += _static_points[i];
    }
    static_mid /= _static_points.size();
    _KDTree tree(_static_points);
    double cost = std::numeric_limits<double>::max();
    int iter = 0;
    while (1)
    {
      cost = 0;
      Eigen::Matrix3d U;
      U.setZero();
      dynamic_mid = _Point3(0.0, 0.0, 0.0);
      for (int i = 0; i < dynamic_points.size(); ++i)
      {
        dynamic_mid += dynamic_points[i];
      }
      dynamic_mid /= dynamic_points.size();
      for (int i = 0; i < __num_random_sample; ++i)
      {
        int randSample = std::rand() % dynamic_points.size();
        _Point3 d = dynamic_points[randSample];
        _Point3 s = _static_points[tree.search_(d)];
        _Point3 ds = s - static_mid;
        _Point3 dd = d - dynamic_mid;
        U(0, 0) = U(0, 0) + ds.x() * dd.x();
        U(0, 1) = U(0, 1) + ds.x() * dd.y();
        U(0, 2) = U(0, 2) + ds.x() * dd.z();
        U(1, 0) = U(1, 0) + ds.y() * dd.x();
        U(1, 1) = U(1, 1) + ds.y() * dd.y();
        U(1, 2) = U(1, 2) + ds.y() * dd.z();
        U(2, 0) = U(2, 0) + ds.z() * dd.x();
        U(2, 1) = U(2, 1) + ds.z() * dd.y();
        U(2, 2) = U(2, 2) + ds.z() * dd.z();
        cost += (s - d).sqlength_();
      }
      Eigen::Matrix4d iRTM;
      iRTM.setIdentity();
      iRTM.block<3, 3>(0, 0) = rotationM;
      iRTM(0, 3) = translationV.x();
      iRTM(1, 3) = translationV.y();
      iRTM(2, 3) = translationV.z();
      Eigen::Matrix4d nRTM = iRTM * RTM;
      RTM = nRTM;
      if (cost < _epsilon || iter == _max_iteration)
        break;
      Eigen::JacobiSVD<Eigen::Matrix3d> svd(U, Eigen::ComputeFullU | Eigen::ComputeFullV);
      Eigen::Matrix3d svdU = svd.matrixU();
      Eigen::Matrix3d svdV = svd.matrixV().transpose();
      rotationM = svdU * svdV;
      translationV = static_mid - dynamic_mid.rotate_(rotationM);
      for (int i = 0; i < dynamic_points.size(); ++i)
      {
        dynamic_points[i] = dynamic_points[i].rotate_(rotationM) + translationV;
      }
      ++iter;
    }
    return RTM;
  }
  Eigen::Matrix4d _ICP::registration_(std::vector<_Point3> dynamic_points, const std::vector<int> &feature_ids) const
  {
    Eigen::Matrix3d rotationM;
    _Point3 translationV(0, 0, 0);
    rotationM.setIdentity();
    Eigen::Matrix4d RTM;
    RTM.setIdentity();
    _Point3 dynamic_mid(0.0, 0.0, 0.0);
    _Point3 static_mid = dynamic_mid;
    for (int i = 0; i < dynamic_points.size(); ++i)
    {
      dynamic_mid += dynamic_points[i];
    }
    dynamic_mid /= dynamic_points.size();
    for (int i = 0; i < _static_points.size(); ++i)
    {
      static_mid += _static_points[i];
    }
    static_mid /= _static_points.size();
    _KDTree tree(_static_points);
    double cost = std::numeric_limits<double>::max();
    int iter = 0;
    while (1)
    {
      cost = 0;
      Eigen::Matrix3d U;
      U.setZero();
      dynamic_mid = _Point3(0.0, 0.0, 0.0);
      for (int i = 0; i < dynamic_points.size(); ++i)
      {
        dynamic_mid += dynamic_points[i];
      }
      dynamic_mid /= dynamic_points.size();
      for (int i = 0; i < feature_ids.size(); ++i)
      {
        _Point3 d = dynamic_points[feature_ids[i]];
        _Point3 s = _static_points[tree.search_(d)];
        _Point3 ds = s - static_mid;
        _Point3 dd = d - dynamic_mid;
        U(0, 0) = U(0, 0) + ds.x() * dd.x();
        U(0, 1) = U(0, 1) + ds.x() * dd.y();
        U(0, 2) = U(0, 2) + ds.x() * dd.z();
        U(1, 0) = U(1, 0) + ds.y() * dd.x();
        U(1, 1) = U(1, 1) + ds.y() * dd.y();
        U(1, 2) = U(1, 2) + ds.y() * dd.z();
        U(2, 0) = U(2, 0) + ds.z() * dd.x();
        U(2, 1) = U(2, 1) + ds.z() * dd.y();
        U(2, 2) = U(2, 2) + ds.z() * dd.z();
        cost += (s - d).sqlength_();
      }
      Eigen::Matrix4d iRTM;
      iRTM.setIdentity();
      iRTM.block<3, 3>(0, 0) = rotationM;
      iRTM(0, 3) = translationV.x();
      iRTM(1, 3) = translationV.y();
      iRTM(2, 3) = translationV.z();
      Eigen::Matrix4d nRTM = iRTM * RTM;
      RTM = nRTM;
      if (cost < _epsilon || iter == _max_iteration)
        break;
      Eigen::JacobiSVD<Eigen::Matrix3d> svd(U, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::Matrix3d svdU = svd.matrixU();
      Eigen::Matrix3d svdV = svd.matrixV().transpose();
      rotationM = svdU * svdV;
      translationV = static_mid - dynamic_mid.rotate_(rotationM);
      for (int i = 0; i < dynamic_points.size(); ++i)
      {
        dynamic_points[i] = dynamic_points[i].rotate_(rotationM) + translationV;
      }
      ++iter;
    }
    return RTM;
  }
} // namespace BGAL

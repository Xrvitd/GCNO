#pragma once
#include <vector>
#include <Eigen/Dense>
#include "BGAL/BaseShape/Point.h"

namespace BGAL
{
	class _ICP
	{
	public:
		_ICP();
		_ICP(const std::vector<_Point3> &static_points);
		Eigen::Matrix4d registration_(std::vector<_Point3> dynamic_points) const;
		Eigen::Matrix4d registration_(std::vector<_Point3> dynamic_points, const std::vector<int> &feature_ids) const;
		void set_max_iteration(const int &max_iteration)
		{
			_max_iteration = max_iteration;
		}
		void set_epsilon(const double &epsilon)
		{
			_epsilon = epsilon;
		}
		void set_num_random_sample(const int &num_random_sample)
		{
			_num_random_sample = num_random_sample;
		}

	private:
		std::vector<_Point3> _static_points;
		int _max_iteration;
		int _num_random_sample;
		double _epsilon;
	};
} // namespace BGAL
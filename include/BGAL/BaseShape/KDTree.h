#pragma once
#include "Point.h"
namespace BGAL 
{
	class _KDTree 
	{
	public:
		_KDTree();
		_KDTree(const std::vector<_Point3>& in_points);
		~_KDTree();
		int search_(const _Point3& in_p) const;
		int search_(const _Point3& in_p, double& min_dist) const;
		std::vector<int> nsearch_(const std::vector<_Point3>& in_ps, const int& k) const;
		std::vector<int> rsearch_(const std::vector<_Point3>& in_p, const double& in_r) const;
	private:
		void build_(const std::vector<_Point3>& in_points);
		void clear_();
		struct _Node 
		{
			int _id;
			_Node* _next[2];
			int _axis;
			_Node() : _id(-1), _axis(-1) 
			{
				_next[0] = _next[1] = nullptr;
			}
		};
	private:
		_Node* _root;
		std::vector<_Point3> _points;
	};
}
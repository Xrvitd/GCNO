#pragma once
#include "Point.h"

namespace BGAL 
{
	class _Segment3 
	{
	private:
		_Point3 _s;
		_Point3 _t;
	public:
		_Segment3();
		_Segment3(const _Point3& in_s, const _Point3& in_t);
		inline const _Point3& sp_() const 
		{
			return _s;
		}
		inline const _Point3& tp_() const 
		{
			return _t;
		}
		inline double length_() const 
		{
			return (_s - _t).length_();
		}
	};
}
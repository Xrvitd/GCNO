#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <cmath>
namespace BGAL
{
	class _BOC
	{
	private:
		static double _precision;
	public:
		enum class _Sign 
		{
			NegativE, ZerO, PositivE, FaileD
		};
		static inline double rand_() 
		{
			return double(rand() / (RAND_MAX + 1.0));
		}
		static inline _Sign sign_(const double _real) 
		{
			if (abs(_real) < _precision)
				return _Sign::ZerO;
			else if (_real > 0)
				return _Sign::PositivE;
			else
				return _Sign::NegativE;
		}
		static inline double PI() 
		{
			return 3.1415926535897932384626433832795;
		}
		static inline void set_precision_(const double in_precision)
		{
			_precision = fabs(in_precision);
		}
		static inline double precision_()
		{
			return _precision;
		}
		//  static int search_files_(const std::string &path, const std::string &ext, std::vector<std::string> &filenames) 
		//{
		//
		//  }
	};
}



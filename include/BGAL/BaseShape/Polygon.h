#pragma once
#include "Point.h"
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Polygon_2.h>
namespace BGAL 
{
	class _Polygon 
	{
	private:
		std::vector<_Point2> _points;
		bool _open_in;
		std::pair<_Point2, _Point2> _bounding_box;
	public:
		struct FaceInfo2 
		{
			FaceInfo2() 
			{

			}
			int nesting_level;
			bool in_domain() 
			{
				return nesting_level % 2 == 1;
			}
		};
		typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
		typedef CGAL::Triangulation_vertex_base_2<K> Vb;
		typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Fbb;
		typedef CGAL::Constrained_triangulation_face_base_2<K, Fbb> Fb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
		typedef CGAL::Exact_predicates_tag Itag;
		typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
		typedef CDT::Point CGAL_Point;
		typedef CGAL::Polygon_2<K> CGAL_Polygon_2;
		typedef CDT::Face_handle CGAL_Face_handle;
		static void mark_domain(CDT& cdt) 
		{
			std::function<void(CDT& ct, CGAL_Face_handle start, int index, std::list<CDT::Edge>& boreder)> mark_domains =
				[](CDT& ct, CGAL_Face_handle start, int index, std::list<CDT::Edge>& border) 
			{
				if (start->info().nesting_level != -1) 
				{
					return;
				}
				std::list<CGAL_Face_handle> queue;
				queue.push_back(start);
				while (!queue.empty()) 
				{
					CGAL_Face_handle fh = queue.front();
					queue.pop_front();
					if (fh->info().nesting_level == -1) 
					{
						fh->info().nesting_level = index;
						for (int i = 0; i < 3; i++) 
						{
							CDT::Edge e(fh, i);
							CGAL_Face_handle n = fh->neighbor(i);
							if (n->info().nesting_level == -1) 
							{
								if (ct.is_constrained(e)) border.push_back(e);
								else queue.push_back(n);
							}
						}
					}
				}
			};
			for (CDT::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it) 
			{
				it->info().nesting_level = -1;
			}
			std::list<CDT::Edge> border;
			mark_domains(cdt, cdt.infinite_face(), 0, border);
			while (!border.empty()) 
			{
				CDT::Edge e = border.front();
				border.pop_front();
				CGAL_Face_handle n = e.first->neighbor(e.second);
				if (n->info().nesting_level == -1) 
				{
					mark_domains(cdt, n, e.first->info().nesting_level + 1, border);
				}
			}
		}

		_Polygon() 
		{

		}
		_Polygon(const std::vector<_Point2>& in_points);
		void start_();
		void insert_(const _Point2& in_p);
		void insert_(const double& in_x, const double& in_y);
		void end_();
		inline int num_() const 
		{
			return _points.size();
		}
		inline const std::pair<_Point2, _Point2>& bounding_box_() const 
		{
			return _bounding_box;
		}
		const _Point2& operator[](const int& in_inx) const;
		_Point2& operator[](const int& in_inx);
		bool is_in_(const _Point2& in_p) const;
		_Point2 nearest_point_(const _Point2& in_p);
		double distance_to_boundary_(const _Point2& in_p);
		int intersection_with_linesegment_(const _Point2& p1, const _Point2& p2, std::vector<std::pair<int, _Point2>>& intersections) const;

		std::vector<_Polygon> constrained_delaunay_triangulation_() const 
		{
			CGAL_Polygon_2 poly;
			for (int i = 0; i < _points.size(); ++i) 
			{
				poly.push_back(CGAL_Point(_points[i].x(), _points[i].y()));
			}
			CDT cdt;
			cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);
			mark_domain(cdt);
			std::vector<_Polygon> triangles;
			for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) 
			{
				if (fit->info().in_domain()) 
				{
					_Polygon tri;
					tri.start_();
					tri.insert_(_Point2(fit->vertex(0)->point().x(), fit->vertex(0)->point().y()));
					tri.insert_(_Point2(fit->vertex(1)->point().x(), fit->vertex(1)->point().y()));
					tri.insert_(_Point2(fit->vertex(2)->point().x(), fit->vertex(2)->point().y()));
					tri.end_();
					triangles.push_back(tri);
				}
			}
			return triangles;
		}
		double area_() const;
		double triangle_area_() const;
		double circumference_() const;
	};
}
#pragma once
#include <set>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "BGAL/BaseShape/Point.h"
#include "BGAL/BaseShape/Polygon.h"
#include "BGAL/BaseShape/Triangle.h"
#include "BGAL/PQP/PQP.h"
// TODO: Path '/' support
namespace BGAL
{
	class _Face_Iterator;
	class _Vertex_Iterator;
	class _FV_Iterator;
	class _Model
	{
		friend class _Face_Iterator;
		friend class _Vertex_Iterator;
		friend class _FV_Iterator;
	public:
		class _MFace : public _Triangle3
		{
		public:
			std::vector<int> _vertices;
			int id;
			_MFace();
			_MFace(const int& id1, const int& id2, const int& id3);
			_MFace(const int& id1, const int& id2, const int& id3, const _Point3& p1, const _Point3& p2, const _Point3& p3);
			int& operator[](int index);
			int operator[](int index) const;
			bool operator<(const _MFace& other) const;
		};
		class _PQP_Query_Resutl
		{
		public:
			int _pos_flag;
			int _triangle_id;
			double _distance;
			_Point3 _nearest_point;
			_PQP_Query_Resutl();
			_PQP_Query_Resutl(const int& in_pos_flag, const int& in_triangle_id, const double& in_distance, const _Point3& in_nearset_point);
		};
		_Model();
		_Model(const std::string& in_file_name);
		void set_name_(const std::string& in_name)
		{
			_name = in_name;
		}
		inline int number_vertices_() const
		{
			return _vertices.size();
		}
		inline int number_faces_() const
		{
			return _faces.size();
		}
		inline const _Point3& vertex_(const int& id) const
		{
			if (id < 0 || id >= _vertices.size())
				throw std::runtime_error("Beyond the index!");
			return _vertices[id];
		}
		inline const _Point3& normal_(const int& id) const
		{
			if (id < 0 || id >= _normals_vertex.size())
				throw std::runtime_error("Beyond the index!");
			return _normals_vertex[id];
		}
		inline const _Point3& normal_vertex_(const int& id) const
		{
			if (id < 0 || id >= _normals_vertex.size())
				throw std::runtime_error("Beyond the index!");
			return _normals_vertex[id];
		}
		inline const _Point3& normal_face_(const int& id) const
		{
			if (id < 0 || id >= _normals_face.size())
				throw std::runtime_error("Beyond the index!");
			return _normals_face[id];
		}
		inline const _Model::_MFace& face_(const int& id) const
		{
			if (id < 0 || id >= _faces.size())
				throw std::runtime_error("Beyond the index!");
			return _faces[id];
		}
		_Face_Iterator face_begin() const;
		inline int face_end() const
		{
			return number_faces_();
		}
		_Vertex_Iterator vertex_begin() const;
		inline int vertex_end() const
		{
			return number_vertices_();
		}
		_FV_Iterator fv_begin(const int& fid) const;
		inline int fv_end(const int& fid) const
		{
			if (fid < 0 || fid >= number_faces_())
			{
				throw std::runtime_error("Beyond the index!");
			}
			return 3;
		}
		void save_obj_file_(const std::string& in_file_name) const;
		void save_scalar_field_obj_file_(const std::string& in_file_name, const std::vector<double>& vals) const;
		void save_pamametrization_obj_file_(const std::string& in_file_name, const std::vector<std::pair<double, double>>& uvs) const;
		inline std::pair<_Point3, _Point3> bounding_box_() const
		{
			return _bounding_box;
		}
		inline std::string name_() const
		{
			return _name;
		}
		void initialization_PQP_();
		std::tuple<_Point3, double, int> nearest_point_(const _Point3& in_point);
		double signed_distance_(const _Point3& in_point);
		double signed_distance_(const _Point3& in_point, _Point3& gradient);
		bool is_in_(const _Point3& in_point);
	protected:
		void read_file_(const std::string& in_file_name);
		void read_obj_file_(const std::string& in_file_name);
		void read_off_file_(const std::string& in_file_name);
		void compute_normal_boundingbox_();
		_PQP_Query_Resutl proximity_query_(const _Point3& in_point);
	protected:
		std::vector<_MFace> _faces;
		std::vector<_Point3> _vertices;
		std::vector<_Point3> _normals_vertex;
		std::vector<_Point3> _normals_face;
		std::set<int> _faces_useless;
		std::string _name;
		std::pair<_Point3, _Point3> _bounding_box;
	private:
		PQP_Model _pqp_model;
	};
}
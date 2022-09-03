#pragma once
#include "Model.h"
#include "BGAL/BaseShape/Line.h"
#include <map>
namespace BGAL {
	class _Edge_Iterator;
	class _FE_Iterator;
	class _FF_Iterator;
	class _VV_Iterator;
	class _VE_Iterator;
	class _VF_Iterator;
	class _ManifoldModel : public _Model 
	{
		friend class _Edge_Iterator;
		friend class _FE_Iterator;
		friend class _FF_Iterator;
		friend class _VV_Iterator;
		friend class _VE_Iterator;
		friend class _VF_Iterator;
	public:
		class _MMEdge : public _Segment3 
		{
		public:
			int _id_left_vertex;
			int _id_right_vertex;
			int _id_opposite_vertex;
			int _id_left_edge;
			int _id_right_edge;
			int _id_reverse_edge;
			int _id_face;
			_MMEdge();
			_MMEdge(const _Point3& in_s, const _Point3& in_t);
		};
		_ManifoldModel();
		_ManifoldModel(const std::string& in_file_name);
		_ManifoldModel(const std::vector<_Point3>& in_vertices, const std::vector<_Model::_MFace>& in_faces);
		_ManifoldModel(const _ManifoldModel& in_mmodel);
		void preprocess_model_();
		inline int number_edges_() const 
		{
			return _edges.size();
		}
		inline _MMEdge edge_(const int& id) const 
		{
			if (id < 0 || id >= _edges.size())
				throw std::runtime_error("Beyond the index!");
			return _edges[id];
		}
		_Edge_Iterator edge_begin() const;
		inline int edge_end() const 
		{
			return number_edges_();
		}
		_FE_Iterator fe_begin(const int& fid) const;
		inline int fe_end(const int& fid) const 
		{
			if (fid < 0 || fid >= number_faces_()) 
			{
				throw std::runtime_error("Beyond the index!");
			}
			return 3;
		}
		_FF_Iterator ff_begin(const int& fid) const;
		inline int ff_end(const int& fid) const 
		{
			if (fid < 0 || fid >= number_faces_()) 
			{
				throw std::runtime_error("Beyond the index!");
			}
			return 3;
		}
		_VV_Iterator vv_begin(const int& vid) const;
		inline int vv_end(const int& vid) const 
		{
			if (vid < 0 || vid >= number_vertices_()) 
			{
				throw std::runtime_error("Beyond the index!");
			}
			if (_degree_of_vertices[vid] == 0) 
			{
				return 0;
			}
			if (edge_(edge_(_neight_edge_of_vertices[vid])._id_reverse_edge)._id_face == -1) 
			{
				return _degree_of_vertices[vid] + 1;
			}
			else 
			{
				return _degree_of_vertices[vid];
			}
		}
		_VE_Iterator ve_begin(const int& vid) const;
		inline int ve_end(const int& vid) const 
		{
			if (vid < 0 || vid >= number_vertices_()) 
			{
				throw std::runtime_error("Beyond the index!");
			}
			if (_degree_of_vertices[vid] == 0) 
			{
				return 0;
			}
			if (edge_(edge_(_neight_edge_of_vertices[vid])._id_reverse_edge)._id_face == -1) 
			{
				return _degree_of_vertices[vid] + 1;
			}
			else 
			{
				return _degree_of_vertices[vid];
			}
		}
		_VF_Iterator vf_begin(const int& vid) const;
		inline int vf_end(const int& vid) const 
		{
			if (vid < 0 || vid >= number_vertices_()) 
			{
				throw std::runtime_error("Beyond the index!");
			}
			return _degree_of_vertices[vid];
		}
		inline int degree_of_vertex_(const int& vid) const 
		{
			if (vid < 0 || vid >= number_vertices_()) 
			{
				throw std::runtime_error("Beyond the index!");
			}
			return _degree_of_vertices[vid];
		}
		inline bool is_boundary_vertex_(const int& vid) const
		{
			if (vid < 0 || vid >= number_vertices_())
			{
				throw std::runtime_error("Beyond the index!");
			}
			if (_neight_edge_of_vertices[vid] != -1 && _edges[_edges[_neight_edge_of_vertices[vid]]._id_reverse_edge]._id_face == -1)
				return true;
			else
				return false;
		}
		inline bool is_boundary_edge_(const int& eid) const
		{
			if (eid < 0 || eid >= number_edges_())
			{
				throw std::runtime_error("Beyond the index!");
			}
			if (_edges[eid]._id_face == -1)
				return true;
			else
				return false;
		}
	protected:
		void creat_edges_from_vertices_faces_();
		void arrange_neighs_of_vertex_face_();
	protected:
		std::vector<_MMEdge> _edges;
		std::vector<int> _neight_edge_of_vertices;
		std::vector<int> _neigh_edge_of_faces;
		std::vector<int> _isolated_vertices;
		std::vector<int> _degree_of_vertices;
	};
}
#pragma once
#include "BGAL/Model/Model.h"
#include "BGAL/Model/Model_Iterator.h"
namespace BGAL
{
  _Model::_MFace::_MFace() : _Triangle3(), id(-1)
  {
    _vertices.resize(3, -1);
  }
  _Model::_MFace::_MFace(const int &id1, const int &id2, const int &id3) : _Triangle3(), id(-1)
  {
    _vertices.resize(3);
    _vertices[0] = id1;
    _vertices[1] = id2;
    _vertices[2] = id3;
  }
  _Model::_MFace::_MFace(const int &id1,
                         const int &id2,
                         const int &id3,
                         const _Point3 &p1,
                         const _Point3 &p2,
                         const _Point3 &p3)
      : _Triangle3(p1, p2, p3), id(-1)
  {
    _vertices.resize(3);
    _vertices[0] = id1;
    _vertices[1] = id2;
    _vertices[2] = id3;
  }
  int &_Model::_MFace::operator[](int index)
  {
    return _vertices[index];
  }
  int _Model::_MFace::operator[](int index) const
  {
    return _vertices[index];
  }
  bool _Model::_MFace::operator<(const _MFace &other) const
  {
    if (_vertices[0] < other[0])
      return true;
    else if (_vertices[0] > other[0])
      return false;
    else
    {
      if (_vertices[1] < other[1])
        return true;
      else if (_vertices[1] > other[1])
        return false;
      else
      {
        if (_vertices[2] < other[2])
          return true;
        else
          return false;
      }
    }
  }
  _Model::_Model() : _name("")
  {
  }
  _Model::_Model(const std::string &in_file_name)
  {
    int folder_loc = in_file_name.rfind("\\") > in_file_name.rfind("/") ? in_file_name.rfind("\\") : in_file_name.rfind("/");
    int dot_loc = in_file_name.rfind('.');
    _name = in_file_name.substr(folder_loc + 1, dot_loc - folder_loc - 1);
    read_file_(in_file_name);
    compute_normal_boundingbox_();
  }
  _Face_Iterator _Model::face_begin() const
  {
    return _Face_Iterator(this, 0);
  }

  _Vertex_Iterator _Model::vertex_begin() const
  {
    return _Vertex_Iterator(this, 0);
  }

  _FV_Iterator _Model::fv_begin(const int &fid) const
  {
    if (fid < 0 || fid >= number_faces_())
    {
      throw std::runtime_error("Beyond the index!");
    }
    return _FV_Iterator(this, fid, 0);
  }
  void _Model::save_obj_file_(const std::string &in_file_name) const
  {
    std::ofstream out(in_file_name);
    int folder_loc = in_file_name.rfind("\\") > in_file_name.rfind("/") ? in_file_name.rfind("\\") : in_file_name.rfind("/");
    int dot_loc = in_file_name.rfind('.');
    out << "g "
        << in_file_name.substr(folder_loc + 1, dot_loc - folder_loc - 1)
        << std::endl;
    for (int i = 0; i < _vertices.size(); ++i)
    {
      out << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << _vertices[i].x() << " "
          << _vertices[i].y() << " " << _vertices[i].z() << std::endl;
    }
    for (int i = 0; i < _faces.size(); ++i)
    {
      if (_faces_useless.find(i) != _faces_useless.end())
      {
        continue;
      }
      out << "f " << _faces[i][0] + 1 << " " << _faces[i][1] + 1 << " " << _faces[i][2] + 1 << std::endl;
    }
    out.close();
  }
  void _Model::save_scalar_field_obj_file_(const std::string &in_file_name, const std::vector<double> &vals) const
  {
    std::ofstream out(in_file_name);
    int folder_loc = in_file_name.rfind("\\") > in_file_name.rfind("/") ? in_file_name.rfind("\\") : in_file_name.rfind("/");
    int dot_loc = in_file_name.rfind('.');
    out << "g " << in_file_name.substr(folder_loc + 1, dot_loc - folder_loc - 1) << std::endl;
    out << "# max value: " << *std::max_element(vals.begin(), vals.end()) << std::endl;
    out << "mtllib BKLineColorBar.mtl" << std::endl;
    out << "usemtl BKLineColorBar" << std::endl;
    for (int i = 0; i < number_vertices_(); ++i)
    {
      out << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << vertex_(i) << std::endl;
    }
    for (int i = 0; i < vals.size(); ++i)
    {
      out << "vt " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << vals[i] << " 0" << std::endl;
    }
    for (int i = 0; i < number_faces_(); ++i)
    {
      if (_faces_useless.find(i) != _faces_useless.end())
      {
        continue;
      }
      out << "f " << _faces[i][0] + 1 << "/" << _faces[i][0] + 1
          << " " << _faces[i][1] + 1 << "/" << _faces[i][1] + 1
          << " " << _faces[i][2] + 1 << "/" << _faces[i][2] + 1 << std::endl;
    }
    out.close();
  }
  void _Model::save_pamametrization_obj_file_(const std::string &in_file_name,
                                              const std::vector<std::pair<double, double>> &uvs) const
  {
    std::ofstream out(in_file_name);
    int folder_loc = in_file_name.rfind("\\") > in_file_name.rfind("/") ? in_file_name.rfind("\\") : in_file_name.rfind("/");
    int dot_loc = in_file_name.rfind('.');
    out << "g " << in_file_name.substr(folder_loc + 1, dot_loc - folder_loc - 1) << std::endl;
    out << "mtllib chessboard.mtl" << std::endl;
    out << "usemtl chessboard" << std::endl;
    for (int i = 0; i < number_vertices_(); ++i)
    {
      out << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << vertex_(i) << std::endl;
    }
    for (int i = 0; i < uvs.size(); ++i)
    {
      out << "vt " << std::setiosflags(std::ios::fixed) << std::setprecision(20) << uvs[i].first << " " << uvs[i].second
          << std::endl;
    }
    for (int i = 0; i < number_faces_(); ++i)
    {
      if (_faces_useless.find(i) != _faces_useless.end())
      {
        continue;
      }
      out << "f " << _faces[i][0] + 1 << "/" << _faces[i][0] + 1
          << " " << _faces[i][1] + 1 << "/" << _faces[i][1] + 1
          << " " << _faces[i][2] + 1 << "/" << _faces[i][2] + 1 << std::endl;
    }
    out.close();
  }
  void _Model::initialization_PQP_()
  {
    _pqp_model.BeginModel();
    PQP_REAL p1[3], p2[3], p3[3];
    for (auto f_it = face_begin(); f_it != face_end(); ++f_it)
    {
      auto fv_it = fv_begin(f_it.id());
      p1[0] = (*fv_it).x();
      p1[1] = (*fv_it).y();
      p1[2] = (*fv_it).z();
      ++fv_it;
      p2[0] = (*fv_it).x();
      p2[1] = (*fv_it).y();
      p2[2] = (*fv_it).z();
      ++fv_it;
      p3[0] = (*fv_it).x();
      p3[1] = (*fv_it).y();
      p3[2] = (*fv_it).z();
      ++fv_it;
      _pqp_model.AddTri(p1, p2, p3, f_it.id());
    }
    _pqp_model.EndModel();
  }
  std::tuple<_Point3, double, int> _Model::nearest_point_(const _Point3 &in_point)
  {
    _PQP_Query_Resutl query_res = proximity_query_(in_point);
    return std::make_tuple(query_res._nearest_point, query_res._distance, query_res._triangle_id);
  }
  double _Model::signed_distance_(const _Point3 &in_point)
  {
    _PQP_Query_Resutl query_res = proximity_query_(in_point);
    double dis1 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[0])).length_();
    double dis2 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[1])).length_();
    double dis3 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[2])).length_();
    _Point3 ave_nom(0, 0, 0);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[0]) * (dis2 + dis3);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[1]) * (dis1 + dis3);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[2]) * (dis1 + dis2);
    ave_nom.normalized_();
    _Point3 v = in_point - query_res._nearest_point;
    v.normalized_();
    if (_BOC::sign_(v.dot_(ave_nom)) == _BOC::_Sign::NegativE)
    {
      return -query_res._distance;
    }
    return query_res._distance;
  }
  double _Model::signed_distance_(const _Point3 &in_point, _Point3 &gradient)
  {
    _PQP_Query_Resutl query_res = proximity_query_(in_point);
    double dis1 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[0])).length_();
    double dis2 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[1])).length_();
    double dis3 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[2])).length_();
    _Point3 ave_nom(0, 0, 0);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[0]) * (dis2 + dis3);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[1]) * (dis1 + dis3);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[2]) * (dis1 + dis2);
    ave_nom.normalized_();
    _Point3 v = in_point - query_res._nearest_point;
    v.normalized_();
    if (_BOC::sign_(v.dot_(ave_nom)) == _BOC::_Sign::NegativE)
    {
      gradient = -v;
      return -query_res._distance;
    }
    gradient = v;
    return query_res._distance;
  }
  bool _Model::is_in_(const _Point3 &in_point)
  {
    _PQP_Query_Resutl query_res = proximity_query_(in_point);
    double dis1 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[0])).length_();
    double dis2 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[1])).length_();
    double dis3 = (query_res._nearest_point - vertex_(face_(query_res._triangle_id)[2])).length_();
    _Point3 ave_nom(0, 0, 0);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[0]) * (dis2 + dis3);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[1]) * (dis1 + dis3);
    ave_nom += normal_vertex_(face_(query_res._triangle_id)[2]) * (dis1 + dis2);
    ave_nom.normalized_();
    _Point3 v = in_point - query_res._nearest_point;
    v.normalized_();
    if (_BOC::sign_(v.dot_(ave_nom)) == _BOC::_Sign::NegativE)
    {
      return true;
    }
    return false;
  }
  void _Model::read_file_(const std::string &in_file_name)
  {
    int dot_loc = (int)in_file_name.rfind('.');
    if (dot_loc == -1)
    {
      throw "File name doesn't contain a dot! " + in_file_name;
    }
    std::string extension = in_file_name.substr(dot_loc + 1);
    if (extension == "obj")
    {
      read_obj_file_(in_file_name);
    }
    else if (extension == "off")
    {
      read_off_file_(in_file_name);
    }
    else
    {
      throw "This format can't be handled! " + extension;
    }
  }
  void _Model::read_obj_file_(const std::string &in_file_name)
  {
    std::ifstream in(in_file_name);
    if (in.fail())
    {
      throw "fail to read file: " + in_file_name;
    }
    _vertices.clear();
    _faces.clear();
    std::string line;
    while (std::getline(in, line))
    {
      std::istringstream sline(line);
      std::string word;
      sline >> word;
      if (word == "v")
      {
        double x, y, z;
        sline >> x >> y >> z;
        _vertices.push_back(_Point3(x, y, z));
      }
      else if (word == "f")
      {
        int id0, id1, idx;
        sline >> id0 >> id1;
        while (sline >> idx)
        {
          _Model::_MFace f(id0 - 1, id1 - 1, idx - 1, _vertices[id0 - 1], _vertices[id1 - 1], _vertices[idx - 1]);
          f.id = _faces.size();
          _faces.push_back(f);
          id1 = idx;
        }
      }
      else
      {
        continue;
      }
    }
    in.close();
  }

  void _Model::read_off_file_(const std::string &in_file_name)
  {
  }
  void _Model::compute_normal_boundingbox_()
  {
    if (_vertices.empty())
      return;
    _normals_vertex.clear();
    _normals_vertex.resize(_vertices.size(), _Point3(0, 0, 0));
    _normals_face.clear();
    _normals_face.resize(_faces.size(), _Point3(0, 0, 0));
    std::vector<bool> visited(_vertices.size(), false);
    for (int i = 0; i < _faces.size(); ++i)
    {
      _Point3 norm =
          (_vertices[_faces[i][1]] - _vertices[_faces[i][0]]).cross_((_vertices[_faces[i][2]] - _vertices[_faces[i][0]]));
      _normals_vertex[_faces[i][0]] += norm;
      _normals_vertex[_faces[i][1]] += norm;
      _normals_vertex[_faces[i][2]] += norm;
      visited[_faces[i][0]] = true;
      visited[_faces[i][1]] = true;
      visited[_faces[i][2]] = true;
      _normals_face[i] = norm.normalize_();
    }
    for (int i = 0; i < _normals_vertex.size(); ++i)
    {
      if (visited[i])
      {
        _normals_vertex[i].normalized_();
      }
    }
    _Point3 ptUp(_vertices[0]);
    _Point3 ptDown(_vertices[0]);
    for (int i = 1; i < _vertices.size(); ++i)
    {
      if (_vertices[i].x() > ptUp.x())
        ptUp[0] = _vertices[i].x();
      else if (_vertices[i].x() < ptDown.x())
        ptDown[0] = _vertices[i].x();
      if (_vertices[i].y() > ptUp.y())
        ptUp[1] = _vertices[i].y();
      else if (_vertices[i].y() < ptDown.y())
        ptDown[1] = _vertices[i].y();
      if (_vertices[i].z() > ptUp.z())
        ptUp[2] = _vertices[i].z();
      else if (_vertices[i].z() < ptDown.z())
        ptDown[2] = _vertices[i].z();
    }
    _bounding_box = std::make_pair(ptDown, ptUp);
  }
  _Model::_PQP_Query_Resutl _Model::proximity_query_(const _Point3 &in_point)
  {
    PQP_DistanceResult dres;
    dres.last_tri = _pqp_model.last_tri;
    PQP_REAL p[3];
    p[0] = in_point.x();
    p[1] = in_point.y();
    p[2] = in_point.z();
    PQP_Distance(&dres, &_pqp_model, p, 0.0, 0.0);
    _PQP_Query_Resutl res;
    res._pos_flag = dres.pos_flag;
    res._triangle_id = dres.last_tri->id;
    res._distance = dres.distance;
    res._nearest_point = _Point3(dres.p1[0], dres.p1[1], dres.p1[2]);
    return res;
  }
  _Model::_PQP_Query_Resutl::_PQP_Query_Resutl()
  {
  }
  _Model::_PQP_Query_Resutl::_PQP_Query_Resutl(const int &in_pos_flag,
                                               const int &in_triangle_id,
                                               const double &in_distance,
                                               const _Point3 &in_nearset_point)
      : _pos_flag(in_pos_flag), _triangle_id(in_triangle_id), _distance(in_distance), _nearest_point(in_nearset_point)
  {
  }
} // namespace BGAL
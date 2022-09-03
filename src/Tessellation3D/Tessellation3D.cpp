#include "BGAL/Tessellation3D/Side3D.h"
#include "BGAL/Tessellation3D/Tessellation3D.h"
namespace BGAL
{
  _Tessellation3D_Skeleton::_Tessellation3D_Skeleton()
  {
    _neights.resize(0);
  }
  _Tessellation3D_Skeleton::_Tessellation3D_Skeleton(Rt &rt, const int &num_vertices)
  {
    _neights.resize(num_vertices);
    for (auto c_it = rt.finite_cells_begin(); c_it != rt.finite_cells_end(); c_it++)
    {
      int i0, i1, i2, i3;
      i0 = c_it->vertex(0)->info();
      i1 = c_it->vertex(1)->info();
      i2 = c_it->vertex(2)->info();
      i3 = c_it->vertex(3)->info();
      if (i0 != -1)
      {
        if (i1 != -1)
        {
          _neights[i0].insert(i1);
          _neights[i1].insert(i0);
        }
        if (i2 != -1)
        {
          _neights[i0].insert(i2);
          _neights[i2].insert(i0);
        }
        if (i3 != -1)
        {
          _neights[i0].insert(i3);
          _neights[i3].insert(i0);
        }
      }
      if (i1 != -1)
      {
        if (i2 != -1)
        {
          _neights[i1].insert(i2);
          _neights[i2].insert(i1);
        }
        if (i3 != -1)
        {
          _neights[i3].insert(i1);
          _neights[i1].insert(i3);
        }
      }
      if (i2 != -1 && i3 != -1)
      {
        _neights[i3].insert(i2);
        _neights[i2].insert(i3);
      }
    }
  }

  _BOC::_Sign _Restricted_Tessellation3D::side_(const int &ip1, const int &ip2, _Symbolic_Point &v)
  {
    const _Point3 &p1 = _sites[ip1];
    const _Point3 &p2 = _sites[ip2];
    _BOC::_Sign res;
    switch (v.flag)
    {
    case 1:
    {
      _Point3 q = _model.vertex_(v.p[0]);
      res = _Side3D::side1_(p1.x(),
                            p1.y(),
                            p1.z(),
                            p2.x(),
                            p2.y(),
                            p2.z(),
                            _weights[ip1],
                            _weights[ip2],
                            q.x(),
                            q.y(),
                            q.z());
      break;
    }
    case 2:
    {
      _Point3 q1 = _model.vertex_(v.p[0]);
      _Point3 q2 = _model.vertex_(v.p[1]);
      const _Point3 &p3 = _sites[v.p[2]];
      res = _Side3D::side2_(p1.x(),
                            p1.y(),
                            p1.z(),
                            p2.x(),
                            p2.y(),
                            p2.z(),
                            p3.x(),
                            p3.y(),
                            p3.z(),
                            _weights[ip1],
                            _weights[ip2],
                            _weights[v.p[2]],
                            q1.x(),
                            q1.y(),
                            q1.z(),
                            q2.x(),
                            q2.y(),
                            q2.z());
      break;
    }
    case 3:
    {
      _Point3 q1 = _model.vertex_(v.p[0]);
      _Point3 q2 = _model.vertex_(v.p[1]);
      _Point3 q3 = _model.vertex_(v.p[2]);
      std::vector<int> ip34;
      for (auto it = v._sym.begin(); it != v._sym.end(); ++it)
      {
        if (*it > 0)
        {
          ip34.push_back(*it);
        }
      }
      const int &ip3 = ip34[0] - 1;
      const int &ip4 = ip34[1] - 1;
      const _Point3 &p3 = _sites[ip3];
      const _Point3 &p4 = _sites[ip4];
      res = _Side3D::side3_(p1.x(),
                            p1.y(),
                            p1.z(),
                            p2.x(),
                            p2.y(),
                            p2.z(),
                            p3.x(),
                            p3.y(),
                            p3.z(),
                            p4.x(),
                            p4.y(),
                            p4.z(),
                            _weights[ip1],
                            _weights[ip2],
                            _weights[ip3],
                            _weights[ip4],
                            q1.x(),
                            q1.y(),
                            q1.z(),
                            q2.x(),
                            q2.y(),
                            q2.z(),
                            q3.x(),
                            q3.y(),
                            q3.z());
      break;
    }
    default:
      throw std::runtime_error("flag error");
      break;
    }
    return res;
  }

  _Restricted_Tessellation3D::_Symbolic_Point _Restricted_Tessellation3D::insec_bisector_(const _Symbolic_Point &p1,
                                                                                          const _Symbolic_Point &p2,
                                                                                          const int &neigh,
                                                                                          const int &center,
                                                                                          const std::vector<int> &current_face_vid) const
  {
    std::set<int> insec = p1.insec_(p2);
    insec.insert(neigh + 1);
    _Symbolic_Point insec_p(insec, center);
    if (insec_p.num_sites_() == 1)
    {
      if (p1.num_sites_() == 0 && p2.num_sites_() == 0)
      {
        insec_p.flag = 2;
        insec_p.p.clear();
        insec_p.p.push_back(p1.p[0]);
        insec_p.p.push_back(p2.p[0]);
        insec_p.p.push_back(neigh);
      }
      else if (p2.num_sites_() == 1)
      {
        insec_p.flag = 2;
        insec_p.p.clear();
        insec_p.p.push_back(p2.p[0]);
        insec_p.p.push_back(p2.p[1]);
        insec_p.p.push_back(neigh);
      }
      else
      {
        insec_p.flag = 2;
        insec_p.p.clear();
        insec_p.p.push_back(p1.p[0]);
        insec_p.p.push_back(p1.p[1]);
        insec_p.p.push_back(neigh);
      }
    }
    else
    {
      insec_p.flag = 3;
      insec_p.p.clear();
      insec_p.p.push_back(current_face_vid[0]);
      insec_p.p.push_back(current_face_vid[1]);
      insec_p.p.push_back(current_face_vid[2]);
    }
    //std::cout << "insec:  ";
    //for (auto tit = insec_p._sym.begin(); tit != insec_p._sym.end(); ++tit)
    //{
    //std::cout << *tit << " ";
    //}
    //std::cout << std::endl;
    return insec_p;
  }

  void _Restricted_Tessellation3D::calculate_()
  {
    std::vector<std::pair<Weighted_point, int>> wps(_num_sites);
    double min_weight = *(std::min_element(_weights.begin(), _weights.end()));
    for (int i = 0; i < _num_sites; ++i)
    {
      wps[i] = std::make_pair(Weighted_point(Point(_sites[i].x(), _sites[i].y(), _sites[i].z()), _weights[i]), i);
    }
    std::pair<_Point3, _Point3> bbox = _model.bounding_box_();
    const _Point3 &min_p = bbox.first;
    const _Point3 &max_p = bbox.second;
    wps.push_back(std::make_pair(Weighted_point(Point(3 * min_p.x() - 2 * max_p.x(),
                                                      3 * min_p.y() - 2 * max_p.y(),
                                                      3 * min_p.z() - 2 * max_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * max_p.x() - 2 * min_p.x(),
                                                      3 * min_p.y() - 2 * max_p.y(),
                                                      3 * min_p.z() - 2 * max_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * max_p.x() - 2 * min_p.x(),
                                                      3 * max_p.y() - 2 * min_p.y(),
                                                      3 * min_p.z() - 2 * max_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * min_p.x() - 2 * max_p.x(),
                                                      3 * max_p.y() - 2 * min_p.y(),
                                                      3 * min_p.z() - 2 * max_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * min_p.x() - 2 * max_p.x(),
                                                      3 * min_p.y() - 2 * max_p.y(),
                                                      3 * max_p.z() - 2 * min_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * max_p.x() - 2 * min_p.x(),
                                                      3 * min_p.y() - 2 * max_p.y(),
                                                      3 * max_p.z() - 2 * min_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * max_p.x() - 2 * min_p.x(),
                                                      3 * max_p.y() - 2 * min_p.y(),
                                                      3 * max_p.z() - 2 * min_p.z()),
                                                min_weight),
                                 -1));
    wps.push_back(std::make_pair(Weighted_point(Point(3 * min_p.x() - 2 * max_p.x(),
                                                      3 * max_p.y() - 2 * min_p.y(),
                                                      3 * max_p.z() - 2 * min_p.z()),
                                                min_weight),
                                 -1));
    Rt rt(wps.begin(), wps.end());
    _skeleton = _Tessellation3D_Skeleton(rt, _num_sites);
    std::vector<bool> face_is_visited(_model.number_faces_(), false);
    std::vector<std::map<int, int>> from_idx_to_locations(_num_sites);
    std::map<int, int> from_locations_to_idx;
    std::vector<std::vector<_Symbolic_Point>> cliped_faces;
    std::queue<std::pair<int, int>> Qfs;
    for (_Face_Iterator f_it = _model.face_begin(); f_it != _model.face_end(); ++f_it)
    {
      if (face_is_visited[f_it.id()])
        continue;
      {
        _Point3 fcenter(0, 0, 0);
        std::vector<int> v_idx;
        std::vector<int> f_idx;
        f_idx.push_back(f_it.id());
        for (auto fe_it = _model.fe_begin(f_it.id()); fe_it != _model.fe_end(f_it.id()); ++fe_it)
        {
          int cv = (*fe_it)._id_left_vertex;
          v_idx.push_back(cv);
          f_idx.push_back(_model.edge_((*fe_it)._id_reverse_edge)._id_face);
          _Point3 cp = _model.vertex_(cv);
          fcenter = fcenter + cp;
        }
        fcenter = fcenter / 3.0;
        int near_site = rt.nearest_power_vertex(K::Point_3(fcenter.x(), fcenter.y(), fcenter.z()))->info();
        std::vector<_Symbolic_Point> fps;
        _Symbolic_Point fp1(-f_idx[0] - 1, -f_idx[1] - 1, -f_idx[3] - 1, near_site);
        fp1.flag = 1;
        fp1.p.push_back(v_idx[0]);
        _Symbolic_Point fp2(-f_idx[0] - 1, -f_idx[1] - 1, -f_idx[2] - 1, near_site);
        fp2.flag = 1;
        fp2.p.push_back(v_idx[1]);
        _Symbolic_Point fp3(-f_idx[0] - 1, -f_idx[2] - 1, -f_idx[3] - 1, near_site);
        fp3.flag = 1;
        fp3.p.push_back(v_idx[2]);
        fps.push_back(fp1);
        fps.push_back(fp2);
        fps.push_back(fp3);
        from_idx_to_locations[near_site][f_idx[0]] = cliped_faces.size();
        from_locations_to_idx[cliped_faces.size()] = f_idx[0];
        cliped_faces.push_back(fps);
        Qfs.push(std::make_pair(near_site, cliped_faces.size() - 1));
      }
      while (!Qfs.empty())
      {
        const int current_site = Qfs.front().first;
        const int current_cliped = Qfs.front().second;
        //std::cout << "current site  " << current_site << std::endl;
        //std::cout << "current cliped  " << current_cliped << std::endl;
        Qfs.pop();
        const int current_face = from_locations_to_idx[current_cliped];
        //std::cout << "current face  " << current_face << std::endl;
        std::vector<_Symbolic_Point> old_cliped = cliped_faces[current_cliped];
        const std::vector<_Symbolic_Point> original_cliped = cliped_faces[current_cliped];
        std::vector<int> current_face_vid;
        current_face_vid.push_back(original_cliped[0].p[0]);
        current_face_vid.push_back(original_cliped[1].p[0]);
        current_face_vid.push_back(original_cliped[2].p[0]);
        std::vector<_Symbolic_Point> update_cliped(0);
        const std::set<int> &planes = _skeleton.neight_(current_site);
        for (auto p = planes.begin(); p != planes.end(); ++p)
        {
          if (old_cliped.size() == 0)
            break;
          //std::cout << "old size  " << old_cliped.size() << std::endl;
          //std::cout << "adj site  " << *p << std::endl;

          _BOC::_Sign pre_state = side_(current_site, *p, old_cliped[0]);
          _BOC::_Sign cur_state = side_(current_site, *p, old_cliped[1]);
          //std::cout << "old cliped[0]:  ";
          //for (auto tit = old_cliped[0]._sym.begin(); tit != old_cliped[0]._sym.end(); ++tit)
          //{
          //  //std::cout << *tit << " ";
          //}
          //std::cout << std::endl << "old cliped[1]:  ";
          //for (auto tit = old_cliped[1]._sym.begin(); tit != old_cliped[1]._sym.end(); ++tit)
          //{
          //std::cout << *tit << " ";
          //}
          //std::cout << std::endl;
          //std::cout << "pre cur: " << (pre_state == _BOC::_Sign::ZerO ? 0 : (pre_state == _BOC::_Sign::PositivE ? 1 : -1)) << " " << (cur_state == _BOC::_Sign::ZerO ? 0 : (cur_state == _BOC::_Sign::PositivE ? 1 : -1)) << std::endl;
          _BOC::_Sign nex_state;
          if (pre_state == _BOC::_Sign::PositivE)
          {
            update_cliped.push_back(old_cliped[0]);
            if (cur_state == _BOC::_Sign::PositivE)
            {
              update_cliped.push_back(old_cliped[1]);
            }
            else if (cur_state == _BOC::_Sign::NegativE)
            {
              update_cliped.push_back(insec_bisector_(old_cliped[0], old_cliped[1], *p, current_site, current_face_vid));
            }
            else if (cur_state == _BOC::_Sign::ZerO)
            {
              _Symbolic_Point updateSym = old_cliped[1];
              updateSym.update_(*p, old_cliped[0]);
              update_cliped.push_back(updateSym);
            }
          }
          else if (pre_state == _BOC::_Sign::NegativE)
          {
            if (cur_state == _BOC::_Sign::PositivE)
            {
              update_cliped.push_back(insec_bisector_(old_cliped[0], old_cliped[1], *p, current_site, current_face_vid));
              update_cliped.push_back(old_cliped[1]);
            }
          }
          else if (pre_state == _BOC::_Sign::ZerO)
          {
            if (cur_state == _BOC::_Sign::PositivE)
            {
              _Symbolic_Point updateSym = old_cliped[0];
              updateSym.update_(*p, old_cliped[1]);
              update_cliped.push_back(updateSym);
              update_cliped.push_back(old_cliped[1]);
            }
            else if (cur_state == _BOC::_Sign::ZerO)
            {
              _BOC::_Sign sp2 = side_(current_site, *p, old_cliped[2]);
              if (sp2 == _BOC::_Sign::PositivE)
              {
                update_cliped.clear();
                continue;
              }
              else if (sp2 == _BOC::_Sign::NegativE)
              {
                old_cliped.clear();
                break;
              }
            }
          }
          //std::cout << "update size:  " << update_cliped.size() << std::endl;
          _BOC::_Sign sp0 = pre_state;
          _BOC::_Sign sp1 = cur_state;
          bool ifbreak = false, ifcontinue = false;
          for (int i = 2; i < old_cliped.size(); ++i)
          {
            //std::cout << "old cliped[i]:  ";
            //for (auto tit = old_cliped[i]._sym.begin(); tit != old_cliped[i]._sym.end(); ++tit)
            //{
            //std::cout << *tit << " ";
            //}
            //std::cout << std::endl;
            nex_state = side_(current_site, *p, old_cliped[i]);
            //std::cout << "nex cur: " << (nex_state == _BOC::_Sign::ZerO ? 0 : (nex_state == _BOC::_Sign::PositivE ? 1 : -1)) << std::endl;
            if (cur_state == _BOC::_Sign::PositivE)
            {
              if (nex_state == _BOC::_Sign::PositivE)
              {
                update_cliped.push_back(old_cliped[i]);
              }
              else if (nex_state == _BOC::_Sign::NegativE)
              {
                update_cliped.push_back(insec_bisector_(old_cliped[i - 1],
                                                        old_cliped[i],
                                                        *p,
                                                        current_site,
                                                        current_face_vid));
              }
              else
              {
                _Symbolic_Point updatesym = old_cliped[i];
                updatesym.update_(*p, old_cliped[i - 1]);
                update_cliped.push_back(updatesym);
              }
            }
            else if (cur_state == _BOC::_Sign::NegativE)
            {
              if (nex_state == _BOC::_Sign::PositivE)
              {
                update_cliped.push_back(insec_bisector_(old_cliped[i - 1],
                                                        old_cliped[i],
                                                        *p,
                                                        current_site,
                                                        current_face_vid));
                update_cliped.push_back(old_cliped[i]);
              }
            }
            else if (cur_state == _BOC::_Sign::ZerO)
            {
              if (nex_state == _BOC::_Sign::PositivE)
              {
                if (pre_state == _BOC::_Sign::NegativE)
                {
                  _Symbolic_Point updateSym = old_cliped[i - 1];
                  updateSym.update_(*p, old_cliped[i]);
                  update_cliped.push_back(updateSym);
                  update_cliped.push_back(old_cliped[i]);
                }
                else if (pre_state == _BOC::_Sign::PositivE)
                {
                  update_cliped.clear();
                  ifcontinue = true;
                  break;
                }
                else
                {
                  throw std::runtime_error("Two consecutive 0s appear in front.");
                }
              }
              else if (nex_state == _BOC::_Sign::ZerO)
              {
                if (pre_state == _BOC::_Sign::PositivE)
                {
                  update_cliped.clear();
                  ifcontinue = true;
                  break;
                }
                else if (pre_state == _BOC::_Sign::NegativE)
                {
                  old_cliped.clear();
                  ifbreak = true;
                  break;
                }
              }
              else if (nex_state == _BOC::_Sign::NegativE)
              {
                if (pre_state == _BOC::_Sign::NegativE)
                {
                  old_cliped.clear();
                  ifbreak = true;
                  break;
                }
              }
            }
            pre_state = cur_state;
            cur_state = nex_state;
            //std::cout << "update size:  " << update_cliped.size() << std::endl;
          }
          if (ifbreak)
            break;
          if (ifcontinue)
            continue;
          if (sp0 == _BOC::_Sign::PositivE)
          {
            if (cur_state == _BOC::_Sign::NegativE)
            {
              update_cliped.push_back(insec_bisector_(old_cliped[old_cliped.size() - 1],
                                                      old_cliped[0],
                                                      *p,
                                                      current_site,
                                                      current_face_vid));
            }
            else if (cur_state == _BOC::_Sign::ZerO)
            {
              if (pre_state == _BOC::_Sign::NegativE)
              {
                _Symbolic_Point updateSym = old_cliped[old_cliped.size() - 1];
                updateSym.update_(*p, old_cliped[0]);
                update_cliped.push_back(updateSym);
              }
              else if (pre_state == _BOC::_Sign::PositivE)
              {
                update_cliped.clear();
                continue;
              }
            }
          }
          else if (sp0 == _BOC::_Sign::NegativE)
          {
            if (cur_state == _BOC::_Sign::PositivE)
            {
              update_cliped.push_back(insec_bisector_(old_cliped[old_cliped.size() - 1],
                                                      old_cliped[0],
                                                      *p,
                                                      current_site,
                                                      current_face_vid));
            }
            else if (cur_state == _BOC::_Sign::ZerO)
            {
              if (pre_state == _BOC::_Sign::NegativE)
              {
                old_cliped.clear();
                break;
              }
            }
          }
          else if (sp0 == _BOC::_Sign::ZerO)
          {
            if (sp1 == _BOC::_Sign::NegativE)
            {
              if (cur_state == _BOC::_Sign::PositivE)
              {
                _Symbolic_Point updateSym = old_cliped[0];
                updateSym.update_(*p, old_cliped[old_cliped.size() - 1]);
                update_cliped.push_back(updateSym);
              }
              else
              {
                old_cliped.clear();
                break;
              }
            }
            else if (sp1 == _BOC::_Sign::PositivE)
            {
              if (cur_state == _BOC::_Sign::PositivE)
              {
                update_cliped.clear();
                continue;
              }
              else if (cur_state == _BOC::_Sign::ZerO)
              {
                update_cliped.clear();
                continue;
              }
            }
          }
          old_cliped = update_cliped;
          update_cliped.clear();
        }
        for (int i = 0; i < old_cliped.size(); ++i)
        {
          for (auto it = old_cliped[i]._sym.begin(); it != old_cliped[i]._sym.end(); ++it)
          {
            if (*it > 0 && from_idx_to_locations[(*it) - 1].find(current_face) == from_idx_to_locations[(*it) - 1].end())
            {
              from_idx_to_locations[(*it) - 1][current_face] = cliped_faces.size();
              from_locations_to_idx[cliped_faces.size()] = current_face;
              std::vector<_Symbolic_Point> copy_cliped = original_cliped;
              copy_cliped[0]._site = (*it) - 1;
              copy_cliped[1]._site = (*it) - 1;
              copy_cliped[2]._site = (*it) - 1;
              cliped_faces.push_back(copy_cliped);
              Qfs.push(std::make_pair((*it) - 1, cliped_faces.size() - 1));
            }
            if (*it < 0 && from_idx_to_locations[current_site].find(-(*it) - 1) == from_idx_to_locations[current_site].end())
            {
              std::vector<int> v_idx;
              std::vector<int> f_idx;
              f_idx.push_back(-(*it) - 1);
              for (auto afe_it = _model.fe_begin(-(*it) - 1); afe_it != _model.fe_end(-(*it) - 1); ++afe_it)
              {
                int cv = (*afe_it)._id_left_vertex;
                v_idx.push_back(cv);
                f_idx.push_back(_model.edge_((*afe_it)._id_reverse_edge)._id_face);
              }
              std::vector<_Symbolic_Point> fps;
              _Symbolic_Point fp1(-f_idx[0] - 1, -f_idx[1] - 1, -f_idx[3] - 1, current_site);
              fp1.flag = 1;
              fp1.p.push_back(v_idx[0]);
              _Symbolic_Point fp2(-f_idx[0] - 1, -f_idx[1] - 1, -f_idx[2] - 1, current_site);
              fp2.flag = 1;
              fp2.p.push_back(v_idx[1]);
              _Symbolic_Point fp3(-f_idx[0] - 1, -f_idx[2] - 1, -f_idx[3] - 1, current_site);
              fp3.flag = 1;
              fp3.p.push_back(v_idx[2]);
              fps.push_back(fp1);
              fps.push_back(fp2);
              fps.push_back(fp3);
              from_idx_to_locations[current_site][-(*it) - 1] = cliped_faces.size();
              from_locations_to_idx[cliped_faces.size()] = -(*it) - 1;
              cliped_faces.push_back(fps);
              Qfs.push(std::make_pair(current_site, cliped_faces.size() - 1));
            }
          }
        }
        cliped_faces[current_cliped] = old_cliped;
        face_is_visited[current_face] = true;
      }
    }
    _vertices.clear();
    _cells.clear();
    _cells.resize(_num_sites);
    _edges.clear();
    _edges.resize(_num_sites);
    std::map<std::set<int>, int> from_Sym_to_vertex;

    for (int i = 0; i < _num_sites; ++i)
    {
      for (auto it = from_idx_to_locations[i].begin(); it != from_idx_to_locations[i].end(); ++it)
      {
        if (cliped_faces[it->second].size() == 0)
          continue;
        if (cliped_faces[it->second].size() < 3)
          throw std::runtime_error("size error");
        int firest_p = -1;
        int second_p = -1;
        if (cliped_faces[it->second][0].flag == 1)
        {
          std::set<int> psym;
          psym.insert(cliped_faces[it->second][0].p[0]);
          if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
          {
            from_Sym_to_vertex[psym] = _vertices.size();
            firest_p = _vertices.size();
            _vertices.push_back(_model.vertex_(cliped_faces[it->second][0].p[0]));
          }
          else
          {
            firest_p = from_Sym_to_vertex[psym];
          }
        }
        else if (cliped_faces[it->second][0].flag == 2)
        {
          std::set<int> psym;
          psym.insert(-cliped_faces[it->second][0].p[0] - 1);
          psym.insert(-cliped_faces[it->second][0].p[1] - 1);
          psym.insert(cliped_faces[it->second][0].p[2] + 1);
          psym.insert(i + 1);
          if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
          {
            from_Sym_to_vertex[psym] = _vertices.size();
            firest_p = _vertices.size();
            int q1 = i;
            int q2 = cliped_faces[it->second][0].p[2];
            _Point3 p1 = _model.vertex_(cliped_faces[it->second][0].p[0]);
            _Point3 p2 = _model.vertex_(cliped_faces[it->second][0].p[1]);
            _Point3 v = (_sites[q2] - _sites[q1]) * 2;
            double d = _weights[q2] - _weights[q1] - 0.5 * v.dot_(_sites[q2] + _sites[q1]);
            double d1 = fabs(p1.dot_(v) + d);
            double d2 = fabs(p2.dot_(v) + d);
            _Point3 cp3 = p1 + (p2 - p1) * d1 / (d1 + d2);
            _vertices.push_back(cp3);
          }
          else
          {
            firest_p = from_Sym_to_vertex[psym];
          }
        }
        else if (cliped_faces[it->second][0].flag == 3)
        {
          std::set<int> psym;
          psym.insert(-cliped_faces[it->second][0].p[0] - 1);
          psym.insert(-cliped_faces[it->second][0].p[1] - 1);
          psym.insert(-cliped_faces[it->second][0].p[2] - 1);
          psym.insert(i + 1);
          std::vector<int> p12;
          int q = -1;
          for (auto tit = cliped_faces[it->second][0]._sym.begin(); tit != cliped_faces[it->second][0]._sym.end();
               ++tit)
          {
            if (*tit > 0)
            {
              psym.insert(*tit);
              p12.push_back(*tit - 1);
            }
            else
              q = -*tit - 1;
          }
          if (psym.size() != 6)
          {
            throw std::runtime_error("psym.size() != 6");
          }
          if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
          {
            _Point3 v1 = (_sites[p12[0]] - _sites[i]) * 2;
            double d1 = _weights[p12[0]] - _weights[i] - 0.5 * v1.dot_(_sites[p12[0]] + _sites[i]);
            _Point3 v2 = (_sites[p12[1]] - _sites[i]) * 2;
            double d2 = _weights[p12[1]] - _weights[i] - 0.5 * v2.dot_(_sites[p12[1]] + _sites[i]);
            _Point3 v0 = _model.normal_face_(q);
            v0.normalized_();
            double d0 = -v0.dot_(_model.vertex_(cliped_faces[it->second][0].p[0]));
            _Point3 cp3 = _Point3::intersection_three_plane(v0, d0, v1, d1, v2, d2);
            from_Sym_to_vertex[psym] = _vertices.size();
            firest_p = _vertices.size();
            _vertices.push_back(cp3);
          }
          else
          {
            firest_p = from_Sym_to_vertex[psym];
          }
        }
        else
        {
          throw std::runtime_error("flag error");
        }
        if (cliped_faces[it->second][1].flag == 1)
        {
          std::set<int> psym;
          psym.insert(cliped_faces[it->second][1].p[0]);
          if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
          {
            from_Sym_to_vertex[psym] = _vertices.size();
            second_p = _vertices.size();
            _vertices.push_back(_model.vertex_(cliped_faces[it->second][1].p[0]));
          }
          else
          {
            second_p = from_Sym_to_vertex[psym];
          }
        }
        else if (cliped_faces[it->second][1].flag == 2)
        {
          std::set<int> psym;
          psym.insert(-cliped_faces[it->second][1].p[0] - 1);
          psym.insert(-cliped_faces[it->second][1].p[1] - 1);
          psym.insert(cliped_faces[it->second][1].p[2] + 1);
          psym.insert(i + 1);
          if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
          {
            from_Sym_to_vertex[psym] = _vertices.size();
            second_p = _vertices.size();
            int q1 = i;
            int q2 = cliped_faces[it->second][1].p[2];
            _Point3 p1 = _model.vertex_(cliped_faces[it->second][1].p[0]);
            _Point3 p2 = _model.vertex_(cliped_faces[it->second][1].p[1]);
            _Point3 v = (_sites[q2] - _sites[q1]) * 2;
            double d = _weights[q2] - _weights[q1] - 0.5 * v.dot_(_sites[q2] + _sites[q1]);
            double d1 = fabs(p1.dot_(v) + d);
            double d2 = fabs(p2.dot_(v) + d);
            _Point3 cp3 = p1 + (p2 - p1) * d1 / (d1 + d2);
            _vertices.push_back(cp3);
          }
          else
          {
            second_p = from_Sym_to_vertex[psym];
          }
        }
        else if (cliped_faces[it->second][1].flag == 3)
        {
          std::set<int> psym;
          psym.insert(-cliped_faces[it->second][1].p[0] - 1);
          psym.insert(-cliped_faces[it->second][1].p[1] - 1);
          psym.insert(-cliped_faces[it->second][1].p[2] - 1);
          psym.insert(i + 1);
          std::vector<int> p12;
          int q = -1;
          for (auto tit = cliped_faces[it->second][1]._sym.begin(); tit != cliped_faces[it->second][1]._sym.end();
               ++tit)
          {
            if (*tit > 0)
            {
              psym.insert(*tit);
              p12.push_back(*tit - 1);
            }
            else
              q = -*tit - 1;
          }
          if (psym.size() != 6)
          {
            throw std::runtime_error("psym.size() != 6");
          }
          if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
          {
            _Point3 v1 = (_sites[p12[0]] - _sites[i]) * 2;
            double d1 = _weights[p12[0]] - _weights[i] - 0.5 * v1.dot_(_sites[p12[0]] + _sites[i]);
            _Point3 v2 = (_sites[p12[1]] - _sites[i]) * 2;
            double d2 = _weights[p12[1]] - _weights[i] - 0.5 * v2.dot_(_sites[p12[1]] + _sites[i]);
            _Point3 v0 = _model.normal_face_(q);
            v0.normalized_();
            double d0 = -v0.dot_(_model.vertex_(cliped_faces[it->second][1].p[0]));
            _Point3 cp3 = _Point3::intersection_three_plane(v0, d0, v1, d1, v2, d2);
            from_Sym_to_vertex[psym] = _vertices.size();
            second_p = _vertices.size();
            _vertices.push_back(cp3);
          }
          else
          {
            second_p = from_Sym_to_vertex[psym];
          }
        }
        else
        {
          throw std::runtime_error("flag error");
        }
        {
          std::set<int> insec_sym = cliped_faces[it->second][0].insec_(cliped_faces[it->second][1]);
          int adj_sites = -1;
          for (auto tit = insec_sym.begin(); tit != insec_sym.end(); ++tit)
          {
            if (*tit > 0)
              adj_sites = *tit - 1;
          }
          if (adj_sites != -1)
          {
            if (_edges[i].find(adj_sites) == _edges[i].end())
            {
              std::vector<std::pair<int, int>> tedges;
              tedges.push_back(std::make_pair(firest_p, second_p));
              _edges[i][adj_sites] = tedges;
            }
            else
            {
              _edges[i][adj_sites].push_back(std::make_pair(firest_p, second_p));
            }
            //_edges[i].push_back(std::make_tuple(firest_p, second_p, adj_sites));
          }
        }

        for (int j = 2; j < cliped_faces[it->second].size(); ++j)
        {
          int third_p = -1;
          if (cliped_faces[it->second][j].flag == 1)
          {
            std::set<int> psym;
            psym.insert(cliped_faces[it->second][j].p[0]);
            if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
            {
              from_Sym_to_vertex[psym] = _vertices.size();
              third_p = _vertices.size();
              _vertices.push_back(_model.vertex_(cliped_faces[it->second][j].p[0]));
            }
            else
            {
              third_p = from_Sym_to_vertex[psym];
            }
          }
          else if (cliped_faces[it->second][j].flag == 2)
          {
            std::set<int> psym;
            psym.insert(-cliped_faces[it->second][j].p[0] - 1);
            psym.insert(-cliped_faces[it->second][j].p[1] - 1);
            psym.insert(cliped_faces[it->second][j].p[2] + 1);
            psym.insert(i + 1);
            if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
            {
              from_Sym_to_vertex[psym] = _vertices.size();
              third_p = _vertices.size();
              int q1 = i;
              int q2 = cliped_faces[it->second][j].p[2];
              _Point3 p1 = _model.vertex_(cliped_faces[it->second][j].p[0]);
              _Point3 p2 = _model.vertex_(cliped_faces[it->second][j].p[1]);
              _Point3 v = (_sites[q2] - _sites[q1]) * 2;
              double d = _weights[q2] - _weights[q1] - 0.5 * v.dot_(_sites[q2] + _sites[q1]);
              double d1 = fabs(p1.dot_(v) + d);
              double d2 = fabs(p2.dot_(v) + d);
              _Point3 cp3 = p1 + (p2 - p1) * d1 / (d1 + d2);
              _vertices.push_back(cp3);
            }
            else
            {
              third_p = from_Sym_to_vertex[psym];
            }
          }
          else if (cliped_faces[it->second][j].flag == 3)
          {
            std::set<int> psym;
            psym.insert(-cliped_faces[it->second][j].p[0] - 1);
            psym.insert(-cliped_faces[it->second][j].p[1] - 1);
            psym.insert(-cliped_faces[it->second][j].p[2] - 1);
            psym.insert(i + 1);
            std::vector<int> p12;
            int q = -1;
            for (auto tit = cliped_faces[it->second][j]._sym.begin(); tit != cliped_faces[it->second][j]._sym.end();
                 ++tit)
            {
              if (*tit > 0)
              {
                psym.insert(*tit);
                p12.push_back(*tit - 1);
              }
              else
                q = -*tit - 1;
            }
            if (psym.size() != 6)
            {
              throw std::runtime_error("psym.size() != 6");
            }
            if (from_Sym_to_vertex.find(psym) == from_Sym_to_vertex.end())
            {
              _Point3 v1 = (_sites[p12[0]] - _sites[i]) * 2;
              double d1 = _weights[p12[0]] - _weights[i] - 0.5 * v1.dot_(_sites[p12[0]] + _sites[i]);
              _Point3 v2 = (_sites[p12[1]] - _sites[i]) * 2;
              double d2 = _weights[p12[1]] - _weights[i] - 0.5 * v2.dot_(_sites[p12[1]] + _sites[i]);
              _Point3 v0 = _model.normal_face_(q);
              v0.normalized_();
              double d0 = -v0.dot_(_model.vertex_(cliped_faces[it->second][j].p[0]));
              _Point3 cp3 = _Point3::intersection_three_plane(v0, d0, v1, d1, v2, d2);
              from_Sym_to_vertex[psym] = _vertices.size();
              third_p = _vertices.size();
              _vertices.push_back(cp3);
            }
            else
            {
              third_p = from_Sym_to_vertex[psym];
            }
          }
          else
          {
            throw std::runtime_error("flag error");
          }
          _cells[i].push_back(std::make_tuple(firest_p, second_p, third_p));
          {
            std::set<int> insec_sym = cliped_faces[it->second][j - 1].insec_(cliped_faces[it->second][j]);
            int adj_sites = -1;
            for (auto tit = insec_sym.begin(); tit != insec_sym.end(); ++tit)
            {
              if (*tit > 0)
                adj_sites = *tit - 1;
            }
            if (adj_sites != -1)
            {
              if (_edges[i].find(adj_sites) == _edges[i].end())
              {
                std::vector<std::pair<int, int>> tedges;
                tedges.push_back(std::make_pair(second_p, third_p));
                _edges[i][adj_sites] = tedges;
              }
              else
              {
                _edges[i][adj_sites].push_back(std::make_pair(second_p, third_p));
              }
              //_edges[i].push_back(std::make_tuple(second_p, third_p, adj_sites));
            }
          }
          second_p = third_p;
        }
        {
          std::set<int> insec_sym = cliped_faces[it->second][0].insec_(cliped_faces[it->second].back());
          int adj_sites = -1;
          for (auto tit = insec_sym.begin(); tit != insec_sym.end(); ++tit)
          {
            if (*tit > 0)
              adj_sites = *tit - 1;
          }
          if (adj_sites != -1)
          {
            if (_edges[i].find(adj_sites) == _edges[i].end())
            {
              std::vector<std::pair<int, int>> tedges;
              tedges.push_back(std::make_pair(second_p, firest_p));
              _edges[i][adj_sites] = tedges;
            }
            else
            {
              _edges[i][adj_sites].push_back(std::make_pair(second_p, firest_p));
            }
            //_edges[i].push_back(std::make_tuple(firest_p, second_p, adj_sites));
          }
        }
      }
    }
  }
  _Restricted_Tessellation3D::_Restricted_Tessellation3D(const _ManifoldModel& in_model)
      :_model(in_model)
  {
      _num_sites = 0;
      _sites.clear();
      _weights.clear();
      _vertices.clear();
      _cells.clear();
  }
  _Restricted_Tessellation3D::_Restricted_Tessellation3D(const _ManifoldModel &in_model,
                                                         const std::vector<_Point3> &in_sites,
                                                         const std::vector<double> &in_weights)
      : _num_sites(in_sites.size()),
        _sites(in_sites),
        _weights(in_weights),
        _model(in_model)
  {
    _vertices.clear();
    _cells.clear();
    _cells.resize(_num_sites);
    calculate_();
  }
  _Restricted_Tessellation3D::_Restricted_Tessellation3D(const _ManifoldModel &in_model,
                                                         const std::vector<_Point3> &in_sites)
      : _num_sites(in_sites.size()),
        _sites(in_sites),
        _model(in_model)
  {
    _vertices.clear();
    _cells.clear();
    _cells.resize(_num_sites);
    _weights.clear();
    _weights.resize(_num_sites, 0);
    calculate_();
  }
  void _Restricted_Tessellation3D::calculate_(const std::vector<_Point3> &in_sites)
  {
    _sites = in_sites;
    _num_sites = _sites.size();
    _weights.clear();
    _weights.resize(_num_sites, 0);
    calculate_();
  }
  void _Restricted_Tessellation3D::calculate_(const std::vector<_Point3>& in_sites, const std::vector<double>& in_weights)
  {
      _sites = in_sites;
      _num_sites = _sites.size();
      _weights = in_weights;
      calculate_();
  }
} // namespace BGAL

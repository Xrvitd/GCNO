#pragma once
// RVD.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <set>
#include <omp.h> 
#include <cstdio>
#include <cstring>
#include <queue>

#include"MyHalfEdgeModel.hpp"
#include"PlaneCut.hpp"
#include"Knn.hpp"
#include"MyPointCloudModel.hpp"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Weighted_point_3.h>
#include <CGAL/bounding_box.h>

typedef CGAL::Simple_cartesian<double> K_T;
typedef K_T::FT FT;
typedef K_T::Point_3 Point_T;
typedef K_T::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K_T> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K_T, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

// for power d
typedef CGAL::Exact_predicates_inexact_constructions_kernel K_P;
typedef CGAL::Regular_triangulation_3<K_P> Regular_triangulation;
typedef K_P::Vector_3                                     Vector;
typedef K_P::Point_3 Point_P;
typedef K_P::Weighted_point_3 Wp;
typedef CGAL::Regular_triangulation_3<K_P>               Rt;


//#define N 1000010


using namespace std;

//vector<int> Repeat;
double gamma = 0.00000000001;

double FeatureWeight = 0.015, noFeatureWeight = 0.005;
//double FeatureWeight = 0.03, noFeatureWeight = 0.01;




struct MyPoint
{
    MyPoint(Eigen::Vector3d a)
    {
        p = a;

    }

    MyPoint(double a, double b, double c)
    {
        p.x() = a;
        p.y() = b;
        p.z() = c;
    }
    Eigen::Vector3d p;

    bool operator<(const MyPoint& a) const
    {



        double dis = (p - a.p).norm();
        if (dis < gamma)
        {
            return false;
        }

        if ((p.x() - a.p.x()) < 0.00000000001 && (p.x() - a.p.x()) > -0.00000000001)
        {
            if ((p.y() - a.p.y()) < 0.00000000001 && (p.y() - a.p.y()) > -0.00000000001)
            {
                return (p.z() < a.p.z());
            }
            return (p.y() < a.p.y());
        }
        return (p.x() < a.p.x());



    }
    bool operator==(const MyPoint& a) const
    {
        if ((p.x() - a.p.x()) < 0.00000000001 && (p.x() - a.p.x()) > -0.00000000001)
        {
            if ((p.y() - a.p.y()) < 0.00000000001 && (p.y() - a.p.y()) > -0.00000000001)
            {
                if ((p.z() - a.p.z()) < 0.00000000001 && (p.z() - a.p.z()) > -0.00000000001)
                {
                    return 1;
                }
            }

        }
        return 0;
    }
};

struct MyFace
{
    MyFace(Eigen::Vector3i a)
    {
        p = a;
    }
    MyFace(int a, int b, int c)
    {
        p.x() = a;
        p.y() = b;
        p.z() = c;
    }
    Eigen::Vector3i p;
    bool operator<(const MyFace& a) const
    {
        if (p.x() == a.p.x())
        {
            if (p.y() == a.p.y())
            {
                return p.z() > a.p.z();
            }
            return p.y() > a.p.y();
        }
        return p.x() > a.p.x();
    }
};



int nowP = 0;



//for power d



vector<Regular_triangulation::Weighted_point> wpoints;
vector<Point_P> points;
vector<double> X, Y, Z;
set<Point_P> Voronoi_vert;
void add_point(double x, double y, double z, double w) {
    wpoints.push_back(Wp(Point_P(x, y, z), w));
    points.push_back(Point_P(x, y, z));
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);

};

void Comput_RPD(string modelName)
{


    vector<string> files;

    //string filePath = "E:\\Dropbox\\SIG-2022-Feature-preserving-recon\\data\\Result_100models\\";
    string ss;
   
    {
        wpoints.clear(); points.clear(); X.clear(); Y.clear(); Z.clear(); Voronoi_vert.clear();
        vector<bool> IsFeature;
        vector<double> Weight;
        //vector<int> Repeat;

        //double FeatureWeight = 0.03, noFeatureWeight = 0.01;
        set<pair<int, int>> L;

        nowP = 0;
        map<pair<MyPoint, MyPoint>, int> RVD;
        map<pair<int, int>, bool> FlagOf2Points;
        map<pair<int, int>, int> linefix;
        map<pair<int, int>, int> Recon;
        vector<Eigen::Vector3d> VersPC_ori, Normal_ori;
        vector<Eigen::Vector3d> VersPC;
        map<MyFace, int> NewFaces;
        map<int, PlaneCut*> PCs;
        map<int, set<MyPoint>> RVDpoints;



        //74488
        //string filePath = "data\\";
        string filePath = "E:\\Dropbox\\MyProjects\\SIG-2022-Feature-preserving-recon\\data\\Results_noise0.0025\\";
        //string filePath = "D:\\SIG-BigExps\\00040123_8fc7d06caf264003a242597a_trimesh_000\\rsmall\\";
       // string modelName = "lucy";
        //string modelName = files[filenum];
        MyHalfEdgeModel* PoissonModel = new MyHalfEdgeModel();
        PoissonModel->ReadObjFile((filePath + modelName + "\\model_poisson.obj").c_str());
        double radis = 0;
        ifstream inRnum(filePath + modelName + "\\radis.txt");
        inRnum >> radis;
        inRnum.close();
        cout << radis << endl;
        //radis = 0.0025;

        FeatureWeight = radis*1.0;
        noFeatureWeight = radis / 3.0;
        //MyPointCloudModel PCmodel;
        //PCmodel.ReadXYZFile(("data\\" + modelName + "\\fandisk50000_6883.xyz").c_str(), true);

        //PCmodel.WriteXYZFile("data\\outtest.xyz", true);
        //VersPC_ori = PCmodel.GetVertices();
        ifstream inFnum(filePath + modelName + "\\FeaturePointNum.txt");

        ifstream inNewPs(filePath + modelName + "\\PoissonPoints_qc.xyz");
        ifstream inOriPs(filePath + modelName + "\\OriPoints_qc.xyz");
        int n = 0;
        double x11, y11, z11;
        while (inNewPs >> x11 >> y11 >> z11)
        {

            n++;
            bool isf;
            inFnum >> isf;
            IsFeature.push_back(isf);
            //Point_T query(VersPC_ori[i].x(), VersPC_ori[i].y(), VersPC_ori[i].z());
            //Point_T closest = tree.closest_point(query);
            Eigen::Vector3d NewP(x11, y11, z11);
            VersPC.push_back(NewP);

            inOriPs >> x11 >> y11 >> z11;

            Eigen::Vector3d NewP2(x11, y11, z11);
            VersPC_ori.push_back(NewP2);

            inOriPs >> x11 >> y11 >> z11;

            Eigen::Vector3d normal(x11, y11, z11);
            Normal_ori.push_back(normal);
        }

        cout << "Read point cloud.\n";

        map<MyPoint, int> Point2ID;
        vector<vector<int>> neighboor;

        //Knn KnnPC(filePath + modelName + "\\knn50_qc.txt", n, 50, true);
        Knn KnnPoisson(filePath + modelName + "\\knn50_poisson.txt", n, 50, false);
        cout << "Read Knn.\n";
        cout << n << endl;





		
        // filePath = "E:\\Dropbox\\SIG-2022-Feature-preserving-recon\\data\\ToRander\\";
        for (int i = 0; i < n; i++)
        {

            if (IsFeature[i])
            {
                add_point(VersPC[i].x(), VersPC[i].y(), VersPC[i].z(), FeatureWeight * FeatureWeight);
                Weight.push_back(FeatureWeight);

            }
            else
            {
                add_point(VersPC[i].x(), VersPC[i].y(), VersPC[i].z(), noFeatureWeight * noFeatureWeight);
                Weight.push_back(noFeatureWeight);
            }
            Point2ID[MyPoint(VersPC[i])] = i;
        }
        Regular_triangulation rt(wpoints.begin(), wpoints.end());
        rt.is_valid();
        cout << "make Regular_triangulation .\n";
        cout << n << endl;
        for (int i = 0; i < n; i++)
        {
            vector<int> tmp;
            neighboor.push_back(tmp);
        }

        for (const Rt::Vertex_handle vh : rt.finite_vertex_handles()) {

            std::vector<Rt::Vertex_handle> f_vertices;

            rt.finite_adjacent_vertices(vh, std::back_inserter(f_vertices));
            vector<int> nb_tmps;
            for (auto nb : f_vertices)
            {
                if (Point2ID.find(MyPoint(nb->point().x(), nb->point().y(), nb->point().z())) == Point2ID.end())
                {
                    cout << "ERROR!~\n\n\n";
                }
                nb_tmps.push_back(Point2ID[MyPoint(nb->point().x(), nb->point().y(), nb->point().z())]);
                //out_P << "v " << nb->point().x() << " " << nb->point().y() << " " << nb->point().z() << "\n";
            }
            if (Point2ID.find(MyPoint(vh->point().x(), vh->point().y(), vh->point().z())) == Point2ID.end())
            {
                cout << "ERROR!~\n\n\n";
            }
            neighboor[Point2ID[MyPoint(vh->point().x(), vh->point().y(), vh->point().z())]] = nb_tmps;
            /*if (Point2ID.find(MyPoint(vh->point().x(), vh->point().y(), vh->point().z())) == Point2ID.end())
            {
                cout << "ERROR!~\n\n\n";
            }*/

        }


        //Polyhedron polyhedron;
        //std::ifstream input("data\\" + modelName + "\\fandisk_poisson2.off");
        //input >> polyhedron;
        ////MyAABBTree* MyTree = new MyAABBTree("data\\cube_poisson.off");
        //Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);





        /*int deletecnt = 0;
        for (int i = 0; i < n; i++)
            Repeat.push_back(-1);

        for (int i = 0; i < n; i++)
        {
            if (KnnPC.neighboor[i].size() < 1)
            {
                Repeat[i] = 1; deletecnt++;
            }
        }
        cout << n - deletecnt << endl;*/
        //Repeat[42291] = 1;
        /*for (int i = 0; i < n; i++)
        {
            if (Repeat[i] != -1)
                continue;

            for (auto p : KnnPC.neighboor[i])
            {
                if (Repeat[p] == -1)
                {
                    double dis = (VersPC[i] - VersPC[p]).norm();
                    if (dis <= 0.001)
                    {
                        Repeat[p] = i;
                    }
                }

            }
        }

        for (int i = 0; i < n; i++)
        {
            for (auto p : KnnPC.neighboor[i])
            {
                if (Weight[i] == FeatureWeight && Weight[p] == noFeatureWeight)
                {
                    double dis = (VersPC[p] - VersPC[i]).norm();
                    if (dis <= 0.001)
                    {
                        Repeat[p] = i;
                    }
                }
            }
        }*/


        auto VersPoisson = PoissonModel->GetVertices();
        auto FacesPoisson = PoissonModel->GetFaces();
        cout << "Read Poisson model.\n";

        //#pragma omp parallel for schedule(dynamic,20)
        ofstream outPts(filePath + modelName + "\\CenterPoints_" + modelName + ".xyz");

        outPts.precision(15);
        outPts.flags(ios::left | ios::fixed);
        outPts.fill('0');







        int OpenmpCnt = 0;
        omp_set_num_threads(16);
#pragma omp parallel for //schedule(dynamic, 20)
        for (int ii = 0; ii < n; ii++)
        {
            // if (ii == 7)
                 //continue;

            /* if (Repeat[ii] != -1)
             {

                 continue;
             }*/

            OpenmpCnt++;
            if (OpenmpCnt % 10000 == 0)
            {
                cout << OpenmpCnt << '\n';
            }

            outPts << VersPC[ii].x() << " " << VersPC[ii].y() << " " << VersPC[ii].z() << "\n";


            PlaneCut* PC = new PlaneCut(VersPC[ii]);
            //ofstream out("data\\knn.xyz");
            //ofstream outf("data\\insideFaces.obj");
            //cout << ii << endl;
            /*Point_T query(VersPC[i].x(), VersPC[i].y(), VersPC[i].z());
            Point_T closest = tree.closest_point(query);
            Eigen::Vector3d NewP(closest.x(), closest.y(), closest.z());*/
            // out << VersPC[i].transpose() << endl;
            nowP = ii;
            //sort(KnnPC.neighboor[ii].begin(), KnnPC.neighboor[ii].end(),cmp);

            map<int, int> opposide;
            map<MyPoint, vector<int>> point2edge;
            // 把knn倒过来，应该会好。

            //for (auto p : KnnPC.neighboor[ii])
            int planeNum = 0;

            //for (int jj = KnnPC.neighboor[ii].size() - 1; jj >= 0; jj--)
            for (int jj = 0; jj < neighboor[ii].size(); jj++)
            {
                //cout << " " << jj << " ";
                auto p = neighboor[ii][jj];
                /*if (Repeat[p] != -1)
                    continue;*/
                if (!IsFeature[p] && IsFeature[ii])
                {
                    //continue;
                }


                //cout << (VersPC[p] - VersPC[ii]).norm() << "\n";
                //cout << p << endl;
                // out << VersPC[p].transpose() << endl;

                auto aa = VersPC[ii];
                auto bb = VersPC[p];
                auto r1 = Weight[ii], r2 = Weight[p];

                //Eigen::Vector3d MidPoint = VersPC[ii] + (VersPC[p] - VersPC[ii])*(Weight[ii]/(Weight[ii]+Weight[p]));  //voronoi
                Eigen::Vector3d MidPoint;
                double lambda = ((r1 * r1 - r2 * r2) / ((bb - aa).norm() * (bb - aa).norm()) + 1.0) / 2.0;
                MidPoint = (1 - lambda) * aa + lambda * bb;

                /*MidPoint.x() = (r1 * r1 - r2 * r2 + bb.x() * bb.x() - aa.x() * aa.x()) / (2 * bb.x() - 2 * aa.x());
                MidPoint.y() = (r1 * r1 - r2 * r2 + bb.y() * bb.y() - aa.y() * aa.y()) / (2 * bb.y() - 2 * aa.y());
                MidPoint.z() = (r1 * r1 - r2 * r2 + bb.z() * bb.z() - aa.z() * aa.z()) / (2 * bb.z() - 2 * aa.z());*///power 

                //out << MidPoint.transpose() << endl;
                K::Point_3 point1(MidPoint.x(), MidPoint.y(), MidPoint.z());
                MidPoint = VersPC[p] - VersPC[ii];

                MidPoint.normalize();
                K::Direction_3 dir(MidPoint.x(), MidPoint.y(), MidPoint.z());

                Plane plane(point1, dir);
                bool ifcut = PC->CutByPlane(plane);


                if (ifcut)
                {


                    opposide[planeNum] = p; planeNum++;
                    //Recon[make_pair(min(p, ii), max(p, ii))] = 1;
                }
                if (IsFeature[ii] && ifcut)
                {
                    linefix[make_pair(ii, p)] = 1;
                }


                /*if (ii == 317)
                {
                    string file = filePath + modelName + "\\debug\\" + to_string(jj) + "BigCube.obj";
                    PC->CuttedMesh->WriteObjFile(file.c_str());
                }*/

            }
            PCs[ii] = PC;
            //continue;
            //cout << endl;
            vector<bool> FlagPlane;
            for (int pl = 0; pl < PC->Planes.size(); pl++)
            {

                if (FlagOf2Points.find(make_pair(min(ii, opposide[pl]), max(ii, opposide[pl]))) == FlagOf2Points.end())
                {
                    FlagOf2Points[make_pair(min(ii, opposide[pl]), max(ii, opposide[pl]))] = true;
                    FlagPlane.push_back(1);
                }
                else
                {
                    FlagPlane.push_back(0);
                }

            }







            set<int> CloseFaces;
            for (auto p : KnnPoisson.neighboor[ii])
            {
                vector<int> PFaces = PoissonModel->GetFacesByPoint(p);
                for (auto f : PFaces)
                {
                    CloseFaces.insert(f);
                }
            }




           


            vector<int> insideF;
            vector<vector<Eigen::Vector3d>> CuttedF;
            vector<bool> aliveF;

            map<MyPoint, int> PointType;
            //vector<int> cuttedFp;
            for (auto f : CloseFaces)
            {

                bool cuted = 0;
                //cout << f << endl;
                int insidePoints = 0;
                auto P1 = VersPoisson[FacesPoisson[f].x()];
                auto P2 = VersPoisson[FacesPoisson[f].y()];
                auto P3 = VersPoisson[FacesPoisson[f].z()];
                K::Point_3 P1_Cgal(P1.x(), P1.y(), P1.z());
                K::Point_3 P2_Cgal(P2.x(), P2.y(), P2.z());
                K::Point_3 P3_Cgal(P3.x(), P3.y(), P3.z());
                MyPoint PP1(P1), PP2(P2), PP3(P3);
                bool vp1 = 0, vp2 = 0, vp3 = 0;
                PointType[PP1] = 0; PointType[PP2] = 0; PointType[PP3] = 0;
                for (int pl = 0; pl < PC->Planes.size(); pl++)
                {

                    if (!FlagPlane[pl])
                    {
                        continue;
                    }

                    auto plane = PC->Planes[pl];
                    if (plane.has_on_positive_side(P1_Cgal))
                    {
                        vp1 = 1;
                        PointType[PP1] = 1;
                    }
                    if (plane.has_on_positive_side(P2_Cgal))
                    {
                        vp2 = 1;
                        PointType[PP2] = 1;
                    }
                    if (plane.has_on_positive_side(P3_Cgal))
                    {
                        vp3 = 1;
                        PointType[PP3] = 1;
                    }
                    double f1, f2, f3;
                    f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                    f2 = plane.a() * P2.x() + plane.b() * P2.y() + plane.c() * P2.z() + plane.d();
                    f3 = plane.a() * P3.x() + plane.b() * P3.y() + plane.c() * P3.z() + plane.d();

                    if ((f1 > 0 && f2 < 0) || (f1 < 0 && f2 > 0))
                    {
                        cuted = 1;
                        //cuttedFp.push_back(opposide[plane]);
                        //Recon[make_pair(min(ii, opposide[pl]), max(ii, opposide[pl]))] = 1;
                        //break;
                    }
                    if ((f2 > 0 && f3 < 0) || (f2 < 0 && f3 > 0))
                    {
                        cuted = 1; //cuttedFp.push_back(opposide[plane]);
                        //Recon[make_pair(min(ii, opposide[pl]), max(ii, opposide[pl]))] = 1;
                        //break;
                    }
                    if ((f1 > 0 && f3 < 0) || (f1 < 0 && f3 > 0))
                    {
                        cuted = 1; //cuttedFp.push_back(opposide[plane]);
                        //Recon[make_pair(min(ii, opposide[pl]), max(ii, opposide[pl]))] = 1;
                        //break;
                    }

                }

                if (cuted)
                {
                    vector<Eigen::Vector3d> thisF;
                    thisF.push_back(VersPoisson[FacesPoisson[f].x()]);
                    thisF.push_back(VersPoisson[FacesPoisson[f].y()]);
                    thisF.push_back(VersPoisson[FacesPoisson[f].z()]);
                    CuttedF.push_back(thisF);

                    //cout << "type " << PointType[thisF[0]] << " " << PointType[thisF[1]] << " " << PointType[thisF[2]] << "\n";
                    aliveF.push_back(1);
                    //CutConnect[f].push_back(ii);
                }

                //if (vp1 == 0)
                //{
                //    insidePoints++;
                //}
                //if (vp2 == 0)
                //{
                //    insidePoints++;
                //}
                //if (vp3 == 0)
                //{
                //    insidePoints++;
                //}
                ////cout << insidePoints << endl;

                //if (insidePoints == 3)
                //{
                //    // in cell
                //    insideF.push_back(f);
                //}
                //if (insidePoints == 0)
                //{
                //    // out of cell

                //}
                //if (insidePoints >0&&insidePoints<3)
                //{
                //    // cut by cell
                //    vector<Eigen::Vector3d> thisF;
                //    thisF.push_back(VersPoisson[FacesPoisson[f].x()]);
                //    thisF.push_back(VersPoisson[FacesPoisson[f].y()]);
                //    thisF.push_back(VersPoisson[FacesPoisson[f].z()]);
                //    CuttedF.push_back(thisF);
                //    //cout << "type " << PointType[thisF[0]] << " " << PointType[thisF[1]] << " " << PointType[thisF[2]] << "\n";
                //    aliveF.push_back(1);
                //    CutConnect[f].push_back(ii);
                //}


            }



            for (int j = 0; j < CuttedF.size(); j++)
            {
                //cout << j << " " << CuttedF.size() << endl;
                auto f = CuttedF[j];
                if (!aliveF[j])
                    continue;

                vector<Eigen::Vector3d> newF;
                bool cgd = 0;
                //newF.push_back(f[0]);
                //cout << "f " << f.size() << endl;
                //cout << "type " << PointType[f[0]] << " " << PointType[f[1]] << " " << PointType[f[2]] << "\n";






                bool fd = 0;
                int fdp = 0;
                for (int pl = 0; pl < PC->Planes.size(); pl++)
                {

                    /*if (!FlagPlane[pl])
                    {
                        continue;
                    }*/

                    auto plane = PC->Planes[pl];
                    if (fd)
                        break;

                    vector<Eigen::Vector3d> newF_tmp;
                    for (int i = 0; i < f.size(); i++)
                    {
                        Eigen::Vector3d P1, P2;
                        P1 = f[i];
                        if (i == f.size() - 1)
                        {
                            P2 = f[0];
                        }
                        else
                        {
                            P2 = f[i + 1];
                        }
                        double f1, f2;
                        f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                        f2 = plane.a() * P2.x() + plane.b() * P2.y() + plane.c() * P2.z() + plane.d();
                        if (fabs(f1) < gamma)
                        {
                            newF_tmp.push_back(P1);
                            MyPoint PP1(P1);
                            point2edge[PP1].push_back(opposide[pl]);
                            continue;
                        }
                        if (fabs(f2) < gamma)
                        {
                            newF_tmp.push_back(P1);
                            MyPoint PP1(P1);
                            point2edge[P2].push_back(opposide[pl]);
                            continue;
                        }
                        if ((f1 > 0 && f2 < 0) || (f1 < 0 && f2 > 0))
                        {
                            //cout << f1 << " " << f2 << endl;

                            Eigen::Vector3d NewPoint;
                            f1 = fabs(f1); f2 = fabs(f2);
                            NewPoint.x() = P1.x() + (P2.x() - P1.x()) * (f1 / (f1 + f2));
                            NewPoint.y() = P1.y() + (P2.y() - P1.y()) * (f1 / (f1 + f2));
                            NewPoint.z() = P1.z() + (P2.z() - P1.z()) * (f1 / (f1 + f2));
                            /*if (PointType.find(NewPoint) != PointType.end())
                            {
                                newF_tmp.push_back(P1);
                                continue;
                            }*/
                            MyPoint PP1(NewPoint);


                            newF_tmp.push_back(P1);
                            newF_tmp.push_back(NewPoint);

                            point2edge[PP1].push_back(opposide[pl]);
                            PointType[NewPoint] = -1;
                            fd = 1;
                        }
                        else
                        {
                            newF_tmp.push_back(P1);
                        }
                    }
                    if (fd)
                    {
                        newF = newF_tmp;
                        fdp = opposide[pl];
                    }
                }







                if (fd == 0)
                {
                    bool isInside = 0;// for test
                    bool onEdge = 1;
                    double maxFo = -99999999.0;
                    for (int i = 0; i < f.size(); i++)
                    {
                        if (PointType[f[i]] == 0)
                        {
                            isInside = 1;
                        }

                        Eigen::Vector3d P1;
                        P1 = f[i];
                        for (int pl = 0; pl < PC->Planes.size(); pl++)
                        {
                            /*if (!FlagPlane[pl])
                            {
                                continue;
                            }*/
                            auto plane = PC->Planes[pl];

                            double f1;
                            f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                            //maxF = max(maxF, f1);
                            if (f1 > maxFo)
                            {
                                maxFo = f1;
                            }
                        }

                    }
                    if (maxFo < gamma)
                    {
                        isInside = 1;
                    }
                    if (isInside)
                    {
                        vector < bool > Nps;
                        vector < Eigen::Vector3d > Newps;
                        for (int i = 0; i < f.size(); i++)
                        {
                            Eigen::Vector3d P1;
                            P1 = f[i];
                            double maxF = -99999999.0;
                            for (int pl = 0; pl < PC->Planes.size(); pl++)
                            {
                                /*if (!FlagPlane[pl])
                                {
                                    continue;
                                }*/
                                auto plane = PC->Planes[pl];

                                double f1, f2;
                                f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                                maxF = max(maxF, f1);
                            }

                            if (fabs(maxF) < gamma)
                            {
                                Nps.push_back(1);
                                Newps.push_back(f[i]);
                            }
                            else
                            {
                                Nps.push_back(0);
                            }

                        }
                        if (Newps.size() < 2)
                        {
                            continue;
                        }


                        /*for (int i = 0; i < Newps.size(); i++)
                        {
                            Eigen::Vector3d P1, P2;
                            P1 = Newps[i];
                            if (i == Newps.size() - 1)
                            {
                                P2 = Newps[0];
                            }
                            else
                            {
                                P2 = Newps[i + 1];
                            }
                            if (P1.x() > P2.x())
                            {
                                RVD[make_pair(P1, P2)] = 1;
                            }
                            else
                            {
                                RVD[make_pair(P2, P1)] = 1;
                            }

                        }
                        for (int i = 0; i < f.size(); i++)
                        {
                            if (PointType[f[i]] == 0)
                            {
                                Eigen::Vector3d P1, P2;

                                int k = i;
                                while (true)
                                {
                                    k++;
                                    if (k == f.size())
                                    {
                                        k = 0;
                                    }
                                    if (Nps[k] == 1)
                                    {
                                        P1 = f[i];
                                        break;
                                    }
                                }
                                k = i;
                                while (true)
                                {
                                    k--;
                                    if (k == -1)
                                    {
                                        k = f.size()-1;
                                    }
                                    if (Nps[k] == 1)
                                    {
                                        P2 = f[i];
                                        break;
                                    }
                                }
                                if (P1.x() > P2.x())
                                {
                                    RVD[make_pair(P1, P2)] = 0;
                                }
                                else
                                {
                                    RVD[make_pair(P2, P1)] = 0;
                                }

                            }
                        }*/




                        vector<Eigen::Vector3d> Drawedpoints;
                        for (int i = 0; i < f.size(); i++)
                        {
                            Eigen::Vector3d P1, P2, Pmid;
                            P1 = f[i];
                            if (i == f.size() - 1)
                            {
                                P2 = f[0];
                            }
                            else
                            {
                                P2 = f[i + 1];
                            }
                            bool V1, V2;
                            V1 = Nps[i];
                            if (i == f.size() - 1)
                            {
                                V2 = Nps[0];
                            }
                            else
                            {
                                V2 = Nps[i + 1];
                            }
                            bool doublecheck = 0;
                            Pmid.x() = P1.x() + (P2.x() - P1.x()) / 2;
                            Pmid.y() = P1.y() + (P2.y() - P1.y()) / 2;
                            Pmid.z() = P1.z() + (P2.z() - P1.z()) / 2;
                            double maxF = -99999999.0;
                            for (int pl = 0; pl < PC->Planes.size(); pl++)
                            {
                                /*if (!FlagPlane[pl])
                                {
                                    continue;
                                }*/
                                auto plane = PC->Planes[pl];

                                double f1, f2;
                                f1 = plane.a() * Pmid.x() + plane.b() * Pmid.y() + plane.c() * Pmid.z() + plane.d();
                                maxF = max(maxF, f1);

                            }
                            if (fabs(maxF) < gamma)
                            {
                                doublecheck = 0;

                            }
                            else
                            {
                                doublecheck = 1;
                            }



                            //if (PointType[P1] != 0&& PointType[P2] != 0)
                            if (V1 && V2 && !doublecheck)
                            {
                                Drawedpoints.push_back(P1);
                                Drawedpoints.push_back(P2);


                                MyPoint p1(P1), p2(P2);
                                RVDpoints[ii].insert(p1);
                                RVDpoints[ii].insert(p2);


                                if (p2 < p1)
                                {

                                    if (RVD.find(make_pair(p1, p2)) == RVD.end())
                                    {
                                        RVD[make_pair(p1, p2)] = ii;
                                    }
                                    else
                                    {
                                        int lst = RVD[make_pair(p1, p2)];
                                        //Recon[make_pair(min(lst, ii), max(lst, ii))] = 1;
                                    }


                                }
                                else
                                {
                                    if (RVD.find(make_pair(p2, p1)) == RVD.end())
                                    {
                                        RVD[make_pair(p2, p1)] = ii;
                                    }
                                    else
                                    {
                                        int lst = RVD[make_pair(p2, p1)];
                                        //Recon[make_pair(min(lst, ii), max(lst, ii))] = 1;
                                    }

                                    //RVD[make_pair(P2, P1)] = 1;
                                }
                            }
                        }
                        map<int, int> cnt;
                        for (auto mp : Drawedpoints)
                        {
                            auto P1 = mp;
                            double maxF = -999999.0;
                            for (int pl = 0; pl < PC->Planes.size(); pl++)
                            {
                                /*if (!FlagPlane[pl])
                                {
                                    continue;
                                }*/
                                auto plane = PC->Planes[pl];
                                double f1;
                                f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                                if (fabs(f1) < gamma)
                                {
                                    int oppoP = opposide[pl];
                                    if (cnt.find(oppoP) == cnt.end())
                                    {
                                        cnt[oppoP] = 1;
                                    }
                                    else
                                    {
                                        cnt[oppoP]++;
                                    }
                                }
                            }

                        }
                        for (auto mp : cnt)
                        {
                            if (mp.second >= 1)
                            {
                                Recon[make_pair(min(mp.first, ii), max(mp.first, ii))] = 1;
                            }
                        }




                    }
                    else
                    {
                        aliveF[j] = 0;

                    }
                }
                else
                {
                    // if(fd)
                    vector<Eigen::Vector3d> newF1, newF2;
                    bool fdd = 0;
                    //cout << "NewF " << newF.size() << endl;
                    for (int i = 0; i < newF.size(); i++)
                    {
                        // cout << "type " << PointType[newF[i]] << "\n";
                        if (PointType[newF[i]] == -1)
                        {
                            PointType[newF[i]] = -2;
                            if (fdd == 0)
                            {
                                newF1.push_back(newF[i]); newF2.push_back(newF[i]);
                                //cout << i << "  111\n"; //cout << i << "  222\n";
                                fdd = 1;

                                continue;
                            }
                            else
                            {
                                newF1.push_back(newF[i]); //cout << i << "  111\n";
                                newF2.push_back(newF[i]); //cout << i << "  222\n";

                                fdd = 0;
                                continue;
                            }



                            /*Eigen::Vector3d P1 = newF[i];
                            double maxF = -99999999.0;
                            for (auto plane : PC->Planes)
                            {
                                double f1, f2;
                                f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                                maxF = max(maxF, f1);

                            }
                            if (fabs(maxF) < 0.0001)
                            {

                            }*/
                        }
                        if (fdd == 0)
                        {
                            newF1.push_back(newF[i]); //cout << i << "  111\n";
                        }
                        else
                        {
                            newF2.push_back(newF[i]); //cout << i << "  222\n";
                        }
                    }
                    aliveF[j] = 0;
                    aliveF.push_back(1);
                    aliveF.push_back(1);
                    CuttedF.push_back(newF1);
                    CuttedF.push_back(newF2);

                    //cuttedFp.push_back(cuttedFp[j]);
                    //cuttedFp.push_back(cuttedFp[j]);




                    //return 0;
                }




            }



            for (auto rvdp : RVDpoints[ii])
            {
                auto P1 = rvdp.p;
                double maxF = -999999.0;
                for (int pl = 0; pl < PC->Planes.size() - 1; pl++)
                {
                    for (int pl2 = pl + 1; pl2 < PC->Planes.size(); pl2++)
                    {
                        auto plane1 = PC->Planes[pl];
                        auto plane2 = PC->Planes[pl2];
                        double f1;
                        f1 = plane1.a() * P1.x() + plane1.b() * P1.y() + plane1.c() * P1.z() + plane1.d();
                        double f2;
                        f2 = plane2.a() * P1.x() + plane2.b() * P1.y() + plane2.c() * P1.z() + plane2.d();
                        if (fabs(f1) < gamma && fabs(f2) < gamma)
                        {
                            int aa, bb, cc;

                            aa = min(ii, min(opposide[pl], opposide[pl2]));
                            cc = max(ii, max(opposide[pl], opposide[pl2]));
                            if (ii != aa && ii != cc)
                                bb = ii;
                            if (opposide[pl] != aa && opposide[pl] != cc)
                                bb = opposide[pl];
                            if (opposide[pl2] != aa && opposide[pl2] != cc)
                                bb = opposide[pl2];

                            MyFace ff(aa, bb, cc);
                            NewFaces[ff] = ii;

                        }



                    }

                }

            }





            //break;
        }



        map<int, bool> ifInsideOtherRVD;
        for (int ii = 0; ii < n; ii++)
        {
            /*if (Repeat[ii] != -1)
                continue;*/
            auto PC = PCs[ii];

            ifInsideOtherRVD[ii] = 0;
            bool flagRvdp = 0;
            for (auto rvdp : RVDpoints[ii])
            {

                bool flagii = 0;
                for (int jj = 0; jj < neighboor[ii].size(); jj++)
                {
                    auto p = neighboor[ii][jj];

                    if (PCs.find(p) == PCs.end())
                    {
                        continue;
                    }
                    if (PCs[p]->Planes.empty())
                    {
                        continue;
                    }
                    auto P1 = rvdp.p;

                    flagii = 0;
                    for (int pl = 0; pl < PCs[p]->Planes.size(); pl++)
                    {
                        auto plane = PCs[p]->Planes[pl];
                        double f1;
                        f1 = plane.a() * P1.x() + plane.b() * P1.y() + plane.c() * P1.z() + plane.d();
                        if (f1 >= 0 || fabs(f1) <= gamma)
                        {
                            flagii = 1;
                            break;
                        }
                    }

                    if (flagii == 0)
                    {
                        break;
                    }
                }
                if (flagii == 1)
                {
                    flagRvdp = 1;
                    break;
                }
            }
            if (!flagRvdp)
            {
                ifInsideOtherRVD[ii] = 1;
            }


        }


        map<MyPoint, int> degree;
        /*for (auto p : RvdPoints2Face)
        {
            degree[p.first] = 0;
        }*/
        ofstream outRVD(filePath + modelName + "\\RVD_" + modelName + ".obj");

        outRVD.precision(15);
        outRVD.flags(ios::left | ios::fixed);
        outRVD.fill('0');

        map<MyPoint, int> Point2IDD;
        int pid = 0;
        for (auto mp : RVD)
        {
            //if (mp.second == 1)
            {
                MyPoint P1 = mp.first.first;
                MyPoint P2 = mp.first.second;

                //MyPoint PP1(P1), PP2(P2);

                if (Point2IDD.find(P1) == Point2IDD.end())
                {
                    pid++;
                    Point2IDD[P1] = pid;
                    outRVD << "v " << P1.p.transpose() << "\n";
                }
                if (Point2IDD.find(P2) == Point2IDD.end())
                {
                    pid++;
                    Point2IDD[P2] = pid;
                    outRVD << "v " << P2.p.transpose() << "\n";
                }
                outRVD << "l " << Point2IDD[P1] << " " << Point2IDD[P2] << "\n";


                if (degree.find(P1) == degree.end())
                {
                    degree[P1] = 0;
                }
                if (degree.find(P2) == degree.end())
                {
                    degree[P2] = 0;
                }
                degree[P1]++;
                degree[P2]++;


            }
        }




        //memset(h, -1, sizeof(h));
        //memset(nxt, -1, sizeof(h));

        //ofstream outl(filePath + modelName + "\\PointCloudConnection.obj");
        //outl.precision(15);
        //outl.flags(ios::left | ios::fixed);
        //outl.fill('0');

        //for (auto p : VersPC_ori)
        //{
        //    /*if (Repeat[ii] != -1)
        //        continue;*/
        //    outl << "v " << p.transpose() << endl;
        //}

        //for (auto mp : Recon)
        //{

        //    add(mp.first.first, mp.first.second);
        //    add(mp.first.second, mp.first.first);


        //    outl << "l " << mp.first.first + 1 << " " << mp.first.second + 1 << "\n";
        //}
        ofstream outRemesh(filePath + modelName + "\\Remesh_" + modelName + ".obj");
        outRemesh.precision(15);
        outRemesh.flags(ios::left | ios::fixed);
        outRemesh.fill('0');


        for (auto p : VersPC_ori)
        {
            /*if (Repeat[ii] != -1)
                continue;*/
            outRemesh << "v " << p.transpose() << endl;
        }

        map<pair<int, int>, int> DegreeOfEdge;



        cout << NewFaces.size() << endl;
        for (auto f : NewFaces)
        {
            if (ifInsideOtherRVD[f.second] == 1)
            {
                continue;
            }


            int aa, bb, cc;
            aa = f.first.p.x();
            bb = f.first.p.y();
            cc = f.first.p.z();

            if (DegreeOfEdge.find(make_pair(min(aa, bb), max(aa, bb))) == DegreeOfEdge.end())
            {
                DegreeOfEdge[make_pair(min(aa, bb), max(aa, bb))] = 0;
            }
            if (DegreeOfEdge.find(make_pair(min(aa, cc), max(aa, cc))) == DegreeOfEdge.end())
            {
                DegreeOfEdge[make_pair(min(aa, cc), max(aa, cc))] = 0;
            }
            if (DegreeOfEdge.find(make_pair(min(cc, bb), max(cc, bb))) == DegreeOfEdge.end())
            {
                DegreeOfEdge[make_pair(min(cc, bb), max(cc, bb))] = 0;
            }
            DegreeOfEdge[make_pair(min(aa, bb), max(aa, bb))]++;
            DegreeOfEdge[make_pair(min(aa, cc), max(aa, cc))]++;
            DegreeOfEdge[make_pair(min(cc, bb), max(cc, bb))]++;




            //outRemesh << "f " << f.first.p.x() + 1 << " " << f.first.p.y() + 1 << " " << f.first.p.z() + 1 << "\n";
        }

        while (true)
        {

            bool vvv = 0;
            for (auto f : NewFaces)
            {
                if (f.second == -1)
                {
                    continue;
                }
                if (ifInsideOtherRVD[f.second] == 1)
                {
                    continue;
                }


                int aa, bb, cc;
                aa = f.first.p.x();
                bb = f.first.p.y();
                cc = f.first.p.z();

                if (DegreeOfEdge[make_pair(min(aa, bb), max(aa, bb))] < 2 || DegreeOfEdge[make_pair(min(aa, cc), max(aa, cc))] < 2 || DegreeOfEdge[make_pair(min(cc, bb), max(cc, bb))] < 2)
                {
                    if (DegreeOfEdge[make_pair(min(aa, bb), max(aa, bb))] > 2 || DegreeOfEdge[make_pair(min(aa, cc), max(aa, cc))] > 2 || DegreeOfEdge[make_pair(min(cc, bb), max(cc, bb))] > 2)
                    {
                        DegreeOfEdge[make_pair(min(aa, bb), max(aa, bb))]--;
                        DegreeOfEdge[make_pair(min(aa, cc), max(aa, cc))]--;
                        DegreeOfEdge[make_pair(min(cc, bb), max(cc, bb))]--;
                        MyFace ff(aa, bb, cc);
                        NewFaces[ff] = -1;
                        vvv = 1;
                    }
                }

            }
            if (!vvv)
            {
                break;
            }
        }
        vector<Eigen::Vector3i> RemeshFs;
        int fid = 0;
        map<pair<int, int>, vector<int>> Edge2Faceid;
        vector<bool> FaceFlag;
        for (auto f : NewFaces)
        {
            if (f.second == -1)
            {
                continue;
            }
            if (ifInsideOtherRVD[f.second] == 1)
            {
                continue;
            }

            Edge2Faceid[make_pair(min(f.first.p.x(), f.first.p.y()), max(f.first.p.x(), f.first.p.y()))].push_back(fid);
            Edge2Faceid[make_pair(min(f.first.p.z(), f.first.p.y()), max(f.first.p.z(), f.first.p.y()))].push_back(fid);
            Edge2Faceid[make_pair(min(f.first.p.x(), f.first.p.z()), max(f.first.p.x(), f.first.p.z()))].push_back(fid);

            FaceFlag.push_back(0);
            RemeshFs.push_back(Eigen::Vector3i(f.first.p.x(), f.first.p.y(), f.first.p.z()));
            fid++;

            //    outRemesh << "f " << f.first.p.x() + 1 << " " << f.first.p.y() + 1 << " " << f.first.p.z() + 1 << "\n";
        }

        // fix normal

        int NormalRightFaceID = -1;
        for (int i = 0; i < RemeshFs.size(); i++)
        {
            auto f = RemeshFs[i];
            if (!IsFeature[f.x()] && !IsFeature[f.y()] && !IsFeature[f.z()])
            {
                NormalRightFaceID = i;

                auto nor1 = (VersPC_ori[f.y()] - VersPC_ori[f.x()]).cross(VersPC_ori[f.z()] - VersPC_ori[f.y()]);
                nor1.normalize();
                auto nor2 = (VersPC_ori[f.x()] - VersPC_ori[f.y()]).cross(VersPC_ori[f.z()] - VersPC_ori[f.x()]);
                nor2.normalize();

                Eigen::Vector3d nor3 = (Normal_ori[f.x()].normalized() + Normal_ori[f.y()].normalized() + Normal_ori[f.z()].normalized()) / 3.0;
                nor3.normalize();
                double dis1 = (nor3 - nor1).norm();
                double dis2 = (nor3 - nor2).norm();

                if (dis2 < dis1)
                {
                    int tmpp = RemeshFs[i].x();
                    RemeshFs[i].x() = RemeshFs[i].y();
                    RemeshFs[i].y() = tmpp;
                }
                FaceFlag[i] = 1;

                break;
            }
        }

        queue<int> FaceQue; FaceQue.push(NormalRightFaceID);
        while (!FaceQue.empty())
        {
            int faceid = FaceQue.front();
            FaceQue.pop();

            auto f = RemeshFs[faceid];
            auto FaceSet1 = Edge2Faceid[make_pair(min(f.x(), f.y()), max(f.x(), f.y()))];
            int aa = f.x(), bb = f.y();
            for (auto ff : FaceSet1)
            {
                if (FaceFlag[ff] == 0)
                {
                    FaceQue.push(ff);
                    FaceFlag[ff] = 1;
                    auto tf = RemeshFs[ff];
                    if (tf.x() == aa)
                    {
                        if (tf.y() == bb)
                        {
                            int tmpf = RemeshFs[ff].y();
                            RemeshFs[ff].y() = RemeshFs[ff].x();
                            RemeshFs[ff].x() = tmpf;
                            continue;
                        }
                    }

                    if (tf.y() == aa)
                    {
                        if (tf.z() == bb)
                        {
                            int tmpf = RemeshFs[ff].y();
                            RemeshFs[ff].y() = RemeshFs[ff].z();
                            RemeshFs[ff].z() = tmpf;
                            continue;
                        }
                    }

                    if (tf.z() == aa)
                    {
                        if (tf.x() == bb)
                        {
                            int tmpf = RemeshFs[ff].z();
                            RemeshFs[ff].z() = RemeshFs[ff].x();
                            RemeshFs[ff].x() = tmpf;
                            continue;
                        }
                    }


                }
            }
            // second edge 
            auto FaceSet2 = Edge2Faceid[make_pair(min(f.z(), f.y()), max(f.z(), f.y()))];
            aa = f.y(); bb = f.z();
            for (auto ff : FaceSet2)
            {
                if (FaceFlag[ff] == 0)
                {
                    FaceFlag[ff] = 1; FaceQue.push(ff);
                    auto tf = RemeshFs[ff];
                    if (tf.x() == aa)
                    {
                        if (tf.y() == bb)
                        {
                            int tmpf = RemeshFs[ff].y();
                            RemeshFs[ff].y() = RemeshFs[ff].x();
                            RemeshFs[ff].x() = tmpf;
                            continue;
                        }
                    }

                    if (tf.y() == aa)
                    {
                        if (tf.z() == bb)
                        {
                            int tmpf = RemeshFs[ff].y();
                            RemeshFs[ff].y() = RemeshFs[ff].z();
                            RemeshFs[ff].z() = tmpf;
                            continue;
                        }
                    }

                    if (tf.z() == aa)
                    {
                        if (tf.x() == bb)
                        {
                            int tmpf = RemeshFs[ff].z();
                            RemeshFs[ff].z() = RemeshFs[ff].x();
                            RemeshFs[ff].x() = tmpf;
                            continue;
                        }
                    }
                }
            }

            // third edge
            auto FaceSet3 = Edge2Faceid[make_pair(min(f.x(), f.z()), max(f.x(), f.z()))];
            aa = f.z(); bb = f.x();
            for (auto ff : FaceSet3)
            {
                if (FaceFlag[ff] == 0)
                {
                    FaceFlag[ff] = 1; FaceQue.push(ff);
                    auto tf = RemeshFs[ff];
                    if (tf.x() == aa)
                    {
                        if (tf.y() == bb)
                        {
                            int tmpf = RemeshFs[ff].y();
                            RemeshFs[ff].y() = RemeshFs[ff].x();
                            RemeshFs[ff].x() = tmpf;
                            continue;
                        }
                    }

                    if (tf.y() == aa)
                    {
                        if (tf.z() == bb)
                        {
                            int tmpf = RemeshFs[ff].y();
                            RemeshFs[ff].y() = RemeshFs[ff].z();
                            RemeshFs[ff].z() = tmpf;
                            continue;
                        }
                    }

                    if (tf.z() == aa)
                    {
                        if (tf.x() == bb)
                        {
                            int tmpf = RemeshFs[ff].z();
                            RemeshFs[ff].z() = RemeshFs[ff].x();
                            RemeshFs[ff].x() = tmpf;
                            continue;
                        }
                    }
                }
            }





        }

        /*for (auto f : RemeshFs)
        {
            outRemesh << "f " << f.x() + 1 << " " << f.y() + 1 << " " << f.z() + 1 << endl;

        }*/


        // delete nonf
        map<pair<int, int>, int> EdgeCnt;
        for (auto f : RemeshFs)
        {
            //outRemesh << "f " << f.x() + 1 << " " << f.y() + 1 << " " << f.z() + 1 << endl;
            int aa = f.x();
            int bb = f.y();
            int cc = f.z();
            if (EdgeCnt.find(make_pair(min(aa, bb), max(aa, bb))) == EdgeCnt.end())
            {
                EdgeCnt[make_pair(min(aa, bb), max(aa, bb))] = 1;
            }
            else
            {
                EdgeCnt[make_pair(min(aa, bb), max(aa, bb))]++;
            }

            if (EdgeCnt.find(make_pair(min(aa, cc), max(aa, cc))) == EdgeCnt.end())
            {
                EdgeCnt[make_pair(min(aa, cc), max(aa, cc))] = 1;
            }
            else
            {
                EdgeCnt[make_pair(min(aa, cc), max(aa, cc))]++;
            }

            if (EdgeCnt.find(make_pair(min(cc, bb), max(cc, bb))) == EdgeCnt.end())
            {
                EdgeCnt[make_pair(min(cc, bb), max(cc, bb))] = 1;
            }
            else
            {
                EdgeCnt[make_pair(min(cc, bb), max(cc, bb))]++;
            }
        }
        for (auto f : RemeshFs)
        {
            int aa = f.x();
            int bb = f.y();
            int cc = f.z();
            /*if (EdgeCnt[make_pair(min(aa, bb), max(aa, bb))] > 2 && (EdgeCnt[make_pair(min(cc, bb), max(cc, bb))]==1 || EdgeCnt[make_pair(min(aa, cc), max(aa, cc))]==1))
            {
                continue;
            }

            if (EdgeCnt[make_pair(min(cc, bb), max(cc, bb))] > 2 && (EdgeCnt[make_pair(min(aa, bb), max(aa, bb))] == 1 || EdgeCnt[make_pair(min(aa, cc), max(aa, cc))] == 1))
            {
                continue;
            }

            if (EdgeCnt[make_pair(min(aa, cc), max(aa, cc))] > 2 && (EdgeCnt[make_pair(min(cc, bb), max(cc, bb))] == 1 || EdgeCnt[make_pair(min(aa, bb), max(aa, bb))] == 1))
            {
                continue;
            }*/

            if (EdgeCnt[make_pair(min(aa, bb), max(aa, bb))] == 1)
            {
                continue;
            }

            if (EdgeCnt[make_pair(min(cc, bb), max(cc, bb))] == 1)
            {
                continue;
            }

            if (EdgeCnt[make_pair(min(aa, cc), max(aa, cc))] == 1)
            {
                continue;
            }


            outRemesh << "f " << f.x() + 1 << " " << f.y() + 1 << " " << f.z() + 1 << endl;
        }

        ofstream outEdge(filePath + modelName + "\\Edges_" + modelName + ".obj");
        outEdge.precision(15);
        outEdge.flags(ios::left | ios::fixed);
        outEdge.fill('0');

        set<pair<int, int>> edges;
        for (auto p : VersPC_ori)
        {
            /*if (Repeat[ii] != -1)
                continue;*/
            outEdge << "v " << p.transpose() << endl;
        }
        for (auto f : RemeshFs)
        {
            //outRemesh << "f " << f.x() + 1 << " " << f.y() + 1 << " " << f.z() + 1 << endl;
            edges.insert(make_pair(min(f.x(), f.y()), max(f.x(), f.y())));
            edges.insert(make_pair(min(f.z(), f.y()), max(f.z(), f.y())));
            edges.insert(make_pair(min(f.x(), f.z()), max(f.x(), f.z())));
        }
        for (auto e : edges)
        {
            outEdge << "l " << e.first + 1 << " " << e.second + 1 << "\n";
        }
        outEdge.close();

        ofstream outFeatureLine(filePath + modelName + "\\FeatureLine_" + modelName + ".obj");
        outFeatureLine.precision(15);
        outFeatureLine.flags(ios::left | ios::fixed);
        outFeatureLine.fill('0');


        for (auto p : VersPC_ori)
        {
            /*if (Repeat[ii] != -1)
                continue;*/
            outFeatureLine << "v " << p.transpose() << endl;
        }

        /*for (auto e : edges)
        {
            if (IsFeature[e.first] && IsFeature[e.second])
                outFeatureLine << "l " << e.first + 1 << " " << e.second + 1 << "\n";
        }*/
        

        MyHalfEdgeModel FinalModel;
        FinalModel.ReadObjFile((filePath + modelName + "\\Remesh_" + modelName + ".obj").c_str());
        auto Fedges = FinalModel.GetEdges();
        auto Ffaces = FinalModel.GetFaces();
        auto Fvecs = FinalModel.GetVertices();
        for (auto e : Fedges)
        {
            if (IsFeature[e.leftVert] && IsFeature[e.rightVert])
            {
                auto f1 = e.indexOfFrontFace;
                auto f2 = Fedges[e.indexOfReverseEdge].indexOfFrontFace;
                Eigen::Vector3d Nf1, Nf2;
				//compute normal of f1
				auto v1 = Fvecs[Ffaces[f1].x()];
				auto v2 = Fvecs[Ffaces[f1].y()];
				auto v3 = Fvecs[Ffaces[f1].z()];
				Nf1 = (v2 - v1).cross(v3 - v1).normalized();
				//compute normal of f2
				v1 = Fvecs[Ffaces[f2].x()];
				v2 = Fvecs[Ffaces[f2].y()];
				v3 = Fvecs[Ffaces[f2].z()];
				Nf2 = (v2 - v1).cross(v3 - v1).normalized();
				//compute angle between f1 and f2
				double angle = acos(Nf1.dot(Nf2));
               /* if (angle > 3.1415926535 / 2)
                {
					angle = 3.1415926535 - angle;
                }*/
				// angle to drgee
				angle = angle * 180 / 3.1415926535;
                if (angle > 40&& angle < 140)
				{
					outFeatureLine << "l " << e.leftVert + 1 << " " << e.rightVert + 1 << "\n";
				}
                
            }
        }
        outFeatureLine.close();
        // output a single feature line model .







    }





    std::cout << "Hello World!\n";





}



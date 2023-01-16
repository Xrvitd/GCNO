#pragma once
// RVD.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <set>
#include <omp.h> 
#include <cstdio>
#include <cstring>
#include <queue>
#include <map>

#include"MyHalfEdgeModel.hpp"
#include"PlaneCut.hpp"
#include"Knn.hpp"
#include"MyPointCloudModel.hpp"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/optimal_bounding_box.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Weighted_point_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Surface_mesh.h>

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

double FeatureWeight = 0.005, noFeatureWeight = 0.005;
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

pair<vector<vector<Eigen::Vector3d>>, pair<vector<double>, map<MyPoint, set<MyPoint>>>> Comput_RPD(string modelpath,string modelName,bool IfNor)
{


    vector<string> files;
    vector<vector<Eigen::Vector3d>> F_VDPs;
    vector<double> Areas;
    map<MyPoint,int> VDPsMap;
    //string filePath = "E:\\Dropbox\\SIG-2022-Feature-preserving-recon\\data\\Result_100models\\";
    string ss;
   
    {
        wpoints.clear(); points.clear(); X.clear(); Y.clear(); Z.clear(); Voronoi_vert.clear();

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
        string filePath = modelpath;
        //string filePath = "D:\\SIG-BigExps\\00040123_8fc7d06caf264003a242597a_trimesh_000\\rsmall\\";
       // string modelName = "lucy";
        //string modelName = files[filenum];
        /*MyHalfEdgeModel* PoissonModel = new MyHalfEdgeModel();
        PoissonModel->ReadObjFile((filePath + "\\model_poisson_"+ modelName +".obj").c_str());*/
        double radis = 0;
        /*ifstream inRnum(filePath  + "\\radis_" + modelName + ".txt");
        inRnum >> radis;
        inRnum.close();
        cout << radis << endl;*/
        //radis = 0.0025;

        FeatureWeight = 1;
        noFeatureWeight = 1;
        //MyPointCloudModel PCmodel;
        //PCmodel.ReadXYZFile(("data\\" + modelName + "\\fandisk50000_6883.xyz").c_str(), true);

        //PCmodel.WriteXYZFile("data\\outtest.xyz", true);
        //VersPC_ori = PCmodel.GetVertices();

        cout << filePath + modelName + ".xyz" << endl;
        ifstream inNewPs(filePath + modelName + ".xyz");
        vector<Point_P> obbPoints;
        int n = 0;
        double x11, y11, z11;
        while (inNewPs >> x11 >> y11 >> z11)
        {
            n++;
            Eigen::Vector3d NewP(x11, y11, z11);
            VersPC.push_back(NewP);
			obbPoints.push_back(Point_P(x11, y11, z11));
            
            inNewPs >> x11 >> y11 >> z11;

            Eigen::Vector3d normal(x11, y11, z11);
            normal.normalize();
            Normal_ori.push_back(normal);
        }
        array<Point_P, 8> obb_points;
        CGAL::oriented_bounding_box(obbPoints, obb_points);
        /*CGAL::Surface_mesh<Point_P>  obb_sm;
        CGAL::make_hexahedron(obb_points[0], obb_points[1], obb_points[2], obb_points[3],
            obb_points[4], obb_points[5], obb_points[6], obb_points[7], obb_sm);*/
        vector<Eigen::Vector3d> OBBmeshPs;
        vector<Eigen::Vector3i> OBBmeshFs;
        
		for (int i = 0; i < 8; i++)
		{
			OBBmeshPs.push_back(Eigen::Vector3d(obb_points[i].x(), obb_points[i].y(), obb_points[i].z()));
		}
        
		OBBmeshFs.push_back(Eigen::Vector3i(1, 0, 2));
		OBBmeshFs.push_back(Eigen::Vector3i(2, 0, 3));
		OBBmeshFs.push_back(Eigen::Vector3i(4, 5, 6));
		OBBmeshFs.push_back(Eigen::Vector3i(4, 6, 7));
		OBBmeshFs.push_back(Eigen::Vector3i(0, 1, 5));
		OBBmeshFs.push_back(Eigen::Vector3i(0, 5, 4));
		OBBmeshFs.push_back(Eigen::Vector3i(2, 3, 7));
		OBBmeshFs.push_back(Eigen::Vector3i(2, 7, 6));
		OBBmeshFs.push_back(Eigen::Vector3i(3, 0, 4));
		OBBmeshFs.push_back(Eigen::Vector3i(7, 3, 4));
		OBBmeshFs.push_back(Eigen::Vector3i(1, 2, 6));
		OBBmeshFs.push_back(Eigen::Vector3i(1, 6, 5));
		//cout << "n:" << n << endl;
        //scale OBBmeshPs
        Eigen::Vector3d AvgP(0, 0, 0);
        for (int i = 0; i < OBBmeshPs.size(); i++)
        {
			AvgP += OBBmeshPs[i];
        }
		AvgP /= OBBmeshPs.size();
        double scaleV = 1.35;
        for (int i = 0; i < OBBmeshPs.size(); i++)
        {
			OBBmeshPs[i] = (OBBmeshPs[i] - AvgP) * scaleV + AvgP;
        }
        
		

        MyBaseModel OBBmesh(OBBmeshPs, OBBmeshFs);
		//OBBmesh.WriteObjFile("..\\..\\data\\OBBmesh.obj");

        cout << "Read point cloud.\n";

        map<MyPoint, int> Point2ID;
        vector<vector<int>> neighboor;



        // filePath = "E:\\Dropbox\\SIG-2022-Feature-preserving-recon\\data\\ToRander\\";
        for (int i = 0; i < n; i++)
        {


            add_point(VersPC[i].x(), VersPC[i].y(), VersPC[i].z(), 0);
            Weight.push_back(0);


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



        /*ofstream outPts(filePath + "\\CenterPoints_" + modelName + ".xyz");

        outPts.precision(15);
        outPts.flags(ios::left | ios::fixed);
        outPts.fill('0');*/





        vector<vector<MyPoint>> VDpoints;
        vector<vector<double>> WindingNum;
        map<MyPoint,set<MyPoint>> neibors;
        int OpenmpCnt = 0;
        //omp_set_num_threads(16);
//#pragma omp parallel for //schedule(dynamic, 20)
        for (int ii = 0; ii < n; ii++)
        {
            // if (ii == 7)
                 //continue;

            /* if (Repeat[ii] != -1)
             {

                 continue;
             }*/
            vector<Eigen::Vector3d> VDPs_tp;
            OpenmpCnt++;
            if (OpenmpCnt % 1 == 0)
            {
                cout << OpenmpCnt <<" "<< neighboor[ii].size()<< '\n';
            }

            //outPts << VersPC[ii].x() << " " << VersPC[ii].y() << " " << VersPC[ii].z() << "\n";

            MyBaseModel VDmesh = OBBmesh;
            PlaneCut* PC = new PlaneCut(&VDmesh);
          //  PlaneCut* PC = new PlaneCut(VersPC[ii]);

            nowP = ii;
            //sort(KnnPC.neighboor[ii].begin(), KnnPC.neighboor[ii].end(),cmp);

            map<int, int> opposide;
            map<MyPoint, vector<int>> point2edge;
            // 把knn倒过来，应该会好。

            //for (auto p : KnnPC.neighboor[ii])
            int planeNum = 0;
            if (neighboor[ii].size() > 70)
            {
                //sort neighboor[ii] by distance
				vector<pair<double, int>> tmp;
				for (auto p : neighboor[ii])
				{
					tmp.push_back(make_pair((VersPC[ii] - VersPC[p]).norm(), p));
				}
				sort(tmp.begin(), tmp.end());
				neighboor[ii].clear();
				for (int ik=0;ik<30;ik++)
				{
					neighboor[ii].push_back(tmp[ik].second);
				}
                
            }

            for (int jj = 0; jj < neighboor[ii].size(); jj++)
            {
                //cout << " " << jj << " ";
                /*if (jj > 50)
                {
                    break;
                }*/

                auto p = neighboor[ii][jj];
                /*if (Repeat[p] != -1)
                    continue;*/



                    //cout << (VersPC[p] - VersPC[ii]).norm() << "\n";
                    //cout << p << endl;
                    // out << VersPC[p].transpose() << endl;

                auto aa = VersPC[ii];
                auto bb = VersPC[p];
                auto r1 = Weight[ii], r2 = Weight[p];

                //Eigen::Vector3d MidPoint = VersPC[ii] + (VersPC[p] - VersPC[ii])*(Weight[ii]/(Weight[ii]+Weight[p]));  //voronoi
                Eigen::Vector3d MidPoint;
                //double lambda = ((r1 * r1 - r2 * r2) / ((bb - aa).norm() * (bb - aa).norm()) + 1.0) / 2.0;
				MidPoint = (aa + bb) / 2.0;

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





            }
            /*string file = filePath + "\\cells\\" + modelName + "_" + to_string(ii) + "cells.obj";
            PC->CuttedMesh->WriteObjFile(file.c_str());*/
            auto VD = PC->CuttedMesh;
            auto VDPs = VD->GetVertices();
			auto VDfaces = VD->GetFaces();
			MyHalfEdgeModel VDHEM(VDPs, VDfaces);
            VDHEM.CreateMeshEdges();
			auto VDHPs = VDHEM.GetVertices();
			auto VDHEs = VDHEM.GetEdges();
			auto VDHFaces = VDHEM.GetFaces();

            vector<MyPoint> vdps_tmp;
            vector<double> windingNum_tmp;
            vector<double> diss_tmp;
            vector<bool> Is3c(VDHPs.size(), 0);
            for (int i =0;i< VDHPs.size();i++)
            {
                auto Fs = VDHEM.GetFacesByPoint(i);
                vector<Eigen::Vector3d> Nors;
                for (auto f : Fs)
                {
                    Eigen::Vector3d nor;
					nor = (VDHPs[VDHFaces[f][0]] - VDHPs[VDHFaces[f][1]]).cross(VDHPs[VDHFaces[f][0]] - VDHPs[VDHFaces[f][2]]);
					nor.normalize();
					Nors.push_back(nor);
                }
                double errorr = 0.0005;
				
                int k = Nors.size();
                vector<int> v1, cnt;
                for (int i = 0; i < k; i++)
                {
                    v1.push_back(i);
                    cnt.push_back(0);
                }
                for (int i = 0; i < k - 1; i++)
                {
                    for (int j = i + 1; j < k; j++)
                    {
                        double dij = min((Nors[i] - (Nors[j])).norm(), (Nors[j] - (Nors[i])).norm());
                        if (dij < errorr)
                        {
                            v1[j] = v1[i];
                        }
                    }
                }
                for (int i = 0; i < k; i++)
                {
                    cnt[v1[i]] = cnt[v1[i]] + 1;
                }
                
                int norNum = 0;
                for (int i = 0; i < k; i++)
                {
                    if (cnt[i] > 0)
                    {
                        norNum++;
                    }
                }
                //cout << norNum << endl;
                if (norNum > 2)
                {
                    //if ((VDHPs[i] - VersPC[ii]).norm() > 0.15)
                    {
                        Is3c[i] = 1;
                        vdps_tmp.push_back(VDHPs[i]);
                        windingNum_tmp.push_back(0);
                        VDPs_tp.push_back(VDHPs[i]);
                        MyPoint mp(VDHPs[i]);
                        if (VDPsMap.find(mp) == VDPsMap.end())
                        {
                            VDPsMap[mp] = 1;
                            //VDPs_tp.push_back(VDHPs[i]);
                            //F_VDPs.second.push_back((VDHPs[i] - VersPC[ii]).norm());
                        }
                    }
                }

            }
            for (int i = 0; i < VDHEs.size(); i++)
            {
                auto nowE = VDHEs[i];
                if (!Is3c[nowE.leftVert] && !Is3c[nowE.rightVert])
                {
                    continue;
                }
                if (Is3c[nowE.leftVert] && Is3c[nowE.rightVert])
                {
                    neibors[VDHPs[nowE.leftVert]].insert(VDHPs[nowE.rightVert]);
                    neibors[VDHPs[nowE.rightVert]].insert(VDHPs[nowE.leftVert]);
                    continue;
                }
                Eigen::Vector3d Nor(0,0,0);
                int Pfrom=-1, Pto = -1;
                if (Is3c[nowE.leftVert])
                {
                    Pfrom = nowE.leftVert;
                    Pto = nowE.rightVert;
                    Nor = VDHPs[nowE.rightVert] - VDHPs[nowE.leftVert];
                    Nor.normalize();
                }
                else
                {
                    Pto = nowE.leftVert;
                    Pfrom = nowE.rightVert;
                    Nor = VDHPs[nowE.leftVert] - VDHPs[nowE.rightVert];
                    Nor.normalize();
                }
                int Pstart = Pfrom;
                bool Final = 0;
                while (1)
                {
                    //VDHEM.GetEdgeByPoint(Pto);
                    bool Fd = 0;
                    auto Fs = VDHEM.GetFacesByPoint(Pto);
                    vector<Eigen::Vector3d> Nors;
                    for (auto f : Fs)
                    {
                        int p[3];
                        p[0] = VDHFaces[f][0];
                        p[1] = VDHFaces[f][1];
                        p[2] = VDHFaces[f][2];
                        for (int j = 0; j < 3; j++)
                        {
                            if (p[j] != Pto && p[j] != Pfrom)
                            {
                                Eigen::Vector3d NewNor = VDHPs[p[j]] - VDHPs[Pto];
								NewNor.normalize();
								if ((NewNor.cross(Nor)).norm() < 0.0001)
								{
									Pfrom = Pto;
									Pto = p[j];
									Nor = NewNor;
									Nor.normalize();
                                    if (Is3c[Pto])
                                    {
                                        Final = 1;
                                    }
                                    Fd = 1;
									break;
								}
                            }
                        }
                        if (Fd)
                        {
							break;
                        }
                    }
					if (!Fd)
					{
						break;
					}
                    if (Final)
                    {
                        break;
                    }

                }
                if (Final)
                {
                    neibors[VDHPs[Pstart]].insert(VDHPs[Pto]);
                    neibors[VDHPs[Pto]].insert(VDHPs[Pstart]);
                }
            }



            VDpoints.push_back(vdps_tmp);
            WindingNum.push_back(windingNum_tmp);

            F_VDPs.push_back(VDPs_tp);
            double maxDis = -99999;
            Eigen::Vector3d FastP(VersPC[ii]);
            for (auto Vp : VDPs_tp)
            {
                double dis = (Vp - VersPC[ii]).norm();
                if (dis > maxDis)
                {
                    maxDis = dis;
                    FastP = Vp;
                }
            }
            K::Point_3 point1(VersPC[ii].x(), VersPC[ii].y(), VersPC[ii].z());
            
            //Eigen::Vector3d dirr = Normal_ori[ii];
            Eigen::Vector3d dirr = FastP - VersPC[ii];
            dirr.normalize();
            
            K::Direction_3 dir(dirr.x(), dirr.y(), dirr.z());

            Plane plane(point1, dir);
            double area = PC->CutByPlaneGetCutArea(plane);
            //cout << area << endl;
            //double area = 1;
            Areas.push_back(area);
            PCs[ii] = PC;
            //continue;
            //cout << endl;
        }
        //缩小box
        
        

        //compute winding number
        //omp_set_num_threads(16);
        // //schedule(dynamic, 20)
        double PI = 3.1415926535;
        double area = 1.0;
        double totarea = 1;
        double maxwd = -9999999;
        double minwd = 9999999;
#pragma omp parallel for
        for (int ii = 0; ii < n; ii++)
        {
            for (int i = 0; i < VDpoints[ii].size(); i++)
            {
                auto queryP = VDpoints[ii][i];
                double wd = 0;
                for (int j = 0; j < n; j++)
                {
                    wd = wd + Areas[j] * (((VersPC[j] - queryP.p).dot(Normal_ori[j])) / (4 * PI * pow(((VersPC[j] - queryP.p).squaredNorm()), 1.5)));
                }
                //wd =  wd;
                WindingNum[ii][i] = wd;
                if (wd > maxwd)
                    maxwd = wd;
                if (wd < minwd)
                    minwd = wd;
				//cout << "Winding Num: " << wd << " " << queryP.p.x() << " " << queryP.p.y() << " " << queryP.p.z() << endl;
            }
        }
		cout << "maxwd:   " << maxwd << endl;
		cout << "minwd:   " << minwd << endl;
        /*ofstream outWwwd(filePath + "\\Out\\" + modelName + "WDDDWDWD.txt");
        int nnn = 500;
        for (int ii = 0; ii < nnn; ii++)
        {
            for (int jj = 0; jj < nnn; jj++)
            {
				Eigen::Vector3d QueryP(-0.6 + ii / 416.0, -0.6 + jj / 416.0, 0);
                double wd = 0;
                for (int j = 0; j < n; j++)
                {
                    wd = wd + Areas[j] * (((VersPC[j] - QueryP).dot(Normal_ori[j])) / (4 * PI * pow(((VersPC[j] - QueryP).squaredNorm()), 1.5)));
                }
				outWwwd << wd << " ";
            }
			outWwwd << endl;
        }
        outWwwd.close();*/





        maxwd = -0.5;
        minwd = 2.0;
		ofstream outWd(filePath +  "\\Out\\" + modelName + "_GTQueryPoints.txt");
        ofstream outWds(filePath + "\\Out\\" + modelName + "_GTWindingNums.txt");
        for (int ii = 0; ii < n; ii++)
        {
            for (int i = 0; i < VDpoints[ii].size(); i++)
            {
                auto queryP = VDpoints[ii][i];
                double dis = (queryP.p - VersPC[ii]).norm();
				//outWd << VDpoints[ii][i].p.x() << " "  << VDpoints[ii][i].p.y() << " " << VDpoints[ii][i].p.z() << " " << (WindingNum[ii][i]-minwd)/(maxwd-minwd)  << " 0 0\n";
                outWds << WindingNum[ii][i] << " " << dis <<"   " << ii << " " << i << endl;

            
                
                if (WindingNum[ii][i] > 1.1)
                {
                    outWd << VDpoints[ii][i].p.x() << " " << VDpoints[ii][i].p.y() << " " << VDpoints[ii][i].p.z() << " " << 0 << " 0 1\n";
                }
                else
                {
                    if (WindingNum[ii][i] < -0.1)
                    {
                        outWd << VDpoints[ii][i].p.x() << " " << VDpoints[ii][i].p.y() << " " << VDpoints[ii][i].p.z() << " " << 0 << " 1 0\n";
                    }
                    else
                    {
                        outWd << VDpoints[ii][i].p.x() << " " << VDpoints[ii][i].p.y() << " " << VDpoints[ii][i].p.z() << " " << (WindingNum[ii][i] - (-0.1)) / (1.1 - (-0.1)) << " 0 0\n";
                    }

                }
            }
        }
        return make_pair(F_VDPs, make_pair(Areas, neibors));
    }
    



    std::cout << "Hello World!\n";


}



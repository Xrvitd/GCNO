/*************************************************************************\

  Copyright 1999 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
  fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             E. Larsen
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:               geom@cs.unc.edu


\**************************************************************************/

//--------------------------------------------------------------------------
// File:   TriDist.cpp
// Author: Eric Larsen
// Description:
// contains SegPoints() for finding closest points on a pair of line
// segments and TriDist() for finding closest points on a pair of triangles
//--------------------------------------------------------------------------

#include "BGAL/PQP/MatVec.h"
#ifdef _WIN32
#include <float.h>
#define isnan _isnan
#endif

//--------------------------------------------------------------------------
// SegPoints() 
//
// Returns closest points between an segment pair.
// Implemented from an algorithm described in
//
// Vladimir J. Lumelsky,
// On fast computation of distance between line segments.
// In Information Processing Letters, no. 21, pages 55-61, 1985.   
//--------------------------------------------------------------------------

void
SegPoints(PQP_REAL VEC[3],
          PQP_REAL X[3], PQP_REAL Y[3],             // closest points
          const PQP_REAL P[3], const PQP_REAL A[3], // seg 1 origin, vector
          const PQP_REAL Q[3], const PQP_REAL B[3]) // seg 2 origin, vector
{
  PQP_REAL T[3], A_dot_A, B_dot_B, A_dot_B, A_dot_T, B_dot_T;
  PQP_REAL TMP[3];

  VmV(T, Q, P);
  A_dot_A = VdotV(A, A);
  B_dot_B = VdotV(B, B);
  A_dot_B = VdotV(A, B);
  A_dot_T = VdotV(A, T);
  B_dot_T = VdotV(B, T);

  // t parameterizes ray P,A 
  // u parameterizes ray Q,B 

  PQP_REAL t, u;

  // compute t for the closest point on ray P,A to
  // ray Q,B

  PQP_REAL denom = A_dot_A * B_dot_B - A_dot_B * A_dot_B;

  t = (A_dot_T * B_dot_B - B_dot_T * A_dot_B) / denom;

  // clamp result so t is on the segment P,A

  if ((t < 0) || isnan(t)) t = 0; else if (t > 1) t = 1;

  // find u for point on ray Q,B closest to point at t

  u = (t * A_dot_B - B_dot_T) / B_dot_B;

  // if u is on segment Q,B, t and u correspond to 
  // closest points, otherwise, clamp u, recompute and
  // clamp t 

  if ((u <= 0) || isnan(u)) {

    VcV(Y, Q);

    t = A_dot_T / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Q, P);
    } else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Q, X);
    } else {
      VpVxS(X, P, A, t);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  } else if (u >= 1) {

    VpV(Y, Q, B);

    t = (A_dot_B + A_dot_T) / A_dot_A;

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VmV(VEC, Y, P);
    } else if (t >= 1) {
      VpV(X, P, A);
      VmV(VEC, Y, X);
    } else {
      VpVxS(X, P, A, t);
      VmV(T, Y, P);
      VcrossV(TMP, T, A);
      VcrossV(VEC, A, TMP);
    }
  } else {

    VpVxS(Y, Q, B, u);

    if ((t <= 0) || isnan(t)) {
      VcV(X, P);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    } else if (t >= 1) {
      VpV(X, P, A);
      VmV(T, Q, X);
      VcrossV(TMP, T, B);
      VcrossV(VEC, B, TMP);
    } else {
      VpVxS(X, P, A, t);
      VcrossV(VEC, A, B);
      if (VdotV(VEC, T) < 0) {
        VxS(VEC, VEC, -1);
      }
    }
  }
}

//--------------------------------------------------------------------------
// TriDist() 
//
// Computes the closest points on two triangles, and returns the 
// distance between them.
// 
// S and T are the triangles, stored tri[point][dimension].
//
// If the triangles are disjoint, P and Q give the closest points of 
// S and T respectively. However, if the triangles overlap, P and Q 
// are basically a random pair of points from the triangles, not 
// coincident points on the intersection of the triangles, as might 
// be expected.
//--------------------------------------------------------------------------

PQP_REAL
TriDist(PQP_REAL P[3], PQP_REAL Q[3],
        const PQP_REAL S[3][3], const PQP_REAL T[3][3]) {
  // Compute vectors along the 6 sides

  PQP_REAL Sv[3][3], Tv[3][3];
  PQP_REAL VEC[3];

  VmV(Sv[0], S[1], S[0]);
  VmV(Sv[1], S[2], S[1]);
  VmV(Sv[2], S[0], S[2]);

  VmV(Tv[0], T[1], T[0]);
  VmV(Tv[1], T[2], T[1]);
  VmV(Tv[2], T[0], T[2]);

  // For each edge pair, the vector connecting the closest points 
  // of the edges defines a slab (parallel planes at head and tail
  // enclose the slab). If we can show that the off-edge vertex of 
  // each triangle is outside of the slab, then the closest points
  // of the edges are the closest points for the triangles.
  // Even if these tests fail, it may be helpful to know the closest
  // points found, and whether the triangles were shown disjoint

  PQP_REAL V[3];
  PQP_REAL Z[3];
  PQP_REAL minP[3], minQ[3], mindd;
  int shown_disjoint = 0;

  mindd = VdistV2(S[0], T[0]) + 1;  // Set first minimum safely high

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      // Find closest points on edges i & j, plus the 
      // vector (and distance squared) between these points

      SegPoints(VEC, P, Q, S[i], Sv[i], T[j], Tv[j]);

      VmV(V, Q, P);
      PQP_REAL dd = VdotV(V, V);

      // Verify this closest point pair only if the distance 
      // squared is less than the minimum found thus far.

      if (dd <= mindd) {
        VcV(minP, P);
        VcV(minQ, Q);
        mindd = dd;

        VmV(Z, S[(i + 2) % 3], P);
        PQP_REAL a = VdotV(Z, VEC);
        VmV(Z, T[(j + 2) % 3], Q);
        PQP_REAL b = VdotV(Z, VEC);

        if ((a <= 0) && (b >= 0)) return sqrt(dd);

        PQP_REAL p = VdotV(V, VEC);

        if (a < 0) a = 0;
        if (b > 0) b = 0;
        if ((p - a + b) > 0) shown_disjoint = 1;
      }
    }
  }

  // No edge pairs contained the closest points.  
  // either:
  // 1. one of the closest points is a vertex, and the
  //    other point is interior to a face.
  // 2. the triangles are overlapping.
  // 3. an edge of one triangle is parallel to the other's face. If
  //    cases 1 and 2 are not true, then the closest points from the 9
  //    edge pairs checks above can be taken as closest points for the
  //    triangles.
  // 4. possibly, the triangles were degenerate.  When the 
  //    triangle points are nearly colinear or coincident, one 
  //    of above tests might fail even though the edges tested
  //    contain the closest points.

  // First check for case 1

  PQP_REAL Sn[3], Snl;
  VcrossV(Sn, Sv[0], Sv[1]); // Compute normal to S triangle
  Snl = VdotV(Sn, Sn);      // Compute square of length of normal

  // If cross product is long enough,

  if (Snl > 1e-15) {
    // Get projection lengths of T points

    PQP_REAL Tp[3];

    VmV(V, S[0], T[0]);
    Tp[0] = VdotV(V, Sn);

    VmV(V, S[0], T[1]);
    Tp[1] = VdotV(V, Sn);

    VmV(V, S[0], T[2]);
    Tp[2] = VdotV(V, Sn);

    // If Sn is a separating direction,
    // find point with smallest projection

    int point = -1;
    if ((Tp[0] > 0) && (Tp[1] > 0) && (Tp[2] > 0)) {
      if (Tp[0] < Tp[1]) point = 0; else point = 1;
      if (Tp[2] < Tp[point]) point = 2;
    } else if ((Tp[0] < 0) && (Tp[1] < 0) && (Tp[2] < 0)) {
      if (Tp[0] > Tp[1]) point = 0; else point = 1;
      if (Tp[2] > Tp[point]) point = 2;
    }

    // If Sn is a separating direction, 

    if (point >= 0) {
      shown_disjoint = 1;

      // Test whether the point found, when projected onto the 
      // other triangle, lies within the face.

      VmV(V, T[point], S[0]);
      VcrossV(Z, Sn, Sv[0]);
      if (VdotV(V, Z) > 0) {
        VmV(V, T[point], S[1]);
        VcrossV(Z, Sn, Sv[1]);
        if (VdotV(V, Z) > 0) {
          VmV(V, T[point], S[2]);
          VcrossV(Z, Sn, Sv[2]);
          if (VdotV(V, Z) > 0) {
            // T[point] passed the test - it's a closest point for 
            // the T triangle; the other point is on the face of S

            VpVxS(P, T[point], Sn, Tp[point] / Snl);
            VcV(Q, T[point]);
            return sqrt(VdistV2(P, Q));
          }
        }
      }
    }
  }

  PQP_REAL Tn[3], Tnl;
  VcrossV(Tn, Tv[0], Tv[1]);
  Tnl = VdotV(Tn, Tn);

  if (Tnl > 1e-15) {
    PQP_REAL Sp[3];

    VmV(V, T[0], S[0]);
    Sp[0] = VdotV(V, Tn);

    VmV(V, T[0], S[1]);
    Sp[1] = VdotV(V, Tn);

    VmV(V, T[0], S[2]);
    Sp[2] = VdotV(V, Tn);

    int point = -1;
    if ((Sp[0] > 0) && (Sp[1] > 0) && (Sp[2] > 0)) {
      if (Sp[0] < Sp[1]) point = 0; else point = 1;
      if (Sp[2] < Sp[point]) point = 2;
    } else if ((Sp[0] < 0) && (Sp[1] < 0) && (Sp[2] < 0)) {
      if (Sp[0] > Sp[1]) point = 0; else point = 1;
      if (Sp[2] > Sp[point]) point = 2;
    }

    if (point >= 0) {
      shown_disjoint = 1;

      VmV(V, S[point], T[0]);
      VcrossV(Z, Tn, Tv[0]);
      if (VdotV(V, Z) > 0) {
        VmV(V, S[point], T[1]);
        VcrossV(Z, Tn, Tv[1]);
        if (VdotV(V, Z) > 0) {
          VmV(V, S[point], T[2]);
          VcrossV(Z, Tn, Tv[2]);
          if (VdotV(V, Z) > 0) {
            VcV(P, S[point]);
            VpVxS(Q, S[point], Tn, Sp[point] / Tnl);
            return sqrt(VdistV2(P, Q));
          }
        }
      }
    }
  }

  // Case 1 can't be shown.
  // If one of these tests showed the triangles disjoint,
  // we assume case 3 or 4, otherwise we conclude case 2, 
  // that the triangles overlap.

  if (shown_disjoint) {
    VcV(P, minP);
    VcV(Q, minQ);
    return sqrt(mindd);
  } else return 0;
}
bool intersectTriangle(float o[], float n[], float *t, float p0[], float p1[], float p2[]) {
  double p[3], e1[3], e2[3];

  e1[0] = p1[0] - p0[0];
  e1[1] = p1[1] - p0[1];
  e1[2] = p1[2] - p0[2];
  e2[0] = p2[0] - p0[0];
  e2[1] = p2[1] - p0[1];
  e2[2] = p2[2] - p0[2];

  p[0] = n[1] * e2[2] - n[2] * e2[1];
  p[1] = n[2] * e2[0] - n[0] * e2[2];
  p[2] = n[0] * e2[1] - n[1] * e2[0];

  double a = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];
  if (fabs(a) < 1.0e-10)
    return false;

  double f = 1.0 / a;
  double s[3];
  s[0] = o[0] - p0[0];
  s[1] = o[1] - p0[1];
  s[2] = o[2] - p0[2];
  double u = f * (s[0] * p[0] + s[1] * p[1] + s[2] * p[2]);
  if ((u < 0.0) || (u > 1.0))
    return false;

  double q[3];
  q[0] = s[1] * e1[2] - s[2] * e1[1];
  q[1] = s[2] * e1[0] - s[0] * e1[2];
  q[2] = s[0] * e1[1] - s[1] * e1[0];

  double v = f * (n[0] * q[0] + n[1] * q[1] + n[2] * q[2]);
  if ((v < 0.0) || ((u + v) > 1.0))
    return false;

  *t = f * (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]);
  return true;
}

bool project3D(float p[3], float q[3], float p0[3], float p1[3]) {
  //  make vector to project onto
  double v0[3];
  v0[0] = p1[0] - p0[1];
  v0[1] = p1[1] - p0[1];
  v0[2] = p1[2] - p0[2];

  //  make vector to project
  double v1[3];
  v1[0] = q[0] - p0[0];
  v1[1] = q[1] - p0[1];
  v1[2] = q[2] - p0[2];

  //  do projection
  double t = v1[0] * v0[0] + v1[1] * v0[1] + v1[2] * v0[2];

  //  test in segment
  if (t < 0) //|| t > len0)
  {
    return false;
  } else {
    double v2[3];
    v2[0] = q[0] - p1[0];
    v2[1] = q[1] - p1[1];
    v2[2] = q[2] - p1[2];
    double t1 = v2[0] * v0[0] + v2[1] * v0[1] + v2[2] * v0[2];
    if (t1 > 0)
      return false;
    {
      //  normalise the vector
      double len0 = (v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
      if (len0 < 1.0e-10)
        t = 0.0;
      else
        t /= len0;
      //  work out point
      p[0] = p0[0] + t * v0[0];
      p[1] = p0[1] + t * v0[1];
      p[2] = p0[2] + t * v0[2];
      return true;
    }
  }
}

double distanceToTriangle(int *pos, float q[], float p[], float v0[], float v1[], float v2[], float *n) {
  float vC[3];
  double e1[3], e2[3];
  e1[0] = v1[0] - v0[0];
  e1[1] = v1[1] - v0[1];
  e1[2] = v1[2] - v0[2];
  e2[0] = v2[0] - v0[0];
  e2[1] = v2[1] - v0[1];
  e2[2] = v2[2] - v0[2];

  //  normal of the triangle plane
  if (n != NULL) {
    vC[0] = n[0];
    vC[1] = n[1];
    vC[2] = n[2];
  } else {
    vC[0] = e1[1] * e2[2] - e1[2] * e2[1];
    vC[1] = e1[2] * e2[0] - e1[0] * e2[2];
    vC[2] = e1[0] * e2[1] - e1[1] * e2[0];

    //double len = sqrt(vC[0]*vC[0]+vC[1]*vC[1]+vC[2]*vC[2]);
    //vC[0]/=len;	vC[1]/=len;	vC[2]/=len;
  }

  double minD = 1.0E+10;

  //  test against triangle
  float t;
  if (intersectTriangle(p, vC, &t, v2, v1, v0)) {
    *pos = 0;
    q[0] = p[0] + t * vC[0];
    q[1] = p[1] + t * vC[1];
    q[2] = p[2] + t * vC[2];
    minD = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
  } else {
    float vv[3][3];
    vv[0][0] = v0[0];
    vv[0][1] = v0[1];
    vv[0][2] = v0[2];
    vv[1][0] = v1[0];
    vv[1][1] = v1[1];
    vv[1][2] = v1[2];
    vv[2][0] = v2[0];
    vv[2][1] = v2[1];
    vv[2][2] = v2[2];
    //  project p1 onto each edge
    int i;
    for (i = 0; i < 3; i++) {
      float pt[3];
      if (project3D(pt, p, vv[i], vv[(i + 1) % 3])) {
        double d = (pt[0] - p[0]) * (pt[0] - p[0]) + (pt[1] - p[1]) * (pt[1] - p[1]) + (pt[2] - p[2]) * (pt[2] - p[2]);
        if (d < minD) {
          minD = d;
          *pos = i + 1;
          q[0] = pt[0];
          q[1] = pt[1];
          q[2] = pt[2];
        }
      }
    }

    //  test against each vertex
    for (i = 0; i < 3; i++) {
      double d = (vv[i][0] - p[0]) * (vv[i][0] - p[0]) + (vv[i][1] - p[1]) * (vv[i][1] - p[1])
          + (vv[i][2] - p[2]) * (vv[i][2] - p[2]);
      if (d < minD) {
        *pos = i + 4;
        q[0] = vv[i][0];
        q[1] = vv[i][1];
        q[2] = vv[i][2];
        minD = d;
      }
    }
  }

  return sqrt(minD);
}

PQP_REAL
PointTriDist(int *posFlag, PQP_REAL q[3], const PQP_REAL p[3], const PQP_REAL tri[3][3]) {
  //float qq[3], pp[3], vv0[3], vv1[3], vv2[3];
  //pp[0] = p[0];	pp[1] = p[1];	pp[2] = p[2];
  //vv0[0] = tri[0][0];	vv0[1] = tri[0][1];	vv0[2] = tri[0][2];
  //vv1[0] = tri[1][0];	vv1[1] = tri[1][1];	vv1[2] = tri[1][2];
  //vv2[0] = tri[2][0];	vv2[1] = tri[2][1];	vv2[2] = tri[2][2];

  //double dis = distanceToTriangle(posFlag, qq, pp, vv0, vv1, vv2, NULL);
  //
  //q[0] = qq[0];	q[1] = qq[1];	q[2] = qq[2];

  //return dis;

  PQP_REAL v[3], e1[3], e2[3];
  v[0] = tri[0][0];
  v[1] = tri[0][1];
  v[2] = tri[0][2];
  e1[0] = tri[1][0] - tri[0][0];
  e1[1] = tri[1][1] - tri[0][1];
  e1[2] = tri[1][2] - tri[0][2];
  e2[0] = tri[2][0] - tri[0][0];
  e2[1] = tri[2][1] - tri[0][1];
  e2[2] = tri[2][2] - tri[0][2];

  PQP_REAL a, b, c, d, e, f;
  PQP_REAL vp[3];
  vp[0] = v[0] - p[0];
  vp[1] = v[1] - p[1];
  vp[2] = v[2] - p[2];

  a = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
  b = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
  c = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
  d = e1[0] * vp[0] + e1[1] * vp[1] + e1[2] * vp[2];
  e = e2[0] * vp[0] + e2[1] * vp[1] + e2[2] * vp[2];
  f = vp[0] * vp[0] + vp[1] * vp[1] + vp[2] * vp[2];

  PQP_REAL det, s, t;
  det = fabs(a * c - b * b);
  s = b * e - c * d;
  t = b * d - a * e;
  if (s + t <= det) {
    if (s < 0) {
      if (t < 0) {
        //region 4
        if (d < 0)    // minimum on edge t=0
        {
          t = 0;
          s = -d >= a ? 1 : -d / a;
        } else    //minimum on edge s=0
        {
          s = 0;
          t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));
        }
//				double tmp0=d;	double tmp1=e;
//				if(tmp1>tmp0)	// minimum on edge s+t=1
//				{
//					t = 0;
////					s = (tmp0<=0?1:(d>=0?0:-d/a));
//					s = (d>=0?0:(-d>=a?1:-d/a));
//				}
//				else	//minimum on edge s=0
//				{
//					s = 0;
////					t = (tmp1<=0?1:(e>=0?0:-e/c));
//					t = (e>=0?0:(-e>=c?1:-e/c));
//				}
      } else {
        //region 3
        s = 0.0;
        t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));
      }
    } else if (t < 0) {
      //region 5
      t = 0.0;
      s = (d >= 0 ? 0 : (-d >= a ? 1 : -d / a));
    } else {
      //region 0
      PQP_REAL invDet = 1.0 / det;
      s *= invDet;
      t *= invDet;
    }
  } else {
    if (s < 0) {
      //region 2
      //s = 0.0;
      //t = 1.0;
      PQP_REAL tmp0 = b + d;
      PQP_REAL tmp1 = c + e;
      if (tmp1 > tmp0)    // minimum on edge s+t=1
      {
        PQP_REAL numer = tmp1 - tmp0;
        PQP_REAL denom = fabs(a - 2.0 * b + c);
        s = (numer >= denom ? 1 : numer / denom);
        t = 1 - s;
      } else    //minimum on edge s=0
      {
        s = 0;
        t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e / c));
      }
    } else if (t < 0) {
      //region 6
      //s = 1.0;
      //t = 0.0;
      PQP_REAL tmp0 = b + e;
      PQP_REAL tmp1 = a + d;
      if (tmp1 > tmp0)    //	minimum on edge s+t=1
      {
        PQP_REAL numer = tmp1 - tmp0;
        PQP_REAL denom = fabs(a - 2.0 * b + c);
        t = (numer >= denom ? 1 : numer / denom);
        s = 1 - t;
      } else //minimum on edge t=0;
      {
        t = 0;
        s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d / a));
      }
    } else {
      //region 1
      PQP_REAL numer = c + e - b - d;
      if (numer <= 0) {
        s = 0.0;
      } else {
        PQP_REAL denom = fabs(a - 2.0 * b + c);
        s = (numer >= denom ? 1 : numer / denom);
      }
      t = 1.0 - s;
    }
  }

  q[0] = v[0] + s * e1[0] + t * e2[0];
  q[1] = v[1] + s * e1[1] + t * e2[1];
  q[2] = v[2] + s * e1[2] + t * e2[2];

  double dist = (p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]);
  dist = sqrt(dist);

  if (s == 0 && t == 0) {
    *posFlag = 4;
    return dist;
  }
  if (s == 1 && t == 0) {
    *posFlag = 5;
    return dist;
  }
  if (s == 0 && t == 1) {
    *posFlag = 6;
    return dist;
  }
  if (s < 1 && t == 0) {
    *posFlag = 1;
    return dist;
  }
  if (s + t == 1 && s < 1 && s > 0 && t < 1 && t > 0) {
    *posFlag = 2;
    return dist;
  }
  if (s == 0 && t < 1) {
    *posFlag = 3;
    return dist;
  }
  if (s + t < 1 && s > 0 && t > 0) {
    *posFlag = 0;
    return dist;
  }

  return dist;
}

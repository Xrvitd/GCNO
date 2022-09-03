#pragma once
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <time.h>
# include <cstring>
# include <vector>
using namespace std;

//****************************************************************************80

void comp_next(int n, int k, int a[], bool* more, int* h, int* t);
int i4_max(int i1, int i2);
int i4_min(int i1, int i2);
int i4_modp(int i, int j);
int i4_wrap(int ival, int ilo, int ihi);
int keast_degree(int rule);
int keast_suborder_num(int rule);
int keast_order_num(int rule);
void keast_rule(int rule, int order_num, double xyz[], double w[]);
int keast_rule_num(void);
int* keast_suborder(int rule, int suborder_num);

void keast_subrule(int rule, int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_01(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_02(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_03(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_04(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_05(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_06(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_07(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_08(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_09(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
void keast_subrule_10(int suborder_num, double suborder_xyzz[],
	double suborder_w[]);
double* monomial_value(int dim_num, int point_num, double x[], int expon[]);
double r8_huge(void);
int r8_nint(double x);
double r8mat_det_4d(double a[]);
double r8vec_dot(int n, double a1[], double a2[]);
int s_len_trim(char* s);
void tetrahedron_reference_to_physical(double t[], int n, double ref[],
	double phy[]);
double tetrahedron_volume(double t[3 * 4]);



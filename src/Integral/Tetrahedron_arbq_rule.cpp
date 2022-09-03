#include "BGAL/Integral/Tetrahedron_arbq_rule.h"

void comp_next(int n, int k, int a[], bool *more, int *h, int *t)

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if (!(*more))
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for (i = 1; i < k; i++)
    {
      a[i] = 0;
    }
  }
  else
  {
    if (1 < *t)
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h - 1];
    a[*h - 1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = (a[k - 1] != n);

  return;
}
//****************************************************************************80

int i4_max(int i1, int i2)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if (i2 < i1)
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min(int i1, int i2)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if (i1 < i2)
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_modp(int i, int j)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if (j == 0)
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit(1);
  }

  value = i % j;

  if (value < 0)
  {
    value = value + abs(j);
  }

  return value;
}
//****************************************************************************80*

int i4_wrap(int ival, int ilo, int ihi)

//****************************************************************************80*
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min(ilo, ihi);
  jhi = i4_max(ilo, ihi);

  wide = jhi + 1 - jlo;

  if (wide == 1)
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp(ival - jlo, wide);
  }

  return value;
}
//****************************************************************************80

int keast_degree(int rule)

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_DEGREE returns the degree of a Keast rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int KEAST_DEGREE, the polynomial degree of exactness of
//    the rule.
//
{
  int degree;

  if (rule == 1)
  {
    degree = 0;
  }
  else if (rule == 2)
  {
    degree = 1;
  }
  else if (rule == 3)
  {
    degree = 2;
  }
  else if (rule == 4)
  {
    degree = 3;
  }
  else if (rule == 5)
  {
    degree = 4;
  }
  else if (rule == 6)
  {
    degree = 4;
  }
  else if (rule == 7)
  {
    degree = 5;
  }
  else if (rule == 8)
  {
    degree = 6;
  }
  else if (rule == 9)
  {
    degree = 7;
  }
  else if (rule == 10)
  {
    degree = 8;
  }
  else
  {
    degree = -1;
    cerr << "\n";
    cerr << "KEAST_DEGREE - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit(1);
  }

  return degree;
}
//****************************************************************************80

int keast_order_num(int rule)

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_ORDER_NUM returns the order of a Keast rule for the triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int KEAST_ORDER_NUM, the order (number of points) of the rule.
//
{
  int order;
  int order_num;
  int *suborder;
  int suborder_num;

  suborder_num = keast_suborder_num(rule);

  suborder = keast_suborder(rule, suborder_num);

  order_num = 0;
  for (order = 0; order < suborder_num; order++)
  {
    order_num = order_num + suborder[order];
  }

  delete[] suborder;

  return order_num;
}
//****************************************************************************80

void keast_rule(int rule, int order_num, double xyz[], double w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_RULE returns the points and weights of a Keast rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int ORDER_NUM, the order (number of points) of the rule.
//
//    Output, double XYZ[3*ORDER_NUM], the points of the rule.
//
//    Output, double W[ORDER_NUM], the weights of the rule.
//
{
  int k;
  int o;
  int s;
  int *suborder;
  int suborder_num;
  double *suborder_w;
  double *suborder_xyzz;
  //
  //  Get the suborder information.
  //
  suborder_num = keast_suborder_num(rule);

  suborder_xyzz = new double[4 * suborder_num];
  suborder_w = new double[suborder_num];

  suborder = keast_suborder(rule, suborder_num);

  keast_subrule(rule, suborder_num, suborder_xyzz, suborder_w);
  //
  //  Expand the suborder information to a full order rule.
  //
  o = 0;

  for (s = 0; s < suborder_num; s++)
  {
    if (suborder[s] == 1)
    {
      xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
      xyz[1 + o * 3] = suborder_xyzz[1 + s * 4];
      xyz[2 + o * 3] = suborder_xyzz[2 + s * 4];
      w[o] = suborder_w[s];
      o = o + 1;
    }
    //
    //  For SUBORDER = 4, we list the coordinates of the generator as
    //
    //    A,B,B,B
    //
    //  and we generate
    //
    //    A, B, B = (1,2,3)
    //    B, B, B = (2,3,4)
    //    B, B, A = (3,4,1)
    //    B, A, B = (4,1,2)
    //
    else if (suborder[s] == 4)
    {
      for (k = 0; k < 4; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[i4_wrap(k, 0, 3) + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[i4_wrap(k + 1, 0, 3) + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[i4_wrap(k + 2, 0, 3) + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    //
    //  For SUBORDER = 6, we list the coordinates of the generator as
    //
    //    A,A,B,B
    //
    //  and we generate
    //
    //    B, A, A = (4,1,2)
    //    A, B, A = (1,4,2)
    //    A, A, B = (1,2,4)
    //
    //    A, B, B = (1,3,4)
    //    B, A, B = (4,2,3)
    //    B, B, A = (4,3,1)
    //
    else if (suborder[s] == 6)
    {
      for (k = 0; k < 3; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[k + o * 3] = suborder_xyzz[2 + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }

      for (k = 0; k < 3; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[2 + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[2 + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[2 + s * 4];
        xyz[k + o * 3] = suborder_xyzz[0 + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    //
    //  For SUBORDER = 12, we list the coordinates of the generator as
    //
    //    A,A,B,C
    //
    //  and we generate
    //
    //    B, A, A
    //    A, B, A
    //    A, A, B
    //
    //    C, A, A
    //    A, C, A
    //    A, A, C
    //
    //    A, B, C
    //    B, C, A
    //    C, A, B
    //    A, C, B
    //    C, B, A
    //    B, A, C
    //
    else if (suborder[s] == 12)
    {
      for (k = 0; k < 3; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[k + o * 3] = suborder_xyzz[2 + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }

      for (k = 0; k < 3; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[0 + s * 4];
        xyz[k + o * 3] = suborder_xyzz[3 + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }

      for (k = 0; k < 3; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[i4_wrap(k + 1, 1, 3) + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[i4_wrap(k + 2, 1, 3) + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[i4_wrap(k + 3, 1, 3) + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }
      for (k = 0; k < 3; k++)
      {
        xyz[0 + o * 3] = suborder_xyzz[i4_wrap(k + 1, 1, 3) + s * 4];
        xyz[1 + o * 3] = suborder_xyzz[i4_wrap(k + 3, 1, 3) + s * 4];
        xyz[2 + o * 3] = suborder_xyzz[i4_wrap(k + 2, 1, 3) + s * 4];
        w[o] = suborder_w[s];
        o = o + 1;
      }
    }
    else
    {
      cerr << "\n";
      cerr << "KEAST_RULE - Fatal error!\n;";
      cerr << "  Illegal SUBORDER(" << s << ") = " << suborder[s] << "\n";
      exit(1);
    }
  }

  delete[] suborder;
  delete[] suborder_xyzz;
  delete[] suborder_w;

  return;
}
//****************************************************************************80

int keast_rule_num(void)

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_RULE_NUM returns the number of Keast rules available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Output, int KEAST_RULE_NUM, the number of rules available.
//
{
  int rule_num;

  rule_num = 10;

  return rule_num;
}
//****************************************************************************80

int *keast_suborder(int rule, int suborder_num)

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBORDER returns the suborders for a Keast rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, int KEAST_SUBORDER[SUBORDER_NUM], the suborders of the rule.
//
{
  int *suborder;

  suborder = new int[suborder_num];

  if (rule == 1)
  {
    suborder[0] = 1;
  }
  else if (rule == 2)
  {
    suborder[0] = 4;
  }
  else if (rule == 3)
  {
    suborder[0] = 1;
    suborder[1] = 4;
  }
  else if (rule == 4)
  {
    suborder[0] = 4;
    suborder[1] = 6;
  }
  else if (rule == 5)
  {
    suborder[0] = 1;
    suborder[1] = 4;
    suborder[2] = 6;
  }
  else if (rule == 6)
  {
    suborder[0] = 6;
    suborder[1] = 4;
    suborder[2] = 4;
  }
  else if (rule == 7)
  {
    suborder[0] = 1;
    suborder[1] = 4;
    suborder[2] = 4;
    suborder[3] = 6;
  }
  else if (rule == 8)
  {
    suborder[0] = 4;
    suborder[1] = 4;
    suborder[2] = 4;
    suborder[3] = 12;
  }
  else if (rule == 9)
  {
    suborder[0] = 1;
    suborder[1] = 4;
    suborder[2] = 4;
    suborder[3] = 4;
    suborder[4] = 6;
    suborder[5] = 12;
  }
  else if (rule == 10)
  {
    suborder[0] = 1;
    suborder[1] = 4;
    suborder[2] = 4;
    suborder[3] = 6;
    suborder[4] = 6;
    suborder[5] = 12;
    suborder[6] = 12;
  }
  else
  {
    cerr << "\n";
    cerr << "KEAST_SUBORDER - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit(1);
  }

  return suborder;
}
//****************************************************************************80

int keast_suborder_num(int rule)

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBORDER_NUM returns the number of suborders for a Keast rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Output, int KEAST_SUBORDER_NUM, the number of suborders of the rule.
//
{
  int suborder_num;

  if (rule == 1)
  {
    suborder_num = 1;
  }
  else if (rule == 2)
  {
    suborder_num = 1;
  }
  else if (rule == 3)
  {
    suborder_num = 2;
  }
  else if (rule == 4)
  {
    suborder_num = 2;
  }
  else if (rule == 5)
  {
    suborder_num = 3;
  }
  else if (rule == 6)
  {
    suborder_num = 3;
  }
  else if (rule == 7)
  {
    suborder_num = 4;
  }
  else if (rule == 8)
  {
    suborder_num = 4;
  }
  else if (rule == 9)
  {
    suborder_num = 6;
  }
  else if (rule == 10)
  {
    suborder_num = 7;
  }
  else
  {
    suborder_num = -1;
    cerr << "\n";
    cerr << "KEAST_SUBORDER_NUM - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit(1);
  }

  return suborder_num;
}
//****************************************************************************80

void keast_subrule(int rule, int suborder_num, double suborder_xyzz[],
                   double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE returns a compressed Keast rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;

  if (rule == 1)
  {
    keast_subrule_01(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 2)
  {
    keast_subrule_02(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 3)
  {
    keast_subrule_03(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 4)
  {
    keast_subrule_04(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 5)
  {
    keast_subrule_05(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 6)
  {
    keast_subrule_06(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 7)
  {
    keast_subrule_07(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 8)
  {
    keast_subrule_08(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 9)
  {
    keast_subrule_09(suborder_num, suborder_xyzz, suborder_w);
  }
  else if (rule == 10)
  {
    keast_subrule_10(suborder_num, suborder_xyzz, suborder_w);
  }
  else
  {
    cerr << "\n";
    cerr << "KEAST_SUBRULE - Fatal error!\n";
    cerr << "  Illegal RULE = " << rule << "\n";
    exit(1);
  }
  //
  //  Renormalize the weights so they sum to 1.
  //
  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = 6.0 * suborder_w[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_01(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_01 returns a compressed Keast rule 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348..
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_01[4 * 1] = {
      0.250000000000000000, 0.250000000000000000,
      0.250000000000000000, 0.250000000000000000};
  double suborder_w_rule_01[1] = {
      0.166666666666666667};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_01[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_01[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_01[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_01[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_01[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_02(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_02 returns a compressed Keast rule 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_02[4 * 1] = {
      0.585410196624968500, 0.138196601125010500,
      0.138196601125010500, 0.138196601125010500};
  double suborder_w_rule_02[1] = {
      0.0416666666666666667};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_02[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_02[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_02[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_02[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_02[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_03(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_03 returns a compressed Keast rule 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_03[4 * 5] = {
      0.250000000000000000, 0.250000000000000000,
      0.250000000000000000, 0.250000000000000000,
      0.500000000000000000, 0.166666666666666667,
      0.166666666666666667, 0.166666666666666667};
  double suborder_w_rule_03[5] = {
      -0.133333333333333333,
      0.075000000000000000};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_03[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_03[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_03[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_03[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_03[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_04(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_04 returns a compressed Keast rule 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_04[4 * 2] = {
      0.568430584196844400, 0.143856471934385200,
      0.143856471934385200, 0.143856471934385200,
      0.500000000000000000, 0.500000000000000000,
      0.000000000000000000, 0.000000000000000000};
  double suborder_w_rule_04[2] = {
      0.0362941783134009000,
      0.00358165890217718333};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_04[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_04[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_04[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_04[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_04[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_05(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_05 returns a compressed Keast rule 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_05[4 * 3] = {
      0.250000000000000000, 0.250000000000000000,
      0.250000000000000000, 0.250000000000000000,
      0.785714285714285714, 0.0714285714285714285,
      0.0714285714285714285, 0.0714285714285714285,
      0.399403576166799219, 0.399403576166799219,
      0.100596423833200785, 0.100596423833200785};
  double suborder_w_rule_05[3] = {
      -0.0131555555555555556,
      0.00762222222222222222,
      0.0248888888888888889};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_05[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_05[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_05[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_05[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_05[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_06(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_06 returns a compressed Keast rule 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_06[4 * 3] = {
      0.500000000000000000, 0.500000000000000000,
      0.000000000000000000, 0.000000000000000000,
      0.698419704324386603, 0.100526765225204467,
      0.100526765225204467, 0.100526765225204467,
      0.0568813795204234229, 0.314372873493192195,
      0.314372873493192195, 0.314372873493192195};
  double suborder_w_rule_06[3] = {
      0.00317460317460317450,
      0.0147649707904967828,
      0.0221397911142651221};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_06[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_06[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_06[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_06[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_06[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_07(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_07 returns a compressed Keast rule 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_07[4 * 4] = {
      0.250000000000000000, 0.250000000000000000,
      0.250000000000000000, 0.250000000000000000,
      0.00000000000000000, 0.333333333333333333,
      0.333333333333333333, 0.333333333333333333,
      0.727272727272727273, 0.0909090909090909091,
      0.0909090909090909091, 0.0909090909090909091,
      0.0665501535736642813, 0.0665501535736642813,
      0.433449846426335728, 0.433449846426335728};
  double suborder_w_rule_07[4] = {
      0.0302836780970891856,
      0.00602678571428571597,
      0.0116452490860289742,
      0.0109491415613864534};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_07[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_07[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_07[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_07[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_07[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_08(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_08 returns a compressed Keast rule 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_08[4 * 4] = {
      0.356191386222544953, 0.214602871259151684,
      0.214602871259151684, 0.214602871259151684,
      0.877978124396165982, 0.0406739585346113397,
      0.0406739585346113397, 0.0406739585346113397,
      0.0329863295731730594, 0.322337890142275646,
      0.322337890142275646, 0.322337890142275646,
      0.0636610018750175299, 0.0636610018750175299,
      0.269672331458315867, 0.603005664791649076};
  double suborder_w_rule_08[4] = {
      0.00665379170969464506,
      0.00167953517588677620,
      0.00922619692394239843,
      0.00803571428571428248};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_08[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_08[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_08[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_08[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_08[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_09(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_08 returns a compressed Keast rule 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_09[4 * 6] = {
      0.250000000000000000, 0.250000000000000000,
      0.250000000000000000, 0.250000000000000000,
      0.765360423009044044, 0.0782131923303186549,
      0.0782131923303186549, 0.0782131923303186549,
      0.634470350008286765, 0.121843216663904411,
      0.121843216663904411, 0.121843216663904411,
      0.00238250666073834549, 0.332539164446420554,
      0.332539164446420554, 0.332539164446420554,
      0.500000000000000000, 0.500000000000000000,
      0.00000000000000000, 0.00000000000000000,
      0.100000000000000000, 0.100000000000000000,
      0.200000000000000000, 0.600000000000000000};
  double suborder_w_rule_09[6] = {
      0.0182642234661087939,
      0.0105999415244141609,
      -0.0625177401143299494,
      0.00489142526307353653,
      0.000970017636684296702,
      0.0275573192239850917};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_09[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_09[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_09[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_09[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_09[s];
  }

  return;
}
//****************************************************************************80

void keast_subrule_10(int suborder_num, double suborder_xyzz[],
                      double suborder_w[])

//****************************************************************************80
//
//  Purpose:
//
//    KEAST_SUBRULE_10 returns a compressed Keast rule 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Patrick Keast,
//    Moderate Degree Tetrahedral Quadrature Formulas,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 55, Number 3, May 1986, pages 339-348.
//
//  Parameters:
//
//    Input, int SUBORDER_NUM, the number of suborders of the rule.
//
//    Output, double SUBORDER_XYZZ[4*SUBORDER_NUM],
//    the barycentric coordinates of the abscissas.
//
//    Output, double SUBORDER_W[SUBORDER_NUM], the suborder weights.
//
{
  int s;
  double suborder_xyzz_rule_10[4 * 7] = {
      0.250000000000000000, 0.250000000000000000,
      0.250000000000000000, 0.250000000000000000,
      0.617587190300082967, 0.127470936566639015,
      0.127470936566639015, 0.127470936566639015,
      0.903763508822103123, 0.0320788303926322960,
      0.0320788303926322960, 0.0320788303926322960,
      0.0497770956432810185, 0.0497770956432810185,
      0.450222904356718978, 0.450222904356718978,
      0.183730447398549945, 0.183730447398549945,
      0.316269552601450060, 0.316269552601450060,
      0.231901089397150906, 0.231901089397150906,
      0.0229177878448171174, 0.513280033360881072,
      0.0379700484718286102, 0.0379700484718286102,
      0.730313427807538396, 0.193746475248804382};
  double suborder_w_rule_10[7] = {
      -0.0393270066412926145,
      0.00408131605934270525,
      0.000658086773304341943,
      0.00438425882512284693,
      0.0138300638425098166,
      0.00424043742468372453,
      0.00223873973961420164};

  for (s = 0; s < suborder_num; s++)
  {
    suborder_xyzz[0 + s * 4] = suborder_xyzz_rule_10[0 + s * 4];
    suborder_xyzz[1 + s * 4] = suborder_xyzz_rule_10[1 + s * 4];
    suborder_xyzz[2 + s * 4] = suborder_xyzz_rule_10[2 + s * 4];
    suborder_xyzz[3 + s * 4] = suborder_xyzz_rule_10[3 + s * 4];
  }

  for (s = 0; s < suborder_num; s++)
  {
    suborder_w[s] = suborder_w_rule_10[s];
  }

  return;
}
//****************************************************************************80

double *monomial_value(int dim_num, int point_num, double x[], int expon[])

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= dim <= dim_num ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points at which the
//    monomial is to be evaluated.
//
//    Input, double X[DIM_NUM*POINT_NUM], the point coordinates.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Output, double MONOMIAL_VALUE[POINT_NUM], the value of the monomial.
//
{
  int dim;
  int point;
  double *value;

  value = new double[point_num];

  for (point = 0; point < point_num; point++)
  {
    value[point] = 1.0;
  }

  for (dim = 0; dim < dim_num; dim++)
  {
    if (0 != expon[dim])
    {
      for (point = 0; point < point_num; point++)
      {
        value[point] = value[point] * pow(x[dim + point * dim_num], expon[dim]);
      }
    }
  }

  return value;
}
//****************************************************************************80

double r8_huge(void)

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number,
//    and is usually defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  return HUGE_VAL;
}
//****************************************************************************80

int r8_nint(double x)

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int s;
  int value;

  if (x < 0.0)
  {
    s = -1;
  }
  else
  {
    s = 1;
  }
  value = s * (int)(fabs(x) + 0.5);

  return value;
}
//****************************************************************************80

double r8mat_det_4d(double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0 + 0 * 4] * (a[1 + 1 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4]) - a[1 + 2 * 4] * (a[2 + 1 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 1 * 4]) + a[1 + 3 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4])) - a[0 + 1 * 4] * (a[1 + 0 * 4] * (a[2 + 2 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 2 * 4]) - a[1 + 2 * 4] * (a[2 + 0 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 0 * 4]) + a[1 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 0 * 4])) + a[0 + 2 * 4] * (a[1 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 1 * 4]) - a[1 + 1 * 4] * (a[2 + 0 * 4] * a[3 + 3 * 4] - a[2 + 3 * 4] * a[3 + 0 * 4]) + a[1 + 3 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4])) - a[0 + 3 * 4] * (a[1 + 0 * 4] * (a[2 + 1 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 1 * 4]) - a[1 + 1 * 4] * (a[2 + 0 * 4] * a[3 + 2 * 4] - a[2 + 2 * 4] * a[3 + 0 * 4]) + a[1 + 2 * 4] * (a[2 + 0 * 4] * a[3 + 1 * 4] - a[2 + 1 * 4] * a[3 + 0 * 4]));

  return det;
}
//****************************************************************************80

double r8vec_dot(int n, double a1[], double a2[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for (i = 0; i < n; i++)
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

int s_len_trim(char *s)

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char *t;

  n = strlen(s);
  t = s + strlen(s) - 1;

  while (0 < n)
  {
    if (*t != ' ')
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

void tetrahedron_reference_to_physical(double t[], int n,
                                       double ref[], double phy[])

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_REFERENCE_TO_PHYSICAL maps reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 4 physical tetrahedron and a point
//    (R,S,T) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y,Z) in physical space.
//
//    This routine will also be correct for an order 10 tetrahedron,
//    if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image tetrahedron are straight, the faces are flat, and
//    the "midside" nodes in the physical tetrahedron are
//    halfway along the edges of the physical tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[3*4], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0,0), (1,0,0),
//    (0,1,0) and (0,0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[3*N], points in the reference tetrahedron
//
//    Output, double PHY[3*N], corresponding points in the
//    physical tetrahedron.
//
{
  int i;
  int j;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < n; j++)
    {
      phy[i + j * 3] = t[i + 0 * 3] * (1.0 - ref[0 + j * 3] - ref[1 + j * 3] - ref[2 + j * 3]) + t[i + 1 * 3] * +ref[0 + j * 3] + t[i + 2 * 3] * +ref[1 + j * 3] + t[i + 3 * 3] * +ref[2 + j * 3];
    }
  }

  return;
}
//****************************************************************************80

double tetrahedron_volume(double tetra[3 * 4])

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
//
{
  double a[4 * 4];
  int i;
  int j;
  double volume;

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 4; j++)
    {
      a[i + j * 4] = tetra[i + j * 3];
    }
  }

  i = 3;
  for (j = 0; j < 4; j++)
  {
    a[i + j * 4] = 1.0;
  }

  volume = fabs(r8mat_det_4d(a)) / 6.0;

  return volume;
}
//****************************************************************************80

//double Tetra_Integration_API(const Geex::vec3 &X0, const Geex::vec3 &C1, const Geex::vec3 &C2, const Geex::vec3 &C3, int rule, std::vector<Geex::vec4> & Points)
//{
//	double node_xyz[3 * 4] = {
//		0.0, 0.0, 0.0,
//		1.0, 0.0, 0.0,
//		0.0, 1.0, 0.0,
//		0.0, 0.0, 1.0 };
//	double node_xyz2[12];
//	node_xyz2[0] = X0.x;
//	node_xyz2[1] = X0.y;
//	node_xyz2[2] = X0.z;
//	node_xyz2[3] = C1.x;
//	node_xyz2[4] = C1.y;
//	node_xyz2[5] = C1.z;
//	node_xyz2[6] = C2.x;
//	node_xyz2[7] = C2.y;
//	node_xyz2[8] = C2.z;
//	node_xyz2[9] = C3.x;
//	node_xyz2[10] = C3.y;
//	node_xyz2[11] = C3.z;
//
//
//	double *w;
//	double *xyz;
//	double *xyz2;
//
//	int order_num = keast_order_num(rule);
//	xyz = new double[3 * order_num];
//	xyz2 = new double[3 * order_num];
//	w = new double[order_num];
//
//	keast_rule(rule, order_num, xyz, w);
//
//	tetrahedron_reference_to_physical(node_xyz2, order_num, xyz, xyz2);
//
//	double volume = tetrahedron_volume(node_xyz2);//��������������
//
//	for (int order = 0; order < order_num; order++)
//	{
//		Points.push_back(Geex::vec4(xyz2[0 + order * 3], xyz2[1 + order * 3], xyz2[2 + order * 3], w[order]));
//	}
//	delete[] w;
//	delete[] xyz;
//	delete[] xyz2;
//
//	return volume;
//}
//
//double Tri_Integration_API(const Geex::vec3 &C1, const Geex::vec3 &C2, const Geex::vec3 &C3, std::vector<Geex::vec4> & Points)
//{
//	Geex::vec3  u = C3 - C1;
//	Geex::vec3  v = C2 - C1;
//
//	Geex::vec3 p1 = C1 + u / 2 + v / 2;
//	Points.push_back(Geex::vec4(p1.x, p1.y, p1.z, 1.0 / 30));
//	Geex::vec3 p2 = C1 + u / 2;
//	Points.push_back(Geex::vec4(p2.x, p2.y, p2.z, 1.0 / 30));
//	Geex::vec3 p3 = C1 + v / 2;
//	Points.push_back(Geex::vec4(p3.x, p3.y, p3.z, 1.0 / 30));
//	Geex::vec3 p4 = C1 + u / 6 + 2 * v / 3;
//	Points.push_back(Geex::vec4(p4.x, p4.y, p4.z, 3.0 / 10));
//	Geex::vec3 p5 = C1 + v / 6 + 2 * u / 3;
//	Points.push_back(Geex::vec4(p5.x, p5.y, p5.z, 3.0 / 10));
//	Geex::vec3 p6 = C1 + v / 6 + u / 6;
//	Points.push_back(Geex::vec4(p6.x, p6.y, p6.z, 3.0 / 10));
//	double area = (cross(u, v).length()) / 2;
//	return area;
//}
#include "mat3d.h"
#include "eig3.h"

// This routine calls a public domain routine (eigen_decomposition), which was 
// downloaded from:
// http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
// According to the header file of eig3.h, the C version in eig3.{h,cpp} 
// was obtained by copying code from the public domain Java library JAMA, 
// and declaring the resulting C code to also be public domain.
// The eig3.{h,cpp} files are included in this package, intact as they were
// downloaded from the Internet.

// This routine written by Jernej Barbic
void eigen_sym(Mat3d & a, Vec3d & eig_val, Vec3d eig_vec[3])
{
  double A[3][3] = { {a[0][0], a[0][1], a[0][2]},
                     {a[1][0], a[1][1], a[1][2]},
                     {a[2][0], a[2][1], a[2][2]} };
  double V[3][3];
  double d[3];
  eigen_decomposition(A, V, d);

  eig_val = Vec3d(d[2],d[1],d[0]);
  eig_vec[0] = Vec3d(V[0][2], V[1][2], V[2][2]);
  eig_vec[1] = Vec3d(V[0][1], V[1][1], V[2][1]);
  eig_vec[2] = Vec3d(V[0][0], V[1][0], V[2][0]);
}


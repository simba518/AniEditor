#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <BasisUpdating.h>
#include <MatrixTools.h>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace LSW_ANI_EDITOR;

BOOST_AUTO_TEST_SUITE(BasisUpdatingTest)

BOOST_AUTO_TEST_CASE(udpateBTest){
  
  const int n = 10*3;
  const int r = 5;
  const MatrixXd uk = MatrixXd::Random(n,r);
  MatrixXd B = MatrixXd::Identity(n,r);
  const MatrixXd Irxr = MatrixXd::Identity(r,r);
  const MatrixXd BtB_old = B.transpose()*B;
  ASSERT_EQ (BtB_old,Irxr);
  
  BasisUpdating::updateB(uk,B,1e-10);
  ASSERT_GT(B.cols(),r);
  const MatrixXd BtB_new = B.transpose()*B;
  const MatrixXd Irxr_new = MatrixXd::Identity(B.cols(),B.cols());
  ASSERT_EQ_SMALL_MAT_TOL (BtB_new,Irxr_new,1e-12);
}

BOOST_AUTO_TEST_CASE(udpateWTest){
  
  const int n = 10*3;
  const int r = 5;
  const MatrixXd pk = MatrixXd::Random(n,r);

  DiagonalMatrix<double,-1> M(n);
  for (int i = 0; i < n; ++i)
    M.diagonal()[i] = i+1.0f;

  MatrixXd W = MatrixXd::Identity(n,r);
  EIGEN3EXT::MGramSchmidt(M,W);
  const MatrixXd Irxr_old = MatrixXd::Identity(W.cols(),W.cols());
  const MatrixXd Wt_M_W_old = W.transpose()*M*W;
  ASSERT_EQ_SMALL_MAT_TOL(Wt_M_W_old,Irxr_old,1e-12);

  BasisUpdating::updateW(pk,M,W,1e-3);
  ASSERT_GT(W.cols(),r);
  const MatrixXd Irxr_new = MatrixXd::Identity(W.cols(),W.cols());
  const MatrixXd Wt_M_W_new = W.transpose()*M*W;
  ASSERT_EQ_SMALL_MAT_TOL (Wt_M_W_new,Irxr_new,1e-12);
}

BOOST_AUTO_TEST_SUITE_END()

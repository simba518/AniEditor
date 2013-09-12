#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <DefGradOperator.h>
#include <NodeRotVec.h>
#include <NodeWarper.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(NodeWarperTest)

BOOST_AUTO_TEST_CASE(testSkewMat){

  Vector3d w;
  w[0] = 1.0f;
  w[1] = 2.0f;
  w[2] = 5.0f;
  const Matrix3d wx = NodeWarper::SkewMat(w);
  ASSERT_EQ(wx,-wx.transpose());
  
  Matrix3d M;
  M <<  0, -5,  2,
	5,  0, -1,
	-2, 1,  0;
  ASSERT_EQ(wx,M);
}

BOOST_AUTO_TEST_CASE(testWarpMat){

  Vector3d w;
  w[0] = 1.0f;
  w[1] = 2.0f;
  w[2] = 5.0f;
  w *= 1e-6;
  const Matrix3d _R = NodeWarper::warpMat(w);
  ASSERT_EQ_SMALL_MAT_TOL (_R, Matrix3d::Identity(3,3),1e-5);

  w *= 1e6;
  const Matrix3d _R2 = NodeWarper::warpMat(w);
  Matrix3d _R3;
  _R3 << -0.0940022,0.0241849,0.209126,
	0.126712,0.0191704,0.366989,
	0.168116,0.387495,0.811379;
  ASSERT_EQ_SMALL_MAT_TOL (_R2, _R3,1e-6);
}

BOOST_AUTO_TEST_CASE(testWarpMatAll){
	  
  VectorV3 w;
  w.push_back(Vector3d(1,2,3));
  w.push_back(Vector3d(10,20,30));
  w.push_back(Vector3d(11,12,13));
  
  SparseMatrix<double> warpedR;
  NodeWarper::warpMat(w,warpedR);

  ASSERT_EQ(warpedR.rows(),(int)w.size()*3);
  ASSERT_EQ(warpedR.cols(),(int)w.size()*3);
  
  const MatrixXd mR = warpedR;
  ASSERT_EQ_SMALL_MAT((mR.block<3,3>)(0,0),NodeWarper::warpMat(w[0]));
  ASSERT_EQ_SMALL_MAT((mR.block<3,3>)(3,3),NodeWarper::warpMat(w[1]));
  ASSERT_EQ_SMALL_MAT((mR.block<3,3>)(6,6),NodeWarper::warpMat(w[2]));

  const VectorXd p = VectorXd::Random(w.size()*3);
  const VectorXd u = warpedR*p;
  ASSERT_EQ_SMALL_VEC_TOL(u.segment(0,3),(NodeWarper::warpMat(w[0])*p.segment(0,3)),3,1e-16);
  ASSERT_EQ_SMALL_VEC_TOL(u.segment(3,3),(NodeWarper::warpMat(w[1])*p.segment(3,3)),3,1e-16);
  ASSERT_EQ_SMALL_VEC_TOL(u.segment(6,3),(NodeWarper::warpMat(w[2])*p.segment(6,3)),3,1e-16);
}

BOOST_AUTO_TEST_CASE(testRotMat){

  // test one rot mat
  const Vector3d w0 = Vector3d::Zero(3);
  const Matrix3d R0 = NodeWarper::rotMat(w0);
  ASSERT_EQ_SMALL_MAT(R0,Matrix3d::Identity(3,3));

  const Vector3d w1 = Vector3d::Random(3);
  const Matrix3d R1 = NodeWarper::rotMat(w1);
  ASSERT_EQ_SMALL_MAT_TOL ((R1*R1.transpose()),Matrix3d::Identity(3,3),1e-8);
  ASSERT_EQ_TOL (R1.determinant(),1,1e-8);

  // test several ws.
  VectorV3 w;
  w.push_back(w0);
  w.push_back(w1);
  SparseMatrix<double> sparseR;
  NodeWarper::rotMat(w,sparseR);
  ASSERT_EQ(sparseR.rows(),6);
  ASSERT_EQ(sparseR.cols(),6);
  ASSERT_EQ(sparseR.nonZeros(),18);
  const MatrixXd R = sparseR;
  ASSERT_EQ_SMALL_MAT (((R.block<3,3>)(0,0)), R0);
  ASSERT_EQ_SMALL_MAT (((R.block<3,3>)(3,3)), R1);
  const MatrixXd RRt = R*R.transpose();
  const MatrixXd I = MatrixXd::Identity(6,6);
  ASSERT_EQ_SMALL_MAT_TOL (RRt,I,1e-8);
  ASSERT_EQ_TOL (R.determinant(),1.0f,1e-8);
}

BOOST_AUTO_TEST_SUITE_END()

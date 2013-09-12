#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <MapMA2RS.h>
#include <RSCoordComp.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(MapMA2RSTest)

BOOST_AUTO_TEST_CASE(testCalculateP){
	  
  const int n = 6; // number of vertices of the tet mesh.
  const int N = 3; // number of tetrahedrons.

  // construct a random G.
  const int rows = 9*N;
  const int cols = 3*n;
  const  MatrixXd MG = MatrixXd::Random(rows, cols);
  SparseMatrix<double> G(rows, cols);
  for (int i = 0; i < rows; ++i){
    for (int j = 0; j < cols; ++j){
	  G.insert(i,j) = MG(i,j);
	}
  }

  // construct a random u.
  const VectorXd u = VectorXd::Random(n*3);

  // compute the RS coordinates using RSCoordComp and MapMA2RS respectively
  VectorXd y1;
  RSCoordComp::constructWithWarp(G,u,y1);
  
  SparseMatrix<double> P;
  MapMA2RS::computeMapMatP(G.rows()/9,P);
  const VectorXd y2 = P*(G*u);

  // compare the results.
  ASSERT_EQ(y1.size(), y2.size());
  ASSERT_EQ_SMALL_VEC(y1,y2,y1.size());
}

BOOST_AUTO_TEST_CASE(testCalculateP9x1){
	  
  // assemble a 9x9 matrix.
  SparseMatrix<double> P;
  MapMA2RS::computeMapMatP(1,P);

  const VectorXd m = VectorXd::Random(9);
  const VectorXd y = P * m;
  VectorXd correct_y(9);

  correct_y[0] = 0.5*(m[7]-m[5]);
  correct_y[1] = 0.5*(m[2]-m[6]);
  correct_y[2] = 0.5*(m[3]-m[1]);

  correct_y[3] = m[0];
  correct_y[4] = 0.5*(m[1]+m[3]);
  correct_y[5] = 0.5*(m[2]+m[6]);

  correct_y[6] = m[4];
  correct_y[7] = 0.5*(m[7]+m[5]);
  correct_y[8] = m[8];

  ASSERT_EQ_SMALL_VEC(y, correct_y, 9);
}

BOOST_AUTO_TEST_SUITE_END()

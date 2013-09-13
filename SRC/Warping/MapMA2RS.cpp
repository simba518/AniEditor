#include <vector>
#include "MapMA2RS.h"
#include <assertext.h>
using namespace LSW_WARPING;
using namespace std;

void MapMA2RS::computeMapMatPGW(const SparseMatrix<double> &G, const MatrixXd &W, MatrixXd &PGW){

  assert_eq (G.rows()%9,0);
  assert_gt (G.rows(),0);
  assert_gt (G.cols(),0);

  SparseMatrix<double> P;
  const int N = G.rows()/9; // number of tetrahedrons.
  computeMapMatP(N,P);

  assert_eq (G.cols(), W.rows());
  assert_eq (P.rows(), 9*N);
  assert_eq (P.cols(), 9*N);

  PGW = P*G*W;
}
	
void MapMA2RS::computeMapMatP(const int N, SparseMatrix<double> &P){
  
  typedef Eigen::Triplet<double> T;
  const int nonzeros = 15*N;

  vector<T> P_triplets;
  P_triplets.reserve(nonzeros);
  
  // for each tetrahedron.
  for (int i = 0; i < N; ++i){

	const int c = 9*i; // position of up-left corner for each tetrahedron.
	P_triplets.push_back( T(c+0, c+5,-0.5) );
	P_triplets.push_back( T(c+0, c+7, 0.5) );
	P_triplets.push_back( T(c+1, c+2, 0.5) );
	P_triplets.push_back( T(c+1, c+6,-0.5) );
	P_triplets.push_back( T(c+2, c+1,-0.5) );

	P_triplets.push_back( T(c+2, c+3, 0.5) );
	P_triplets.push_back( T(c+3, c+0, 1.0) );
	P_triplets.push_back( T(c+4, c+1, 0.5) );
	P_triplets.push_back( T(c+4, c+3, 0.5) );
	P_triplets.push_back( T(c+5, c+2, 0.5) );

	P_triplets.push_back( T(c+5, c+6, 0.5) );
	P_triplets.push_back( T(c+6, c+4, 1.0) );
	P_triplets.push_back( T(c+7, c+5, 0.5) );
	P_triplets.push_back( T(c+7, c+7, 0.5) );
	P_triplets.push_back( T(c+8, c+8, 1.0) );
  }

  const int rows = 9*N;
  P.resize( rows, rows );
  P.reserve( nonzeros );
  P.setFromTriplets( P_triplets.begin(),P_triplets.end() );
  P.makeCompressed();
}

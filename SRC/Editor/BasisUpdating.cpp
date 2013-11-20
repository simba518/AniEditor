#include <BasisUpdating.h>
using namespace LSW_ANI_EDITOR;

// update basis
void BasisUpdating::update(const MatrixXd&Uk,MatrixXd&B,MatrixXd&W,MatrixXd&subK,bool useOptSK){
  updateB(Uk,B,tolForUpdateB);
  update(Uk,W,subK,useOptSK);
}

void BasisUpdating::update(const MatrixXd&Uk,MatrixXd &W,MatrixXd &subK,bool useOptSubK){

  MatrixXd Pk;
  computeLinearP(Uk,Pk);
  updateW(Pk,M,W,tolForUpdateW);
  if(useOptSubK){
	const MatrixXd tmpK = subK;
	subK = W.transpose()*K*W;
	subK.topLeftCorner(tmpK.rows(),tmpK.cols()) = tmpK;
  }else{
	subK = W.transpose()*K*W;
  }
}

// compute reduced keyframes
void BasisUpdating::computeZk(const MatrixXd&Uk,const MatrixXd &B,const MatrixXd &W,MatrixXd &Zk){
  const MatrixXd U = B*(B.transpose()*Uk);
  computeZk(U,W,Zk);
}

void BasisUpdating::computeZk(const MatrixXd&Uk,const MatrixXd &W,MatrixXd &Zk){
  MatrixXd Pk;
  computeLinearP(Uk,Pk);
  Eigen::LDLT<MatrixXd> solver(W.transpose()*W);
  Zk = solver.solve(W.transpose()*Pk);
}

int BasisUpdating::updateB(const MatrixXd&Uk,MatrixXd &B,const double tolB){

  assert_eq(Uk.rows(),B.rows());
  const int oldCols = B.cols();
  MatrixXd Bext(B.rows(),B.cols()+Uk.cols());
  Bext.leftCols(B.cols()) = B;
  int totalBasis = B.cols();
  for (int i = 0; i < Uk.cols(); ++i){
	const double u_rms = RMS_norm(Uk.col(i));
	if (u_rms > 0.0f){
	  const int c = totalBasis;
	  Bext.col(c) = Uk.col(i);
	  for (int j = 0; j < c; ++j)
		Bext.col(c) -= (Uk.col(i).dot(Bext.col(j)))*Bext.col(j);
	  if(RMS_norm(Bext.col(c))/u_rms >= tolB){
		Bext.col(c).normalize();
		totalBasis ++;
	  }
	}
  }
  if (totalBasis > oldCols)
	B = Bext.leftCols(totalBasis);
  return (totalBasis-oldCols);
}

int BasisUpdating::updateW(const MatrixXd&Pk,const DiagonalMatrix<double,-1>&M,
							MatrixXd &W,const double tolW){

  assert_eq(Pk.rows(),W.rows());
  assert_eq(Pk.rows(),M.rows());
  const int oldCols = W.cols();
  MatrixXd Wext(W.rows(),W.cols()+Pk.cols());
  Wext.leftCols(W.cols()) = W;
  int totalBasis = W.cols();

  for (int i = 0; i < Pk.cols(); ++i){
	const double p_rms = RMS_norm2_M(Pk.col(i),M);
	if(p_rms > 0.0f){
	  const int c = totalBasis;
	  Wext.col(c) = Pk.col(i);
	  for (int j = 0; j < c; ++j)
		Wext.col(c) -= (Pk.col(i).dot(M*Wext.col(j)))*Wext.col(j);
	  const double alpha = Wext.col(c).dot(M*Wext.col(c));
	  if (alpha/(Wext.rows()*p_rms) > tolW){
		Wext.col(c) /= sqrt(alpha);
		totalBasis ++;
	  }
	}
  }
  if(totalBasis > oldCols)
	W = Wext.leftCols(totalBasis);
  return (totalBasis-oldCols);
}

void BasisUpdating::precompute(){

  // PG
  SparseMatrix<double> P;
  DefGradOperator::compute(tetmesh,G);
  MapMA2RS::computeMapMatP(G.rows()/9,P);
  PG = P*G;

  // K,M
  ElasticForceTetFullStVK ela(tetmesh);
  VectorXd x0;
  tetmesh->nodes(x0);
  K = ela.K(x0)*(-1.0f);
  MassMatrix mass;
  mass.compute(M,*tetmesh);

  // rm_fixed_nodes
  // genReshapeMatrix(tetmesh->nodes().size()*3,3,fixed_nodes,rm_fixed_nodes);
  /// @todo remove fixed nodes
  // K = rm_fixed_nodes*(K*rm_fixed_nodes.transpose());
  // M = rm_fixed_nodes*(diagM*rm_fixed_nodes.transpose());

  // initialize solver for computing linear displacement
  SparseMatrix<double> con_Matrix;
  computeConM(fixed_nodes,con_Matrix,tetmesh->nodes().size());
  vector<Triplet<double> > tri;
  const SparseMatrix<double> L = PG.transpose()*PG;
  const SparseMatrix<double> ct = con_Matrix.transpose();
  EIGEN3EXT::addToBlock(tri,L,0,0);
  EIGEN3EXT::addToBlock(tri,con_Matrix,L.rows(),0);
  EIGEN3EXT::addToBlock(tri,ct,0,L.cols());
  A.resize(L.rows()+con_Matrix.rows(),L.cols()+con_Matrix.rows());
  A.setFromTriplets(tri.begin(),tri.end());
  A.makeCompressed();
  solver.compute(A);
}

void BasisUpdating::computeLinearP(const MatrixXd &U,MatrixXd &LinearP)const{

  VectorXd b(A.rows());
  b.setZero();
  LinearP.resize(U.rows(),U.cols());
  VectorXd y;
  for (int k = 0; k < U.cols(); ++k){
	const VectorXd &u = U.col(k);
	RSCoordComp::constructWithoutWarp(G,u,y);
	b.head(PG.cols()) = PG.transpose()*y;
	LinearP.col(k) = solver.solve(b).head(LinearP.rows());
  }
}

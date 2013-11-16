#include <eigen3/Eigen/Dense>
#include <MatrixIO.h>
#include <assertext.h>
#include <Euler2ReducedCoord.h>
#include <ConMatrixTools.h>
#include <RS2Euler.h>
#include <SparseMatrixTools.h>
#include <ElasticForceTetFullStVK.h>
#include <MassMatrix.h>
#include <MatrixTools.h>
#include <SparseGenEigenSolver.h>
#include <string>
#include <set>
using namespace std;
using namespace Eigen;
using namespace UTILITY;
using namespace LSW_WARPING;
using namespace EIGEN3EXT;

const string d = "/home/simba/Workspace/AnimationEditor/Data/beam/";
SparseMatrix<double> P,G,PG,con_Matrix, rm_fixed_nodes;
SparseMatrix<double> K,M;
pTetMesh tetmesh;
MatrixXd W,U,zk,LinearU,warpedU;
VectorXd Lambda;
RS2Euler rs2euler;
const int old_red_dim = 3;

/// init U
void loadU(){

  TRACE_FUN();
  VectorXi keyf(5);
  keyf << 0,20,40*2,50*2,70*2;
  MatrixXd all_U;
  bool succ = EIGEN3EXT::load(d+"W80_B40_C97_FixCen/U400_swing_right.b",all_U);
  assert(succ);
  U.resize(all_U.rows(),keyf.size());
  for (size_t i = 0; i < keyf.size(); ++i){
	const int j = keyf[i];
	if (j < 0)
	  U.col(i) = -1.0f*all_U.col(-j);
	else
	  U.col(i) = all_U.col(j);
  }
}

/// init tetmesh, lambda, G and hatW.
void initModelMtl(){
  
  TRACE_FUN();

  tetmesh = pTetMesh(new TetMesh);
  bool succ = tetmesh->load(d+"/mesh/beam.abq"); assert(succ);

  succ = DefGradOperator::compute(tetmesh,G); assert(succ);
  MapMA2RS::computeMapMatP(G.rows()/9,P);
  PG = P*G;
  
  vector<int> fixed_nodes;  
  succ = loadVec(d+"/W80_B40_C97_FixEnd/con_nodes.bou",fixed_nodes,UTILITY::TEXT); 
  assert(succ);
  computeConM(fixed_nodes,con_Matrix,tetmesh->nodes().size());

  set<int> fixed_nodes_set;
  for (int i = 0; i < fixed_nodes.size(); ++i)
    fixed_nodes_set.insert(fixed_nodes[i]);
  genReshapeMatrix(tetmesh->nodes().size()*3,3,fixed_nodes_set,rm_fixed_nodes);

  rs2euler.setTetMesh(tetmesh);
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();
}

/// compute K,M,W, and Lambda
void computeKMW(){

  TRACE_FUN();
  
  tetmesh->setSingleMaterial(1000.0f,2e6,0.45);
  ElasticForceTetFullStVK ela(tetmesh);
  VectorXd x0;
  tetmesh->nodes(x0);
  K = ela.K(x0)*(-1.0f);
  MassMatrix mass;
  DiagonalMatrix<double,-1> diagM;
  mass.compute(diagM,*tetmesh);

  K = rm_fixed_nodes*(K*rm_fixed_nodes.transpose());
  const SparseMatrix<double> KLower = EIGEN3EXT::getLower(K);
  M = rm_fixed_nodes*(diagM*rm_fixed_nodes.transpose());

  bool succ = EIGEN3EXT::EigenSparseGenEigenSolver::solve(KLower,M,W,Lambda,old_red_dim);
  assert(succ);
}

/// compute LinearU
void computeLinearU(){
  
  TRACE_FUN();

  vector<Triplet<double> > tri;
  const SparseMatrix<double> L = PG.transpose()*PG;
  const SparseMatrix<double> ct = con_Matrix.transpose();
  EIGEN3EXT::addToBlock(tri,L,0,0);
  EIGEN3EXT::addToBlock(tri,con_Matrix,L.rows(),0);
  EIGEN3EXT::addToBlock(tri,ct,0,L.cols());
  SparseMatrix<double> A(L.rows()+con_Matrix.rows(),L.cols()+con_Matrix.rows());
  A.setFromTriplets(tri.begin(),tri.end());
  A.makeCompressed();
  SparseQR<SparseMatrix<double>,COLAMDOrdering<int> > solver(A);
  solver.compute(A);

  VectorXd b(A.rows());
  b.setZero();
  LinearU.resize(U.rows(),U.cols());
  for (int k = 0; k < U.cols(); ++k){
    const VectorXd &u = U.col(k);
	VectorXd y;
	RSCoordComp::constructWithoutWarp(G,u,y);
	b.head(PG.cols()) = PG.transpose()*y;
	LinearU.col(k) = solver.solve(b).head(LinearU.rows());
  }
  
  LinearU = rm_fixed_nodes*LinearU;
}

/// compute zk, and new W.
void compute(){

  TRACE_FUN();

  // expand W using mass-PCA
  const MatrixXd tmpW1 = W;
  W.resize(tmpW1.rows(),tmpW1.cols()+LinearU.cols()-1);
  W.leftCols(tmpW1.cols()) = tmpW1;
  W.rightCols(LinearU.cols()-1) = LinearU.rightCols(LinearU.cols()-1);
  EIGEN3EXT::MGramSchmidt(M,W);

  // update W, lambda
  const MatrixXd subK = W.transpose()*K*W;
  SelfAdjointEigenSolver<MatrixXd> eigenK(subK);
  W *= eigenK.eigenvectors();
  Lambda = eigenK.eigenvalues();

  // const MatrixXd tmpW2 = W;
  // W.resize(tmpW1.rows(),tmpW1.cols()+1);
  // W = tmpW2.leftCols(W.cols());

  // insert zeros to W,LinearU
  W = rm_fixed_nodes.transpose()*W;
  LinearU = rm_fixed_nodes.transpose()*LinearU;

  Eigen::LDLT<MatrixXd> solver(W.transpose()*W);
  zk = solver.solve(W.transpose()*LinearU);
}

/// compute warpedU
void computeWarpedU(){
  
  TRACE_FUN();

  warpedU.resize(LinearU.rows(),LinearU.cols());
  for (int k = 0; k < LinearU.cols(); ++k){
    const VectorXd &p = LinearU.col(k);
  	VectorXd y;
  	RSCoordComp::constructWithWarp(G,p,y);
	VectorXd u;
	rs2euler.reconstruct(y,u);
	warpedU.col(k) = u;
  }
}

/// save results
void save(){

  TRACE_FUN();
  bool succ = EIGEN3EXT::write(d+"KeyW/key_z.b",zk); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW/lambda.b",Lambda); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW/W.b",W); assert(succ);
  succ = tetmesh->writeVTK(d+"KeyW/tempt/old_key_u",U); assert(succ);
  succ = tetmesh->writeVTK(d+"KeyW/tempt/old_key_warpedU",warpedU); assert(succ);
  succ = tetmesh->writeVTK(d+"KeyW/tempt/old_key_linearU",LinearU); assert(succ);
  const MatrixXd lw = (W*2000.0f);
  succ = tetmesh->writeVTK(d+"KeyW/tempt/new_key_linearW",lw); assert(succ);
}

int main(int argc, char *argv[]){

  loadU();
  initModelMtl();
  computeKMW();
  computeLinearU();
  compute();
  computeWarpedU();
  save();
  return 0;
}

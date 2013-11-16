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
MatrixXd B,W,unwarpB, W_ma;
VectorXd Lambda;
RS2Euler rs2euler;
const int r = 30;

/// init tetmesh, lambda, G and U
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

  succ = load(d+"/W80_B40_C97_FixEnd/nlmode80.b",B); assert(succ);
}

/// compute K,M, W_ma, remove fixed nodes
void computeKM(){

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

  const int r = B.cols();
  bool succ = EIGEN3EXT::EigenSparseGenEigenSolver::solve(KLower,M,W_ma,Lambda,r);
  W_ma = rm_fixed_nodes.transpose()*W_ma;
  assert(succ);
}

/// compute unwarpB, remove fixed nodes.
void computeUnwarpB(){
  
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
  unwarpB.resize(B.rows(),B.cols());
  for (int k = 0; k < B.cols(); ++k){
    const VectorXd &u = B.col(k);
	VectorXd y;
	RSCoordComp::constructWithoutWarp(G,u,y);
	b.head(PG.cols()) = PG.transpose()*y;
	unwarpB.col(k) = solver.solve(b).head(unwarpB.rows());
  }
  
  unwarpB = rm_fixed_nodes*unwarpB;
}

/// compute W, Lambda using B
void computeW(){

  TRACE_FUN();
  
  W = unwarpB;
  EIGEN3EXT::MGramSchmidt(M,W);

  // update W, lambda
  const MatrixXd subK = W.transpose()*K*W;
  SelfAdjointEigenSolver<MatrixXd> eigenK(subK);
  W *= eigenK.eigenvectors();
  Lambda = eigenK.eigenvalues();

  // insert zeros to W,LinearU
  W = rm_fixed_nodes.transpose()*W;
  unwarpB = rm_fixed_nodes.transpose()*unwarpB;
}

void save(){
  
  TRACE_FUN();
  bool succ = EIGEN3EXT::write(d+"RedRS/lambda.b",Lambda); assert(succ);
  succ = EIGEN3EXT::write(d+"RedRS/W.b",W); assert(succ);
  succ = EIGEN3EXT::write(d+"RedRS/B.b",B); assert(succ);

  B *= 10.0f;
  succ = tetmesh->writeVTK(d+"RedRS/tempt/B",B); assert(succ);
  unwarpB *= 10.0f;
  succ = tetmesh->writeVTK(d+"RedRS/tempt/unwarpB",unwarpB); assert(succ);
  W *= 1000.0f;
  succ = tetmesh->writeVTK(d+"RedRS/tempt/W",W); assert(succ);

  W_ma *= 1000.0f;
  succ = tetmesh->writeVTK(d+"RedRS/tempt/W_ma",W_ma); assert(succ);
}

int main(int argc, char *argv[]){

  initModelMtl();
  computeKM();
  computeUnwarpB();
  computeW();
  save();
  return 0;
}

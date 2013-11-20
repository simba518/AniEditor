#include <eigen3/Eigen/Dense>
#include <MatrixIO.h>
#include <AuxTools.h>
#include <BasisUpdating.h>
using namespace LSW_ANI_EDITOR;

const string d = "/home/simba/Workspace/AnimationEditor/Data/beam/";
pTetMesh tetmesh;
MatrixXd B,W,U,zk;
MatrixXd un_update_B,un_update_W,un_update_zk;
VectorXd Lambda,un_update_Lambda;
pBasisUpdating basisUpdate;

/// init U
void loadU(){

  TRACE_FUN();
  VectorXi keyf(5);
  keyf << 0,20,40*2,50*2,70*2;
  MatrixXd all_U;
  bool succ = EIGEN3EXT::load(d+"W80_B40_C97_FixCen/U400_swing_right.b",all_U);

  const int rm[] = {0,1,2,3,4,5,6,7,8,9,17,18,19,20,21,22,23,30,31,32,33,34,35,36,43,44,45,46,60,61,62,64,65,67,68,69,70,72,75,76,77,78,79,80,81,82,83,84,85,92,93,94,95,96,97,98,105,106,107,108,109,110,111,112,119,120,121,122,123,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,157,165,166,167,168,169,179,180,181,182,183,185,186,189,190,191,192,193,194,195,196,197,198,199,200,201,202,214,216,217,218,219,220,221,222,223,224,225,239,240,241,243,244,245,246,247,248,249,250,251,252,253,264,265,266,267,268,273,278,279,287,288,291,293,302,303,304,305,306,307,308,309,312,313,314,316,317,318,319,320,322,323,324,325,330,331,332,333,334,335,336,339,340,343,344,345,347,348,349,351,353,354,356,357};
  const int len = 189;
  for (int i = 0; i < len; ++i){
	all_U.row(rm[i]*3+0).setZero();
	all_U.row(rm[i]*3+1).setZero();
	all_U.row(rm[i]*3+2).setZero();
  }

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

/// init tetmesh, basisUpdate, B, W.
void initModelMtl(){
  
  TRACE_FUN();

  tetmesh = pTetMesh(new TetMesh);
  bool succ = tetmesh->load(d+"/mesh/beam.abq"); assert(succ);
  tetmesh->setSingleMaterial(1000.0f,2e6,0.45);

  set<int> fixed_nodes;  
  succ = load(d+"/W80_B40_C97_FixEnd/con_nodes.bou",fixed_nodes,UTILITY::TEXT);
  assert(succ);

  basisUpdate = pBasisUpdating(new BasisUpdating(tetmesh,fixed_nodes,1e-3,1e-3));

  succ = EIGEN3EXT::load(d+"RedRS_new/W3.b",W); assert(succ);
  succ = EIGEN3EXT::load(d+"RedRS_new/B20.b",B); assert(succ);
  succ = EIGEN3EXT::load(d+"RedRS_new/lambda3.b",un_update_Lambda); assert(succ);
  cout << "un update Lambda: " << un_update_Lambda.transpose() << endl;

  un_update_B = B;
  un_update_W = W;
}

/// compute B,W,lambda,zk
void updateBasis(){

  TRACE_FUN();

  MatrixXd subK;
  basisUpdate->update(U,B,W,subK,false);
  SelfAdjointEigenSolver<MatrixXd> eigenK(subK);
  W *= eigenK.eigenvectors();
  Lambda = eigenK.eigenvalues();
  cout << "Lambda: " << Lambda.transpose() << endl;

  basisUpdate->computeZk(U,B,W,zk);
  basisUpdate->computeZk(U,un_update_B,un_update_W,un_update_zk);
}

/// save results
void save(){

  const string r = TOSTR(un_update_W.cols())+"_"+TOSTR(W.cols())+".b";

  TRACE_FUN();
  bool succ = EIGEN3EXT::write(d+"KeyW_new/key_z"+r,zk); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW_new/lambda"+r,Lambda); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW_new/W"+r,W); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW_new/B"+r,B); assert(succ);

  succ = EIGEN3EXT::write(d+"KeyW_new/un_update_zk"+r,un_update_zk); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW_new/un_update_Lambda"+r,un_update_Lambda); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW_new/un_update_W"+r,un_update_W); assert(succ);
  succ = EIGEN3EXT::write(d+"KeyW_new/un_update_B"+r,un_update_B); assert(succ);

  succ = tetmesh->writeVTK(d+"KeyW_new/tempt/old_key_u"+r,U); assert(succ);
  const MatrixXd lw = (W*2000.0f);
  succ = tetmesh->writeVTK(d+"KeyW_new/tempt/new_key_linearW"+r,lw); assert(succ);
}

int main(int argc, char *argv[]){

  loadU();
  initModelMtl();
  updateBasis();
  save();
  return 0;
}

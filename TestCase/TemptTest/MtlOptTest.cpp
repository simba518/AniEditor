#include <iostream>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <AuxTools.h>
#include <MatrixIO.h>
#include <MASimulatorAD.h>
#include <TetMeshEmbeding.h>
#include <RS2Euler.h>
#include <RSCoordComp.h>
#include <DefGradOperator.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include "MtlOpt.h"
using namespace std;
using namespace EIGEN3EXT;
using namespace UTILITY;
using namespace LSW_WARPING;

typedef struct _MtlOptModel{

  _MtlOptModel(){

	dataDir = "/home/simba/Workspace/AnimationEditor/Data/beam/";
	T = 200;
	h = 0.1f;
	alphaK = 0.001f;
	alphaM = 0.001f;
	Kid.resize(5);
	Kid << 0,50,100,150,199;
	loadLambda(dataDir+"eigen_values80.b");
	z0 = VectorXd::Zero(lambda.size());
	z0[0] = 100;
	z0[6] = 100;
	initVolObj();
  }
  void initVolObj(){

	const string evect = dataDir+"/eigen_vectors80.b";
	TEST_ASSERT( load(evect,W) );

	TEST_ASSERT( volobj.loadObjMesh(dataDir+"beam.obj"));
	TEST_ASSERT( volobj.loadTetMesh(dataDir+"beam.abq"));
	TEST_ASSERT( volobj.loadWeights(dataDir+"beamW.txt"));
	vector<int> fixed_nodes;
	TEST_ASSERT( loadVec(dataDir+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
	rs2euler.setTetMesh(volobj.getTetMesh());
	rs2euler.setFixedNodes(fixed_nodes);
	rs2euler.precompute();
	TEST_ASSERT( DefGradOperator::compute(volobj.getTetMesh(),G) );
	MapMA2RS::computeMapMatPGW(G,W,PGW);	
  }
  bool loadLambda(const string fname){
	const bool succ = load(fname,lambda);
	assert(succ);
	return succ;
  }
  void produceSimRlst(){

	const int r = lambda.size();
	LSW_ANI_EDITOR::MASimulatorAD sim;
	sim.setTimeStep(h);
	sim.setEigenValues(lambda);
	sim.setStiffnessDamping(alphaK);
	sim.setMassDamping(alphaM);
	const VectorXd zero = VectorXd::Zero(r);
	sim.setIntialStatus(zero,z0);

	Z.resize(r,T);
	for (int i = 0; i < T; ++i){
	  Z.col(i) = sim.getEigenZ();
	  sim.forward(zero);
	}

	Kz.resize(r,Kid.size());
	for (int i = 0; i < Kid.size(); ++i){
	  const int f = Kid[i];
	  assert_in( f,0,Z.cols() );
	  Kz.col(i) = Z.col(f);
	}
  }
  void initMtlOpt(RedSpaceTimeEnergyAD &ad){
	ad.setT(T);
	ad.setTimestep(h);
	ad.setDamping(alphaK,alphaM);
	ad.setK(lambda);
	ad.setKeyframes(Kz,Kid);
  }
  void initSolver(const SXMatrix &E, const VSX &x){
	fun = CasADi::SXFunction(x,E);
	solver = CasADi::IpoptSolver(fun);
	solver.setOption("generate_hessian",true);
	solver.init();
	vector<double> x0(x.size(),0);
	solver.setInput(x0,CasADi::NLP_X_INIT);
  }
  void solve(){
	solver.solve();
	double cost;
	solver.getOutput(cost,CasADi::NLP_COST);
	cout << "optimal cost: " << cost << endl;
  }
  void getZfromSolver(MatrixXd &Z){
	
  }
  void saveRlst(){
	
	saveMesh(Z,"input");
	saveMesh(Kz,"key");
	MatrixXd newZ;
	getZfromSolver(newZ);
	saveMesh(newZ,"input");
  }
  void saveMesh(const MatrixXd &Z,const string fname){

	for (int i = 0; i < Z.cols(); ++i){
	  const string ffname=dataDir+"tempt/meshes/"+fname+TOSTR(i)+".vtk";
	  saveMeshOneZ(Z.col(i),ffname);
	}
  }
  void saveMeshOneZ(const VectorXd &z,const string fname){

	VectorXd y,u;
	const VectorXd p = W*z;
	RSCoordComp::constructWithWarp(G,p,y);
	TEST_ASSERT ( rs2euler.reconstruct(y,u) );
	volobj.interpolate(u);
	TEST_ASSERT (LSW_SIM::MeshVtkIO::write(volobj.getObjMesh(),fname));
  }

  string dataDir;
  MatrixXd W;
  SparseMatrix<double> G;
  MatrixXd PGW;
  TetMeshEmbeding volobj;
  RS2Euler rs2euler;

  int T;
  double h;
  double alphaK;
  double alphaM;
  VectorXd lambda;
  VectorXd z0;
  MatrixXd Z;
  MatrixXd Kz;
  VectorXi Kid;

  CasADi::SXFunction fun;
  CasADi::IpoptSolver solver;

}MtlOptModel;

BOOST_AUTO_TEST_SUITE(MtlOptTest)

BOOST_AUTO_TEST_CASE(checkRedSpaceTimeEnergyAD){
  
  RedSpaceTimeEnergyAD energy;
  VectorXd lambda(2);
  lambda << 1,2;
  energy.setT(6);
  energy.setTimestep(0.1);
  energy.setDamping(0.2,0.3);
  energy.setK(lambda);
  VectorXi Kid(2);
  Kid << 0,5;
  MatrixXd Kz(2,2);
  Kz << 1,1,2,2;
  energy.setKeyframes(Kz,Kid);
  energy.assembleEnergy();

  const VSX &z = energy.getVarZ();
  ASSERT_EQ (z.size(),8);
  ASSERT_EQ (energy.reducedDim(),2);
}

BOOST_AUTO_TEST_CASE(OptZ){
 
  MtlOptModel model;
  model.produceSimRlst();
  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);
  ad.assembleEnergy();
  model.initSolver(ad.getEnergy(),ad.getVarZ());
  model.solve();
  model.saveRlst();
}

BOOST_AUTO_TEST_SUITE_END()

#include "MtlOptModel.h"
#include "HarmonicOscillator.h"
#include <JsonFilePaser.h>
using namespace UTILITY;

MtlOptModel::_MtlOptModel(){

  dataDir = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string jsf = dataDir+"/mtlopt.ini";
  JsonFilePaser jsonf;
  TEST_ASSERT ( jsonf.open(jsf) );
  TEST_ASSERT ( jsonf.read("T",T) );
  TEST_ASSERT ( jsonf.read("h",h) );
  TEST_ASSERT ( jsonf.read("alphaK",alphaK) );
  TEST_ASSERT ( jsonf.read("alphaM",alphaM) );
  vector<int> keyid, modes;
  TEST_ASSERT ( jsonf.read("keyId", keyid) );
  TEST_ASSERT ( jsonf.read("modes", modes) );
  vector<double> initZ;
  TEST_ASSERT ( jsonf.read("z0", initZ) );
  assert_eq (initZ.size(), modes.size());

  Kid.resize(keyid.size());
  for (int i = 0; i < Kid.size(); ++i){
	assert_in(keyid[i],0,T-1);
	Kid[i] = keyid[i];
  }

  testModeId.resize(modes.size());
  z0.resize(testModeId.size());
  for (int i = 0; i < testModeId.size(); ++i){
    testModeId[i] = modes[i];
	z0[i] = initZ[i];
  }

  loadLambda(dataDir+"eigen_values80.b");
  const VectorXd tlam = lambda;
  lambda.resize(testModeId.size());
  for (int i = 0; i < testModeId.size(); ++i)
	lambda[i] = tlam[testModeId[i]];
  initVolObj();
}
bool MtlOptModel::loadLambda(const string fname){
  const bool succ = load(fname,lambda);
  assert(succ);
  return succ;
}
void MtlOptModel::initVolObj(){

  const string evect = dataDir+"/eigen_vectors80.b";
  MatrixXd tW;
  TEST_ASSERT( load(evect,tW) );
  W.resize(tW.rows(),testModeId.size());
  for (int i = 0; i < testModeId.size(); ++i){
	W.col(i) = tW.col(testModeId[i]);
  }	

  TEST_ASSERT( volobj.loadObjMesh(dataDir+"beam.obj"));
  TEST_ASSERT( volobj.loadTetMesh(dataDir+"beam.abq"));
  TEST_ASSERT( volobj.loadWeights(dataDir+"interp-weights.txt"));
  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(dataDir+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  rs2euler.setTetMesh(volobj.getTetMesh());
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();
  TEST_ASSERT( DefGradOperator::compute(volobj.getTetMesh(),G) );
  MapMA2RS::computeMapMatPGW(G,W,PGW);	
}
void MtlOptModel::produceSimRlst(){

  const int r = redDim();
  const VectorXd zero = VectorXd::Zero(r);
  HarmonicOscillatorSet<double> sim(lambda,alphaK,alphaM,z0,zero);
  Z = sim.generateSequence<MatrixXd>(0,h,T);
  assert_eq(Z.rows(),r);
  assert_eq(Z.cols(),T);
  
  CorrectZ.resize(r,(T-Kid.size()));
  int ci = 0;
  for (int i = 0; i < T; ++i){
	if( RedSpaceTimeEnergyAD::isKeyframe(Kid,i) < 0){
	  CorrectZ.col(ci) = Z.col(i);
	  ci ++;
	}
  }

  Kz.resize(r,Kid.size());
  for (int i = 0; i < Kid.size(); ++i){
	const int f = Kid[i];
	assert_in( f,0,Z.cols() );
	Kz.col(i) = Z.col(f);
  }
}
void MtlOptModel::extrangeKeyframes(){
  ASSERT_EQ(Kid.size(),9);
  Kid << 0,50,20,30,40,12,85,140,199;
}
void MtlOptModel::initMtlOpt(RedSpaceTimeEnergyAD &ad)const{

  ad.setT(T);
  ad.setTimestep(h);
  ad.setDamping(alphaK,alphaM,lambda);
  ad.setK(lambda);
  ad.setKeyframes(Kz,Kid);
}
void MtlOptModel::initMtlData(MtlDataModel &model){

  model.setT(T);
  model.setTimestep(h);
  model.setDamping(alphaK,alphaM,lambda);
  model.setLambda(lambda);
  model.setKeyframes(Kz,Kid);
  model.setSubZ(CorrectZ);
}
void MtlOptModel::initSolver(const SXMatrix &E, const VSX &x){

  allVars = x;
  fun = CasADi::SXFunction(x,E);
  solver = CasADi::IpoptSolver(fun);
  solver.setOption("generate_hessian",true);
  solver.init();
	
  const int r = redDim();
  const int lenZ = (T-Kid.size())*r;
  vector<double> x0(x.size(),0);
  vector<double> lowerB(x0.size(),-std::numeric_limits<double>::infinity());
  for (int i = lenZ; i < x0.size(); ++i){
	x0[i] = 0.01f;
	lowerB[i] = 0.0f;
  }
  solver.setInput(x0,CasADi::NLP_X_INIT);
  solver.setInput(lowerB,CasADi::NLP_LBX);
}
void MtlOptModel::solve(){

  solver.solve();
  double cost;
  solver.getOutput(cost,CasADi::NLP_COST);
  cout<<"input data: " << endl;
  cout<<"lambda: " << lambda.transpose() << endl << endl;
  cout<< endl << "optimal cost: " << cost << endl << endl;
  const VectorXd xx = getOutput();
  const int lenZ = redDim()*(T-Kid.size());
  cout<< "opt mtl: " << xx.tail(xx.size()-lenZ).transpose() << endl << endl;
}
VectorXd MtlOptModel::getOutput(){

  VectorXd x;
  const int r = redDim();
  vector<double> vx(allVars.size());
  solver.getOutput(vx,CasADi::NLP_X_OPT);
  x.resize(vx.size());
  for (size_t i = 0; i < vx.size(); ++i)
	x[i] = vx[i];
  ASSERT_GE(x.size(),r*(T-Kid.size()));
  return x;
}
void MtlOptModel::getZfromSolver(MatrixXd &Z){

  const VectorXd x = getOutput();
  const VectorXd kz = Map<VectorXd>(&Kz(0,0),Kz.size());
  Z = assembleFullZ(x,kz,Kid,redDim());
}
void MtlOptModel::saveRlst(){
	
  saveMesh(Z,"input");
  saveMesh(Kz,"key");
  MatrixXd newZ;
  getZfromSolver(newZ);
  saveMesh(newZ,"new");
}
void MtlOptModel::saveMesh(const MatrixXd &Z,const string fname){

  for (int i = 0; i < Z.cols(); ++i){
	const string ffname=dataDir+"tempt/meshes/"+fname+TOSTR(i)+".vtk";
	saveMeshOneZ(Z.col(i),ffname);
  }
}
void MtlOptModel::saveMeshOneZ(const VectorXd &z,const string fname){

  VectorXd y,u;
  const VectorXd p = W*z;
  RSCoordComp::constructWithWarp(G,p,y);
  TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  volobj.interpolate(u);
  TEST_ASSERT (LSW_SIM::MeshVtkIO::write(volobj.getObjMesh(),fname));
}
void MtlOptModel::computeEnergy(const VectorXd &X){
  VectorXd rlst;
  CASADI::evaluate(fun,X,rlst);
  cout<< "E(X) = " << rlst.transpose() << endl;
}
int MtlOptModel::redDim()const{
  return lambda.size();
}
MatrixXd MtlOptModel::assembleFullZ(const VectorXd&subZ,const VectorXd&keyZ,const VectorXi&Kid, const int r){

  const VectorXd &x = subZ;
  const int T = (subZ.size()+keyZ.size())/r;
  assert_eq(T*r,subZ.size()+keyZ.size());
  MatrixXd Z(r,T);
  int xpos = 0;

  for (size_t i = 0; i < T; ++i){
	const int k = RedSpaceTimeEnergyAD::isKeyframe(Kid,i);
	if(k >= 0){
	  assert_le ((k+1)*r,keyZ.size());
	  Z.col(i) = keyZ.segment(k*r,r);
	}else{
	  assert_le (xpos,x.size()-r);
	  Z.col(i) = x.segment(xpos,r);
	  xpos = xpos+r;
	}
  }

  return Z;
}

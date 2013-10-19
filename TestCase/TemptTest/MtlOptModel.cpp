#include "MtlOptModel.h"

MtlOptModel::_MtlOptModel(){

  dataDir = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  T = 200;
  h = 0.3f;
  alphaK = 0.001f;
  alphaM = 0.001f;
  Kid.resize(6);
  Kid << 0,12,50,85,140,199;

  testModeId.resize(3);
  testModeId<<0,1,2;
  z0.resize(testModeId.size());
  z0[0] = 10000;
  z0[1] = 7000;
  z0[2] = 7000;

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
  LSW_ANI_EDITOR::MASimulatorAD sim;
  sim.setTimeStep(h);
  sim.setEigenValues(lambda);
  sim.setStiffnessDamping(alphaK);
  sim.setMassDamping(alphaM);
  const VectorXd zero = VectorXd::Zero(r);
  sim.setIntialStatus(zero,z0);

  Z.resize(r,T);
  CorrectZ.resize(r,(T-Kid.size()));
  int ci = 0;
  for (int i = 0; i < T; ++i){
	Z.col(i) = sim.getEigenZ();
	if( RedSpaceTimeEnergyAD::isKeyframe(Kid,i) < 0){
	  CorrectZ.col(ci) = sim.getEigenZ();
	  ci ++;
	}
	sim.forward(zero);
  }

  Kz.resize(r,Kid.size());
  for (int i = 0; i < Kid.size(); ++i){
	const int f = Kid[i];
	assert_in( f,0,Z.cols() );
	Kz.col(i) = Z.col(f);
  }
}
void MtlOptModel::extrangeKeyframes(){
  ASSERT_EQ(Kid.size(),6);
  Kid << 0,50,12,85,140,199;
}
void MtlOptModel::initMtlOpt(RedSpaceTimeEnergyAD &ad)const{

  ad.setT(T);
  ad.setTimestep(h);
  ad.setDamping(alphaK,alphaM,lambda);
  ad.setK(lambda);
  ad.setKeyframes(Kz,Kid);
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

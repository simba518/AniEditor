#include "MtlOptModel.h"
#include <HarmonicOscillator.h>
#include <ConNodesOfFrame.h>
#include <VTKWriter.h>
#include <JsonFilePaser.h>
using namespace UTILITY;
using namespace LSW_SIM;

MtlOptModel::_MtlOptModel(const string initf){

  JsonFilePaser jsonf;
  TEST_ASSERT ( jsonf.open(initf) );
  TEST_ASSERT ( jsonf.read("T",T) );
  TEST_ASSERT ( jsonf.read("h",h) );
  TEST_ASSERT ( jsonf.read("alphaK",alphaK) );
  TEST_ASSERT ( jsonf.read("alphaM",alphaM) );

  vector<int> keyid, modes;
  TEST_ASSERT ( jsonf.read("keyId", keyid) );
  TEST_ASSERT ( jsonf.read("modes", modes) );
  if(modes.size() <= 0){
	int number_modes = 0;
	TEST_ASSERT ( jsonf.read("number_modes", number_modes) );
	assert_gt(number_modes,0);
	modes.resize(number_modes);
	for (int i = 0; i < number_modes; ++i){
	  modes[i] = i;
	}
  }

  // partial constraints
  penaltyCon = 10.0f;
  string partial_con_str;
  if (jsonf.readFilePath("partial_constraints",partial_con_str)){ 
	
	TEST_ASSERT ( jsonf.read("penaltyCon",penaltyCon) );
	ConNodesOfFrameSet AllConNodes;
	TEST_ASSERT ( AllConNodes.load(partial_con_str) );
	const set<pConNodesOfFrame> &con_groups=AllConNodes.getConNodeGroups();
	BOOST_FOREACH(pConNodesOfFrame pc, con_groups){
	  if ( pc != NULL && !pc->isEmpty() ) {
		const int frame_id = pc->getFrameId();
		const VectorXd bary_uc = pc->getBarycenterUc();
		const vector<set<int> > &groups = pc->getConNodesSet();
		vector<int> v;
		for (int i = 0; i < groups.size(); ++i)
		  v.push_back(*(groups[i].begin()));
		assert_eq (v.size()*3,bary_uc.size());

		conFrames.push_back(frame_id);
		uc.push_back(bary_uc);
		conNodes.push_back(v);
	  }
	}
  }

  vector<double> initZ;
  TEST_ASSERT ( jsonf.read("z0", initZ) );
  assert_gt(initZ.size(),0);
  if(initZ.size() != modes.size()){
	initZ.resize(modes.size());
	for (int i = 1; i< modes.size(); ++i){
	  initZ[i] = initZ[0];
	}
  }
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

  loadLambda(initf);
  const VectorXd tlam = lambda;
  lambda.resize(testModeId.size());
  for (int i = 0; i < testModeId.size(); ++i)
	lambda[i] = tlam[testModeId[i]];

  MatrixXd kz;
  TEST_ASSERT( jsonf.readMatFile("keyZ",kz) );
  Kz = kz.topLeftCorner(redDim(),Kid.size());
  CorrectZ = MatrixXd(redDim(),T);
  CorrectZ.setZero();

  initVolObj(initf);
  initWarper(jsonf);
}
bool MtlOptModel::loadLambda(const string initf){

  JsonFilePaser jsonf;
  TEST_ASSERT ( jsonf.open(initf) );
  string fname;
  TEST_ASSERT ( jsonf.readFilePath("eigen_values",fname) );
  const bool succ = load(fname,lambda);
  assert(succ);
  return succ;
}
void MtlOptModel::initVolObj(const string initf){

  JsonFilePaser jsonf;
  TEST_ASSERT ( jsonf.open(initf) );
  string eigenfname, objfname, volfname, weightfname, connodefname;
  TEST_ASSERT ( jsonf.readFilePath("eigen_vectors",eigenfname) );
  TEST_ASSERT ( jsonf.readFilePath("obj_mesh_file",objfname) );
  TEST_ASSERT ( jsonf.readFilePath("vol_filename",volfname) );
  TEST_ASSERT ( jsonf.readFilePath("vol2obj_weights_filename",weightfname) );
  TEST_ASSERT ( jsonf.readFilePath("con_nodes_file",connodefname) );

  MatrixXd tW;
  TEST_ASSERT( load(eigenfname,tW) );
  W.resize(tW.rows(),testModeId.size());
  for (int i = 0; i < testModeId.size(); ++i){
	W.col(i) = tW.col(testModeId[i]);
  }	

  TEST_ASSERT( volobj.loadObjMesh(objfname));
  TEST_ASSERT( volobj.loadTetMesh(volfname));
  TEST_ASSERT( volobj.loadWeights(weightfname));
  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(connodefname,fixed_nodes,UTILITY::TEXT) );
  rs2euler.setTetMesh(volobj.getTetMesh());
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();
  TEST_ASSERT( DefGradOperator::compute(volobj.getTetMesh(),G) );
  MapMA2RS::computeMapMatPGW(G,W,PGW);	
}
void MtlOptModel::initMtlData(MtlDataModel &model){

  model.setT(T);
  model.setTimestep(h);
  model.setDamping(alphaK,alphaM,lambda);
  model.setLambda(lambda);
  model.setKeyframes(Kz,Kid);
  model.setSubZ(CorrectZ);

  model.setPartialCon(conFrames,conNodes,uc);
  model.setPenaltyCon(penaltyCon);
  model.setWarper(warper);

}
void MtlOptModel::saveMesh(const MatrixXd &Z,const string fname){

  for (int i = 0; i < Z.cols(); ++i){
	const string ffname=fname+TOSTR(i)+".vtk";
	saveMeshOneZ(Z.col(i),ffname);
  }
}
void MtlOptModel::saveMeshOneZ(const VectorXd &z,const string fname){

  VectorXd y,u;
  const VectorXd p = W*z;
  RSCoordComp::constructWithWarp(G,p,y);
  TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  volobj.interpolate(u);
  TEST_ASSERT (volobj.getObjMesh()->writeVTK(fname));
}
void MtlOptModel::saveUc(const string fname,const VectorXd &uc,const vector<int> &nid)const{

  VVec3d v(uc.size()/3);
  const VVec3d &nodes = volobj.getTetMesh()->nodes();
  for (size_t i = 0; i < v.size(); ++i){
	v[i][0] = uc[i*3+0];
	v[i][1] = uc[i*3+1];
	v[i][2] = uc[i*3+2];
	v[i] += nodes[nid[i]];
  }
  VTKWriter writer;
  writer.addPoints(v);
  TEST_ASSERT( writer.write(fname) );
}
void MtlOptModel::saveUc(const string fname)const{
  
  for (int i = 0; i < uc.size(); ++i){
	const string ffname=fname+TOSTR(i)+".vtk";
	saveUc(ffname,uc[i],conNodes[i]);
  }
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
void MtlOptModel::initWarper(JsonFilePaser &inf){

  MatrixXd B;
  vector<int> cubP;
  vector<double> cubW;
  TEST_ASSERT(inf.readMatFile("nonlinear_basis",B));
  TEST_ASSERT(inf.readVecFile("cubaturePoints",cubP,UTILITY::TEXT));
  TEST_ASSERT(inf.readVecFile("cubatureWeights",cubW));
  warper = pRedRSWarperAD( new RedRSWarperAD(rs2euler,B,W,cubP,cubW) );
}

void MtlOptModel::print()const{
  
  cout << "initial values: " << lambda.transpose() << endl;
  cout << "number of modes: " << lambda.size() << endl;
  cout << "number of nodes: " << volobj.getTetMesh()->nodes().size() << endl;
  cout << "number of tetrahedrons: " << volobj.getTetMesh()->tets().size() << endl;
}

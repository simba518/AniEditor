#include "MtlOptModel.h"
#include <ConNodesOfFrame.h>
#include <VTKWriter.h>
using namespace UTILITY;
using namespace LSW_SIM;

MtlOptModel::MtlOptModel(const string initf){

  // load mtl
  JsonFilePaser jsonf;
  TEST_ASSERT ( jsonf.open(initf) );
  TEST_ASSERT ( jsonf.read("T",T) );
  TEST_ASSERT ( jsonf.read("h",h) );
  TEST_ASSERT ( jsonf.read("alphaK",alphaK) );
  TEST_ASSERT ( jsonf.read("alphaM",alphaM) );

  // read modes, and lambda  
  VectorXd t_lambda;
  TEST_ASSERT ( jsonf.readVecFile("eigen_values",t_lambda) );
  TEST_ASSERT ( jsonf.read("modes", testModeId) );
  if(testModeId.size() <= 0){
  	int number_modes = 0;
  	TEST_ASSERT ( jsonf.read("number_modes", number_modes) );
	if (number_modes <= 0)
	  number_modes = t_lambda.size();
  	testModeId.resize(number_modes);
  	for (int i = 0; i < number_modes; ++i)
  	  testModeId[i] = i;
  }

  lambda.resize(testModeId.size());
  ASSERT_LE(testModeId.size(),t_lambda.size());
  for (int i = 0; i < testModeId.size(); ++i)
  	lambda[i] = t_lambda[testModeId[i]];

  // keyframe constraints
  TEST_ASSERT ( jsonf.read("keyId", Kid) );
  MatrixXd kz;
  TEST_ASSERT( jsonf.readMatFile("keyZ",kz) );
  this->Kz = kz.topLeftCorner(redDim(),Kid.size());

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

  // load tet mesh
  string tetf;
  tetmesh = pTetMesh(new TetMesh);
  TEST_ASSERT( jsonf.readFilePath("vol_filename",tetf) );
  TEST_ASSERT( tetmesh->load(tetf) );

  // load or compute hatW
  if(!jsonf.readMatFile("rs_basis",hatW)){
	MatrixXd tempt_W;
	TEST_ASSERT( jsonf.readMatFile("eigen_vectors",tempt_W) );
	ASSERT_GE(tempt_W.cols(),testModeId.size());
	W.resize(tempt_W.rows(),testModeId.size());
	for (int i = 0; i < testModeId.size(); ++i)
	  W.col(i) = tempt_W.col(testModeId[i]);
	
	SparseMatrix<double> G;
	TEST_ASSERT( DefGradOperator::compute(tetmesh,G) );
	MapMA2RS::computeMapMatPGW(G,W,hatW);
  }
  ASSERT_GE(hatW.cols(),testModeId.size());
  if (hatW.cols() != testModeId.size()){
	const MatrixXd tempt_hatW = hatW;
	hatW.resize(hatW.rows(),testModeId.size());
	for (int i = 0; i < hatW.cols(); ++i)
	  hatW.col(i) = tempt_hatW.col(testModeId[i]);
  }

  // init warper
  initWarper(jsonf);
}

void MtlOptModel::initMtlData(MtlDataModel &model){

  model.setT(T);
  model.setTimestep(h);
  model.setDamping(alphaK,alphaM,lambda);
  model.setLambda(lambda);
  model.setKeyframes(Kz,Kid);
  MatrixXd z = MatrixXd(redDim(),T);
  z.setZero();
  model.setSubZ(z);

  model.setPartialCon(conFrames,conNodes,uc);
  model.setPenaltyCon(penaltyCon);
  model.setWarper(warper);
}

void MtlOptModel::saveUc(const string fname,const VectorXd&uc,const vector<int>&nid)const{

  VVec3d v(uc.size()/3);
  const VVec3d &nodes = tetmesh->nodes();
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

void MtlOptModel::initWarper(JsonFilePaser &inf){

  vector<int> fixed_nodes;
  TEST_ASSERT( inf.readVecFile("con_nodes_file",fixed_nodes,UTILITY::TEXT) );
  rs2euler.setTetMesh(tetmesh);
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();

  MatrixXd B;
  vector<int> cubP;
  vector<double> cubW;
  if (inf.readMatFile("nonlinear_basis",B) &&
	  inf.readVecFile("cubaturePoints",cubP,UTILITY::TEXT)&&
	  inf.readVecFile("cubatureWeights",cubW)
	  && W.cols() == redDim()){
	warper = pRedRSWarperAD( new RedRSWarperAD(rs2euler,B,W,cubP,cubW) );
  }
}

void MtlOptModel::print()const{
  cout << "initial values: " << lambda.transpose() << endl;
  cout << "number of modes: " << lambda.size() << endl;
  cout << "number of nodes: " << tetmesh->nodes().size() << endl;
  cout << "number of tetrahedrons: " << tetmesh->tets().size() << endl;
  cout << "time step: " << h << endl;
  cout << "alphak: " << alphaK << endl;
  cout << "alpham: " << alphaM << endl;
}

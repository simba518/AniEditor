#include "MtlOptModel.h"
#include <PartialConstraints.h>
#include <VTKWriter.h>
using namespace UTILITY;

MtlOptModel::MtlOptModel(const string initf){

  // load mtl
  JsonFilePaser jsonf;
  TEST_ASSERT ( jsonf.open(initf) );
  TEST_ASSERT ( jsonf.read("T",T) );
  TEST_ASSERT ( jsonf.read("h",h) );
  TEST_ASSERT ( jsonf.read("alpha_k",alphaK) );
  TEST_ASSERT ( jsonf.read("alpha_m",alphaM) );
  
  // read modes, and lambda  
  VectorXd t_lambda;
  TEST_ASSERT ( jsonf.readVecFile("eigenvalues",t_lambda) );
  double lambda0_scale = 1.0f;
  if ( jsonf.read("lambda0_scale", lambda0_scale) ){
	assert_gt(lambda0_scale,0);
	t_lambda *= lambda0_scale;
  }

  int number_modes = 0;
  TEST_ASSERT ( jsonf.read("rw", number_modes) );
  if (number_modes <= 0)
	number_modes = t_lambda.size();
  testModeId.resize(number_modes);
  for (int i = 0; i < number_modes; ++i)
	testModeId[i] = i;

  lambda.resize(testModeId.size());
  ASSERT_LE(testModeId.size(),t_lambda.size());
  for (int i = 0; i < testModeId.size(); ++i)
  	lambda[i] = t_lambda[testModeId[i]];

  if(!jsonf.read("rs",r_s)){
	r_s = lambda.size();
  }

  // keyframe constraints
  int fix_head = 0, fix_end = 0;
  TEST_ASSERT ( jsonf.read("fix_head", fix_head) );
  TEST_ASSERT ( jsonf.read("fix_end", fix_end) );
  Kid.clear();
  for (int i = 0; i < fix_head; ++i){
    Kid.push_back(i);
  }
  for (int i = T-fix_end; i < T; ++i){
    Kid.push_back(i);
  }
  this->Kz.resize(selectedRedDim(),Kid.size());
  this->Kz.setZero();

  // partial constraints
  partialConPenalty = 10.0f;
  string partial_con_str;
  if (jsonf.readFilePath("partial_con",partial_con_str)){
	
  	TEST_ASSERT ( jsonf.read("con_penalty",partialConPenalty) );
  	PartialConstraintsSet AllConNodes;
  	TEST_ASSERT ( AllConNodes.load(partial_con_str) );
  	const set<pPartialConstraints> &con_groups=AllConNodes.getPartialConSet();
  	BOOST_FOREACH(pPartialConstraints pc, con_groups){
  	  if ( pc != NULL && !pc->isEmpty() ) {
  		vector<int> conNod;
		VectorXd conPos;
		pc->getPartialCon(conNod,conPos);
  		conFrames.push_back(pc->getFrameId());
  		conNodes.push_back(conNod);
  		uc.push_back(conPos);
  	  }
  	}
  }

  if (!jsonf.read("s_penalty", sPenalty)){
	sPenalty = 0.0f;
  }

  // load tet mesh
  string tetf;
  tetmesh = pTetMesh(new TetMesh);
  TEST_ASSERT( jsonf.readFilePath("vol_file",tetf) );
  TEST_ASSERT( tetmesh->load(tetf) );

  // compute hatW
  MatrixXd tempt_W;
  TEST_ASSERT( jsonf.readMatFile("eigenvectors",tempt_W) );
  ASSERT_GE(tempt_W.cols(),testModeId.size());
  W.resize(tempt_W.rows(),testModeId.size());
  for (int i = 0; i < testModeId.size(); ++i)
	W.col(i) = tempt_W.col(testModeId[i]);
  for (int i = 0; i < W.cols(); ++i){
	W.col(i) *= 1.0f/sqrt(lambda[i]);
  }
	
  SparseMatrix<double> G;
  TEST_ASSERT( DefGradOperator::compute(tetmesh,G) );
  MapMA2RS::computeMapMatPGW(G,W,hatW);
  ASSERT_EQ(hatW.cols(),testModeId.size());

  if(!jsonf.readMatFile("Z_initial",Z_initial)){
	Z_initial.resize(0,0);
	Z_initial.setZero();
  }

  // init warper
  initWarper(jsonf);
}

void MtlOptModel::initMtlData(MtlDataModel &model){

  const VectorXd subLambda = lambda.head(r_s);
  model.setT(T);
  model.setTimestep(h);
  model.setDamping(alphaK,alphaM,subLambda);
  model.setLambda(subLambda);
  model.setFixframes(Kz,Kid);
  model.resetWSelectMatrix(lambda.size(),subLambda.size());
  model.setInitialLambda(lambda);

  if(Z_initial.rows() != r_s || Z_initial.size() != model.getZdim()){
	const int subT = model.getZdim()/r_s;
	Z_initial.resize(r_s,subT);
	Z_initial.setZero();
  }
  model.setSubZ(Z_initial);

  model.setPartialCon(conFrames,conNodes,uc);
  model.setConPenalty(partialConPenalty,partialConPenalty);
  model.setSPenalty(sPenalty);
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

  TEST_ASSERT( inf.readVecFile("fixed_nodes",fixed_nodes,UTILITY::TEXT) );
  rs2euler.setTetMesh(tetmesh);
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();

  if (inf.readMatFile("nonlinear_basis",B) &&
	  inf.readVecFile("cubature_points",cubP,UTILITY::TEXT)&&
	  inf.readVecFile("cubature_weights",cubW)
	  && W.cols() == redDim()){
	warper = pRedRSWarperAD( new RedRSWarperAD(rs2euler,B,W,cubP,cubW) );
  }
}

void MtlOptModel::print()const{
  
  cout << "initial values: " << lambda.head(r_s).transpose() << endl;
  cout << "number of selected modes: " << r_s << endl;
  cout << "number of modes: " << lambda.size() << endl;
  cout << "number of nodes: " << tetmesh->nodes().size() << endl;
  cout << "number of tetrahedrons: " << tetmesh->tets().size() << endl;
  cout << "time step: " << h << endl;
  cout << "alphak: " << alphaK << endl;
  cout << "alpham: " << alphaM << endl;

  cout << "keyframe id: ";
  for (int i = 0; i < Kid.size(); ++i) cout << Kid[i] << " ";
  cout<< endl << "partial penalty: " << partialConPenalty << endl;
  
}

void MtlOptModel::saveMeshOneZ(const VectorXd &z,const string fname){
  
  cout << "SAVE TO: " << fname << endl;
  VectorXd u;
  if (warper){
	INFO_LOG("use reduced rs method");
	warper->warp(z,u);
  }else{
	INFO_LOG("full rs method");
	TEST_ASSERT ( rs2euler.reconstruct(hatW*z,u) );
  }
  TEST_ASSERT ( tetmesh->writeVTK(fname,u) );
}

void MtlOptModel::updateW(const MatrixXd &W){

  this->W = W;
  warper = pRedRSWarperAD( new RedRSWarperAD(rs2euler,B,W,cubP,cubW) );
}

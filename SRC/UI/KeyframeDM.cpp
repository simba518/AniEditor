#include "KeyframeDM.h"
#include "VolApproxTri.h"
#include <JsonFilePaser.h>
#include <MatrixIO.h>
using namespace EIGEN3EXT;
using namespace UTILITY;
using namespace ANI_EDIT_UI;

KeyframeDM::KeyframeDM(pTetMeshEmbeding vo, pAniDataModel ani, pAniEditDM dm){

  vol_obj = vo;
  data_model = dm;
  animation = ani;
  applyConOnAll = true;
  begin_keyframe = 1;
  max_keyframe = 20;
  keyframe_step = 1;
  max_keynodes = 50;

}

bool KeyframeDM::initialize(const string initfile){
  
  string cub_elements_str,cub_weight_str;
  JsonFilePaser jsonf;
  if ( !jsonf.open(initfile) ){
	ERROR_LOG("failed to open "<< initfile);
	return false;
  }

  bool succ = jsonf.readFilePath("cubature_points", cub_elements_str);
  succ &= jsonf.readFilePath("cubature_weights", cub_weight_str);
  if (succ){
	succ = loadCubatures(cub_elements_str, cub_weight_str);
  }
  succ &= jsonf.read("begin_keyframe", begin_keyframe);
  succ &= jsonf.read("max_keyframe", max_keyframe);
  succ &= jsonf.read("keyframe_step", keyframe_step);
  succ &= jsonf.read("max_keynodes", max_keynodes);

  string keyframe_u;
  if (jsonf.readFilePath("keyframe_u",keyframe_u)){
	loadVolKeyframes(keyframe_u);
  }

  vector<string> keyframe_vol_files;
  jsonf.readFilePath("keyframe_vol", keyframe_vol_files,false);
  DEBUG_LOG("keyframe_vol_files.size() = " << keyframe_vol_files.size());
  vector<int> frame_id;
  if(jsonf.read("keyframe_vol_f", frame_id)){
	succ = true;
	assert_eq( keyframe_vol_files.size(), frame_id.size() );
	for (int i = 0; i < frame_id.size(); ++i){
	  DEBUG_LOG("keyframe id : " << frame_id[i]);
	  assert_ge(frame_id[i],0);
	  succ &= loadKeyVolMesh(keyframe_vol_files[i], frame_id[i]);
	}
  }

  return succ;
  
}

bool KeyframeDM::loadCubatures(const string cub_elements_str,const string cub_weight_str){
  
  vector<int> cub_elements;
  vector<double> cub_weights;
  bool succ = loadVec(cub_elements_str,cub_elements,UTILITY::TEXT);
  if (succ && cub_elements.size() > 0){
	succ = loadVec(cub_weight_str,cub_weights);
  }
  if(succ){
	assert_eq(cub_elements.size(),cub_weights.size());
	for (size_t i =0; i < cub_weights.size(); ++i)
	  cub_weights_elements.push_back(make_pair(cub_weights[i],cub_elements[i]));
  }
  sort(cub_weights_elements.begin(),cub_weights_elements.end());

  return succ;
}

bool KeyframeDM::loadObjKeyframe(const string filename){

  pObjmesh objMeshKey = pObjmesh(new Objmesh());
  if ( !objMeshKey->load(filename) ){
	return false;
  }

  VolApproxTri volapprox;
  volapprox.setEmbedingMesh(vol_obj);
  volapprox.setElasticPenalty(0.000001);
  if(!volapprox.generateKeyVolMesh(objMeshKey)){
	ERROR_LOG("failed to generate the volume keyframe mesh.");
	return false;
  }

  setKeyframe(volapprox.getVolKeyU(), currentFrameNum());
  return true;
}

void KeyframeDM::setKeyframe(const VectorXd &volKeyU,const int frame_id){
  
  vol_key_U.push_back(volKeyU);
  keyframe_id.push_back(frame_id);

  // select constraint according the partial constraints
  const int desired_points = max_keynodes;
  set<int> con_nodes_set;
  vector<int> con_nodes;
  const VVec4i &tets = vol_obj->getTetMesh()->tets();

  for (int i = 0; i < cub_weights_elements.size(); ++i){

	const int elem = cub_weights_elements[i].second;
	assert_in(elem,0,tets.size()-1);
	con_nodes_set.insert(tets[elem][0]);
	con_nodes_set.insert(tets[elem][1]);
	con_nodes_set.insert(tets[elem][2]);
	con_nodes_set.insert(tets[elem][3]);
	if (con_nodes_set.size() >= desired_points)
	  break;
  }
  
  BOOST_FOREACH(const int n, con_nodes_set)	con_nodes.push_back(n);

  // convert to partial constraints
  PartialConstraintsSet &par_con = partialCon;
  MatrixXd key_points(con_nodes.size()*3,volKeyU.cols());

  const VVec3d &nodes = vol_obj->getTetMesh()->nodes();

  Matrix<double,3,-1> pos_con(3,con_nodes.size());
  const VectorXd &uf = volKeyU;
  for (size_t i = 0; i < con_nodes.size(); ++i){

	const int n = con_nodes[i];
	assert_in(n*3,0,uf.size()-3);

	pos_con(0,i) = uf(n*3+0);
	pos_con(1,i) = uf(n*3+1);
	pos_con(2,i) = uf(n*3+2);

	key_points(i*3+0,0) = nodes[n][0]+pos_con(0,i);
	key_points(i*3+1,0) = nodes[n][1]+pos_con(1,i);
	key_points(i*3+2,0) = nodes[n][2]+pos_con(2,i);
  }

  assert_eq(con_nodes.size()*3, pos_con.size());
  par_con.addConNodes(con_nodes,frame_id);
  par_con.updatePc(pos_con,frame_id);

}

void KeyframeDM::addConNodes(const vector<int> &group){

  if ( group.size() <=0 ) return;

  if (!applyConOnAll){
	const int frame_id = currentFrameNum();
	assert(isKeyframe(frame_id) >= 0);
	partialCon.addConNodes(group,frame_id);
	resetPartialCon(frame_id);
  }else{
	const set<pPartialConstraints> par = partialCon.getPartialConSet();
	BOOST_FOREACH(pPartialConstraints ele, par){
	  const int frame_id = ele->getFrameId();
	  partialCon.addConNodes(group,frame_id);
	  resetPartialCon(frame_id);
	}
  }
}

void KeyframeDM::rmConNodes(const vector<int> &group){

  if (group.size() <=0 ) return;

  if (!applyConOnAll){
	const int frame_id = currentFrameNum();
	assert(isKeyframe(frame_id) >= 0);
	partialCon.rmConNodes(group,frame_id);
	resetPartialCon(frame_id);
  }else{
	const set<pPartialConstraints> par = partialCon.getPartialConSet();
	BOOST_FOREACH(pPartialConstraints ele, par){
	  const int frame_id = ele->getFrameId();
	  partialCon.rmConNodes(group,frame_id);
	  resetPartialCon(frame_id);
	}
  }
}

void KeyframeDM::resetPartialCon(const int frame_id){
  
  assert(isKeyframe(frame_id) >= 0);
  const VectorXd &u = getKeyU(frame_id);
  pPartialConstraints_const par = partialCon.getPartialCon(frame_id);
  Matrix<double,3,-1> pc;
  getSubUc(par->getConNodesSet(),u,pc);
  partialCon.updatePc(pc,frame_id);
}

void KeyframeDM::getSubUc(const vector<set<int> > &groups,const VectorXd &full_u,Matrix<double,3,-1> &sub_u)const{

  int nodes = 0;
  BOOST_FOREACH(const set<int>& s, groups)
	nodes += s.size();

  sub_u.resize(3,nodes);
  int index = 0;
  BOOST_FOREACH(const set<int>& s, groups){
	BOOST_FOREACH(const int i, s){
	  assert_in(i*3,0,full_u.size()-3);
	  sub_u.col(index) = full_u.segment<3>(i*3);
	  index++;
	}
  }
}

bool KeyframeDM::loadKeyVolMesh(const string filename,const int frame_id){

  TetMesh tetmesh;
  bool succ = tetmesh.load(filename);
  if (succ){
	VectorXd u, r;
	vol_obj->getTetMesh()->nodes(r);
	tetmesh.nodes(u);
	assert_eq(u.size(), r.size());
	u -= r;
	this->setKeyframe(u,frame_id);
  }
  return succ;
}

bool KeyframeDM::loadVolKeyframes(const string filename){
  
  MatrixXd U;
  const bool succ = load(filename,U);
  if (succ) setKeyframeSequence(U);
  return succ;
}

bool KeyframeDM::loadKeyVolMesh(const string filename){
  return loadKeyVolMesh(filename, currentFrameNum());;
}

void KeyframeDM::setKeyframeSequence(const MatrixXd &U){

  assert_gt(max_keyframe,0);
  assert_gt(keyframe_step,0);
  assert_ge(begin_keyframe,0);
  clear();
  for (int i = 0; i < U.cols() && i < max_keyframe; i+=keyframe_step){
    setKeyframe(U.col(i+begin_keyframe), i);
  }
}

bool KeyframeDM::saveKeyVolMeshes(const string filename)const{

  pTetMesh_const tetmesh = vol_obj->getTetMesh();
  bool succ = false;
  if (tetmesh && vol_key_U.size() > 0){
	
	MatrixXd U(vol_key_U[0].size(), vol_key_U.size());
	for (int i = 0; i < vol_key_U.size(); ++i){
	  U.col(i) = vol_key_U[i];
	}
	
	succ = tetmesh->write(filename,U);
	succ &= tetmesh->writeVTK(filename,U);
  }
  return succ;
}

bool KeyframeDM::saveKeyObjMesh(const string filename)const{

  TRACE_FUN();
  pObjmesh_const objmesh = vol_obj->getObjMesh();
  bool succ = false;
  if (objmesh && isKeyframe()>=0){
	double x,y,z;
	getSceneTranslate(x,y,z);
	VectorXd u = getKeyU();
	for (int i = 0; i < u.size()/3; ++i){
	  u[i*3+0] += x;
	  u[i*3+1] += y;
	  u[i*3+2] += z;
	}
	vol_obj->interpolate(u);
	succ = objmesh->write(filename);
	assert(succ);
  }
  return succ;
}

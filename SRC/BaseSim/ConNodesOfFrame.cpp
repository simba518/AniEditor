#include <algorithm> 
#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>
#include <assertext.h>
#include <Log.h>
#include "ConNodesOfFrame.h"
using namespace std;
using namespace LSW_SIM;

void ConNodesOfFrame::setConNodes
(const vector<set<int> > &con_nodes_set, const VectorXd &b_uc, 
 const VectorXd &rest_barycenters,const int frame_id){
  
  assert_ge (frame_id,0);
  assert_eq (b_uc.size(), (int)(con_nodes_set.size()*3));
  assert_eq (rest_barycenters.size(), b_uc.size());
  
  this->con_nodes_set.setGroup( con_nodes_set );
  this->barycenter_uc = b_uc;
  this->barycenter_rest = rest_barycenters;
  this->frame_id = frame_id;
}

VectorXd ConNodesOfFrame::getBarycenter()const{

  return barycenter_rest + barycenter_uc;
}

bool ConNodesOfFrame::equal(const ConNodesOfFrame &other,const double tol)const{

  const vector<set<int> > &other_con_nodes_set = other.getConNodesSet();
  const VectorXd &other_barycenter_uc = other.getBarycenterUc();
  const VectorXd &other_barycenter = other.getBarycenter();
  const VectorXd this_barycenter = this->getBarycenter();
  const int other_frame_id = other.getFrameId();

  const bool is_equal =  (other_frame_id == frame_id)&&
	((other_barycenter-this_barycenter).norm()<tol)&&
	((other_barycenter_uc-barycenter_uc).norm()<tol)&&
	(other_con_nodes_set == this->getConNodesSet());

  return is_equal;
}

ConNodesOfFrameSet::ConNodesOfFrameSet(const ConNodesOfFrameSet &other){
  
  copy(other);
}

ConNodesOfFrameSet&ConNodesOfFrameSet::operator=(const ConNodesOfFrameSet &o){

  return copy(o);
}

ConNodesOfFrameSet &ConNodesOfFrameSet::copy(const ConNodesOfFrameSet &other){

  this->clear();
  const set<pConNodesOfFrame> &other_node_set=other.getConNodeGroups();
  set<pConNodesOfFrame>::const_iterator it=other_node_set.begin();
  for (; it != other_node_set.end(); ++it){
	node_set.insert(pConNodesOfFrame(new ConNodesOfFrame(*(*it))));
  }
  return (*this);
}

bool ConNodesOfFrameSet::load(const string filename){

  ifstream in_file;
  in_file.open (filename.c_str());
  // failed to open file 
  if (!in_file.is_open()){
	ERROR_LOG("failed to open file: "<<filename);
	return false;
  }

  this->clear();
  int node_set_size = 0;
  int frame_id = 0;
  int barycenter_uc_len = 0;
  int con_nodes_set_size = 0;
  int con_nodes_set_i_size = 0;
  int node_id = 0;
  VectorXd barycenter_uc, barycenter_rest;
  vector<set<int> > con_nodes_set;
  set<int> con_nodes_set_i;

  // length of node_set
  in_file >> node_set_size;

  for (int f_num=0; f_num < node_set_size; ++f_num){

	con_nodes_set.clear();

	// frame_id
	in_file >> frame_id;

	// barycenter_uc,  barycenter_rest
	in_file >> barycenter_uc_len;
	barycenter_uc.resize(barycenter_uc_len);
	barycenter_rest.resize(barycenter_uc_len);
	for (int i = 0; i < barycenter_uc_len; ++i){
	  in_file >> barycenter_uc[i];
	}
	for (int i = 0; i < barycenter_uc_len; ++i){
	  in_file >> barycenter_rest[i];
	}

	// number of con_nodes_set
	in_file >> con_nodes_set_size;

	for (int k = 0; k < con_nodes_set_size; ++k){

	  // length of con_nodes_set[i]
	  in_file >> con_nodes_set_i_size;
	  con_nodes_set_i.clear();
	  // con_nodes_set[i]
	  for (int i = 0; i < con_nodes_set_i_size; ++i){
		in_file >> node_id;
		con_nodes_set_i.insert(node_id);
	  }
	  con_nodes_set.push_back(con_nodes_set_i);
	}

	// record for one frame
	pConNodesOfFrame cp = pConNodesOfFrame(new ConNodesOfFrame());
	cp->setConNodes(con_nodes_set, barycenter_uc,barycenter_rest,frame_id);
	node_set.insert(cp);
  }
  return true;
}

/**
 * format of the file is:
 * length of node_set
 * frame_id
 * length of barycenter_uc
 * barycenter_uc
 * barycenter_rest
 * number of con_nodes_set
 * length of con_nodes_set[i]
 * con_nodes_set[i]
 */
bool ConNodesOfFrameSet::save(const string filename)const{

  ofstream out_file;
  out_file.open (filename.c_str());
  // failed to open file 
  if (!out_file.is_open()){
	ERROR_LOG("failed to open file: "<<filename);
	return false;
  }

  // no constrained nodes
  if (node_set.size() <= 0){
	out_file.close();
	return true;
  }

  // length of node_set
  out_file << node_set.size() << endl;

  BOOST_FOREACH(pConNodesOfFrame groups, node_set){

	// frame_id
	out_file << groups->getFrameId() << endl;	

	// barycenter_uc,  barycenter_rest
	const VectorXd &uc = groups->getBarycenterUc();
	const VectorXd &rc = groups->getBarycenterRest();
	out_file << uc.size() << endl;
	for (int i = 0; i < uc.size(); ++i){
	  out_file<< setprecision(10) << uc[i] << " ";
	}
	out_file << endl;
	for (int i = 0; i < rc.size(); ++i){
	  out_file<< setprecision(10) << rc[i] << " ";
	}
	out_file << endl;
	
	// number of con_nodes_set
	const vector<set<int> > &nodes = groups->getConNodesSet();

	out_file<< nodes.size() << endl;

	BOOST_FOREACH(const set<int> &con_set, nodes){

	  // length of con_nodes_set[i]
	  out_file<< con_set.size() << endl;
	  // con_nodes_set[i]
	  BOOST_FOREACH(int ele, con_set){
		out_file << ele << " ";
	  }
	  out_file << endl;
	}
  }
  out_file.close();
  return true;
}

/**
 * save the target constrained positions, which is usefull to export the path
 * for control.
 */
bool ConNodesOfFrameSet::saveConPositions(const string filename)const{

  ofstream out_file;
  out_file.open (filename.c_str());
  // failed to open file 
  if (!out_file.is_open()){
	ERROR_LOG("failed to open file: "<<filename);
	return false;
  }

  // no constrained nodes
  if (node_set.size() <= 0){
	out_file.close();
	return true;
  }

  // frame_id
  out_file << "frames" << endl;
  BOOST_FOREACH(pConNodesOfFrame groups, node_set){
	out_file << groups->getFrameId() << endl;
  }
  
  // barycenter, 
  out_file << "positions" << endl;
  BOOST_FOREACH(pConNodesOfFrame groups, node_set){
	const VectorXd bc = groups->getBarycenter();
	for (int i = 0; i < bc.size(); ++i){
	  out_file << bc[i] << " ";
	}
	out_file << endl;
  }
  out_file.close();
  return true;
}

void ConNodesOfFrameSet::addConNodeGroup(pConNodesOfFrame node_group){

  if (node_group->isEmpty()){
	return ;
  }

  set<pConNodesOfFrame>::iterator it = node_set.begin();
  for (; it != node_set.end(); ++it){
    if ((*it)->getFrameId() == node_group->getFrameId()){
	  node_set.erase(it);
	  break;
	}
  }
  node_set.insert (node_group);
}

void ConNodesOfFrameSet::addConNodeGroup(const ConNodesOfFrame &node_group){
  
  pConNodesOfFrame p_node_group = pConNodesOfFrame
	(new ConNodesOfFrame(node_group.getConNodesSet(),
						 node_group.getBarycenterUc(),
						 node_group.getBarycenterRest(),
						 node_group.getFrameId()));
  addConNodeGroup(p_node_group);
}

bool ConNodesOfFrameSet::updateConNodeGroup(const VectorXd &uc,const int f){

  set<pConNodesOfFrame>::iterator it = node_set.begin();
  for (; it != node_set.end(); ++it){
    if ((*it)->getFrameId() == f){

	  pConNodesOfFrame con_nodes = (*it);
	  con_nodes->setUc (uc);
	  node_set.erase(it);
	  node_set.insert(con_nodes);
	  break;
	}
  }
  return it != node_set.end();
}

pConNodesOfFrame_const ConNodesOfFrameSet::getConNodeGroup(const int f)const{

  pConNodesOfFrame_const pc;
  set<pConNodesOfFrame>::iterator it=node_set.begin();
  for (; it != node_set.end(); ++it){
    if ( (*it)->getFrameId() == f){
  	  pc = *it;
  	  break;
  	}
  }
  return pc;
}

pConNodesOfFrame ConNodesOfFrameSet::getConNodeGroup(const int frame_id){
  
  pConNodesOfFrame pc;
  set<pConNodesOfFrame>::iterator it=node_set.begin();
  for (; it != node_set.end(); ++it){
    if ( (*it)->getFrameId() == frame_id){
  	  pc = *it;
  	  break;
  	}
  }
  return pc;
}

bool ConNodesOfFrameSet::equal(const ConNodesOfFrameSet &other,const double tol)const{

  const set<pConNodesOfFrame> &other_node_set = other.getConNodeGroups();
  if (other_node_set.size() != node_set.size()){
	return false;
  }
  
  set<pConNodesOfFrame>::const_iterator it=node_set.begin();
  set<pConNodesOfFrame>::const_iterator o_it=other_node_set.begin();
  for (; it != node_set.end(); ++it, ++o_it){

	if (!( (**it).equal( (**o_it), tol) )){
	  return false;
	}
  }
  return true;
}

void ConNodesOfFrameSet::addConNodeGroup(const vector<int> &con_nodes,const int frame_id){

  pConNodesOfFrame pc = getConNodeGroup(frame_id);
  if (pc != NULL){
	pc->addConNodeGroup(con_nodes);
  }else{
	pc = pConNodesOfFrame(new ConNodesOfFrame(frame_id) );
	pc->addConNodeGroup(con_nodes);
	this->addConNodeGroup(pc);
  }
}

void ConNodesOfFrameSet::rmConNodeGroup( const vector<int> &con_nodes,const int frame_id ){
  
  pConNodesOfFrame pc = getConNodeGroup(frame_id);
  if (pc != NULL){
	pc->rmConNodeGroup(con_nodes);
  }
}

void ConNodesOfFrameSet::setConPos(const VectorXd &uc, const VectorXd &rest,const int frame_id){
  
  pConNodesOfFrame pc = getConNodeGroup(frame_id);
  if (pc != NULL){
	pc->updateConPos(uc,rest);
  }
}

vector<ConNodesOfFrame> ConNodesOfFrameSet::getSortedConNodeGroups()const{
  
  vector<ConNodesOfFrame> conVec;
  BOOST_FOREACH(pConNodesOfFrame con, node_set){
	if (!con->isEmpty()){
	  conVec.push_back(*con);
	}
  }
  std::sort (conVec.begin(), conVec.end());
  return conVec;
}

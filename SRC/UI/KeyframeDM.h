#ifndef _KEYFRAMEDM_H_
#define _KEYFRAMEDM_H_

#include <QObject>
#include <QtGui/QMainWindow>
#include <FileDialog.h>
#include <AniDataModel.h>
#include <PartialConstraints.h>
#include <TetMeshEmbeding.h>
#include "AniEditDM.h"
using namespace QGLVEXT;
using namespace UTILITY;

namespace ANI_EDIT_UI{
  
  /**
   * @class KeyframeDM keyframe data model.
   * 
   */
  class KeyframeDM{
	
  public:
	KeyframeDM(pTetMeshEmbeding vol_obj, pAniDataModel animation, pAniEditDM dm);

	bool initialize(const string initfile);
	bool loadCubatures(const string cub_elements_str,const string cub_weights_str);
	bool loadObjKeyframe(const string filename);
	bool loadVolKeyframes(const string filename);
	bool loadKeyVolMesh(const string filename);
	bool loadKeyVolMesh(const string filename,const int frame_id);
	bool saveKeyVolMeshes(const string filename)const;
	bool saveKeyObjMesh(const string filename)const;

	void setKeyframe(const VectorXd &volKeyU,const int frame_id);
	void setKeyframeSequence(const MatrixXd &U);

	void addConNodes(const vector<int> &group);
	void rmConNodes(const vector<int> &group);
	void clear(){
	  keyframe_id.clear();
	  vol_key_U.clear();
	}

	int currentFrameNum()const{
	  if(animation)
		return animation->currentFrameNum();
	  return 0;
	}
	const VectorXd &getKeyU()const{
	  return getKeyU(currentFrameNum());
	}
	const VectorXd &getKeyU(const int fid)const{
	  const int index = isKeyframe(fid);
	  assert_in(index,0, vol_key_U.size()-1);
	  return vol_key_U[index];
	}
	const vector<VectorXd> &getAllKeyU()const{
	  return vol_key_U;
	}
	const vector<int> &getKeyframesId()const{
	  return keyframe_id;
	}
	int isKeyframe()const{
	  return isKeyframe(currentFrameNum());
	}
	int isKeyframe(const int fid)const{
	  for (int i = 0; i < keyframe_id.size(); ++i){
		if (fid == keyframe_id[i])
		  return i;
	  }
	  return -1;
	}
	pTetMesh_const getVolMesh()const{
	  return vol_obj->getTetMesh();
	}
	pObjmesh_const getObjMesh()const{
	  vol_obj->interpolate(getKeyU());
	  return vol_obj->getObjMesh();
	}
	pObjmesh_const getObjMesh(const int keyf_index)const{
	  assert_in(keyf_index,0, vol_key_U.size()-1);
	  vol_obj->interpolate(vol_key_U[keyf_index]);
	  return vol_obj->getObjMesh();
	}
	pPartialConstraints_const getPartialCon()const{
	  return partialCon.getPartialCon(currentFrameNum());
	}
	pPartialConstraints getPartialCon(){
	  return partialCon.getPartialCon(currentFrameNum());
	}
	const PartialConstraintsSet &getPartialConSet()const{
	  return partialCon;
	}
	void getSceneTranslate(double &x,double &y,double &z)const{
	  data_model->getSceneTranslate(x,y,z);
	}

	bool saveParitalCon(const string filename)const{
	  return partialCon.write(filename);
	}

  protected:
	void getSubUc(const vector<set<int> > &groups,const VectorXd &full_u,Matrix<double,3,-1> &sub_u)const;
	void resetPartialCon(const int frame_id);

  protected:
	vector<std::pair<double, int> > cub_weights_elements;
	pAniDataModel animation;
	pAniEditDM data_model;
	pTetMeshEmbeding vol_obj; 
	PartialConstraintsSet partialCon;
	vector<int> keyframe_id;
	vector<VectorXd> vol_key_U;
	bool applyConOnAll;
	int begin_keyframe;
	int max_keyframe;
	int keyframe_step;
	int max_keynodes;
  };
  typedef boost::shared_ptr<KeyframeDM> pKeyframeDM;

  class KeyframeDM_UI:public QObject,public KeyframeDM{
	
	Q_OBJECT
	
  public:
    KeyframeDM_UI(QMainWindow*win,pTetMeshEmbeding vo,pAniDataModel ani,pAniEditDM dm):
	  KeyframeDM(vo,ani,dm), data_model(dm){
	  this->file_dialog = pFileDialog(new FileDialog(win));
	}
	bool initialize(const string initfile){
	  const bool succ = KeyframeDM::initialize(initfile);
	  if (succ)	applyKeyframes();
	  return succ;
	}

  public slots:
	void loadObjKeyframe(){
	  const string fname = file_dialog->load("obj","load keyframe");
	  if(fname.size() >0) file_dialog->warning(KeyframeDM::loadObjKeyframe(fname));
	}
	void loadVolKeyframes(){
	  const string fname = file_dialog->load("b","load vol keyframe sequence");
	  if(fname.size() >0) file_dialog->warning(KeyframeDM::loadVolKeyframes(fname));
	}
	void loadKeyVolMesh(){
	  const string fname = file_dialog->load("vol","load vol keyframe mesh");
	  if(fname.size() >0) file_dialog->warning(KeyframeDM::loadKeyVolMesh(fname));
	}
	void saveKeyVolMeshes()const{
	  const string fname = file_dialog->save("abq","save vol keyframe meshes");
	  if(fname.size()>0)file_dialog->warning(KeyframeDM::saveKeyVolMeshes(fname));
	}
	void saveKeyObjMesh()const{
	  const string fname = file_dialog->save("obj","save obj keyframe mesh");
	  if(fname.size()>0)file_dialog->warning(KeyframeDM::saveKeyObjMesh(fname));
	}
	void addConNodes(const vector<int> &group){
	  KeyframeDM::addConNodes(group);
	}
	void rmConNodes(const vector<int> &group){
	  KeyframeDM::addConNodes(group);
	}
	void applyKeyframes(){
	  if(data_model)
		data_model->setKeyframes(KeyframeDM::getAllKeyU(), KeyframeDM::getKeyframesId());
	  // if(data_model)
	  // 	data_model->addPartialCon(getPartialConSet());
	}
	void toogleApplyConOnAll(){
	  applyConOnAll = applyConOnAll ? false:true;
	}
	void recordCurrentFrameAsKeyframe(){
	  assert(data_model);
	  KeyframeDM::setKeyframe(data_model->getVolFullU(),currentFrameNum());
	}
	void clear(){
	  KeyframeDM::clear();
	}
	
  private:
	pFileDialog file_dialog;
	pAniEditDM data_model;
  };
  typedef boost::shared_ptr<KeyframeDM_UI> pKeyframeDM_UI;
  
}//end of namespace

#endif /*_KEYFRAMEDM_H_*/

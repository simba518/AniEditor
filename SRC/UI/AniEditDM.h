#ifndef _ANIEDITDM_H_
#define _ANIEDITDM_H_

#include <boost/shared_ptr.hpp>
#include <AniDataModel.h>
#include <PartialConstraints.h>
#include <BaseInterpolator.h>
#include <TetMeshEmbeding.h>
#include <MatrixIO.h>
using namespace QGLVEXT;
using namespace LSW_ANI_EDITOR;
using namespace UTILITY;

namespace ANI_EDIT_UI{
  
  /**
   * @class AniEditDM animation edit data model.
   * 
   */
  class AniEditDM{
	
  public: 
	AniEditDM(pTetMeshEmbeding vol_obj, pAniDataModel animation);
	AniEditDM(pTetMeshEmbeding vol_obj,pAniDataModel animation,pBaseInterpolator interp);

	bool initialize(const string filename,const bool create_interpolator=true);

	int totalFrameNum()const;
	int reducedDim()const;
	int fullDim()const;
	int currentFrameNum()const{
	  if(_animation)
		return _animation->currentFrameNum();
	  return 0;
	}
	double getMaxRadius()const{
	  if (vol_obj){
		return vol_obj->getMaxRadius();
	  }
	  return 0;
	}

	pObjmesh_const getOutputObjMesh(const bool translate=false)const{
	  return getOutputObjMesh(currentFrameNum(), translate);
	}
	pObjmesh_const getOutputObjMesh(const int f,const bool translate=false)const;
	pObjmesh_const getInputObjMesh(const bool translate=false)const{
	  return getInputObjMesh( currentFrameNum(), translate );
	}
	pObjmesh_const getInputObjMesh(const int f,const bool translate=false)const;
	pObjmesh_const getKeyframeObjMesh()const;
	pTetMeshEmbeding_const getVolObjMesh()const{
	  return vol_obj;
	}
	pTetMesh_const getVolMesh()const{
	  return vol_obj->getTetMesh();
	}
	pTetMesh getVolMesh(){
	  return vol_obj->getTetMesh();
	}
	pObjmesh_const getScene()const{
	  return scene;
	}

	const VectorXd &getVolFullU()const{
	  return getVolFullU(currentFrameNum());
	}
	const VectorXd &getVolFullU(const int frame, const bool translate=false)const;
	const VectorXd &getUnwarpedU(const int frame, const bool translate=false)const;
	const VectorXd &getUnwarpedU()const{
	  return getUnwarpedU(currentFrameNum());
	}
	const VectorXd &getUforConstraint(const int frame)const;
	const VectorXd &getUforConstraint()const{
	  return getUforConstraint(currentFrameNum());
	}
	const VectorXd &getInputU()const{
	  return getInputU(currentFrameNum());
	}
	const VectorXd &getInputU(const int frame, const bool translate=false)const;

	pBaseInterpolator_const getInterpolator()const{
	  return interpolator;
	}

	// partial constraints
	const Matrix<double,3,-1> getXc(const int group)const;
	void updateXc(const Matrix<double,3,-1> &xc,const int group_id);
	pPartialConstraints_const getPartialCon()const{
	  return _partialCon.getPartialCon(currentFrameNum());
	}
	pPartialConstraints getPartialCon(){
	  return _partialCon.getPartialCon(currentFrameNum());
	}
	pPartialConstraints_const getPartialConX0()const{
	  return _partialCon_x0.getPartialCon(currentFrameNum());
	}
	pPartialConstraints getPartialConX0(){
	  return _partialCon_x0.getPartialCon(currentFrameNum());
	}
	void addConNodes(const vector<int> &groups);
	void rmConNodes(const vector<int> &groups);
	void updateConPos(const Matrix<double,3,-1> &uc,const int group_id);
	void updateConPos(const Matrix<double,3,-1> &uc);
	vector<set<int> > getConNodes()const;
	const Matrix<double,3,-1> getUc(const int group)const;
	const PartialConstraintsSet &getAllConNodes()const{return _partialCon;}
	void addPartialCon(const PartialConstraintsSet&parcons);
	void removeAllPosCon();
	bool loadConPath(const string filename);
	void getSceneTranslate(double &tx, double &ty, double &tz)const{
	  getSceneTranslate(tx,ty,tz,currentFrameNum());
	}
	void getSceneTranslate(double &tx, double &ty, double &tz, const int frame)const;
	
	// keyframes
	void setKeyframes(const vector<VectorXd> &U){
	  if (U.size() > 0){
		assert(interpolator);
		MatrixXd Uk(U[0].size(), U.size());
		for (int i = 0; i < U.size(); ++i)
		  Uk.col(i) = U[i];
		interpolator->setKeyframes(Uk);
	  }
	}
	void setKeyframes(const vector<VectorXd> &U, const vector<int> &kid){
	  assert_eq(U.size(), kid.size());
	  if (U.size() > 0){
		assert(interpolator);
		MatrixXd Uk(U[0].size(), U.size());
		for (int i = 0; i < U.size(); ++i)
		  Uk.col(i) = U[i];
		interpolator->setKeyframes(Uk,kid);
	  }
	}

	// warping.
	void useWarp(const bool use_warp){
	  if (interpolator)
		interpolator->useWarp(use_warp);
	}
	bool isUseWarp()const{
	  if (interpolator)
		return interpolator->isUseWarp();
	  return false;
	}

	// interpolate
	bool interpolate();
	bool toggleInterpolate(){
	  enableInterpolate = enableInterpolate ? false:true;
	  return enableInterpolate;
	}

	// print information
	void print()const;

	// file io
	bool saveParitalCon(const string filename)const;
	bool loadParitalCon(const string filename);
	void printConPositions()const;
	bool saveOutputVolMeshesVTK(const string filename)const;
	bool saveOutputObjMeshesVTK(const string filename)const;
	bool saveOutputMeshes(const string filename)const;
	bool saveInputMeshes(const string filename)const;
	bool saveCurrentOutputMesh(const string filename)const;
	bool saveCurrentInputMesh(const string filename)const;
	bool saveVolFullU(const string filename)const;
	bool saveCurrentOutputVolMesh(const string filename);
	bool loadCurrentReducedEdits(const string filename);
	bool saveCurrentReducedEdits(const string filename)const;
	bool loadAllReducedEdits(const string filename);
	bool saveAllReducedEdits(const string filename)const;
	bool saveCurrentVolFullU(const string filename)const;
	bool saveSceneSequence(const string filename, const bool VTK=false)const;
	bool saveSceneSequenceVTK(const string filename)const{
	  return saveSceneSequence(filename, true);
	}

	bool editable()const;

	// additional animation for comparison
	bool loadAnimation(const string filename){
	  return EIGEN3EXT::load(filename, additionalVolU);
	}
	bool hasAdditionalAni()const{
	  const int f = currentFrameNum();
	  return (f>=0 && f < additionalVolU.size());
	}
	pObjmesh_const getAdditionalAniObjMesh(){
	  assert(hasAdditionalAni());
	  const VectorXd &vol_u = additionalVolU[currentFrameNum()];
	  vol_obj->interpolate(vol_u);
	  return vol_obj->getObjMesh();
	}
	void clearAdditionalAni(){
	  additionalVolU.clear();
	}
	bool saveAdditionalAniObjMeshesVTK(const string filename)const;
	bool saveAdditionalAniObjMeshes(const string filename)const;
	bool savePartialConBalls(const string filename)const;

  protected:
	VectorXd barycenOfRestShape(const vector<set<int> >&g)const; 
	VectorXd barycenOfVolU(const vector<set<int> >&g)const;
	string getZeroStr(const int frame, const int T)const;
	void resetPartialCon(const int frame, const bool interp = true);
	void getSubUc(const vector<set<int> > &groups,const VectorXd &full_u,Matrix<double,3,-1> &sub_u)const;
	void translate(VectorXd &u, const int frame)const;
	
  private:
	pAniDataModel _animation;
	pBaseInterpolator interpolator;
	pTetMeshEmbeding vol_obj;
	PartialConstraintsSet _partialCon;
	PartialConstraintsSet _partialCon_x0;
	bool enableInterpolate;
	vector<VectorXd> additionalVolU; // animation for comparison.
	VectorXd rest_shape;

	// scene and translation
	pObjmesh scene;
	vector<double> scene_translation; // len%3 = 0
  };
  
  typedef boost::shared_ptr<AniEditDM> pAniEditDM;
  typedef boost::shared_ptr<const AniEditDM> pAniEditDM_const;
  
}//end of namespace

#endif /*_ANIEDITDM_H_*/

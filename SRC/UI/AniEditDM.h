#ifndef _ANIEDITDM_H_
#define _ANIEDITDM_H_

#include <boost/shared_ptr.hpp>
#include <AniDataModel.h>
#include <ConNodesOfFrame.h>
#include <BaseInterpolator.h>
#include <TetMeshEmbeding.h>
#include "DragTrajectoryRecord.h"
using namespace QGLVEXT;
using namespace LSW_ANI_EDITOR;
using namespace LSW_SIM;

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

	pObjmesh_const getOutputObjMesh()const{
	  return getOutputObjMesh(currentFrameNum());
	}
	pObjmesh_const getOutputObjMesh(const int f)const;
	pObjmesh_const getInputObjMesh()const{
	  return getInputObjMesh( currentFrameNum() );
	}
	pObjmesh_const getInputObjMesh(const int f)const;
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

	const VectorXd &getVolFullU()const{
	  return getVolFullU(currentFrameNum());
	}
	const VectorXd &getVolFullU(const int frame)const;
	const VectorXd &getUnwarpedU()const;
	const VectorXd &getUforConstraint()const;
	const VectorXd &getInputU()const{
	  return getInputU(currentFrameNum());
	}
	const VectorXd &getInputU(const int frame)const;

	pBaseInterpolator_const getInterpolator()const{
	  return interpolator;
	}

	// constraints
	void addConNodes(const vector<int> &groups);
	void rmConNodes(const vector<int> &groups);
	void updateConPos(const VectorXd &uc);
	vector<set<int> > getConNodes()const;
	const ConNodesOfFrameSet &getAllConNodes()const{return con_nodes_for_edit;}
	void removeAllPosCon();

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

	// print information
	void print()const;

	// file io
	bool saveParitalCon(const string filename)const;
	bool loadParitalCon(const string filename);
	bool loadConPath(const string filename);
	bool saveOutputVolMeshesVTK(const string filename)const;
	bool saveOutputMeshes(const string filename)const;
	bool saveInputMeshes(const string filename)const;
	bool saveCurrentOutputMesh(const string filename)const;
	bool saveCurrentInputMesh(const string filename)const;
	bool saveVolFullU(const string filename)const;
	bool saveCurrentOutputVolMesh(const string filename);
	bool loadCurrentReducedEdits(const string filename);
	bool saveCurrentReducedEdits(const string filename)const;
	bool loadAllReducedEdits(const string filename)const;
	bool saveAllReducedEdits(const string filename)const;
	bool saveCurrentVolFullU(const string filename)const;

	// record drag operations
	void setRecord(const bool recrod){
	  this->record_drag = recrod;
	}
	bool isRecordDrag()const{
	  return record_drag;
	}
	bool saveDragRecord(const string filename)const{
	  return drag_record.save(filename);
	}
	bool loadDragRecord(const string filename){
	  return drag_record.load(filename);
	}
	bool playRecodDrag();

	bool editable()const;
	void computeConNodeTrajectory();
	const vector<VectorXd> &getConNodeTraj()const{
	  return con_node_traj;
	}

  protected:
	VectorXd barycenOfRestShape(const vector<set<int> >&g)const; 
	VectorXd barycenOfVolU(const vector<set<int> >&g)const;
	string getZeroStr(const int frame, const int T)const;
	void recordDragOp(const int frame_num, const VectorXd &uc){
	  if (isRecordDrag()){
		drag_record.record(frame_num, uc);
	  }
	}
	
  private:
	pAniDataModel _animation;
	pBaseInterpolator interpolator; // perform the edit operation.
	ConNodesOfFrameSet con_nodes_for_warping; // con nodes for warping.
	ConNodesOfFrameSet con_nodes_for_edit; // con nodes under edit.
	pTetMeshEmbeding vol_obj; // store both vol and obj mesh.
	DragTrajectoryRecord drag_record;
	bool record_drag; // if true then record the drag operation.
	vector<VectorXd> con_node_traj; // trajectory of contraint nodes.
  };
  
  typedef boost::shared_ptr<AniEditDM> pAniEditDM;
  typedef boost::shared_ptr<const AniEditDM> pAniEditDM_const;
  
}//end of namespace

#endif /*_ANIEDITDM_H_*/

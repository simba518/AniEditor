#ifndef _VOLOBJMESHCTRL_H_
#define _VOLOBJMESHCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <QtGui/QMainWindow>
#include <VolObjMesh.h>
#include <FileDialog.h>
using namespace QGLVEXT;
using namespace LSW_SIM;

namespace LSW_BASE_UI{
  
  /**
   * @class VolObjMeshCtrl manager the VolObjMesh object: load obj mesh, volume
   * mesh, interpolate weights, and precompute interpolate weights.
   * 
   */
  class VolObjMeshCtrl: public QObject{

	Q_OBJECT

  public: 
	VolObjMeshCtrl(QMainWindow *main_win, pVolObjMesh volobjmesh);
	bool initialize(const string filename);
	pVolObjMesh getVolObjMesh(){
	  return volobjmesh;
	}
	pVolObjMesh_const getVolObjMesh()const{
	  return volobjmesh;
	} 
	pVolumetricMesh_const getVolMesh()const{
	  assert (volobjmesh != NULL);
	  return volobjmesh->getVolMesh();
	}
	pObjRenderMesh_const getVolSurface()const{
	  assert (volobjmesh != NULL);
	  return volobjmesh->getVolSurface();
	}
	pObjRenderMesh_const getObjMesh()const{
	  assert (volobjmesh != NULL);
	  return volobjmesh->getObjMesh();
	}
	
  public slots:
	void loadObjMesh();
	void loadVolMesh();
	void loadVolSurface();
	void loadInterpWeights();
	void saveInterpWeights()const;
	void togglePhoneShading();

	void precomputeInterpWeights();
	void sendMsgResetScene()const;
 
  signals:
	// (x0,y0,z0) is the corner with the smallest coordinates of the bounding box.
	// (x1,y1,z1) is the corner with the largest coordinates of the bounding box.
	void resetSceneMsg(double x0,double y0,double z0,double x1,double y1,double z1)const;
	
	void update();
	
  private:
	pFileDialog file_dialog;
	pVolObjMesh volobjmesh;
  };
  
  typedef boost::shared_ptr<VolObjMeshCtrl> pVolObjMeshCtrl;
  
}//end of namespace

#endif /*_VOLOBJMESHCTRL_H_*/

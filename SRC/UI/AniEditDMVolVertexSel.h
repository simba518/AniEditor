#ifndef _ANIEDITDMVOLVERTEXSEL_H_
#define _ANIEDITDMVOLVERTEXSEL_H_

#include <boost/shared_ptr.hpp>
#include <Selectable.h>
#include <SelectCtrl.h>
#include "AniEditDM.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class AniEditDMVolVertexSel implement the drawWithNames() function to
   * select the vertex of the volume mesh of current frame of AniEditDM data
   * model.
   */
  class AniEditDMVolVertexSel:public Selectable, public SelectObserver{
	
  public: 
	AniEditDMVolVertexSel(pAniEditDM data_model);

	int totalEleNum ()const;
	void drawWithNames ()const;

	void addSelection(const vector<int> &sel_ids){

	  if (data_model != NULL)
		data_model->addConNodes(sel_ids);
	}
	void removeSelection(const vector<int> &sel_ids){

	  if (data_model != NULL)
		data_model->rmConNodes(sel_ids);
	}

  protected:
	bool hasVolMesh()const{

	  bool has = false;
	  if(data_model != NULL && data_model->getVolObjMesh() != NULL)
		has = (data_model->getVolObjMesh()->getVolMesh() != NULL);
	  return has;
	}
	pVolumetricMesh_const getVolMesh()const{
	  return data_model->getVolObjMesh()->getVolMesh();
	}

	void drawVertice()const;
	
  private:
	pAniEditDM data_model;
	
  };
  
  typedef boost::shared_ptr<AniEditDMVolVertexSel> pAniEditDMVolVertexSel;
  
}//end of namespace

#endif /*_ANIEDITDMVOLVERTEXSEL_H_*/

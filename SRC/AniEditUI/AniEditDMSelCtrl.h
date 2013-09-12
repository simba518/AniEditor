#ifndef _ANIEDITDMSELCTRL_H_
#define _ANIEDITDMSELCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <SelectCtrl.h>
#include "AniEditDMVolVertexSel.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class AniEditDMSelCtrl control the select operation of AniEditDM object.
   * 
   */
  class AniEditDMSelCtrl: public QObject{

	Q_OBJECT
	
  public:
	AniEditDMSelCtrl(pQGLViewerExt viewer, pAniEditDM data_model){

	  vol_vert_sel=pAniEditDMVolVertexSel(new AniEditDMVolVertexSel(data_model));
	  sel_ctrl = pSelectCtrl(new SelectCtrl(viewer,vol_vert_sel));
	  sel_ctrl->setObserver(vol_vert_sel);
	}
	
  public slots:
	void togglePrintSelEle(){ sel_ctrl->togglePrintSelEle(); }
	void enableSelOp(){	sel_ctrl->enableSelOp(); }
	void disableSelOp(){ sel_ctrl->disableSelOp(); }
	bool toggleSelOp(){	return sel_ctrl->toggleSelOp();	}

  private:
	pSelectCtrl sel_ctrl;
	pAniEditDMVolVertexSel vol_vert_sel;
  };
  
  typedef boost::shared_ptr<AniEditDMSelCtrl> pAniEditDMSelCtrl;
  
}//end of namespace

#endif /*_ANIEDITDMSELCTRL_H_*/

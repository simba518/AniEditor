#include <MatrixIO.h>
#include <MatrixTools.h>
#include <ConMatrixTools.h>
#include <JsonFilePaser.h>
#include <AuxTools.h>
#include <DrawCurves.h>
#include <Log.h>
#include "MtlOptInterpolator.h"
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

MtlOptInterpolator::MtlOptInterpolator():WarpInterpolator(){

  _ctrlF = pCtrlForceEnergy(new CtrlForceEnergy());
  _mtlOpt = pMtlOptEnergy(new MtlOptEnergy());
  _ctrlFSolver = pNoConIpoptSolver(new NoConIpoptSolver(_ctrlF));
  _mtlOptSolver = pNoConIpoptSolver(new NoConIpoptSolver(_mtlOpt));
  _optMtl = false;
}

bool MtlOptInterpolator::init(const string init_filename){

  TRACE_FUN();
  bool succ = WarpInterpolator::init(init_filename);
  JsonFilePaser json_f;
  succ &= json_f.open(init_filename);

  if (succ){
	double h, alpha_k, alpha_m, penaltyCon,penaltyKey;
	succ &= json_f.read("h",h);
	succ &= json_f.read("alpha_k",alpha_k);
	succ &= json_f.read("alpha_m",alpha_m);
	succ &= json_f.read("penaltyCon",penaltyCon);
	succ &= json_f.read("penaltyKey",penaltyKey);
	if (succ){

	  _ctrlF->setTimestep(h);
	  _ctrlF->setTotalFrames(getT());
	  _ctrlF->setPenaltyCon(penaltyCon,penaltyKey);
	  const VectorXd D = alpha_m*VectorXd::Ones(Lambda.size())+Lambda*alpha_k;
	  _ctrlF->setMtl(Lambda,D);

	  _mtlOpt->setTimestep(h);
	  _mtlOpt->setTotalFrames(getT());
	  _mtlOpt->setMtl(Lambda,alpha_k,alpha_m);
	}
  }
  if (succ){
	_ctrlF->setRedWarper(nodeWarper);	
	if(!json_f.read("optimize_mtl",_optMtl)) _optMtl = false;
  }
  succ &= _ctrlFSolver->initialize();
  succ &= _mtlOptSolver->initialize();
  if (!succ) ERROR_LOG("failed to initialize MtlOptInterpolator.");
  return succ;
}

bool MtlOptInterpolator::interpolate(){

  if (modalDisplayer.isShowModalModes()){
	return modalDisplayer.showModalModes(_ctrlF->getLambda(),delta_z);
  }

  _ctrlF->setPartialCon(con_frame_id,con_nodes,uc);

  const int MAX_IT = 100;
  const double TOL = 1e-3;
  _ctrlFSolver->setPrintLevel(5);
  _mtlOptSolver->setPrintLevel(5);
  _ctrlFSolver->setTol(0.01);
  _mtlOptSolver->setTol(0.01);
  _ctrlFSolver->setMaxIt(60000);
  _mtlOptSolver->setMaxIt(60000);

  bool succ = _ctrlFSolver->solve();
  double objValue = _ctrlF->getObjValue();

  for (int it = 0; (it < MAX_IT) && _optMtl; ++it){

	_mtlOpt->setZ( _ctrlF->getZ() );
	succ = _mtlOptSolver->solve();
	_ctrlF->setKD(_mtlOpt->getK(),_mtlOpt->getD());
	succ = _ctrlFSolver->solve();
	if( fabs(_ctrlF->getObjValue() - objValue) <= TOL )
	  break;
	objValue = _ctrlF->getObjValue();
	INFO_LOG("----------------OUTTER ITERATION: "<< it);
  }

  EIGEN3EXT::convert(_ctrlF->getU()*_ctrlF->getZ(),delta_z);

  return succ;
}

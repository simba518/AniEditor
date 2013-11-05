#include <MatrixIO.h>
#include <MatrixTools.h>
#include <ConMatrixTools.h>
#include <JsonFilePaser.h>
#include <AuxTools.h>
#include <DrawCurves.h>
#include <Log.h>
#include "MtlOptInterpolatorAdj.h"
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

MtlOptInterpolatorAdj::MtlOptInterpolatorAdj():WarpInterpolator(){

  ctrlF = pCtrlForceEnergy(new CtrlForceEnergy());
  mtlOpt = pMtlOptEnergy(new MtlOptEnergy());
  ctrlFSolver = pNoConIpoptSolver(new NoConIpoptSolver(ctrlF));
  mtlOptSolver = pNoConIpoptSolver(new NoConIpoptSolver(mtlOpt));
  _optMtl = false;
}

bool MtlOptInterpolatorAdj::init(const string init_filename){

  TRACE_FUN();
  bool succ = WarpInterpolator::init(init_filename);
  JsonFilePaser json_f;
  if (succ){
	succ = json_f.open(init_filename);
  }
  if (succ){
	double h, alpha_k, alpha_m, penaltyCon,penaltyKey;
	succ &= json_f.read("h",h);
	succ &= json_f.read("alpha_k",alpha_k);
	succ &= json_f.read("alpha_m",alpha_m);
	succ &= json_f.read("penaltyCon",penaltyCon);
	succ &= json_f.read("penaltyKey",penaltyKey);
	if (succ){

	  ctrlF->setTimestep(h);
	  ctrlF->setTotalFrames(getT());
	  ctrlF->setPenaltyCon(penaltyCon,penaltyKey);
	  const VectorXd D = alpha_m*VectorXd::Ones(Lambda.size())+Lambda*alpha_k;
	  ctrlF->setMtl(Lambda,D);

	  mtlOpt->setTimestep(h);
	  mtlOpt->setTotalFrames(getT());
	  mtlOpt->setMtl(Lambda,alpha_k,alpha_m);
	}
  }
  if (succ){
	succ = initWarper(json_f);
	if (succ)  ctrlF->setRedWarper(nodeWarper);
  }
  if (succ){
	if(!json_f.read("optimize_mtl",_optMtl)) _optMtl = false;
  }
  succ &= ctrlFSolver->initialize();
  succ &= mtlOptSolver->initialize();
  if (!succ) ERROR_LOG("failed to initialize MtlOptInterpolatorAdj.");
  return succ;
}

bool MtlOptInterpolatorAdj::interpolate(){

  if (modalDisplayer.isShowModalModes()){
	return modalDisplayer.showModalModes(ctrlF->getLambda(),delta_z);
  }

  ctrlF->setPartialCon(con_frame_id,con_nodes,uc);
  ctrlF->setZ0(delta_z[0]);
  ctrlF->setV0(VectorXd::Ones(reducedDim())*0.0f);
  const int MAX_IT = 200;
  const double TOL = 1e-3;
  ctrlFSolver->setPrintLevel(5);
  mtlOptSolver->setPrintLevel(5);
  ctrlFSolver->setTol(0.01);
  mtlOptSolver->setTol(0.01);
  ctrlFSolver->setMaxIt(500);
  mtlOptSolver->setMaxIt(500);
  bool succ = ctrlFSolver->solve();
  double objValue = ctrlF->getObjValue();

  static MatrixXd V,Z;
  for (int it = 0; (it < MAX_IT) && _optMtl; ++it){
	ctrlF->forward(V,Z);
	mtlOpt->setVZ(V,Z);
	succ = mtlOptSolver->solve();

	ctrlF->setKD(mtlOpt->getK(),mtlOpt->getD());
	succ = ctrlFSolver->solve();
	if( fabs(ctrlF->getObjValue() - objValue) <= TOL )
	  break;
	objValue = ctrlF->getObjValue();
	INFO_LOG("----------------OUTTER ITERATION: "<< it);
  }
  ctrlF->forward(V,Z);
  Z = ctrlF->getU()*Z;
  assert_eq(Z.rows(),reducedDim());
  assert_eq(Z.cols(),getT());
  EIGEN3EXT::convert(Z,delta_z);

  return succ;
}

#include <MatrixTools.h>
#include <ConMatrixTools.h>
#include <Log.h>
#include "RedRSInterpolator.h"
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

bool RedRSInterpolator::init (const string init_filename){
  
  TRACE_FUN();
  bool succ = WarpInterpolator::init(init_filename);
  succ &= anieditor.initialize(init_filename);
  anieditor.setTotalFrame(u_ref.size());
  if (succ){
	anieditor.precompute();
  }
  ERROR_LOG_COND("failed to initialize RedRSInterpolator.",succ);
  return succ;
}

bool RedRSInterpolator::interpolate (){

  if (modalDisplayer.isShowModalModes()){
	return modalDisplayer.showModalModes(anieditor.getEigenvalues(),delta_z);
  }

  vector<MatrixXd > Cs(con_frame_id.size());
  const int r = reducedDim();
  const int iterNum = 6; ///@todo
  if (con_frame_id.size() > 0){
	for (int it = 0; it < iterNum; ++it){

	  // get partial constraints
	  VectorXd new_uc;
	  for (size_t f = 0; f < con_frame_id.size(); ++f){

		Cs[f].resize(con_nodes[f].size()*3,r);
		const VectorXd &z = delta_z[ con_frame_id[f] ];
		for (size_t i = 0; i < con_nodes[f].size(); ++i){

		  const int con_node = con_nodes[f][i];
		  const Vector3d u3 = nodeWarper->warp(z,con_frame_id[f],con_node);
		  Eigen::Matrix<double,3,-1> C;
		  nodeWarper->jacobian(z,f,con_node,C);
		  const Vector3d c3 = uc[f].segment(i*3,3)+C*z-u3;
		  new_uc.conservativeResize(new_uc.size()+3);
		  new_uc.segment(new_uc.size()-3,3) = c3;
		  Cs[f].block(3*i,0,3,r) = C;
		}
	  }

	  // solve
	  VectorXd zc;
	  anieditor.computeConstrainedFrames(Cs,new_uc,zc);

	  assert_eq(zc.size(), con_frame_id.size()*r);
	  for (size_t i = 0; i < con_frame_id.size(); ++i){
		delta_z[con_frame_id[i]] = zc.segment(r*i,r);
	  }

	  /// @todo test
	  static VectorXd lastZ0;
	  assert_gt (con_frame_id.size(),0);
	  if ( lastZ0.size() == delta_z[con_frame_id[0]].size() ){
		cout<< "z diff: " << (lastZ0 - delta_z[con_frame_id[0]]).norm() << endl;
	  }
	  lastZ0 = delta_z[con_frame_id[0]];
	}
  }

  // get the results
  anieditor.updateZ();
  const int T = getT();
  const VectorXd &z = anieditor.getEditResult();
  assert_eq(z.size(),(T-anieditor.numBoundaryFrames())*r);
  delta_z.resize(T);
  int count = 0;
  for (int i = 0; i < T; ++i){
  	if ( !anieditor.isBoundaryFrame(i) ){
  	  delta_z[i] = z.segment(count*r,r);
  	  count ++;
  	}else{
  	  delta_z[i].setZero();
  	}
  }
  return true;
}

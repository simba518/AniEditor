#ifndef _MODALWARPINGEXT_H_
#define _MODALWARPINGEXT_H_

#include <interfaces/ipopt/ipopt_solver.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <Timer.h>
#include <AniSeqWarper.h>
#include <NodeRotVec.h>
#include <NodeWarper.h>
#include <NodeUnwarper.h>
#include <UnwarperAD.h>
using UTILITY::Timer;
using CasADi::IpoptSolver;

namespace LSW_WARPING{
  
  /**
   * @class ModalWarpingExt
   * 
   */
  class ModalWarpingExt:public AniSeqWarper{
	
  public:
	void setRefU(const vector<VectorXd> &ref_u){
	  this->ref_p = ref_u;
	  ref_w.resize(ref_u.size());
	}
	bool precompute(){

	  // compute G
	  bool succ = true;
	  assert (tetmesh != NULL);
	  assert_gt (tetmesh->nodes().size(),3);
	  DefGradOperator::compute(tetmesh,G);
	  assert_eq(G.rows(),tetmesh->tets().size()*9);
	  assert_eq(G.cols(),tetmesh->nodes().size()*3);

	  // initialize UnwarperAD
	  // Timer timer;
	  // timer.start();
	  // energy = pUnwarpEnergyAD( new UnwarpEnergyAD(tetmesh,0.5f,0.3f,0.2f) );
	  // timer.stop("**************init unwarp energy: ");
	  
	  // compute ref_w and ref_p
	  assert_eq(ref_p.size(), ref_w.size());
	  for (size_t i = 0; i < ref_w.size(); ++i){
		// computeWP(i);
		computeWPWithoutSolve(i);
	  }
	  return succ;
	}
	void warp(const VectorXd &p,int f_id,VectorXd &u){

	  assert_in(f_id,0,ref_p.size()-1);
	  VectorV3 new_w;
	  NodeRotVec::compute(tetmesh,G,p,new_w);
	  for (size_t i = 0; i < new_w.size(); ++i){
		new_w[i] += ref_w[f_id][i];
	  }
	  const VectorXd new_p = ref_p[f_id] + p;
	  NodeWarper::warp(new_w, new_p, u);
	}

	const VectorV3 &getInputRotVec(int f_id)const{
	  assert_in(f_id,0,ref_w.size()-1);
	  return ref_w[f_id];
	}
	const VectorXd &getInputLocalDisp(int f_id)const{
	  assert_in(f_id,0,ref_p.size()-1);
	  return ref_p[f_id];
	}
	void getWarpedR(const VectorXd &p, int f_id, SparseMatrix<double> &warpedR)const{

	  assert_in(f_id,0,ref_p.size()-1);
	  assert_in(f_id,0,ref_w.size()-1);
	  VectorV3 w;
	  NodeRotVec::compute(tetmesh,G,p,w);
	  assert_eq(w.size(),ref_w[f_id].size());
	  for (size_t i = 0; i < w.size(); ++i){
		w[i] += ref_w[f_id][i];
	  }
	  NodeWarper::warpMat(w,warpedR);
	}
	void getR(const VectorXd &p, int f_id, SparseMatrix<double> &R){

	  assert_in(f_id,0,ref_p.size()-1);
	  assert_in(f_id,0,ref_w.size()-1);
	  VectorV3 w;
	  NodeRotVec::compute(tetmesh,G,p,w);
	  assert_eq(w.size(),ref_w[f_id].size());
	  for (size_t i = 0; i < w.size(); ++i){
		w[i] += ref_w[f_id][i];
	  }
	  NodeWarper::rotMat(w,R);
	}

  protected:
	void getInitialValue(int f_id,vector<double>&wp)const{

	  VectorV3 w0;
	  VectorXd p0;
	  NodeUnwarper::compute(tetmesh,G,ref_p[f_id],w0,p0);
	  assert_eq(w0.size()*3,p0.size());

	  wp.resize(p0.size()*2);
	  for(size_t i=0; i < w0.size(); ++i){

		wp[i*3+0] = w0[i][0];
		wp[i*3+1] = w0[i][1];
		wp[i*3+2] = w0[i][2];
		if (w0[i].norm() < 1e-8){
		  wp[i*3+0] = 1e-8;
		}
	  }
	  for(int i = 0; i < p0.size(); ++i){
		wp[i+w0.size()*3] = p0[i];
	  }
	}
	void computeWP(int f_id){

	  // init energy function
	  assert_eq(ref_p[f_id].size(), tetmesh->nodes().size()*3);
	  energy->genEnergyFun(ref_p[f_id]);
	  
	  // init solver
	  IpoptSolver solver(energy->getEnergyFun());
	  solver.setOption("generate_hessian",false);
	  solver.setOption("tol",1.0f);
	  solver.setOption("max_iter",100);
	  solver.init();

	  // set intial value
	  vector<double> wp;
	  getInitialValue(f_id,wp);
	  solver.setInput(wp,CasADi::NLP_X_INIT);

	  // solve and get result
	  solver.solve();
	  vector<double> uopt(energy->getDim());
	  solver.getOutput(uopt,CasADi::NLP_X_OPT);
	  ref_w[f_id].resize(ref_p[f_id].size()/3);
	  for (size_t i = 0; i < ref_w[f_id].size(); ++i){
		ref_w[f_id][i][0] = uopt[i*3+0];
		ref_w[f_id][i][1] = uopt[i*3+1];
		ref_w[f_id][i][2] = uopt[i*3+2];
	  }
	  for (int i = 0; i < ref_p[f_id].size(); ++i){
		ref_p[f_id][i] = uopt[ref_p[f_id].size()+i];
	  }
	}
	void computeWPWithoutSolve(int f_id){
	  const VectorXd u = ref_p[f_id];
	  NodeUnwarper::compute(tetmesh,G,u,ref_w[f_id],ref_p[f_id]);
	}

  private:
	SparseMatrix<double> G;
	pUnwarpEnergyAD energy;
	vector<VectorV3> ref_w; // input rotation vectors.
	vector<VectorXd> ref_p; // input local displacements.
  };
  
  typedef boost::shared_ptr<ModalWarpingExt> pModalWarpingExt;
  
}//end of namespace

#endif /*_MODALWARPINGEXT_H_*/

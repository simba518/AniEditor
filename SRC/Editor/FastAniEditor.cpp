#include <algorithm>
#include <stdexcept>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/UmfPackSupport>
#include <JsonFilePaser.h>
#include "SpaceTimeHessianAutoDiff.h"
#include "FastAniEditor.h"
using namespace std;
using namespace IEDS;

FastAniEditor::FastAniEditor(){

  hessian_comp = pBaseSpaceTimeHessian(new SpaceTimeHessianAutoDiff());
  boundary_frames.push_back(0);
  boundary_frames.push_back(1);
}

bool FastAniEditor::initialize(const string &ini_file){

  UTILITY::JsonFilePaser json_f;
  bool succ = json_f.open(ini_file);
  if (succ){
	double h,alpha_k,alpha_m;
	VectorXd eigen_values;
	vector<int> b_frames;

	succ &= json_f.read("h",h);
	succ &= json_f.read("alpha_k",alpha_k);
	succ &= json_f.read("alpha_m",alpha_m);
	succ &= json_f.readVecFile("eigen_values",eigen_values);
	json_f.read("boundary_frames",b_frames);
	
	if(succ){
	  setTimeStep(h);
	  setStiffnessDamping(alpha_k);
	  setMassDamping(alpha_m);
	  setEigenValues(eigen_values);
	  if (b_frames.size() > 0){
		setBoundayFrames(b_frames);
	  }
	}
  }
  return succ;
}

// calculate the Hessian Matrix, and remove the rows and columns in
// boundary_frames, which are used as boudary conditions
bool FastAniEditor::precompute(){

  bool succ = false;
  VecT temp_H_triplet;
  if(hessian_comp != NULL){
	succ = hessian_comp->compute_HTriplet(temp_H_triplet);
  }
  if(succ){
	
	H_triplet.reserve(temp_H_triplet.size());
	for (int i = 0; i < (int)temp_H_triplet.size(); ++i){
	
	  const Eigen::Triplet<double> &Tri = temp_H_triplet[i];
	  const int row = (int)Tri.row();
	  const int col = (int)Tri.col();
	  const int frame_id0 = row/reducedDim();
	  const int frame_id1 = col/reducedDim();
	  if (!isBoundaryFrame(frame_id0) && !isBoundaryFrame(frame_id1)){
		const int fc0 = removedFrameNumBefore(frame_id0);
		const int fc1 = removedFrameNumBefore(frame_id1);
		assert_ge(fc0,0);
		assert_ge(fc1,0);
		const int rm0 = fc0*reducedDim();
		const int rm1 = fc1*reducedDim();
		assert_in(row-rm0,0,totalDim());
		assert_in(col-rm1,0,totalDim());
		H_triplet.push_back(Eigen::Triplet<double>(row-rm0,col-rm1,Tri.value()));
	  }
	}
  }else{
	ERROR_LOG("failed to compute the Hessian Matrix H");
  }
  return succ;
}

/**
 * Assign C to A.bolck(r0,c0, C.rows(),C.cols()).
 */
void blockAssign(VecT &A_triplet,const MatrixXd &C,const int r0,const int c0){

  assert_ge(r0,0);
  assert_ge(c0,0);
  for (int i = 0; i < C.rows(); ++i){
    for (int j = 0; j < C.cols(); ++j){
	  A_triplet.push_back( Eigen::Triplet<double>(r0+i,c0+j,C(i,j)) );
	}
  }
}

void FastAniEditor::setPosCon(const vector<int> &con_frame_id,
							  const vector<MatrixXd> &C,
							  const vector<VectorXd> &delta_uc){

  // get dimensions
  const int r = reducedDim();
  int n = totalDim();
  int nz = H_triplet.size();
  for (int i = 0; i < (int)C.size(); ++i){
    n += C[i].rows();
  	nz += C.size();
  }

  // set triplets of A
  VecT A_triplet;
  A_triplet.reserve(nz);

  // assign H_triplet to A_triplet
  for (int i = 0; i < (int)H_triplet.size(); ++i){
	
	const Eigen::Triplet<double> &Tri = H_triplet[i];
	const int row = (int)Tri.row();
	const int col = (int)Tri.col();
	A_triplet.push_back( Eigen::Triplet<double>(row,col,Tri.value() ));
  }

  // assign C and C^t to A_triplet
  int r0 = totalDim();
  assert_ge(r0,0);
  for (int i = 0; i < (int)C.size(); ++i){

  	assert_in (con_frame_id[i],0,totalFrames()-1);
	assert ( !isRemoved(con_frame_id[i]) );
	const int fc = removedFrameNumBefore(con_frame_id[i]);
	if ( fc >= 0 ){
	  const int c0 = (con_frame_id[i]-fc)*r;
	  blockAssign (A_triplet, C[i],r0,c0 );
	  blockAssign (A_triplet, C[i].transpose(),c0,r0 );
	  r0 += C[i].rows();
	}
  }

  // initialize A from A_triplet
  A.resize (n,n);
  A.reserve (nz);
  A.setFromTriplets( A_triplet.begin(), A_triplet.end() );

  // resize and set P
  P.resize(n);
  P.setZero();
  int p_c0 = totalDim();
  for (int i = 0; i < (int)delta_uc.size(); ++i){

    P.segment(p_c0,delta_uc[i].size()) = delta_uc[i];
  	p_c0 += delta_uc[i].size();
  }
}

void FastAniEditor::setTotalFrame(const int T){

  assert ( hessian_comp != NULL);
  assert_ge ( T,5 );
  hessian_comp->setTotalFrame(T);
  correctBoundaryFrames();
}

void FastAniEditor::setTimeStep(const double h){
  
  assert ( hessian_comp != NULL);
  hessian_comp->setTimeStep(h);
}

void FastAniEditor::setStiffnessDamping(const double s){

  assert ( hessian_comp != NULL);
  hessian_comp->setStiffnessDamping(s);
}

void FastAniEditor::setMassDamping(const double s){

  assert ( hessian_comp != NULL);
  hessian_comp->setMassDamping(s);
}

void FastAniEditor::setEigenValues(const VectorXd &eigen_values){
  
  assert ( hessian_comp != NULL);
  hessian_comp->setEigenValues(eigen_values);
}

int FastAniEditor::reducedDim()const{
  
  int r = 0;
  if(hessian_comp != NULL){
	r = hessian_comp->reducedDim();
  }
  return r;
}

int FastAniEditor::totalFrames()const{
  
  int T = 0;
  if(hessian_comp!=NULL){
	T = hessian_comp->totalFrames();
  }
  return T;
}

int FastAniEditor::totalDim()const{
  
  const int r = reducedDim();
  const int T = totalFrames();
  assert_ge(r,0);
  assert_ge(T-numBoundaryFrames(),0);
  return (T-numBoundaryFrames())*r;
}

/// Solve A X = P for X, where X = [Z , lambda]^t, and P = [O , P^m]^t.
bool FastAniEditor::applyEdits(){

  assert_eq(A.rows(), P.size());
  assert_eq(A.rows(), A.cols());
  bool succ = true;

  UmfPackLU<SparseMatrix<double> > solver;
  solver.compute(A);
  if(solver.info()!=Eigen::Success) {
  	ERROR_LOG("failed to decompse A for solving A X = P");
	succ = false;
  }

  VectorXd X(A.rows());
  X.setZero();
  if (succ){

	X = solver.solve(P);
	if(solver.info()!=Eigen::Success) {
	  succ = false;
	  ERROR_LOG("failed to solve for A X = P.");
	}
  }

  const int r = reducedDim();
  const int T = totalFrames();
  const int z_len = r*T;

  Z.resize(z_len);
  Z.setZero();
  for (int frame_id = 0; frame_id < totalFrames(); ++frame_id){
	
	const int fc = removedFrameNumBefore(frame_id);
	if (fc >= 0){
	  Z.segment( frame_id*r,r ) = X.segment( (frame_id-fc)*r,r );
	}
  }
  return succ;
}

bool FastAniEditor::isRemoved(const int frame_id)const{
  
  for (int i = 0; i < (int)boundary_frames.size(); ++i){
    if (boundary_frames[i] == frame_id){
	  return true;
	}
  }
  return false;
}

int FastAniEditor::removedFrameNumBefore(const int frame_id)const{

  int i = 0;
  for (; i < (int)boundary_frames.size(); ++i){
    if (boundary_frames[i] == frame_id){
	  i = -1;
	  break;
	}else if (boundary_frames[i] > frame_id){
	  break;
	}
  }
  return i;
}

void FastAniEditor::setBoundayFrames(const vector<int> &b_frames){

  boundary_frames.clear();
  for (size_t i =0; i < b_frames.size(); ++i){
	assert_ge(b_frames[i],0);
	boundary_frames.push_back(b_frames[i]);
  }
  sort( boundary_frames.begin(), boundary_frames.end() );
  if (boundary_frames.size() <= 1){

	string msg = "boundary frames should be >= 2,";
	msg += "and we will use frame 0 and 1 as boundary frames instead";
	WARN_LOG(msg);

	boundary_frames.clear();
	boundary_frames.push_back(0);
	boundary_frames.push_back(1);
  }
}

int FastAniEditor::numBoundaryFrames()const{
  
  return (int)boundary_frames.size();
}

bool FastAniEditor::editable(const int frame_id)const{
  
  if (frame_id >=0 && frame_id < totalFrames()){
	return !( this->isBoundaryFrame(frame_id) );
  }
  return false;
}

bool FastAniEditor::isBoundaryFrame(const int frame_id)const{
  
  for (size_t i = 0; i < boundary_frames.size(); ++i){
    if (frame_id == boundary_frames[i]){
	  return true;
	}
  }
  return false;
}

void FastAniEditor::correctBoundaryFrames(){

  vector<int> b_frames = boundary_frames;
  boundary_frames.clear();
  for (size_t i = 0; i < b_frames.size(); ++i){
	if (b_frames[i] >= 0 && b_frames[i] < totalFrames()){
	  boundary_frames.push_back(b_frames[i]);
	}
  }
  assert_ge(boundary_frames.size(),2);
}

#include <JsonFilePaser.h>
#include <SpaceTimeHessianAutoDiff.h>
#include <FastIntuitiveAniEditor.h>
#include <SpaceTimeHessianInverse.h>
using namespace IEDS;
using namespace UTILITY;

// read: h, alpha_k, alpha_m, r,T,eigen_values, boundary_frames.
bool FastIntuitiveAniEditor::initialize(const string &ini_file){

  UTILITY::JsonFilePaser json_f;
  bool succ = json_f.open(ini_file);
  if (succ){

	succ &= json_f.read("h",h);
	succ &= json_f.read("alpha_k",alpha_k);
	succ &= json_f.read("alpha_m",alpha_m);
	succ &= json_f.readVecFile("eigen_values",eigen_values);
	vector<int> b_frames;
	json_f.read("boundary_frames",b_frames);
	if (b_frames.size() >= 2){
	  setBoundayFrames(b_frames);
	}
  }
  return succ;
}
	
// assemble H and calculate inv(H)
void FastIntuitiveAniEditor::precompute(){

  // compute H
  SpaceTimeHessianAutoDiff H_comp;
  H_comp.setTotalFrame(T);
  H_comp.setTimeStep(h);
  H_comp.setStiffnessDamping(alpha_k);
  H_comp.setMassDamping(alpha_m);
  H_comp.setEigenValues(eigen_values);
  SparseMatrix<double> H;
  Timer timer;
  timer.start();
  H_comp.compute_H(H);
  timer.stop("compute H: ");

  // remove boundary frames from H
  std::set<int> remove_rows_set;
  for (size_t i = 0; i < boundary_frames.size(); ++i){
    remove_rows_set.insert(boundary_frames[i]);
  }
  SparseMatrix<double> P;
  timer.start();
  genReshapeMatrix(H.rows(),reducedDim(),remove_rows_set,P);
  H = P*H*P.transpose();
  timer.stop("remove boundaries of H: ");
  
  // compute inv(H)
  SpaceTimeHessianInverse::inverse(H,reducedDim(),invH);

  // cout << "(H*invH).norm = " << (H*invH).norm() << endl << endl;
  cout << "condition(H): "<< H.norm()*invH.norm() << endl;
  const SparseMatrix<double> I = eye(H.rows(),1.0f);
  cout << "(H*invH-I).norm(): " << (H*invH-I).norm() << endl;
}

// manipulate
void FastIntuitiveAniEditor::setConstrainedFrames(const vector<int> &con_frame_id){

  const int r = reducedDim();
  assert_gt(r,0);
  if(con_frame_id.size() > 0){

	std::set<int> reserve_frames_set;
	for (size_t i = 0; i < con_frame_id.size(); ++i){
	  if (i < con_frame_id.size()-1){
		// check: con_frame_id should be sorted.
		assert_gt(con_frame_id[i+1],con_frame_id[i]);
	  }
	  const int frame = con_frame_id[i];
	  assert_in(frame,0,totalFrame()-1);
	  const int rm = removedFrameNumBefore(frame);
	  assert_in(rm,0,frame);
	  const int rm_frame = frame-rm;
	  reserve_frames_set.insert(rm_frame);
	}
	SparseMatrix<double> P;
	genReshapeMatrix(invH.rows(),r,reserve_frames_set,P,false);
	invHcZ = invH*P.transpose();
	invHc = P*invHcZ;
	const int n = con_frame_id.size();
	invHcSet.resize(n*(1+n)/2);
	int count = 0;
	for (size_t i = 0; i < con_frame_id.size(); ++i){
	  for (size_t j = i; j < con_frame_id.size(); ++j){
		invHcSet[count] = block(invHc,i*r,j*r,r,r);
		count ++;
	  }
	}
  }else{
	clearConstraints();
  }
}

void FastIntuitiveAniEditor::computeConstrainedFrames(const vector<MatrixXd> &Cs,
													  const VectorXd &uc,
													  VectorXd &zc){
  
  assert(hasConstraints());
  MatrixXd A(uc.size(),uc.size());
  int r0 = 0;
  int c0 = 0;
  int count = 0;
  for (size_t i = 0; i < Cs.size(); ++i){
  	const int r = Cs[i].rows();
  	const int tc = c0;
  	for (size_t j = i; j < Cs.size(); ++j){
  	  const int c = Cs[j].rows();
  	  A.block(r0,c0,r,c) = Cs[i]*invHcSet[count]*(Cs[j].transpose());
  	  count ++;
  	  c0 += c;
  	}
  	c0 = tc + Cs[i].rows();
  	r0 += r;
  }
  for (int i = 0; i < A.rows(); ++i){
    for (int j = 0; j < i; ++j){
  	  A(i,j) = A(j,i);
  	}
  }
  const VectorXd _lambda = A.llt().solve(uc);
  const int r = reducedDim();
  C_lambda.resize(r*Cs.size());
  count = 0;
  for (size_t i = 0; i < Cs.size(); ++i){
  	C_lambda.segment(i*r,r) = Cs[i].transpose()*_lambda.segment(count,Cs[i].rows());
  	count += Cs[i].rows();
  }
  zc = invHc*C_lambda;
}

void FastIntuitiveAniEditor::updateZ(){
  if(hasConstraints()){
	assert_eq(invHcZ.cols(),C_lambda.size());
	z = invHcZ*C_lambda;
  }else{
	z.resize((totalFrame()-numBoundaryFrames())*reducedDim());
	z.setZero();
  }
}

int FastIntuitiveAniEditor::removedFrameNumBefore(const int frame_id)const{

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

void FastIntuitiveAniEditor::correctBoundaryFrames(){

  vector<int> b_frames = boundary_frames;
  boundary_frames.clear();
  for (size_t i = 0; i < b_frames.size(); ++i){
	if (b_frames[i] >= 0 && b_frames[i] < totalFrame()){
	  boundary_frames.push_back(b_frames[i]);
	}
  }
  assert_ge(boundary_frames.size(),2);
}

void FastIntuitiveAniEditor::clearConstraints(){
  
  invHcZ.resize(0,0);
  invHc.resize(0,0);
  invHcSet.resize(0);
}

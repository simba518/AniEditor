#include <MatrixIO.h>
#include <MatrixTools.h>
#include <Log.h>
#include <RSWarperExt.h>
#include <RedRSWarperExt.h>
#include <Timer.h>
#include <PythonFFT.h>
#include "MtlOptInterpolator.h"
#include "MtlOptSolverSQP.h"
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

MtlOptInterpolator::MtlOptInterpolator(){

  dataModel = pMtlOptDM(new MtlOptDM());
  reducedRS = pReducedRS(new ReducedRS());
  mtlopt = pMtlOptSolver(new MtlOptSolverSQP(dataModel, reducedRS));
  mtlopt->setMaxFrequency(1000.0f);
  simulator = pMtlOptSimulator(new MtlOptSimulatorAnalytic(dataModel));
}

bool MtlOptInterpolator::init(const string init_filename){

  TRACE_FUN();
  JsonFilePaser json_f;
  bool succ = json_f.open(init_filename);

  if (succ){

	json_f.readFilePath("save_S_to", save_S_to, false);

	// init warper
	succ &= loadUref(json_f,u_ref);
	succ &= initWarper(json_f);

	// init data model
	int T = 0, simulate_len = 0;
	succ &= json_f.read("T",T);
	assert_gt(T,0);

	int fixHead = 0, fixEnd = 0;
	double h, ak, am;
	succ &= json_f.read("h",h);
	succ &= json_f.read("alpha_k",ak);
	succ &= json_f.read("alpha_m",am);
	succ &= json_f.read("fix_head", fixHead);
	succ &= json_f.read("fix_end", fixEnd);	

	dataModel->setTotalFrames(T);
	dataModel->setDamping(ak,am);
	dataModel->fixHeadFrames(fixHead);
	dataModel->fixTailFrames(fixEnd);

	VectorXd lambda;
	succ &= json_f.readVecFile("eigenvalues",lambda);
	double lambda0_scale = 1.0f;
	if (json_f.read("lambda0_scale", lambda0_scale)){
	  lambda *= lambda0_scale;
	}
	INFO_LOG("the scaled initial lambda is: "<<lambda.transpose());
	
	int S_rows = 0;
	if( json_f.read("rw", S_rows) && S_rows > 0){
	  assert_in(S_rows,1,lambda.size());
	  const VectorXd t = lambda;
	  lambda = t.head(S_rows);
	}

	int S_columns_min = 0;
	int S_columns_max = 0;
	succ &= json_f.read("rs_min", S_columns_min);
	succ &= json_f.read("rs_max", S_columns_max);

	assert_in(S_columns_max,0, lambda.size());
	assert_in(S_columns_min,0, lambda.size());
	dataModel->initVariables(lambda, S_columns_min);

	// init energies
	double penaltyCon = 1.0f, penaltyS = 1.0f;
	succ &= json_f.read("con_penalty",penaltyCon);
	succ &= json_f.read("s_penalty",penaltyS);

	mtlopt->getEcS()->setPenalty(penaltyCon);
	mtlopt->getEcZ()->setPenalty(penaltyCon);

	// @todo use penalty for keyframe respectively.
	double penaltyKey = 1.0f;
	if(!json_f.read("keyf_penalty",penaltyKey)){
	  penaltyKey = 1.0f;
	}
	mtlopt->getEkS()->setPenalty(penaltyKey);
	mtlopt->getEkZ()->setPenalty(penaltyKey);

	mtlopt->getEsS()->setPenalty(penaltyS);
	if (mtlopt->getEmZ()){
	  double penalty_v = 0.0, penalty_z = 0.0f;
	  if (!json_f.read("v_penalty",penalty_v))
		penalty_v = 0.0f;
	  if (!json_f.read("z_penalty",penalty_z))
		penalty_z = 0.0f;
	  mtlopt->getEmZ()->setPenalty(penalty_v, penalty_z);
	  int mechanic_con_frames = 2;
	  if (!json_f.read("mechanic_con_frames",mechanic_con_frames))
		mechanic_con_frames = 2;
	  dataModel->setMechanicConFrames(mechanic_con_frames);
	}
	
	// init solver
	mtlopt->setR(S_columns_min, S_columns_max);
	int maxIt, maxIt_s, maxIt_z;
	double tol, tol_s_f, tol_s_g, tol_z_f, tol_z_g;
	succ &= json_f.read("opt_mtl_maxIt",maxIt);
	succ &= json_f.read("opt_mtl_tol",tol);
	succ &= json_f.read("opt_z_tol_g",tol_z_g);
	succ &= json_f.read("opt_z_tol_f",tol_z_f);
	succ &= json_f.read("opt_z_maxit",maxIt_z);
	succ &= json_f.read("opt_s_tol_g",tol_s_g);
	succ &= json_f.read("opt_s_tol_f",tol_s_f);
	succ &= json_f.read("opt_s_maxit",maxIt_s);
	mtlopt->setMaxIt(maxIt, maxIt_z, maxIt_s);
	mtlopt->setTol(tol, tol_z_f, tol_z_g, tol_s_f, tol_s_g);

	bool useMtlopt = true;
	succ &= json_f.read("opt_mtl",useMtlopt);
	mtlopt->setMtlOpt(useMtlopt);

	// init simulator
	succ &= json_f.read("simulate_len",simulate_len);
	assert_ge(simulate_len,0);
	simulator->setTotalFrames(simulate_len);

  }

  ERROR_LOG_COND("failed to initialize WarpInterpolator.", succ);

  return succ;
}

bool MtlOptInterpolator::interpolate(){
  
  Timer timer;
  timer.start();
  bool succ = mtlopt->solve();
  simulator->simulate();
  const double t = timer.stop("",false);
  cout << "interpolate time: " << t << endl;

  // if (save_S_to.size() > 0){
  // 	const bool succ = EIGEN3EXT::write(save_S_to, dataModel->S);
  // 	ERROR_LOG_COND("failed to save S to" << save_S_to, succ);
  // }

  { //save some results: S, lambda
	// bool succ = EIGEN3EXT::write("./tempt/opt_S.b", dataModel->S); assert(succ);
	// succ = EIGEN3EXT::write("./tempt/opt_lambda.b", dataModel->Lambda); assert(succ);
	// succ = EIGEN3EXT::write("./tempt/opt_scaled_W.b", scaled_W); assert(succ);
	
	// MatrixXd Z, Ef;
	// getAllRedDisp(Z);
	// getAllRedCtrlF(Ef);
	// succ = EIGEN3EXT::write("./tempt/opt_Z.b", Z); assert(succ);
	// succ = EIGEN3EXT::write("./tempt/opt_F.b", Ef); assert(succ)
  }

  return succ;
}

void MtlOptInterpolator::setKeyframes(const MatrixXd &Uk){
  
  MatrixXd zk;
  unwarp.convert(Uk,zk);
  vector<int> keyf;
  for (int i = 0; i < zk.cols(); ++i){
	keyf.push_back(i);
  }
  dataModel->setKeyfames(keyf, zk);
  mtlopt->getEkZ()->updateConstraints();
  
  { /// @todo tempt code
	PythonFFT pfft;
	for (int i = 0; i < zk.rows(); ++i){
	  const VectorXd zi_curve = zk.row(i).transpose();
	  pfft.add(string("z"+TOSTR(i)),zi_curve,1.0f);
	  pfft.write("fft(zk)",string("./tempt/figures/fft"));
	}
  }

}

void MtlOptInterpolator::setKeyframes(const MatrixXd &Uk,const vector<int>& keyf){

  assert_eq(Uk.cols(), keyf.size());
  MatrixXd zk;
  unwarp.convert(Uk,zk);
  dataModel->setKeyfames(keyf, zk);
  mtlopt->getEkZ()->updateConstraints();
}

const VectorXd& MtlOptInterpolator::getInterpU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  assert(reducedRS);
  assert(dataModel);

  const MatrixXd &S = dataModel->S;
  const MatrixXd &Z = dataModel->Z;
  assert_eq(Z.cols()+simulator->totalFrames(), getT());
  assert_eq(Z.rows(), S.cols());
  assert_eq(S.rows(), (reducedRS->getHatW().cols()));

  VectorXd zi;
  if (frame_id < Z.cols()){
  	zi = Z.col(frame_id);
  }else{
  	zi = simulator->getZ(frame_id-Z.cols());
  }

  { // display mode
	// const VectorXd tmp = zi;
  	// zi.setZero();
	// const int m = 2;
	// zi[m] = tmp[m];
  }

  zi = S*zi;
  reducedRS->warp( zi, frame_id, full_u_i );
  return full_u_i;
}

void MtlOptInterpolator::setConGroups(const int frame_id,const vector<set<int> >&group, 
									  const Eigen::Matrix<double,3,-1> &pc){

  if(!editable(frame_id)){
	ERROR_LOG("frame "<<frame_id<<" can not constrained.");
	return ;
  }

  assert_ge(pc.size(),3);
  const VectorXd uc = Map<VectorXd>(const_cast<double*>(&pc(0,0)),pc.size());

  vector<int> g;
  BOOST_FOREACH(const set<int> & s, group){
	BOOST_FOREACH(const int nodeId, s){
	  g.push_back(nodeId);
	}
  }
  assert_eq(g.size()*3, uc.size());
  dataModel->setConGroups(frame_id, g, uc);
  mtlopt->updateConstraints();

}

void MtlOptInterpolator::setUc(const int frame_id, const Eigen::Matrix<double,3,-1> &pc){
 
  if(!editable(frame_id)){
  	ERROR_LOG("frame "<<frame_id<<" can not constrained.");
  	return ;
  }
  assert_ge(pc.size(),3);
  const VectorXd &uc = Map<VectorXd>(const_cast<double*>(&pc(0,0)),pc.size());

  dataModel->setUc(frame_id,uc);
  mtlopt->updateConstraints();
}

void MtlOptInterpolator::sortCubPoints(vector<int> &cubPoints,vector<double>&cubWeights){

  assert_eq(cubPoints.size(), cubWeights.size());
  set<pair<int,double> > s;
  for (int i = 0; i < cubPoints.size(); ++i){
    s.insert(pair<int,double>(cubPoints[i], cubWeights[i]));
  }

  cubPoints.clear();
  cubWeights.clear();
  set<pair<int,double> >::iterator it = s.begin();
  for (; it != s.end(); ++it){
	DEBUG_LOG(it->first <<", " << it->second);
    cubPoints.push_back(it->first);
	cubWeights.push_back(it->second);
  }
}

bool MtlOptInterpolator::initWarper(JsonFilePaser &inf){

  pRedRSWarperExt nodeWarper;
  bool succ = true;

  vector<int> cubP;
  vector<double> cubW;
  pRSWarperExt warper = pRSWarperExt(new RSWarperExt());

  { // init nodeWarper
	string vol_filename;
	bool succ = inf.readFilePath("vol_file",vol_filename);
	pTetMesh tet_mesh = pTetMesh(new TetMesh());
	if (succ) succ = tet_mesh->load(vol_filename);
	if (succ){ warper->setTetMesh(tet_mesh); }

	vector<int> fixed_nodes;
	if(inf.readVecFile("fixed_nodes",fixed_nodes,UTILITY::TEXT))
	  warper->setFixedNodes(fixed_nodes);

	warper->setRefU(u_ref);
	if(succ) succ = warper->precompute();

	MatrixXd B,W;
	succ &= inf.readMatFile("eigenvectors",W);

	VectorXd lambda;
	succ &= inf.readVecFile("eigenvalues",lambda);
	double lambda0_scale = 1.0f;
	if (inf.read("lambda0_scale", lambda0_scale)){
	  lambda *= lambda0_scale;
	}
	for (int i = 0; i < W.cols(); ++i){
	  W.col(i) *= 1.0f/sqrt(lambda[i]);
	}
	succ &= inf.readMatFile("nonlinear_basis",B);

	int S_rows = 0;
	if( inf.read("rw", S_rows) && S_rows > 0){
	  assert_in(S_rows,1,W.cols());
	  const MatrixXd t = W;
	  W = t.leftCols(S_rows);
	}
	DEBUG_LOG("W.rows(): " << W.rows() << endl);
	DEBUG_LOG("W.cols(): " << W.cols() << endl);
	DEBUG_LOG("W.norm(): " << W.norm() << endl);

	scaled_W = W;
	if (succ && inf.readVecFile("cubature_points",cubP,UTILITY::TEXT)
		&& inf.readVecFile("cubature_weights",cubW)){
	  sortCubPoints(cubP, cubW);
	  nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W,cubP,cubW));
	}else{
	  nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W));
	}
  }

  { // init reducedRS using nodeWarper
	if (succ){

	  reducedRS->setBP(nodeWarper->getB(), nodeWarper->getP());
	  reducedRS->setHatW(nodeWarper->getHatW());
	  reducedRS->setSqrtV(nodeWarper->getSqrtV(),false);
	  reducedRS->setRefY(nodeWarper->getRefY());
	  reducedRS->checkDimensions();
	}
  }

  { // init ReducedRSUnwarp using node nodeWarper
	if (succ){
	  unwarp.setHatW(nodeWarper->getHatW());
	  const SparseMatrix<double> &G = warper->get_G();
	  set<int> cubSet;
	  for (int i = 0; i < cubP.size(); ++i){
		if (i > 0){
		  assert_gt(cubP[i],cubP[i-1]);
		}
		cubSet.insert(cubP[i]);
	  }
     const SparseMatrix<double>S=EIGEN3EXT::genReshapeMatrix<double>(G.rows(),9,cubSet,false);
	 const SparseMatrix<double> subG = S*G;
	 unwarp.setSubG(subG);
	}
  }
  
  return succ;
}

bool MtlOptInterpolator::loadUref(JsonFilePaser &json_f,vector<VectorXd> &u_ref)const{

  int T = 0;
  if( !json_f.read("T",T) )
	return false;
  assert_ge(T,3);
  int simulate_len = 0;
  if ( json_f.read("simulate_len",simulate_len) ){
	assert_ge(simulate_len,0);
	T += simulate_len;
  }

  MatrixXd W;
  if ( !json_f.readMatFile("eigenvectors",W) )
	return false;

  const int nx3 = W.rows();
  assert_gt(nx3,0);
  assert_eq(nx3%3,0);

  MatrixXd Uin;
  if ( json_f.readMatFile("input_animation",Uin) ){
	assert_eq(Uin.rows(), nx3);
	assert_ge(Uin.cols(), T);
	u_ref.resize(T);
	for (int i = 0; i < T; ++i) 
	  u_ref[i] = Uin.col(i);
  }else{
	u_ref.resize(T);
	for (int i = 0; i < T; ++i) {
	  u_ref[i].resize(nx3);
	  u_ref[i].setZero();
	}
  }

  return true;
}

void MtlOptInterpolator::printSolveInf()const{

  INFO_LOG("lambda0: "<<dataModel->Lambda0.transpose());
  dataModel->print();

  MatrixXd Z;
  getAllRedDisp(Z);

  const MatrixXd newZ = Z.transpose();
  cout << "curve_new_z: "<<newZ.rows()<<" "<<newZ.cols()<<" ";
  cout << Map<VectorXd >(const_cast<double*>(&newZ(0,0)),newZ.size()).transpose() << endl;

  cout << "norms of zi: ";
  for (int i = 0; i < dataModel->Z.rows(); ++i){
    cout << " " << dataModel->Z.row(i).norm();
  }
  cout << "\n\nS :" << dataModel->S << endl;
  cout<< endl << "final optimal value: " << mtlopt->funValue() << endl;

  cout << "\n\n\nSTxS: " << dataModel->S.transpose()*dataModel->S << endl;
  
  const MatrixXd StS = dataModel->S.transpose()*dataModel->S;
  Eigen::SelfAdjointEigenSolver<MatrixXd> eigenK(StS);
  cout << "\n\n\n eign(S)" << eigenK.eigenvalues() << endl;
  
  double minS = dataModel->S(0,0);
  for (int i = 0; i < dataModel->S.rows(); i++){
    for (int j = 0; j < dataModel->S.cols(); ++j){
	  if (dataModel->S(i,j) < minS)
		minS = dataModel->S(i,j);
	}
  }
  cout << "mins = " << minS << endl;

  MatrixXd Ef(reducedDim(), getT()-2);
  { // print Ef
	cout << "Ef = : ";
	for (int i = 1; i <= getT()-2; ++i){
	  const VectorXd &z0 = Z.col(i-1);
	  const VectorXd &z1 = Z.col(i);
	  const VectorXd &z2 = Z.col(i+1);
	  const VectorXd &d = dataModel->Damping;
	  const VectorXd &la = dataModel->Lambda;
	  const VectorXd df = (z2-2*z1+z0)+d.asDiagonal()*(z2-z1)+la.asDiagonal()*z1;
	  Ef.col(i-1) = df;
	  const double dd = df.norm();
	  cout << dd*dd << " ";
	}
  }

  const MatrixXd new_Ef = Ef.transpose();
  cout << "\n\ncurve_new_ef: "<<new_Ef.rows()<<" "<<new_Ef.cols()<<" ";
  cout << Map<VectorXd >(const_cast<double*>(&new_Ef(0,0)),new_Ef.size()).transpose() << endl;

  { // print Ek
  	cout << "\n\nEk = : ";
	if (dataModel->keyframes.size() == dataModel->Z.cols()){
	  const vector<int> &keyframes = dataModel->keyframes;
	  const MatrixXd &Zk = dataModel->Zk;
	  assert_eq(dataModel->S.rows(), Zk.rows());
	  assert_eq(dataModel->S.cols(), dataModel->Z.rows());
	  const MatrixXd Z = dataModel->S*dataModel->Z;
	  assert_eq(Zk.cols(), keyframes.size());
	  for (int i = 0; i < keyframes.size(); ++i){
		const int f = keyframes[i];
		if (i > 0){
		  assert_gt(keyframes[i], keyframes[i-1]);
		}
		assert_in(f,0,Z.cols()-1);
		assert_eq(Zk.rows(), Z.rows());
		const double dd = (Zk.col(i)-Z.col(f)).norm();
		cout << dd*dd << " ";
	  }
	}
  }

  if (mtlopt->optimizeMtl()){
	cout<< endl << "optimize mtl: " << 1 << endl;
  }
}

void MtlOptInterpolator::getAllRedDisp(MatrixXd &Z)const{
  
  Z.resize(reducedDim(), getT());
  Z.leftCols(dataModel->Z.cols()) = dataModel->Z;
  for (int i = 0; i < simulator->totalFrames(); ++i){
    Z.col(i+dataModel->Z.cols()) = simulator->getZ(i);
  }
}

void MtlOptInterpolator::getAllRedCtrlF(MatrixXd &Ef)const{

  MatrixXd Z;
  getAllRedDisp(Z);

  Ef.resize(reducedDim(), getT()-2);
  for (int i = 1; i <= getT()-2; ++i){
	const VectorXd &z0 = Z.col(i-1);
	const VectorXd &z1 = Z.col(i);
	const VectorXd &z2 = Z.col(i+1);
	const VectorXd &d = dataModel->Damping;
	const VectorXd &la = dataModel->Lambda;
	Ef.col(i-1) = (z2-2*z1+z0)+d.asDiagonal()*(z2-z1)+la.asDiagonal()*z1;
  }
}

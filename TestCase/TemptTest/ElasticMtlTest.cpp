#include <string>
#include <eigen3/Eigen/Dense>
#include <Eigen/SVD>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <CASADITools.h>
#include <MatrixIO.h>
#include <MatrixTools.h>
#include <DefGradOperator.h>
#include <RSCoordComp.h>
#include <MapMA2RS.h>
#include <volumetricMesh.h>
#include <volumetricMeshLoader.h>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;
using namespace CASADI;  

class ElasticMtlOpt{
  
public:
  ElasticMtlOpt(double h=0.1f,int PCADim=80, double ak=0.0f,double am=0.0f):
	_h(h),_alphaK(ak),_alphaM(am),_PCADim(PCADim){}
  bool loadVolMesh(const string filename){
	_volMesh = pVolumetricMesh(VolumetricMeshLoader::load(filename.c_str()));
	return (_volMesh != NULL);
  }
  bool loadAniSeq(const string filename){
	return EIGEN3EXT::load(filename,_U);
  }
  void compute(){
	computeRS();
	PCAonRS();
	computeK();
	decomposeK();
  }
  
  void computeRS(){

	// given _volMesh, _U, compute _Urs;
	TRACE_FUN();
	SparseMatrix<double> G;
	const bool succ= LSW_WARPING::DefGradOperator::compute(_volMesh,G);
	assert(succ);
	assert_ge(_U.cols(),3);
	_Urs.resize(G.rows(),_U.cols());
	VectorXd y;
	for (int i = 0; i < _U.cols(); ++i){
	  if(_U.col(i).norm () > 1e-8){
		LSW_WARPING::RSCoordComp::constructWithoutWarp(G,_U.col(i),y);
		assert_eq(y.size(),_Urs.rows());
		if(y.norm() != y.norm()){
		  WARN_LOG("RS coordinates containes nan number.");
		  for (int j = 0; j < y.size(); ++j){
			if(y[j] != y[j])
			  y[j] = 0.0f;
		  }
		}
		_Urs.col(i) = y;
	  }else{
		_Urs.col(i).setZero();
		WARN_LOG("column "<<i<<"of U is too small");
	  }
	}
  }
  void PCAonRS(){
	// given _Urs, make PCA to compute _Wrs and _zrs 
	TRACE_FUN();
	assert(!(_Urs.norm() != _Urs.norm()));
	assert_gt(_Urs.norm(),0);

	JacobiSVD<MatrixXd> svd(_Urs, ComputeThinU | ComputeThinV);
	const int nsigv = svd.singularValues().size();
	const int reserved = _PCADim<= nsigv? _PCADim:nsigv;
	const VectorXd sigv = svd.singularValues().segment(0,reserved);
	_Wrs = svd.matrixU().leftCols(reserved);
	_zrs = svd.matrixV().leftCols(reserved).transpose();
	// cout << sigv << endl;
	for (int i = 0; i < reserved; ++i)
	  _zrs.row(i) *= sigv[i];
	INFO_LOG("Wz-U" << ((_Wrs*_zrs) - _Urs).norm());
  }
  void computeK(){
	// given _zrs,h, compute _K
	TRACE_FUN();
	assert_gt(_h,0.0f);
	assert_ge(_zrs.cols(),3);
	const int T = _zrs.cols();
	const int r = _zrs.rows();
	CasADi::SXMatrix K;
	vector<CasADi::SX> s;
	produceSymetricMat("k",r,s,K);
	assert_eq(K.size1(),r);
	assert_eq(K.size2(),r);
	
	vector<CasADi::SXMatrix> z(T);
	for (int i = 0; i < T; ++i){
	  convert((VectorXd)_zrs.col(i),z[i]);
	  assert_eq(z[i].size1(),r);
	  assert_eq(z[i].size2(),1);
	}

	CasADi::SXMatrix E = 0;
	for (int i = 1; i < T-1; ++i){
	  const CasADi::SXMatrix za = (z[i+1]-z[i]*2.0f+z[i-1])/(_h*_h);
	  const CasADi::SXMatrix d = za+K.mul(z[i]);
	  for (int j = 0; j < r; ++j)
		E += d.elem(j,0)*d.elem(j,0);
	}
	CasADi::SXFunction E_fun(s,E);
	E_fun.init();

	const CasADi::SXMatrix G_SX = E_fun.grad();
	CasADi::SXFunction G_fun(s,G_SX);
	G_fun.init();
	const VectorXd x0 = VectorXd::Zero(s.size());
	VectorXd b;
	evaluate(G_fun,x0,b);
	// cout<< "b: " << b << endl << endl;
	assert_gt(b.norm(),0);

	const CasADi::SXMatrix H_SX = E_fun.hess();
	const MatrixXd H = convert<double>(H_SX);
	assert_eq(H.rows(),b.rows());
	// cout<< "H: " << H << endl;

	JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
	const VectorXd k = svd.solve(-b);
	cout<< "Hx+b: " << (H*k + b).norm() << endl;

	VectorXd diff;
	evaluate(E_fun,k,diff);
	cout << "diff: " << diff << endl;

	_K.resize(r,r);
	for (int i = 0; i < r; ++i)
	  for (int j = 0; j < r; ++j)
		_K(i,j) = k[symIndex(i,j)];
	assert_eq(_K,(_K.transpose()));
	// cout<< "_K: " << _K << endl;

	/////////////////// test
	const VectorXd x = VectorXd::Random(s.size());
	VectorXd out,zeroOut;
	evaluate(E_fun,x,out);
	evaluate(E_fun,x0,zeroOut);
	ASSERT_EQ(out.size(),1);
	ASSERT_EQ(zeroOut.size(),1);
	const double val = (x.transpose()*H*x)(0,0)*0.5f+b.transpose()*x+zeroOut[0];
	ASSERT_EQ_TOL( val,(out[0]),(1e-10*val));
	//////////////////////
  }
  void decomposeK(){
	// given _Wrs, _K, compute _Lambda, _B, and _W
	TRACE_FUN();
	SelfAdjointEigenSolver<MatrixXd> eigensolver(_K);
	assert_eq(eigensolver.info(),Success);
	const MatrixXd &W = eigensolver.eigenvectors();
	const VectorXd &La = eigensolver.eigenvalues();
	// cout<< "La: " << La << endl;

	int start = 0;
	for (;start<La.size() && La[start]<=0;++start);
	assert_lt(start,La.size());

	_Lambda = La.segment(start,La.size()-start);
	_B = W.block(0,start,W.rows(),W.cols()-start);
	_W = _Wrs*_B;
  }
  static void produceSymetricMat(const string name,const int dim,vector<CasADi::SX> &s,CasADi::SXMatrix &SM){

	assert_gt(dim,0);
	assert_gt(name.size(),0);
	s = makeSymbolic(dim*(1+dim)/2,name);

	SM.makeDense(dim, dim, 0.0f);
	for (int i = 0; i < SM.size1(); ++i){
	  for (int j = 0; j < SM.size2(); ++j){
		const int n = symIndex(i,j);
		assert_in(n,0,s.size()-1);
		SM.elem(i,j) = s[n];
	  }
	}
  }
  static int symIndex(int r,int c){
	// 0
	// 1,2
	// 3,4,5
	// 6,7,8,9
	const int rr = r>=c?r:c;
	const int cc = r>=c?c:r;
	const int n = (rr+1)*(rr+2)/2-1-(rr-cc);
	return n;
  }
  
public:
  double _h;
  double _alphaK;
  double _alphaM;
  int _PCADim;
  pVolumetricMesh _volMesh;
  MatrixXd _U;

  MatrixXd _Urs;
  MatrixXd _Wrs;
  MatrixXd _zrs;
  MatrixXd _K;
  MatrixXd _B;
  VectorXd _Lambda;
  MatrixXd _W;
};

BOOST_AUTO_TEST_SUITE(ElasticMtlTest)

BOOST_AUTO_TEST_CASE(symIndexTest){

  MatrixXd M(3,3);
  M << 1,2,3,  2,4,5, 3,5,6;
  ASSERT_EQ(M,M.transpose());
  
  VectorXd v(6);
  v << 1,2,4,3,5,6;
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
	  const int n = ElasticMtlOpt::symIndex(i,j);
	  ASSERT_GE(n,0);
	  ASSERT_LT(n,v.size());
	  ASSERT_EQ(v[n],M(i,j));
	}
  }  
}

BOOST_AUTO_TEST_CASE(produceSymetricMatTest){

  vector<CasADi::SX> s;
  CasADi::SXMatrix K;
  int dim = 2;
  ElasticMtlOpt::produceSymetricMat("x",2,s,K);
  ASSERT_EQ(K.elem(0,0).toString(),std::string("x_0"));
  ASSERT_EQ(K.elem(0,1).toString(),std::string("x_1"));
  ASSERT_EQ(K.elem(1,0).toString(),std::string("x_1"));
  ASSERT_EQ(K.elem(1,1).toString(),std::string("x_2"));
}

BOOST_AUTO_TEST_CASE(produceCorrectData){
  
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string volstr = data+"/sim-mesh.hinp";
  pVolumetricMesh volMesh;
  volMesh=pVolumetricMesh(VolumetricMeshLoader::load(volstr.c_str()));
  ASSERT(volMesh != NULL);

  const string wstr = data+"eigen_vectors80.b";
  MatrixXd W_t;
  ASSERT( load(wstr,W_t));
  const int r = 20;
  const MatrixXd W = W_t.leftCols(r);

  SparseMatrix<double> G;
  ASSERT( LSW_WARPING::DefGradOperator::compute(volMesh,G) );

  MatrixXd PGW;
  LSW_WARPING::MapMA2RS::computeMapMatPGW(G,W,PGW);
  const string pgw = data+"/tempt/PGW.b";
  ASSERT( write(pgw,PGW) );
  
  const MatrixXd invPGW = PseudoInverse(PGW);
  const string invpgw = data+"/tempt/invPGW.b";
  ASSERT( write(invpgw,invPGW));
}

BOOST_AUTO_TEST_CASE(computeTest){
  
  ElasticMtlOpt mtlOpt(0.01f,20);
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  ASSERT(mtlOpt.loadVolMesh(data+"/sim-mesh.hinp"));
  ASSERT(mtlOpt.loadAniSeq(data+"mtlOptU.b"));
  const MatrixXd U = mtlOpt._U.leftCols(50);
  mtlOpt._U = U;
  cout<< "mtlOptU.norm(): " << mtlOpt._U.norm () << endl;
  mtlOpt.compute();
  // cout<< mtlOpt._zrs.cols() << endl;
  // cout<< mtlOpt._zrs.rows() << endl;
}

BOOST_AUTO_TEST_SUITE_END()

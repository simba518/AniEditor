#include <MtlOptEnergy.h>
using namespace LSW_ANI_EDITOR;

void MtlOptEnergy::setZ(const MatrixXd &Z){

  const int T = getT();
  const int r = reducedDim();
  assert_eq(Z.rows(),r);
  assert_eq(Z.cols(),T);
  const MatrixXd &Z0 = Z.block(0,0,r,T-2);
  const MatrixXd &Z2 = Z.block(0,2,r,T-2);

  _Z1 = Z.block(0,1,r,T-2);
  _Z2v = (Z2-_Z1)/_h;
  _Z2a = (Z2-2*_Z1+Z0)/(_h*_h);
}

double MtlOptEnergy::fun(const double *x){

  computeKD(x,_K,_D);
  const double objValue = (_Z2a+_D*_Z2v+_K*_Z1).norm();
  return objValue*objValue*0.5f;
}

void MtlOptEnergy::grad(const double *x,double *g){
  
  const int r = reducedDim();
  const double &ak = x[0];
  const double &am = x[1];
  const MatrixXd &A = Map<MatrixXd>(const_cast<double*>(&x[2]),r,r);
  _K = A.transpose()*A;

  {// grad ak, am
	const MatrixXd F1 = _K*_Z1+_Z2a;
	MatrixXd Fm = ak*_K*_Z2v+F1;
	MatrixXd Gm = _K*_Z2v;
	MatrixXd Hm = am*_Z2v+F1;
	const VectorXd &E = Map<VectorXd>(&_Z2v(0,0),_Z2v.size());
	const VectorXd &F = Map<VectorXd>(&Fm(0,0),Fm.size());
	const VectorXd &G = Map<VectorXd>(&Gm(0,0),Gm.size());
	const VectorXd &H = Map<VectorXd>(&Hm(0,0),Hm.size());
	g[1] = (E.transpose()*(am*E+F))(0,0); // grad am
	g[0] = (G.transpose()*(ak*G+H))(0,0); // grad ak
  }

  {// grad A.
	const MatrixXd &E = _Z2v;
	const MatrixXd &At = A.transpose();
	const double &J = ak;
	const MatrixXd &L = _Z1;
	const MatrixXd M = am*_Z2v+_Z2a;
	const MatrixXd &Et = E.transpose();
	const double &Jt = J;
	const MatrixXd &Lt = L.transpose();
	const MatrixXd &Mt = M.transpose();
	const MatrixXd gA = (At*A*L*Lt*At).transpose()+(A*L*Et*At*A*Jt)+(A*L*Mt)+
	  (A*At*A*L*Lt)+(L*Et*At*A*Jt*At).transpose()+(A*M*Lt)+
	  (A*Jt*At*A*L*Et)+(E*Et*At*A*Jt*J*At).transpose()+(A*Jt*M*Et)+
	  (Jt*At*A*L*Et*At).transpose()+(A*E*Et*At*A*Jt*J)+(A*E*Mt*J);
	memcpy(&g[2],&gA(0,0),sizeof(double)*r*r);
  }
}

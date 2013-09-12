#ifndef _REDRSWARPEREXT_H_
#define _REDRSWARPEREXT_H_

#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <Timer.h>
#include <RSWarperExt.h>
#include <ComputeBj.h>
#include <SparseMatrixTools.h>
#include <RSCoordComp.h>
#include <MapMA2RS.h>
using namespace IEDS;

namespace LSW_WARPING{
  
  /**
   * @class RedRSWarperGrad compute the gradient of each tetrahedron with
   * respect to z using the code generated by maxima, thus very fast.
   * 
   */
  class RedRSWarperGrad{
	
  public:
	inline void jacobian(const VectorXd &y,
						 const vector<double> &Sqrt_V,
						 const int nodeId,
						 const MatrixXd &hatW,
						 const MatrixXd &BP,
						 Eigen::Matrix<double,3,-1> &Jz){

	  const int elemNum = y.size()/9;
	  static Eigen::Matrix<double,9,9> J;
	  dBdz.resize(y.size(),hatW.cols());

	  // i = 0
	  grad(Sqrt_V[0],&y[0],&J(0,0));
	  Jz = (BP.block(nodeId*3,0,3,9)*J)*hatW.block(0,0,9,hatW.cols());

	  for (int i = 1; i < elemNum; ++i){
		grad(Sqrt_V[i],&y[i*9],&J(0,0));
		Jz += (BP.block(nodeId*3,i*9,3,9)*J)*hatW.block(9*i,0,9,hatW.cols());
		// dBdz.block(i*9,0,9,hatW.cols()) = J*
	  }
	  // Jz = BP.block(nodeId*3,0,3,BP.cols())*dBdz;
	}

	// compute the gradient of one b(y), generated by maxima.
	inline void grad(double V,const double *s,double *out){
  
	  //declare temporary variables
	  double x1=(s[0]*s[0]*s[0]*s[0]);
	  double x2=(s[0]*s[0]);
	  double x3=(s[1]*s[1]);
	  double x4=(s[1]*s[1]*s[1]*s[1]);
	  double x5=2*x2;
	  double x6=2*x3;
	  double x7=x6+x5;
	  double x8=(s[2]*s[2]);
	  double x9=(s[2]*s[2]*s[2]*s[2]);
	  double x10=1/(x9+x7*x8+x4+2*x2*x3+x1);
	  double x11=2*s[0]*x3;
	  double x12=2*s[0]*x8;
	  double x13=-2*s[0]*x8;
	  double x14=x8+x3+x2;
	  double x15=sqrt(x14);
	  double x16=cos(x15);
	  double x17=(x13-2*s[0]*x3)*x16;
	  double x18=x17+x12+x11;
	  double x19=-x2*s[1];
	  double x20=(s[1]*s[1]*s[1]);
	  double x21=s[1]*x8;
	  double x22=x2*s[1];
	  double x23=-x20;
	  double x24=(s[0]*s[0]*s[0]);
	  double x25=-x24;
	  double x26=-s[0]*x3;
	  double x27=x26+x25;
	  double x28=x27*s[2];
	  double x29=-s[1]*x8;
	  double x30=(s[2]*s[2]*s[2]);
	  double x31=-s[0]*x30;
	  double x32=(x31+x29+x28+x23+x22)*x16;
	  double x33=x32+x21+x20+x19;
	  double x34=-x2;
	  double x35=(x3+x34)*s[2];
	  double x36=x24*s[1];
	  double x37=s[0]*x20;
	  double x38=-x3;
	  double x39=(x38+x2)*s[2];
	  double x40=s[0]*s[1]*x8;
	  double x41=-x30;
	  double x42=(x41+x40+x39+x37+x36)*x16;
	  double x43=x42+x30+x35;
	  double x44=-s[0]*x8;
	  double x45=x44+x26;
	  double x46=sin(x15);
	  double x47=s[0]*s[2];
	  double x48=x47+x22;
	  double x49=-s[0]*s[1];
	  double x50=x2*s[2]+x49;
	  double x51=1/x15;
	  double x52=(x31+x21+x28+x20+x19)*x16;
	  double x53=x52+x29+x23+x22;
	  double x54=-2*s[0]*x3*x16;
	  double x55=x54+x11;
	  double x56=2*s[0]*s[1]*s[2];
	  double x57=x2*x3;
	  double x58=-2*s[0]*s[1]*s[2];
	  double x59=x2*x8;
	  double x60=(x59+x58+x57+x1)*x16;
	  double x61=x60+x56;
	  double x62=-x1*s[1];
	  double x63=-x2*x20;
	  double x64=s[0]*x3;
	  double x65=x64+x24;
	  double x66=x65*s[2];
	  double x67=s[0]*x30;
	  double x68=x67-x2*s[1]*x8+x66+x63+x62;
	  double x69=pow(s[0],5);
	  double x70=x24*x3;
	  double x71=s[0]*x9;
	  double x72=x71+(x64+2*x24)*x8+x70+x69;
	  double x73=-x24*s[1];
	  double x74=-s[0]*x20;
	  double x75=(x74+x73)*s[2];
	  double x76=(x6+x2)*x8;
	  double x77=-s[0]*s[1]*x30;
	  double x78=x9+x77+x76+x75+x4+x57;
	  double x79=-s[0]*s[1]*x8;
	  double x80=(x41+x79+x39+x74+x73)*x16;
	  double x81=x80+x30+x35;
	  double x82=(x59+x56+x57+x1)*x16;
	  double x83=x82+x58;
	  double x84=2*s[0]*x8*x16;
	  double x85=x84+x13;
	  double x86=(x57+x1)*s[2];
	  double x87=x2*x30+x40+x86+x37+x36;
	  double x88=(x37+x36)*s[2];
	  double x89=s[0]*s[1]*x30;
	  double x90=x9+x89+x76+x88+x4+x57;
	  double x91=x27*x8-s[0]*x4-2*x24*x3-x69;
	  double x92=-2*x2*s[1];
	  double x93=2*x2*s[1]*x16;
	  double x94=x93+x92;
	  double x95=s[0]*x8;
	  double x96=x23+x19;
	  double x97=x96*s[2];
	  double x98=-s[1]*x30;
	  double x99=(x98+x44+x97+x64+x25)*x16;
	  double x100=x99+x95+x26+x24;
	  double x101=(x3*x8+x56+x4+x57)*x16;
	  double x102=x101+x58;
	  double x103=-pow(s[1],5);
	  double x104=-s[1]*x9;
	  double x105=x104+(x19-2*x20)*x8+x103+x63;
	  double x106=(x20+x22)*s[2];
	  double x107=s[1]*x30;
	  double x108=x107+s[0]*x3*x8+x106+s[0]*x4+x70;
	  double x109=x3+x5;
	  double x110=x9+x89+x109*x8+x88+x57+x1;
	  double x111=(x107+x44+x106+x64+x25)*x16;
	  double x112=x111+x95+x26+x24;
	  double x113=2*x2*s[1];
	  double x114=2*s[1]*x8;
	  double x115=-2*s[1]*x8;
	  double x116=(x115+x92)*x16;
	  double x117=x116+x114+x113;
	  double x118=(x41+x79+x35+x74+x73)*x16;
	  double x119=x118+x30+x39;
	  double x120=-s[1]*s[2];
	  double x121=x120+x64;
	  double x122=x29+x19;
	  double x123=x3*s[2]+s[0]*s[1];
	  double x124=-x2*x3;
	  double x125=-x4;
	  double x126=(-x3*x8+x56+x125+x124)*x16;
	  double x127=x126+x58;
	  double x128=(x41+x40+x35+x37+x36)*x16;
	  double x129=x128+x30+x39;
	  double x130=2*s[1]*x8*x16;
	  double x131=x130+x115;
	  double x132=-x1;
	  double x133=-2*x2;
	  double x134=-x9;
	  double x135=x134+x89+(x38+x133)*x8+x88+x124+x132;
	  double x136=x3*x30+x79+(x4+x57)*s[2]+x74+x73;
	  double x137=x96*x8;
	  double x138=x137+x103-2*x2*x20+x62;
	  double x139=-2*x2*s[2];
	  double x140=2*x2*s[2]*x16;
	  double x141=x140+x139;
	  double x142=x38+x34;
	  double x143=x142*x8;
	  double x144=(x134+x143+x56)*x16;
	  double x145=x144+x58;
	  double x146=(x107+x95+x106+x26+x25)*x16;
	  double x147=x146+x44+x64+x24;
	  double x148=-2*x3;
	  double x149=pow(s[2],5);
	  double x150=-x149+(x148+x34)*x30+(x125+x124)*s[2];
	  double x151=-2*x2*x3;
	  double x152=x89+x143+x88+x125+x151+x132;
	  double x153=x71+x98+x65*x8+x97;
	  double x154=(x134+x143+x58)*x16;
	  double x155=x154+x56;
	  double x156=2*x3*s[2];
	  double x157=-2*x3*s[2]*x16;
	  double x158=x157+x156;
	  double x159=(x67+x29+x66+x20+x22)*x16;
	  double x160=x159+x21+x23+x19;
	  double x161=x77+x143+x75+x125+x151+x132;
	  double x162=x149+x109*x30+x86;
	  double x163=x104+x31+x137+x28;
	  double x164=(x107+x44+x106+x64+x24)*x16;
	  double x165=x164+x95+x26+x25;
	  double x166=(x31+x29+x28+x20+x22)*x16;
	  double x167=x166+x21+x23+x19;
	  double x168=(x148+x133)*s[2];
	  double x169=x7*s[2]*x16;
	  double x170=x169+x168;
	  double x171=x44+x120;
	  double x172=x29+x47;
	  double x173=x3+x2;
	  double x174=1/x14;
	  double x175=x174*((x8+x3)*x16+x2)*V;
	  double x176=1/(x15*x15*x15);
	  double x177=x15*(s[0]*s[1]*x16+x49)*V;
	  double x178=-x176*((x41+x142*s[2])*x46*V+x177);
	  double x179=x15*(s[0]*s[2]*x16-s[0]*s[2])*V;
	  double x180=-x176*((x21+x20+x22)*x46*V+x179);
	  double x181=-x176*((x30+x173*s[2])*x46*V+x177);
	  double x182=x174*((x8+x2)*x16+x3)*V;
	  double x183=x15*(s[1]*s[2]*x16+x120)*V;
	  double x184=-x176*((x44+x26+x25)*x46*V+x183);
	  double x185=-x176*((x29+x23+x19)*x46*V+x179);
	  double x186=-x176*((x95+x64+x24)*x46*V+x183);
	  double x187=x174*(x173*x16+x8)*V;

	  //outputs
	  out[0]=x10*(x15*(x50*x46*s[5]+x48*x46*s[4]+x45*x46*s[3]+x45*x46)*V+(x43*s[5]+x33*s[4]+x18*s[3]+x17+x12+x11)*V);
	  out[1]=x10*(x15*(x50*x46*s[7]+x48*x46*s[6]+x45*x46*s[4]+x48*x46)*V+(x43*s[7]+x33*s[6]+x18*s[4]+x32+x21+x20+x19)*V);
	  out[2]=x10*(x15*(x50*x46*s[8]+x48*x46*s[7]+x45*x46*s[5]+x50*x46)*V+(x43*s[8]+x33*s[7]+x18*s[5]+x42+x30+x35)*V);
	  out[3]=-x51*x10*((x78*x46*s[5]+x72*x46*s[4]+x68*x46*s[3]+x68*x46)*V+x15*(x61*s[5]+x55*s[4]+x53*s[3]+x52+x29+x23+x22)*V);
	  out[4]=-x51*x10*((x78*x46*s[7]+x72*x46*s[6]+x68*x46*s[4]+x72*x46)*V+x15*(x61*s[7]+x55*s[6]+x53*s[4]+x54+x11)*V);
	  out[5]=-x51*x10*((x78*x46*s[8]+x72*x46*s[7]+x68*x46*s[5]+x78*x46)*V+x15*(x61*s[8]+x55*s[7]+x53*s[5]+x60+x56)*V);
	  out[6]=x51*x10*((x91*x46*s[5]+x90*x46*s[4]+x87*x46*s[3]+x87*x46)*V+x15*(x85*s[5]+x83*s[4]+x81*s[3]+x80+x30+x35)*V);
	  out[7]=x51*x10*((x91*x46*s[7]+x90*x46*s[6]+x87*x46*s[4]+x90*x46)*V+x15*(x85*s[7]+x83*s[6]+x81*s[4]+x82+x58)*V);
	  out[8]=x51*x10*((x91*x46*s[8]+x90*x46*s[7]+x87*x46*s[5]+x91*x46)*V+x15*(x85*s[8]+x83*s[7]+x81*s[5]+x84+x13)*V);
	  out[9]=x51*x10*((x110*x46*s[5]+x108*x46*s[4]+x105*x46*s[3]+x105*x46)*V+x15*(x102*s[5]+x100*s[4]+x94*s[3]+x93+x92)*V);
	  out[10]=x51*x10*((x110*x46*s[7]+x108*x46*s[6]+x105*x46*s[4]+x108*x46)*V+x15*(x102*s[7]+x100*s[6]+x94*s[4]+x99+x95+x26+x24)*V);
	  out[11]=x51*x10*((x110*x46*s[8]+x108*x46*s[7]+x105*x46*s[5]+x110*x46)*V+x15*(x102*s[8]+x100*s[7]+x94*s[5]+x101+x58)*V);
	  out[12]=x10*(x15*(x123*x46*s[5]+x122*x46*s[4]+x121*x46*s[3]+x121*x46)*V+(x119*s[5]+x117*s[4]+x112*s[3]+x111+x95+x26+x24)*V);
	  out[13]=x10*(x15*(x123*x46*s[7]+x122*x46*s[6]+x121*x46*s[4]+x122*x46)*V+(x119*s[7]+x117*s[6]+x112*s[4]+x116+x114+x113)*V);
	  out[14]=x10*(x15*(x123*x46*s[8]+x122*x46*s[7]+x121*x46*s[5]+x123*x46)*V+(x119*s[8]+x117*s[7]+x112*s[5]+x118+x30+x39)*V);
	  out[15]=x51*x10*((x138*x46*s[5]+x136*x46*s[4]+x135*x46*s[3]+x135*x46)*V+x15*(x131*s[5]+x129*s[4]+x127*s[3]+x126+x58)*V);
	  out[16]=x51*x10*((x138*x46*s[7]+x136*x46*s[6]+x135*x46*s[4]+x136*x46)*V+x15*(x131*s[7]+x129*s[6]+x127*s[4]+x128+x30+x39)*V);
	  out[17]=x51*x10*((x138*x46*s[8]+x136*x46*s[7]+x135*x46*s[5]+x138*x46)*V+x15*(x131*s[8]+x129*s[7]+x127*s[5]+x130+x115)*V);
	  out[18]=x51*x10*((x153*x46*s[5]+x152*x46*s[4]+x150*x46*s[3]+x150*x46)*V+x15*(x147*s[5]+x145*s[4]+x141*s[3]+x140+x139)*V);
	  out[19]=x51*x10*((x153*x46*s[7]+x152*x46*s[6]+x150*x46*s[4]+x152*x46)*V+x15*(x147*s[7]+x145*s[6]+x141*s[4]+x144+x58)*V);
	  out[20]=x51*x10*((x153*x46*s[8]+x152*x46*s[7]+x150*x46*s[5]+x153*x46)*V+x15*(x147*s[8]+x145*s[7]+x141*s[5]+x146+x44+x64+x24)*V);
	  out[21]=-x51*x10*((x163*x46*s[5]+x162*x46*s[4]+x161*x46*s[3]+x161*x46)*V+x15*(x160*s[5]+x158*s[4]+x155*s[3]+x154+x56)*V);
	  out[22]=-x51*x10*((x163*x46*s[7]+x162*x46*s[6]+x161*x46*s[4]+x162*x46)*V+x15*(x160*s[7]+x158*s[6]+x155*s[4]+x157+x156)*V);
	  out[23]=-x51*x10*((x163*x46*s[8]+x162*x46*s[7]+x161*x46*s[5]+x163*x46)*V+x15*(x160*s[8]+x158*s[7]+x155*s[5]+x159+x21+x23+x19)*V);
	  out[24]=-x10*(x15*(x173*s[2]*x46*s[5]+x172*x46*s[4]+x171*x46*s[3]+x171*x46)*V+(x170*s[5]+x167*s[4]+x165*s[3]+x164+x95+x26+x25)*V);
	  out[25]=-x10*(x15*(x173*s[2]*x46*s[7]+x172*x46*s[6]+x171*x46*s[4]+x172*x46)*V+(x170*s[7]+x167*s[6]+x165*s[4]+x166+x21+x23+x19)*V);
	  out[26]=-x10*(x15*(x173*s[2]*x46*s[8]+x172*x46*s[7]+x171*x46*s[5]+x173*s[2]*x46)*V+(x170*s[8]+x167*s[7]+x165*s[5]+x169+x168)*V);
	  out[27]=x175;
	  out[28]=0;
	  out[29]=0;
	  out[30]=x178;
	  out[31]=0;
	  out[32]=0;
	  out[33]=x180;
	  out[34]=0;
	  out[35]=0;
	  out[36]=x181;
	  out[37]=x175;
	  out[38]=0;
	  out[39]=x182;
	  out[40]=x178;
	  out[41]=0;
	  out[42]=x184;
	  out[43]=x180;
	  out[44]=0;
	  out[45]=x185;
	  out[46]=0;
	  out[47]=x175;
	  out[48]=x186;
	  out[49]=0;
	  out[50]=x178;
	  out[51]=x187;
	  out[52]=0;
	  out[53]=x180;
	  out[54]=0;
	  out[55]=x181;
	  out[56]=0;
	  out[57]=0;
	  out[58]=x182;
	  out[59]=0;
	  out[60]=0;
	  out[61]=x184;
	  out[62]=0;
	  out[63]=0;
	  out[64]=x185;
	  out[65]=x181;
	  out[66]=0;
	  out[67]=x186;
	  out[68]=x182;
	  out[69]=0;
	  out[70]=x187;
	  out[71]=x184;
	  out[72]=0;
	  out[73]=0;
	  out[74]=x185;
	  out[75]=0;
	  out[76]=0;
	  out[77]=x186;
	  out[78]=0;
	  out[79]=0;
	  out[80]=x187;
	}

  private:
	MatrixXd dBdz;
  };

  /**
   * @class RedRSWarperExt reduced RS warper.
   * 
   */
  class RedRSWarperExt{
	
  public:
	RedRSWarperExt(pRSWarperExt full_warper,const MatrixXd &NLBasis,const MatrixXd &LBasis){
	  init(full_warper,NLBasis,LBasis);
	  BP = B*P;
	}
	RedRSWarperExt(pRSWarperExt full_warper,const MatrixXd &NLBasis,const MatrixXd &LBasis,
				   const vector<int> &selPoints, const vector<double> &cubWeights){
	  
	  init(full_warper,NLBasis,LBasis);
	  setCubature(selPoints, cubWeights);
	}

	// warp all nodes using z of frame f.
	void warp(const VectorXd &z,int frame_id, VectorXd &u){

	  static VectorXd b;
	  computeB(z,frame_id,b);
	  assert_eq(P.cols(),b.size());
	  const VectorXd q = P*b;
	  u = B*q;
	}

	// warp one node using z of frame f.
	Vector3d warp(const VectorXd &z,int frame_id, int node){

	  static VectorXd b;
	  computeB(z,frame_id,b);
	  assert_eq(BP.cols(),b.size());
	  const Vector3d u3 = BP.block(node*3,0,3,BP.cols())*b;
	  return u3;
	}

	// return the jacobian of one node with respect to z of frame f
	void jacobian(const VectorXd &z,int frame_id,int node,Eigen::Matrix<double,3,-1> &J){

	  VectorXd y = hatW*z + ref_y[frame_id];
	  for (int i = 0; i < y.size(); i+=9){
		y[i] += 1e-8;
	  }
	  grad.jacobian(y,Sqrt_V,node,hatW,BP,J);
	}

  protected:
	void computeB(const VectorXd &z,int frame_id, VectorXd &b){

	  VectorXd y = hatW*z + ref_y[frame_id];
	  b.resize(y.size());
	  for (int i = 0; i < y.size(); i+=9){
		y[i] += 1e-10; //@todo
		ComputeBj::compute(&y[i],Sqrt_V[i/9],&b[i]);
	  }
	}
	void setCubature(const vector<int> &selPoints, const vector<double> &cubWeights){
	  
	  // re-shape matrix
	  SparseMatrix<double> S;
	  set<int> selPointsSet;
	  for (size_t i = 0; i < selPoints.size(); ++i){
		assert_in(selPoints[i],0,P.cols()/9);
		selPointsSet.insert(selPoints[i]);
	  }
	  EIGEN3EXT::genReshapeMatrix(P.cols(),9,selPointsSet,S,false);
	  
	  // Pj, BPj, hatWj
	  assert_eq(selPoints.size(),cubWeights.size());
	  for (size_t i = 0; i < selPointsSet.size(); ++i){
		assert_gt(cubWeights[i],0);
		P.block(0,selPoints[i]*9,P.rows(),9) *= cubWeights[i];
	  }
	  P = P*S.transpose();
	  BP = B*P;
	  hatW = S*hatW;

	  // Sqrt_Vj
	  vector<int> sorted_selPoints = selPoints;
	  std::sort(sorted_selPoints.begin(), sorted_selPoints.end());
	  vector<double> tempt(sorted_selPoints.size());
	  for (size_t i = 0; i < sorted_selPoints.size(); ++i){
		assert_in(sorted_selPoints[i],0,(int)Sqrt_V.size());
		tempt[i] = Sqrt_V[sorted_selPoints[i]];
	  }
	  Sqrt_V = tempt;	  

	  // ref_y
	  VectorXd yf(sorted_selPoints.size()*9);
	  for (size_t f = 0; f < ref_y.size(); ++f){
		for (size_t i = 0; i < sorted_selPoints.size(); ++i){
		  yf.segment(i*9,9) = ref_y[f].segment(sorted_selPoints[i]*9,9);
		}
		ref_y[f] = yf;
	  }
	  
	}
	void init(pRSWarperExt full_warper,const MatrixXd &NLBasis,const MatrixXd &LBasis){

	  B = NLBasis;
	  assert(full_warper);
	  const LSW_WARPING::RS2Euler &rs2euler = full_warper->getRS2Euler();
	  const SparseMatrix<double> &A = rs2euler.get_L();
	  const MatrixXd T = (B.transpose()*(A*B));
	  P = T.inverse();
	  cout << "cond(P) = " << P.norm()*T.norm() << endl;
	  P = P*(B.transpose()*rs2euler.get_VG_t());
	  MapMA2RS::computeMapMatPGW(rs2euler.get_G(),LBasis,hatW);

	  Sqrt_V = rs2euler.get_Sqrt_V();
	  ref_y = full_warper->getInputRSCoord();
	}
	
  private:
	MatrixXd B; // nonlinear basis
	MatrixXd P; // (B^t*A*B)^{-1}B^T*((VG)^t)
	MatrixXd BP; // B*P.
	MatrixXd hatW;
	RedRSWarperGrad grad;
	vector<double> Sqrt_V;
	vector<VectorXd> ref_y;
  };

  typedef boost::shared_ptr<RedRSWarperExt> pRedRSWarperExt;

}//end of namespace

#endif /*_REDRSWARPEREXT_H_*/

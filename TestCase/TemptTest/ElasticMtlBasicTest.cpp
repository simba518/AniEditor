#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <UnitTestAssert.h>
#include <DrawCurves.h>
#include <AuxTools.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include <VolObjMesh.h>
#include <MASimulatorAD.h>
#include "ElasticMtlOpt.h"
#include "HarmonicOscillator.h"
#include "DFT.h"
#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace LSW_SIM;
using namespace ANI_EDIT;
using namespace UTILITY;
using namespace Eigen;

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
  const int dim = 2;
  ElasticMtlOpt::produceSymetricMat("x",dim,s,K);
  ASSERT_EQ(K.elem(0,0).toString(),std::string("x_0"));
  ASSERT_EQ(K.elem(0,1).toString(),std::string("x_1"));
  ASSERT_EQ(K.elem(1,0).toString(),std::string("x_1"));
  ASSERT_EQ(K.elem(1,1).toString(),std::string("x_2"));

  const int dim1 = 1;
  ElasticMtlOpt::produceSymetricMat("x",dim1,s,K);
  ASSERT_EQ(K.size1(),1);
  ASSERT_EQ(K.size2(),1);
  ASSERT_EQ(K.elem(0,0).toString(),std::string("x_0"));
}

BOOST_AUTO_TEST_SUITE_END()

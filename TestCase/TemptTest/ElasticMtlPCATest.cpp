#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <UnitTestAssert.h>
#include "ElasticMtlOpt.h"
using namespace ANI_EDIT;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(ElasticMtlTest)

BOOST_AUTO_TEST_CASE(checkPCAuseOneMode){

  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string hatWstr = data+"/tempt/PGW.b";
  MatrixXd hatW;
  TEST_ASSERT( load(hatWstr,hatW) );

  const int T = 100;
  const int mid = 6;
  
  ElasticMtlOpt mtlotp(0.1,1);
  mtlotp._Urs.resize(hatW.rows(),T);
  VectorXd z(T);
  for (int i = 0; i < T; ++i){
	z[i] = 0.02f*i;
	mtlotp._Urs.col(i) = hatW.col(mid)*z[i];
  }
  mtlotp.PCAonRS();
  cout << (mtlotp._zrs.row(0).transpose()-z*hatW.col(mid).norm()).norm() << endl;
  cout << (hatW.col(mid)/hatW.col(mid).norm()-mtlotp._Wrs.col(0)).norm() << endl;
}

BOOST_AUTO_TEST_CASE(checkPCAuseTwoModes){

  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string hatWstr = data+"/tempt/PGW.b";
  MatrixXd hatW;
  TEST_ASSERT( load(hatWstr,hatW) );

  const int T = 100;
  const int mid0 = 6;
  const int mid1 = 7;
  
  ElasticMtlOpt mtlotp(0.1,2);
  mtlotp._Urs.resize(hatW.rows(),T);
  MatrixXd z(2,T);
  for (int i = 0; i < T; ++i){
	z(0,i) = 0.02f*i;
	z(1,i) = 0.005f*i;
	mtlotp._Urs.col(i) = hatW.col(mid0)*z(0,i)+hatW.col(mid1)*z(1,i);
  }
  mtlotp.PCAonRS();

  cout<<mtlotp._Wrs.col(0).transpose()*mtlotp._Wrs.col(1)<<endl;
  cout<<hatW.col(mid0).normalized().transpose()*hatW.col(mid1).normalized()<<endl;

  cout<<(mtlotp._zrs.row(0)-z.row(0)*hatW.col(mid0).norm()).norm() << endl;
  cout<<(mtlotp._zrs.row(1)-z.row(1)*hatW.col(mid1).norm()).norm() << endl;

  cout<<(hatW.col(mid0).normalized()-mtlotp._Wrs.col(0)).norm()<<endl;
  cout<<(hatW.col(mid1).normalized()-mtlotp._Wrs.col(1)).norm()<<endl;
}

BOOST_AUTO_TEST_SUITE_END()

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

BOOST_AUTO_TEST_CASE(computeKuseRSSimulation){

  LSW_ANI_EDITOR::MASimulatorAD simulator;
  const double h = 0.5f;
  const int r = 1;
  VectorXd lambda(r);
  for (int i = 0; i < r; ++i){
	lambda[i] = 0.5f*(i+1.0f);
  }
  simulator.setTimeStep(h);
  simulator.setEigenValues(lambda);
  simulator.setStiffnessDamping(0.0f);
  simulator.setMassDamping(0.0f);
  const VectorXd v0 = VectorXd::Zero(r);
  VectorXd z0(r);
  for (int i = 0; i < r; ++i){
	z0[i] = 1.0f/(i+1.0f);
  }
  simulator.setIntialStatus(v0,z0);

  const int T = 50;
  ElasticMtlOpt mtlopt(h,r);
  mtlopt._zrs.resize(r,T);
  for (int f = 0; f < T; ++f){
	mtlopt._zrs.col(f) = simulator.getEigenZ();
    simulator.forward(VectorXd::Zero(r));
  }
  mtlopt.computeK();

  PythonScriptDraw2DCurves<VectorXd> curves;
  const VectorXd z = mtlopt._zrs.row(0).transpose();
  curves.add("z_old",z,h,0.0f,"--");

  HarmonicOscillator<double> fz(mtlopt._K(0,0),mtlopt._D(0,0),z0[0],v0[0]);
  const VectorXd z2 = fz.generateSequence<VectorXd>(0.0f,h,z.size());
  curves.add("z_new",z2,h,0.0f,"*");
  TEST_ASSERT( curves.write("a.py") );
  
  cout << mtlopt._K(0,0) << endl << endl;
  VectorXd outreal,outimag;
  computeDFT(z,outreal,outimag);
  curves.clear();
  curves.add("real",outreal,h,0,"o");
  curves.add("imag",outimag,h,0,"*");
  TEST_ASSERT( curves.write("b.py") );

}

BOOST_AUTO_TEST_SUITE_END()

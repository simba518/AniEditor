#ifndef _DFT_H_
#define _DFT_H_

#include <math.h>
#include <assertext.h>
#include <AuxTools.h>
using std::string;
using std::stringstream;

/* 
 * Computes the discrete Fourier transform (DFT) of the given vector.
 * All the array arguments must have the same length.
 * @see http://nayuki.eigenstate.org/page/how-to-implement-the-discrete-fourier-transform
 * 
 * @note the original implementation in the site is not correct, and the result
 * of this version is equal to the result of maxima.
 */
template <class VECTOR>
inline void computeDFT(const VECTOR &inreal, const VECTOR &inimag,VECTOR &outreal,VECTOR &outimag){

  assert_eq(inreal.size(), inimag.size());
  const double n = inreal.size();
  outreal.resize(n);
  outimag.resize(n);

  for (double k = 0.0f; k < n; k+=1.0f) {  /* For each output element */
	double sumreal = 0.0f;
	double sumimag = 0.0f;
	for (double t = 0.0f; t < n; t += 1.0f) {  /* For each input element */
	  sumreal+=(inreal[t]*cos(2*M_PI*t*k/n) - inimag[t]*sin(2*M_PI*t*k/n));
	  sumimag+=(inreal[t]*sin(2*M_PI*t*k/n) + inimag[t]*cos(2*M_PI*t*k/n));
	}
	outreal[k] = sumreal/n;
	outimag[k] = sumimag/n;
  }
}

template <class VECTOR>
inline void computeDFT(const VECTOR &inreal, VECTOR &outreal,VECTOR &outimag){

  VECTOR inimag;
  inimag.resize(inreal.size());
  for (int i = 0; i < inreal.size(); ++i)
    inimag[i] = 0;
  computeDFT(inreal,inimag,outreal,outimag);
}

class PythonFFT{
  
public:
  template <class VECTOR>
  void add(const string name,const VECTOR& s,const double h){
	stringstream data;
	data << name<<"=[";
	for (int i = 0; i < (int)s.size(); ++i){
	  data << std::setprecision(14) << s[i];
	  if((int)s.size()-1!=i) data<<", ";
	}
	data << "]\n";
	data << name<<"_h="<<h<<"\n\n";
	signalData += data.str();
	if(signalList.size()>0){
	  signalList = signalList+","+name;
	  stepList = stepList+","+name+"_h";
	}else{
	  stepList = name+"_h";
	  signalList = name;
	}
  }
  void write(const string fname,const string saveFigTo=""){
	OUTFILE(file,string(fname+".py").c_str());
	file << head();
	file << signalData;
	file << "\nsignalList=[" << signalList <<"]\n";
	file << "\nstepList=[" << stepList <<"]\n";
	file << end(saveFigTo);
  }
  void clear(){
	signalData = "";
	signalList = "";
	stepList = "";
  }
  
protected:
  static string head(){
	string h = "import numpy as np\n";
	h += "import scipy\n";
	h += "import scipy.fftpack\n";
	h += "import matplotlib\n";
	h += "import matplotlib.pyplot as plt\n";
	h += "import sys\n\n";
	return h;
  }
  static string end(const string saveFigTo){
	string e = "\nfor i in range(0,len(signalList)):\n";
    e += "\tsignal = signalList[i]\n";
    e += "\tfourier = np.fft.fft(signal)\n";
    e += "\tfreq = np.fft.fftfreq(len(signal),stepList[i])\n";
    e += "\tplt.stem(freq,np.absolute(fourier),\"-\")\n";
    e += "\tplt.grid(True)\n";
	if (saveFigTo.size()>0){
	  e += string("\tplt.savefig(\"")+saveFigTo+"_\"+str(i)+\".png\")\n";
	}else{
	  e += "\tplt.show()\n";
	}
    e += "\tplt.clf()";
	return e;
  }
  
private:
  string signalData;
  string signalList;
  string stepList;
};

#endif /* _DFT_H_ */

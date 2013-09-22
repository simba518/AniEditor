#ifndef _DFT_H_
#define _DFT_H_

#include <math.h>
#include <assertext.h>

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

#endif /* _DFT_H_ */

#ifndef _TENSORUTILS_H_
#define _TENSORUTILS_H_
#include "NuGamma/TensorInc.h"

namespace TensorUtils {

  // -----------------------------------------------------
  // Common tensors
  // -----------------------------------------------------

  TensorOrder2* LeviCivitaOrder2();
  TensorOrder3* LeviCivitaOrder3();
  TensorOrder4* LeviCivitaOrder4();
  TensorOrder2* Metric();
  TensorOrder2* UnityOrder2(unsigned int i = 4);
  TensorOrder3* UnityOrder3(unsigned int i = 4);
  TensorOrder4* UnityOrder4(unsigned int i = 4);
  TensorOrder2* DiracMatrix(unsigned int i);
  TensorOrder2* DiracTensor(unsigned int a, unsigned int b);

  // -----------------------------------------------------
  // Common utilities for tensors
  // -----------------------------------------------------

  TComplex   DotProd(const TensorOrder1 &t1, const TensorOrder1 &t2);
  TComplex   DotProdMetric(const TensorOrder1 &t1, const TensorOrder1 &t2);

  TensorOrder1* CrossProd(const TensorOrder1 &t1, const TensorOrder1 &t2);
  TensorOrder2* SlashDirac(const TensorOrder1 &t1);
}

#endif

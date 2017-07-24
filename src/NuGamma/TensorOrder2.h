#ifndef _TENSORORDER2_H_
#define _TENSORORDER2_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <complex>

#include "NuGamma/TensorInc.h"

namespace TensorUtils{
  

  class TensorOrder2: public Tensor{
    
  private:

    bool         CheckIndex(const unsigned int i, const unsigned int j) const;
    unsigned int GetGlobalIndex(const unsigned int i, const unsigned int j) const;

  public:

    TensorOrder2(const TensorOrder2& t1);
    TensorOrder2(const char * Name, unsigned int Dim1, unsigned int Dim2);
    TensorOrder2(const char * Name = "none", unsigned int Dim = 4);

    void Print() const;
    
    TComplex At(const unsigned int i, const unsigned int j) const;
    void Set(const unsigned int i, const unsigned int j, TComplex c);
    TComplex& operator()(const unsigned int i1, const unsigned int i2);
    TComplex  operator()(const unsigned int i1, const unsigned int i2) const;
  
    using Tensor::operator=;
    using Tensor::SumOver;
    using Tensor::Add;
    using Tensor::Subtract;
    using Tensor::Multiply;
    using Tensor::Divide;

    using Tensor::IsEqual;
    using Tensor::IsDifferent;
    using Tensor::IsBiggerOrEqual;
    using Tensor::IsSmallerOrEqual;
    using Tensor::IsBigger;
    using Tensor::IsSmaller;

    using Tensor::Conjugate;
    using Tensor::Real;
    using Tensor::Imaginary;
    using Tensor::Abs;
    using Tensor::Minus;
    using Tensor::OneOverElementWise;
    using Tensor::MultiplyElementWise;
    using Tensor::DivideElementWise;

    using Tensor::AssertGoodSize;
    using Tensor::GoodSize;

    friend TensorOrder2 Real(const TensorOrder2 &t1);
    friend TensorOrder2 Imaginary(const TensorOrder2 &t1);
    friend TensorOrder2 Abs(const TensorOrder2 &t1);
    friend TensorOrder2 operator-(const TensorOrder2 &t1);

    friend TensorOrder2 operator+(const int c, const TensorOrder2 &t1);
    friend TensorOrder2 operator-(const int c, const TensorOrder2 &t1);
    friend TensorOrder2 operator*(const int c, const TensorOrder2 &t1);
    friend TensorOrder2 operator/(const int c, const TensorOrder2 &t1);

    friend TensorOrder2 operator+(const TensorOrder2 &t1, const int c);
    friend TensorOrder2 operator-(const TensorOrder2 &t1, const int c);
    friend TensorOrder2 operator*(const TensorOrder2 &t1, const int c);
    friend TensorOrder2 operator/(const TensorOrder2 &t1, const int c);

    friend TensorOrder2 operator+(TComplex c, const TensorOrder2 &t1);
    friend TensorOrder2 operator-(TComplex c, const TensorOrder2 &t1);
    friend TensorOrder2 operator*(TComplex c, const TensorOrder2 &t1);
    friend TensorOrder2 operator/(TComplex c, const TensorOrder2 &t1);

    friend TensorOrder2 operator+(const TensorOrder2 &t1, TComplex c);
    friend TensorOrder2 operator-(const TensorOrder2 &t1, TComplex c);
    friend TensorOrder2 operator*(const TensorOrder2 &t1, TComplex c);
    friend TensorOrder2 operator/(const TensorOrder2 &t1, TComplex c);

    friend TensorOrder2 operator+(double d, const TensorOrder2 &t1);
    friend TensorOrder2 operator-(double d, const TensorOrder2 &t1);
    friend TensorOrder2 operator*(double d, const TensorOrder2 &t1);
    friend TensorOrder2 operator/(double d, const TensorOrder2 &t1);

    friend TensorOrder2 operator+(const TensorOrder2 &t1, double d);
    friend TensorOrder2 operator-(const TensorOrder2 &t1, double d);
    friend TensorOrder2 operator*(const TensorOrder2 &t1, double d);
    friend TensorOrder2 operator/(const TensorOrder2 &t1, double d);

    friend TensorOrder2 operator+(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend TensorOrder2 operator-(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend TensorOrder2 operator*(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend TensorOrder1 operator*(const TensorOrder2 &t1, const TensorOrder1 &t2);

    friend bool operator==(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend bool operator!=(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend bool operator>=(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend bool operator<=(const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend bool operator> (const TensorOrder2 &t1, const TensorOrder2 &t2);
    friend bool operator< (const TensorOrder2 &t1, const TensorOrder2 &t2);

    friend bool operator==(double d, const TensorOrder2 &t1);
    friend bool operator!=(double d, const TensorOrder2 &t1);
    friend bool operator>=(double d, const TensorOrder2 &t1);
    friend bool operator<=(double d, const TensorOrder2 &t1);
    friend bool operator> (double d, const TensorOrder2 &t1);
    friend bool operator< (double d, const TensorOrder2 &t1);

    friend bool operator==(const TensorOrder2 &t1, double d);
    friend bool operator!=(const TensorOrder2 &t1, double d);
    friend bool operator>=(const TensorOrder2 &t1, double d);
    friend bool operator<=(const TensorOrder2 &t1, double d);
    friend bool operator> (const TensorOrder2 &t1, double d);
    friend bool operator< (const TensorOrder2 &t1, double d);

    friend bool operator==(TComplex d, const TensorOrder2 &t1);
    friend bool operator!=(TComplex d, const TensorOrder2 &t1);
    friend bool operator>=(TComplex d, const TensorOrder2 &t1);
    friend bool operator<=(TComplex d, const TensorOrder2 &t1);
    friend bool operator> (TComplex d, const TensorOrder2 &t1);
    friend bool operator< (TComplex d, const TensorOrder2 &t1);

    friend bool operator==(const TensorOrder2 &t1, TComplex d);
    friend bool operator!=(const TensorOrder2 &t1, TComplex d);
    friend bool operator>=(const TensorOrder2 &t1, TComplex d);
    friend bool operator<=(const TensorOrder2 &t1, TComplex d);
    friend bool operator> (const TensorOrder2 &t1, TComplex d);
    friend bool operator< (const TensorOrder2 &t1, TComplex d);

    ClassDef(TensorOrder2,0);
  };
}
#endif

#ifndef _TENSORORDER1_H_
#define _TENSORORDER1_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <complex>

#include "NuGamma/TensorInc.h"

namespace TensorUtils{
 
  class TensorOrder1: public Tensor{

  private:

    bool CheckIndex(const unsigned int i) const;
    unsigned int  GetGlobalIndex(const unsigned int i) const;

  public:

    //TensorOrder1(const Tensor& t1);
    TensorOrder1(const TensorOrder1& t1);
    TensorOrder1(const char * Name = "none", unsigned int Dim = 4);

    void Print() const;
    
    TComplex At(const unsigned int i) const;
    void Set(const unsigned int i, TComplex c);
    TComplex& operator()(const unsigned int i);
    TComplex  operator()(const unsigned int i) const;

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

    friend TensorOrder1 Real(const TensorOrder1 &t1);
    friend TensorOrder1 Imaginary(const TensorOrder1 &t1);
    friend TensorOrder1 Abs(const TensorOrder1 &t1);
    friend TensorOrder1 operator-(const TensorOrder1 &t1);

    friend TensorOrder1 operator+(const int c, const TensorOrder1 &t1);
    friend TensorOrder1 operator-(const int c, const TensorOrder1 &t1);
    friend TensorOrder1 operator*(const int c, const TensorOrder1 &t1);
    friend TensorOrder1 operator/(const int c, const TensorOrder1 &t1);

    friend TensorOrder1 operator+(const TensorOrder1 &t1, const int c);
    friend TensorOrder1 operator-(const TensorOrder1 &t1, const int c);
    friend TensorOrder1 operator*(const TensorOrder1 &t1, const int c);
    friend TensorOrder1 operator/(const TensorOrder1 &t1, const int c);

    friend TensorOrder1 operator+(TComplex c, const TensorOrder1 &t1);
    friend TensorOrder1 operator-(TComplex c, const TensorOrder1 &t1);
    friend TensorOrder1 operator*(TComplex c, const TensorOrder1 &t1);
    friend TensorOrder1 operator/(TComplex c, const TensorOrder1 &t1);

    friend TensorOrder1 operator+(const TensorOrder1 &t1, TComplex c);
    friend TensorOrder1 operator-(const TensorOrder1 &t1, TComplex c);
    friend TensorOrder1 operator*(const TensorOrder1 &t1, TComplex c);
    friend TensorOrder1 operator/(const TensorOrder1 &t1, TComplex c);

    friend TensorOrder1 operator+(double d, const TensorOrder1 &t1);
    friend TensorOrder1 operator-(double d, const TensorOrder1 &t1);
    friend TensorOrder1 operator*(double d, const TensorOrder1 &t1);
    friend TensorOrder1 operator/(double d, const TensorOrder1 &t1);

    friend TensorOrder1 operator+(const TensorOrder1 &t1, double d);
    friend TensorOrder1 operator-(const TensorOrder1 &t1, double d);
    friend TensorOrder1 operator*(const TensorOrder1 &t1, double d);
    friend TensorOrder1 operator/(const TensorOrder1 &t1, double d);

    friend TensorOrder1 operator+(const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend TensorOrder1 operator-(const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend TensorOrder1 operator*(const TensorOrder2 &td2, const TensorOrder1 &td1);

    friend bool operator==(const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend bool operator!=(const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend bool operator>=(const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend bool operator<=(const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend bool operator> (const TensorOrder1 &t1, const TensorOrder1 &t2);
    friend bool operator< (const TensorOrder1 &t1, const TensorOrder1 &t2);

    friend bool operator==(double d, const TensorOrder1 &t1);
    friend bool operator!=(double d, const TensorOrder1 &t1);
    friend bool operator>=(double d, const TensorOrder1 &t1);
    friend bool operator<=(double d, const TensorOrder1 &t1);
    friend bool operator> (double d, const TensorOrder1 &t1);
    friend bool operator< (double d, const TensorOrder1 &t1);

    friend bool operator==(const TensorOrder1 &t1, double d);
    friend bool operator!=(const TensorOrder1 &t1, double d);
    friend bool operator>=(const TensorOrder1 &t1, double d);
    friend bool operator<=(const TensorOrder1 &t1, double d);
    friend bool operator> (const TensorOrder1 &t1, double d);
    friend bool operator< (const TensorOrder1 &t1, double d);

    friend bool operator==(TComplex d, const TensorOrder1 &t1);
    friend bool operator!=(TComplex d, const TensorOrder1 &t1);
    friend bool operator>=(TComplex d, const TensorOrder1 &t1);
    friend bool operator<=(TComplex d, const TensorOrder1 &t1);
    friend bool operator> (TComplex d, const TensorOrder1 &t1);
    friend bool operator< (TComplex d, const TensorOrder1 &t1);

    friend bool operator==(const TensorOrder1 &t1, TComplex d);
    friend bool operator!=(const TensorOrder1 &t1, TComplex d);
    friend bool operator>=(const TensorOrder1 &t1, TComplex d);
    friend bool operator<=(const TensorOrder1 &t1, TComplex d);
    friend bool operator> (const TensorOrder1 &t1, TComplex d);
    friend bool operator< (const TensorOrder1 &t1, TComplex d);

    ClassDef(TensorOrder1,0);
  };
}
#endif

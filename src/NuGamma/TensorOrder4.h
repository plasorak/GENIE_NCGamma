#ifndef _TENSORORDER4_H_
#define _TENSORORDER4_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <complex>

#include "NuGamma/Tensor.h"

namespace TensorUtils{

  class TensorOrder4: public Tensor{
    
  private:

    bool CheckIndex(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l) const;
    unsigned int GetGlobalIndex(const unsigned int i1, const unsigned int i2, const unsigned int i3, const unsigned int i4) const;

  public:
    TensorOrder4(const TensorOrder4& t1);
    TensorOrder4(const char * Name, unsigned int Dim1, unsigned int Dim2, unsigned int Dim3, unsigned int Dim4);
    TensorOrder4(const char * Name = "none", unsigned int Dim = 4);
    TComplex At(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l) const;
    void Set(unsigned int i, unsigned int j, unsigned int k, unsigned int l, TComplex c);
    TComplex& operator()(const unsigned int i1, const unsigned int i2, const unsigned int i3, const unsigned int i4);
    TComplex  operator()(const unsigned int i1, const unsigned int i2, const unsigned int i3, const unsigned int i4) const;
    void Print() const;
  
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

    friend TensorOrder4 Real(const TensorOrder4 &t1);
    friend TensorOrder4 Imaginary(const TensorOrder4 &t1);
    friend TensorOrder4 Abs(const TensorOrder4 &t1);
    friend TensorOrder4 operator-(const TensorOrder4 &t1);

    friend TensorOrder4 operator+(const int c, const TensorOrder4 &t1);
    friend TensorOrder4 operator-(const int c, const TensorOrder4 &t1);
    friend TensorOrder4 operator*(const int c, const TensorOrder4 &t1);
    friend TensorOrder4 operator/(const int c, const TensorOrder4 &t1);

    friend TensorOrder4 operator+(const TensorOrder4 &t1, const int c);
    friend TensorOrder4 operator-(const TensorOrder4 &t1, const int c);
    friend TensorOrder4 operator*(const TensorOrder4 &t1, const int c);
    friend TensorOrder4 operator/(const TensorOrder4 &t1, const int c);

    friend TensorOrder4 operator+(TComplex c, const TensorOrder4 &t1);
    friend TensorOrder4 operator-(TComplex c, const TensorOrder4 &t1);
    friend TensorOrder4 operator*(TComplex c, const TensorOrder4 &t1);
    friend TensorOrder4 operator/(TComplex c, const TensorOrder4 &t1);

    friend TensorOrder4 operator+(const TensorOrder4 &t1, TComplex c);
    friend TensorOrder4 operator-(const TensorOrder4 &t1, TComplex c);
    friend TensorOrder4 operator*(const TensorOrder4 &t1, TComplex c);
    friend TensorOrder4 operator/(const TensorOrder4 &t1, TComplex c);

    friend TensorOrder4 operator+(double d, const TensorOrder4 &t1);
    friend TensorOrder4 operator-(double d, const TensorOrder4 &t1);
    friend TensorOrder4 operator*(double d, const TensorOrder4 &t1);
    friend TensorOrder4 operator/(double d, const TensorOrder4 &t1);

    friend TensorOrder4 operator+(const TensorOrder4 &t1, double d);
    friend TensorOrder4 operator-(const TensorOrder4 &t1, double d);
    friend TensorOrder4 operator*(const TensorOrder4 &t1, double d);
    friend TensorOrder4 operator/(const TensorOrder4 &t1, double d);

    friend TensorOrder4 operator+(const TensorOrder4 &t1, const TensorOrder4 &t2);
    friend TensorOrder4 operator-(const TensorOrder4 &t1, const TensorOrder4 &t2);

    friend bool operator==(const TensorOrder4 &t1, const TensorOrder4 &t2);
    friend bool operator!=(const TensorOrder4 &t1, const TensorOrder4 &t2);
    friend bool operator>=(const TensorOrder4 &t1, const TensorOrder4 &t2);
    friend bool operator<=(const TensorOrder4 &t1, const TensorOrder4 &t2);
    friend bool operator> (const TensorOrder4 &t1, const TensorOrder4 &t2);
    friend bool operator< (const TensorOrder4 &t1, const TensorOrder4 &t2);

    friend bool operator==(double d, const TensorOrder4 &t1);
    friend bool operator!=(double d, const TensorOrder4 &t1);
    friend bool operator>=(double d, const TensorOrder4 &t1);
    friend bool operator<=(double d, const TensorOrder4 &t1);
    friend bool operator> (double d, const TensorOrder4 &t1);
    friend bool operator< (double d, const TensorOrder4 &t1);

    friend bool operator==(const TensorOrder4 &t1, double d);
    friend bool operator!=(const TensorOrder4 &t1, double d);
    friend bool operator>=(const TensorOrder4 &t1, double d);
    friend bool operator<=(const TensorOrder4 &t1, double d);
    friend bool operator> (const TensorOrder4 &t1, double d);
    friend bool operator< (const TensorOrder4 &t1, double d);

    friend bool operator==(TComplex d, const TensorOrder4 &t1);
    friend bool operator!=(TComplex d, const TensorOrder4 &t1);
    friend bool operator>=(TComplex d, const TensorOrder4 &t1);
    friend bool operator<=(TComplex d, const TensorOrder4 &t1);
    friend bool operator> (TComplex d, const TensorOrder4 &t1);
    friend bool operator< (TComplex d, const TensorOrder4 &t1);

    friend bool operator==(const TensorOrder4 &t1, TComplex d);
    friend bool operator!=(const TensorOrder4 &t1, TComplex d);
    friend bool operator>=(const TensorOrder4 &t1, TComplex d);
    friend bool operator<=(const TensorOrder4 &t1, TComplex d);
    friend bool operator> (const TensorOrder4 &t1, TComplex d);
    friend bool operator< (const TensorOrder4 &t1, TComplex d);

    ClassDef(TensorOrder4,0);
  };

} 


#endif 

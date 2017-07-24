#ifndef _TENSORORDER3_H_
#define _TENSORORDER3_H_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>
#include <complex>

#include "NuGamma/TensorInc.h"

namespace TensorUtils{
  

  class TensorOrder3: public Tensor{
    
  private:
    
    bool CheckIndex(const unsigned int i, const unsigned int j, const unsigned int k)const;
    unsigned int GetGlobalIndex(const unsigned int i1, const unsigned int i2, const unsigned int i3)const;
  
  public:

    TensorOrder3(const TensorOrder3& t1);
    TensorOrder3(const char * Name, unsigned int Dim1, unsigned int Dim2, unsigned int Dim3);
    TensorOrder3(const char * Name = "none", unsigned int Dim = 4);
    
    TComplex At(const unsigned int i, const unsigned int j, const unsigned int k)const;
    void Set(const unsigned int i, const unsigned int j, const unsigned int k, TComplex c);
    TComplex& operator()(const unsigned int i1, const unsigned int i2, const unsigned int i3);
    TComplex operator()(const unsigned int i1, const unsigned int i2, const unsigned int i3) const;
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

    friend TensorOrder3 Real(const TensorOrder3 &t1);
    friend TensorOrder3 Imaginary(const TensorOrder3 &t1);
    friend TensorOrder3 Abs(const TensorOrder3 &t1);
    friend TensorOrder3 operator-(const TensorOrder3 &t1);

    friend TensorOrder3 operator+(const int c, const TensorOrder3 &t1);
    friend TensorOrder3 operator-(const int c, const TensorOrder3 &t1);
    friend TensorOrder3 operator*(const int c, const TensorOrder3 &t1);
    friend TensorOrder3 operator/(const int c, const TensorOrder3 &t1);

    friend TensorOrder3 operator+(const TensorOrder3 &t1, const int c);
    friend TensorOrder3 operator-(const TensorOrder3 &t1, const int c);
    friend TensorOrder3 operator*(const TensorOrder3 &t1, const int c);
    friend TensorOrder3 operator/(const TensorOrder3 &t1, const int c);

    friend TensorOrder3 operator+(TComplex c, const TensorOrder3 &t1);
    friend TensorOrder3 operator-(TComplex c, const TensorOrder3 &t1);
    friend TensorOrder3 operator*(TComplex c, const TensorOrder3 &t1);
    friend TensorOrder3 operator/(TComplex c, const TensorOrder3 &t1);

    friend TensorOrder3 operator+(const TensorOrder3 &t1, TComplex c);
    friend TensorOrder3 operator-(const TensorOrder3 &t1, TComplex c);
    friend TensorOrder3 operator*(const TensorOrder3 &t1, TComplex c);
    friend TensorOrder3 operator/(const TensorOrder3 &t1, TComplex c);

    friend TensorOrder3 operator+(double d, const TensorOrder3 &t1);
    friend TensorOrder3 operator-(double d, const TensorOrder3 &t1);
    friend TensorOrder3 operator*(double d, const TensorOrder3 &t1);
    friend TensorOrder3 operator/(double d, const TensorOrder3 &t1);

    friend TensorOrder3 operator+(const TensorOrder3 &t1, double d);
    friend TensorOrder3 operator-(const TensorOrder3 &t1, double d);
    friend TensorOrder3 operator*(const TensorOrder3 &t1, double d);
    friend TensorOrder3 operator/(const TensorOrder3 &t1, double d);

    friend TensorOrder3 operator+(const TensorOrder3 &t1, const TensorOrder3 &t2);
    friend TensorOrder3 operator-(const TensorOrder3 &t1, const TensorOrder3 &t2);

    friend bool operator==(const TensorOrder3 &t1, const TensorOrder3 &t2);
    friend bool operator!=(const TensorOrder3 &t1, const TensorOrder3 &t2);
    friend bool operator>=(const TensorOrder3 &t1, const TensorOrder3 &t2);
    friend bool operator<=(const TensorOrder3 &t1, const TensorOrder3 &t2);
    friend bool operator> (const TensorOrder3 &t1, const TensorOrder3 &t2);
    friend bool operator< (const TensorOrder3 &t1, const TensorOrder3 &t2);

    friend bool operator==(double d, const TensorOrder3 &t1);
    friend bool operator!=(double d, const TensorOrder3 &t1);
    friend bool operator>=(double d, const TensorOrder3 &t1);
    friend bool operator<=(double d, const TensorOrder3 &t1);
    friend bool operator> (double d, const TensorOrder3 &t1);
    friend bool operator< (double d, const TensorOrder3 &t1);

    friend bool operator==(const TensorOrder3 &t1, double d);
    friend bool operator!=(const TensorOrder3 &t1, double d);
    friend bool operator>=(const TensorOrder3 &t1, double d);
    friend bool operator<=(const TensorOrder3 &t1, double d);
    friend bool operator> (const TensorOrder3 &t1, double d);
    friend bool operator< (const TensorOrder3 &t1, double d);

    friend bool operator==(TComplex d, const TensorOrder3 &t1);
    friend bool operator!=(TComplex d, const TensorOrder3 &t1);
    friend bool operator>=(TComplex d, const TensorOrder3 &t1);
    friend bool operator<=(TComplex d, const TensorOrder3 &t1);
    friend bool operator> (TComplex d, const TensorOrder3 &t1);
    friend bool operator< (TComplex d, const TensorOrder3 &t1);

    friend bool operator==(const TensorOrder3 &t1, TComplex d);
    friend bool operator!=(const TensorOrder3 &t1, TComplex d);
    friend bool operator>=(const TensorOrder3 &t1, TComplex d);
    friend bool operator<=(const TensorOrder3 &t1, TComplex d);
    friend bool operator> (const TensorOrder3 &t1, TComplex d);
    friend bool operator< (const TensorOrder3 &t1, TComplex d);
 
    ClassDef(TensorOrder3,0);
  };
}
#endif 

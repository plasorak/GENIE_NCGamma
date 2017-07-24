#include "NuGamma/TensorInc.h"
#include "TString.h"

namespace TensorUtils{

  bool TensorOrder2::CheckIndex(const unsigned int i, const unsigned int j) const{
    return (i < Dim[0] && j < Dim[1]);
  };

  unsigned int TensorOrder2::GetGlobalIndex(const unsigned int i,
					  const unsigned int j) const{
    return i * Dim[0] + j;
  };

  TensorOrder2::TensorOrder2(const char * Name, unsigned int SizeDim1, unsigned int SizeDim2){
    this->Name = Name;
    this->Dim.reserve(2);
    this->Dim.push_back(SizeDim1);
    this->Dim.push_back(SizeDim2);
    this->Order = Dim.size();
    this->Element.reserve(Dim[0]*Dim[1]);

    TComplex Zero(0.,0.);

    for(unsigned int i1 = 0; i1 < this->Dim[0]; i1++)
      for(unsigned int i2 = 0; i2 < this->Dim[1]; i2++)
	this->Element.push_back(Zero);
      
    
  };
  TensorOrder2::TensorOrder2(const TensorOrder2& t1){
    t1.AssertGoodSize(2);
    this->Name = t1.Name;
    this->Order = t1.Order;
    this->Dim.reserve(2);
    this->Dim.push_back(t1.Dim[0]);
    this->Dim.push_back(t1.Dim[1]);
    this->Element.reserve(Dim[0]*Dim[1]);

    for(unsigned int i = 0; i < t1.Element.size(); i++)
      this->Element.push_back(t1.Element[i]);
  };

  TensorOrder2::TensorOrder2(const char * Name, unsigned int SizeDim){
    this->Name = Name;
    this->Dim.reserve(2);
    this->Dim.push_back(SizeDim);
    this->Dim.push_back(SizeDim);
    this->Order = Dim.size();
    this->Element.reserve(Dim[0]*Dim[1]);

    TComplex Zero(0.,0.);

    for(unsigned int i1 = 0; i1 < this->Dim[0]; i1++)
      for(unsigned int i2 = 0; i2 < this->Dim[1]; i2++)
	this->Element.push_back(Zero);
    
    
  };

  TComplex TensorOrder2::At(const unsigned int i, const unsigned int j) const{
    if(!CheckIndex(i, j))
      exit(1);
    return Element[GetGlobalIndex(i, j)];

  };

  void TensorOrder2::Set(const unsigned int i, const unsigned int j, TComplex c){
    if(!CheckIndex(i, j))
      exit(1);
    
    Element[GetGlobalIndex(i, j)] = c;
      
  };

  void TensorOrder2::Print() const{
    std::cout << "______________________________________________________________________" << std::endl;
    std::cout << "TensorOrder" << this->Dim.size() << " with name '" << this->Name << "'" << std::endl;
    for(unsigned int i = 0; i < Order; i++){
      std::cout << "Dimension(" << i << ") size: " << this->Dim[i] << std::endl;
    }
    std::cout << "Element(<position0>, <position1> = <globalposition>)" << std::endl;
    for(unsigned int i = 0; i < this->Dim[0]; i++){
      for(unsigned int j = 0; j < this->Dim[1]; j++){
	std::cout << " Element(" << i << ", " << j << " = " << this->GetGlobalIndex(i, j) << ") = " << this->At(i, j) << std::endl;
      }
    }
    std::cout << "______________________________________________________________________" << std::endl;

  };
  
  TComplex& TensorOrder2::operator()(const unsigned int i1, const unsigned int i2){
    if(!CheckIndex(i1, i2))
      exit(1);
    return Element[GetGlobalIndex(i1, i2)];
  };

  TComplex  TensorOrder2::operator()(const unsigned int i1, const unsigned int i2) const{
    if(!CheckIndex(i1, i2))
      exit(1);
    return Element[GetGlobalIndex(i1, i2)];
  };


  TensorOrder2 Real(const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Real();
    return t;
  };

  TensorOrder2 Imaginary(const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Imaginary();
    return t;
  };

  TensorOrder2 Abs(const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Abs();
    return t;
  };
  
  TensorOrder2 operator-(const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Minus();
    return t;
  };
  
  TensorOrder2 operator+(const int c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder2 operator-(const int c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder2 operator*(const int c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder2 operator/(const int c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder2 operator+(const TensorOrder2 &t1, const int c){
    TensorOrder2 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder2 operator-(const TensorOrder2 &t1, const int c){
    TensorOrder2 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder2 operator*(const TensorOrder2 &t1, const int c){
    TensorOrder2 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder2 operator/(const TensorOrder2 &t1, const int c){
    TensorOrder2 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder2 operator+(TComplex c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder2 operator-(TComplex c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder2 operator*(TComplex c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder2 operator/(TComplex c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder2 operator+(const TensorOrder2 &t1, TComplex c){
    TensorOrder2 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder2 operator-(const TensorOrder2 &t1, TComplex c){
    TensorOrder2 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder2 operator*(const TensorOrder2 &t1, TComplex c){
    TensorOrder2 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder2 operator/(const TensorOrder2 &t1, TComplex c){
    TensorOrder2 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder2 operator+(double c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder2 operator-(double c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder2 operator*(double c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder2 operator/(double c, const TensorOrder2 &t1){
    TensorOrder2 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder2 operator+(const TensorOrder2 &t1, double c){
    TensorOrder2 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder2 operator-(const TensorOrder2 &t1, double c){
    TensorOrder2 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder2 operator*(const TensorOrder2 &t1, double c){
    TensorOrder2 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder2 operator/(const TensorOrder2 &t1, double c){
    TensorOrder2 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder2 operator+(const TensorOrder2 &t1, const TensorOrder2 &t2){
    TensorOrder2 t(t1);
    t.Add(t2);
    return t;
  };

  TensorOrder2 operator-(const TensorOrder2 &t1, const TensorOrder2 &t2){
    TensorOrder2 t(t1);
    t.Subtract(t2);
    return t;
  };

  TensorOrder2 operator*(const TensorOrder2 &t1, const TensorOrder2 &t2){

    // First check the dimension of the tensors
    t1.AssertGoodSize(2);
    t2.AssertGoodSize(2);
    // First check the dimension of the matrix
    // (i * j) (j * k) = (i * k)
    t1.AssertGoodSize(1, t2.GetOrderDim(1));

    TensorOrder2 result(Form("(%sX%s)",t1.GetName(), t2.GetName()), t1.GetOrderDim(0), t2.GetOrderDim(1));

    for(unsigned int i = 0; i < t1.GetOrderDim(0); i++)
      for(unsigned int k = 0; k < t2.GetOrderDim(1); k++)
	for(unsigned int j = 0; j < t1.GetOrderDim(1); j++)
	  result(i,k) += t1(i,j)*t2(j,k);
    
    return result;
  };

  bool operator==(const TensorOrder2 &t1, const TensorOrder2 &t2){ return t1.IsEqual(t2);          };
  bool operator!=(const TensorOrder2 &t1, const TensorOrder2 &t2){ return t1.IsDifferent(t2);      };
  bool operator>=(const TensorOrder2 &t1, const TensorOrder2 &t2){ return t1.IsBiggerOrEqual(t2);  };
  bool operator<=(const TensorOrder2 &t1, const TensorOrder2 &t2){ return t1.IsSmallerOrEqual(t2); };
  bool operator> (const TensorOrder2 &t1, const TensorOrder2 &t2){ return t1.IsBigger(t2);         };
  bool operator< (const TensorOrder2 &t1, const TensorOrder2 &t2){ return t1.IsSmaller(t2);        };

  bool operator==(double d, const TensorOrder2 &t1){ return t1.IsEqual(d);          };
  bool operator!=(double d, const TensorOrder2 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(double d, const TensorOrder2 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(double d, const TensorOrder2 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (double d, const TensorOrder2 &t1){ return t1.IsSmaller(d);        };
  bool operator< (double d, const TensorOrder2 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder2 &t1, double d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder2 &t1, double d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder2 &t1, double d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder2 &t1, double d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder2 &t1, double d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder2 &t1, double d){ return t1.IsSmaller(d);        };

  bool operator==(TComplex d, const TensorOrder2 &t1){ return t1.IsEqual(d);          };
  bool operator!=(TComplex d, const TensorOrder2 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(TComplex d, const TensorOrder2 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(TComplex d, const TensorOrder2 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (TComplex d, const TensorOrder2 &t1){ return t1.IsSmaller(d);        };
  bool operator< (TComplex d, const TensorOrder2 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder2 &t1, TComplex d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder2 &t1, TComplex d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder2 &t1, TComplex d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder2 &t1, TComplex d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder2 &t1, TComplex d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder2 &t1, TComplex d){ return t1.IsSmaller(d);        };


}

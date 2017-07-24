#include "NuGamma/TensorInc.h"
#include "TString.h"

namespace TensorUtils{

  bool TensorOrder1::CheckIndex(const unsigned int i) const{
    return (i < Dim[0]);
  };

  unsigned int TensorOrder1::GetGlobalIndex(const unsigned int i) const{
    return i;
  };
  
  TensorOrder1::TensorOrder1(const char *Name, unsigned int SizeDim){
    this->Name = Name;
    Dim.clear();
    Dim.reserve(1);
    Dim.push_back(SizeDim);
    Order = Dim.size();
    Element.clear();
    Element.reserve(Dim[0]);

    TComplex Zero(0.,0.);
    for(unsigned int i = 0; i < Dim[0]; i++){
      Element.push_back(Zero);
    }       
  };


  TensorOrder1::TensorOrder1(const TensorOrder1& t1){
    t1.AssertGoodSize(1);
    this->Name = t1.Name;
    this->Order = t1.Order;
    for(unsigned int i = 0; i < t1.Element.size(); i++)
      this->Element.push_back(t1.Element[i]);
    for(unsigned int i = 0; i < t1.Dim.size(); i++)
      this->Dim.push_back(t1.Dim[i]);
  };

  TComplex TensorOrder1::At(const unsigned int i) const{
    if(!CheckIndex(i))
      exit(1);
    return Element[GetGlobalIndex(i)];

  };

  void TensorOrder1::Set(const unsigned int i, TComplex c){
    if(!CheckIndex(i))
      exit(1);
    Element[GetGlobalIndex(i)] = c;

  };

  void TensorOrder1::Print() const{
    std::cout << "______________________________________________________________________" << std::endl;
    std::cout << "TensorOrder"<< this->Dim.size() << " with name '" << Name << "'" << std::endl;
    std::cout << "Dimension(" << 0 << ") size: " << Dim[0] << std::endl;
    std::cout << "Element(<position> = <globalposition>)" << std::endl;
    for(unsigned int j = 0; j < Dim[0]; j++){
      std::cout << " Element(" << j << " = " << this->GetGlobalIndex(j) << ") = " << this->At(j) << std::endl;
    }
    std::cout << "______________________________________________________________________" << std::endl;
  };

  TComplex& TensorOrder1::operator()(const unsigned int i){
    if(!CheckIndex(i))
      exit(1);
    return Element[GetGlobalIndex(i)];
  };

  TComplex TensorOrder1::operator()(const unsigned int i) const{
    if(!CheckIndex(i))
      exit(1);
    return Element[GetGlobalIndex(i)];
  };

  TensorOrder1 Real(const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Real();
    return t;
  };

  TensorOrder1 Imaginary(const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Imaginary();
    return t;
  };

  TensorOrder1 Abs(const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Abs();
    return t;
  };
  
  TensorOrder1 operator-(const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Minus();
    return t;
  };
  
  
  TensorOrder1 operator+(const int c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder1 operator-(const int c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder1 operator*(const int c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder1 operator/(const int c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder1 operator+(const TensorOrder1 &t1, const int c){
    TensorOrder1 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder1 operator-(const TensorOrder1 &t1, const int c){
    TensorOrder1 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder1 operator*(const TensorOrder1 &t1, const int c){
    TensorOrder1 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder1 operator/(const TensorOrder1 &t1, const int c){
    TensorOrder1 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder1 operator+(TComplex c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder1 operator-(TComplex c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder1 operator*(TComplex c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder1 operator/(TComplex c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder1 operator+(const TensorOrder1 &t1, TComplex c){
    TensorOrder1 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder1 operator-(const TensorOrder1 &t1, TComplex c){
    TensorOrder1 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder1 operator*(const TensorOrder1 &t1, TComplex c){
    TensorOrder1 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder1 operator/(const TensorOrder1 &t1, TComplex c){
    TensorOrder1 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder1 operator+(double c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder1 operator-(double c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder1 operator*(double c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder1 operator/(double c, const TensorOrder1 &t1){
    TensorOrder1 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder1 operator+(const TensorOrder1 &t1, double c){
    TensorOrder1 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder1 operator-(const TensorOrder1 &t1, double c){
    TensorOrder1 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder1 operator*(const TensorOrder1 &t1, double c){
    TensorOrder1 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder1 operator/(const TensorOrder1 &t1, double c){
    TensorOrder1 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder1 operator+(const TensorOrder1 &t1, const TensorOrder1 &t2){
    TensorOrder1 t(t1);
    t.Add(t2);
    return t;
  };

  TensorOrder1 operator-(const TensorOrder1 &t1, const TensorOrder1 &t2){
    TensorOrder1 t(t1);
    t.Subtract(t2);
    return t;
  };

  TensorOrder1 operator*(const TensorOrder2 &t1, const TensorOrder1 &t2){

    // First check the dimension of the tensors
    t1.AssertGoodSize(2);
    t2.AssertGoodSize(1);

    // Then check that size are good for multiplication
    // (i * j) (j) = (i)
    t1.AssertGoodSize(1,t2.GetOrderDim(0));

    TensorOrder1 result(Form("(%sX%s)",t1.GetName(), t2.GetName()),t1.GetOrderDim(0));
    
    for(unsigned int i = 0; i < t1.GetOrderDim(0); i++)
      for(unsigned int j = 0; j < t2.GetOrderDim(0); j++){
	result(i) += t1(i,j)*t2(j);
      }
    return result;
  };

  bool operator==(const TensorOrder1 &t1, const TensorOrder1 &t2){ return t1.IsEqual(t2);          };
  bool operator!=(const TensorOrder1 &t1, const TensorOrder1 &t2){ return t1.IsDifferent(t2);      };
  bool operator>=(const TensorOrder1 &t1, const TensorOrder1 &t2){ return t1.IsBiggerOrEqual(t2);  };
  bool operator<=(const TensorOrder1 &t1, const TensorOrder1 &t2){ return t1.IsSmallerOrEqual(t2); };
  bool operator> (const TensorOrder1 &t1, const TensorOrder1 &t2){ return t1.IsBigger(t2);         };
  bool operator< (const TensorOrder1 &t1, const TensorOrder1 &t2){ return t1.IsSmaller(t2);        };

  bool operator==(double d, const TensorOrder1 &t1){ return t1.IsEqual(d);          };
  bool operator!=(double d, const TensorOrder1 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(double d, const TensorOrder1 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(double d, const TensorOrder1 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (double d, const TensorOrder1 &t1){ return t1.IsSmaller(d);        };
  bool operator< (double d, const TensorOrder1 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder1 &t1, double d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder1 &t1, double d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder1 &t1, double d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder1 &t1, double d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder1 &t1, double d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder1 &t1, double d){ return t1.IsSmaller(d);        };

  bool operator==(TComplex d, const TensorOrder1 &t1){ return t1.IsEqual(d);          };
  bool operator!=(TComplex d, const TensorOrder1 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(TComplex d, const TensorOrder1 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(TComplex d, const TensorOrder1 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (TComplex d, const TensorOrder1 &t1){ return t1.IsSmaller(d);        };
  bool operator< (TComplex d, const TensorOrder1 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder1 &t1, TComplex d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder1 &t1, TComplex d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder1 &t1, TComplex d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder1 &t1, TComplex d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder1 &t1, TComplex d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder1 &t1, TComplex d){ return t1.IsSmaller(d);        };



}

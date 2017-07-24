#include "NuGamma/TensorInc.h"

namespace TensorUtils{

  bool TensorOrder3::CheckIndex(const unsigned int i, const unsigned int j, const unsigned int k) const{
    return (i < Dim[0] && j < Dim[1] && k < Dim[2]);
  };

  unsigned int TensorOrder3::GetGlobalIndex(const unsigned int i1,
				 const unsigned int i2,
				 const unsigned int i3) const{
    return Dim[0] * Dim[1] * i1 + Dim[1] * i2 + i3;
  };

  TensorOrder3::TensorOrder3(const TensorOrder3& t1){
    t1.AssertGoodSize(3);
    this->Name = t1.Name;
    this->Order = t1.Order;
    for(unsigned int i = 0; i < t1.Element.size(); i++)
      this->Element.push_back(t1.Element[i]);
    for(unsigned int i = 0; i < t1.Dim.size(); i++)
      this->Dim.push_back(t1.Dim[i]);
  };

  TensorOrder3::TensorOrder3(const char * Name, unsigned int SizeDim1, unsigned int SizeDim2, unsigned int SizeDim3){
    this->Name = Name;
    Dim.clear();
    Dim.reserve(3);
    Dim.push_back(SizeDim1);
    Dim.push_back(SizeDim2);
    Dim.push_back(SizeDim3);
    Order = Dim.size();
    Element.clear();
    Element.reserve(Dim[0]*Dim[1]*Dim[2]);

    TComplex Zero(0.,0.);
    for(unsigned int i1 = 0; i1 < Dim[0]; i1++){
      for(unsigned int i2 = 0; i2 < Dim[1]; i2++){
	for(unsigned int i3 = 0; i3 < Dim[2]; i3++){
	  Element.push_back(Zero);
	}
      }
    }
  };

  TensorOrder3::TensorOrder3(const char * Name, unsigned int SizeDim){
    this->Name = Name;
    Dim.clear();
    Dim.reserve(3);
    Dim.push_back(SizeDim);
    Dim.push_back(SizeDim);
    Dim.push_back(SizeDim);
    Order = Dim.size();
    Element.clear();
    Element.reserve(Dim[0]*Dim[1]*Dim[2]);

    TComplex Zero(0.,0.);
    for(unsigned int i1 = 0; i1 < Dim[0]; i1++){
      for(unsigned int i2 = 0; i2 < Dim[1]; i2++){
	for(unsigned int i3 = 0; i3 < Dim[2]; i3++){
	  Element.push_back(Zero);
	}
      }
    }
  };

  TComplex TensorOrder3::At(const unsigned int i, const unsigned int j, const unsigned int k)const{
    if(!CheckIndex(i, j, k))
      exit(1);

    return Element[GetGlobalIndex(i, j, k)];

  };

  void TensorOrder3::Set(const unsigned int i, const unsigned int j, const unsigned int k, TComplex c){
    if(!CheckIndex(i, j, k))
      exit(1);
    
    Element[GetGlobalIndex(i, j, k)] = c;
    
  };

  void TensorOrder3::Print() const{
    std::cout << "______________________________________________________________________" << std::endl;
    std::cout << "Tensor Order" << this->Dim.size() << " with name '" << this->Name << "'" << std::endl;
    for(unsigned int i = 0; i < Order; i++){
      std::cout << "Dimension(" << i << ") size: " << Dim[i] << std::endl;
    }
    std::cout << "Element(<position0>, <position1>, <position2> = <globalposition>)" << std::endl;
    for(unsigned int i = 0; i < Dim[0]; i++){
      for(unsigned int j = 0; j < Dim[1]; j++){
	for(unsigned int k = 0; k < Dim[2]; k++){
	  std::cout << "Element(" << i << ", " << j << ", " << k << " = " << this->GetGlobalIndex(i,j,k) << ") = " << this->At(i, j, k) << std::endl;
	}
      }
    }
    std::cout << "______________________________________________________________________" << std::endl;

  };

  TComplex& TensorOrder3::operator()(const unsigned int i1,
				   const unsigned int i2,
				   const unsigned int i3){
    if(!CheckIndex(i1, i2, i3))
      exit(1);
    return Element[GetGlobalIndex(i1, i2, i3)];
  };

  TComplex TensorOrder3::operator()(const unsigned int i1,
				  const unsigned int i2,
				  const unsigned int i3) const{
    if(!CheckIndex(i1, i2, i3))
      exit(1);
    return Element[GetGlobalIndex(i1, i2, i3)];
  };

  TensorOrder3 Real(const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Real();
    return t;
  };

  TensorOrder3 Imaginary(const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Imaginary();
    return t;
  };

  TensorOrder3 Abs(const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Abs();
    return t;
  };
  TensorOrder3 operator-(const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Minus();
    return t;
  };
  
  TensorOrder3 operator+(const int c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder3 operator-(const int c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder3 operator*(const int c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder3 operator/(const int c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder3 operator+(const TensorOrder3 &t1, const int c){
    TensorOrder3 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder3 operator-(const TensorOrder3 &t1, const int c){
    TensorOrder3 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder3 operator*(const TensorOrder3 &t1, const int c){
    TensorOrder3 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder3 operator/(const TensorOrder3 &t1, const int c){
    TensorOrder3 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder3 operator+(TComplex c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder3 operator-(TComplex c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder3 operator*(TComplex c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder3 operator/(TComplex c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder3 operator+(const TensorOrder3 &t1, TComplex c){
    TensorOrder3 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder3 operator-(const TensorOrder3 &t1, TComplex c){
    TensorOrder3 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder3 operator*(const TensorOrder3 &t1, TComplex c){
    TensorOrder3 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder3 operator/(const TensorOrder3 &t1, TComplex c){
    TensorOrder3 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder3 operator+(double c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder3 operator-(double c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder3 operator*(double c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder3 operator/(double c, const TensorOrder3 &t1){
    TensorOrder3 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder3 operator+(const TensorOrder3 &t1, double c){
    TensorOrder3 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder3 operator-(const TensorOrder3 &t1, double c){
    TensorOrder3 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder3 operator*(const TensorOrder3 &t1, double c){
    TensorOrder3 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder3 operator/(const TensorOrder3 &t1, double c){
    TensorOrder3 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder3 operator+(const TensorOrder3 &t1, const TensorOrder3 &t2){
    TensorOrder3 t(t1);
    t.Add(t2);
    return t;
  };

  TensorOrder3 operator-(const TensorOrder3 &t1, const TensorOrder3 &t2){
    TensorOrder3 t(t1);
    t.Subtract(t2);
    return t;
  };

  bool operator==(const TensorOrder3 &t1, const TensorOrder3 &t2){ return t1.IsEqual(t2);          };
  bool operator!=(const TensorOrder3 &t1, const TensorOrder3 &t2){ return t1.IsDifferent(t2);      };
  bool operator>=(const TensorOrder3 &t1, const TensorOrder3 &t2){ return t1.IsBiggerOrEqual(t2);  };
  bool operator<=(const TensorOrder3 &t1, const TensorOrder3 &t2){ return t1.IsSmallerOrEqual(t2); };
  bool operator> (const TensorOrder3 &t1, const TensorOrder3 &t2){ return t1.IsBigger(t2);         };
  bool operator< (const TensorOrder3 &t1, const TensorOrder3 &t2){ return t1.IsSmaller(t2);        };

  bool operator==(double d, const TensorOrder3 &t1){ return t1.IsEqual(d);          };
  bool operator!=(double d, const TensorOrder3 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(double d, const TensorOrder3 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(double d, const TensorOrder3 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (double d, const TensorOrder3 &t1){ return t1.IsSmaller(d);        };
  bool operator< (double d, const TensorOrder3 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder3 &t1, double d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder3 &t1, double d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder3 &t1, double d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder3 &t1, double d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder3 &t1, double d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder3 &t1, double d){ return t1.IsSmaller(d);        };

  bool operator==(TComplex d, const TensorOrder3 &t1){ return t1.IsEqual(d);          };
  bool operator!=(TComplex d, const TensorOrder3 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(TComplex d, const TensorOrder3 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(TComplex d, const TensorOrder3 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (TComplex d, const TensorOrder3 &t1){ return t1.IsSmaller(d);        };
  bool operator< (TComplex d, const TensorOrder3 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder3 &t1, TComplex d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder3 &t1, TComplex d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder3 &t1, TComplex d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder3 &t1, TComplex d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder3 &t1, TComplex d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder3 &t1, TComplex d){ return t1.IsSmaller(d);        };


}

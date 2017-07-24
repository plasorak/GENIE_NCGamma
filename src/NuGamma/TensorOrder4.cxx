#include "NuGamma/TensorInc.h"

namespace TensorUtils{

  bool TensorOrder4::CheckIndex(const unsigned int i, const unsigned int j, const unsigned int k, const unsigned int l) const{
    if(i < Dim[0] && j < Dim[1] && k < Dim[2] && l < Dim[3])
      return true;
    return false;
  };

  unsigned int TensorOrder4::GetGlobalIndex(const unsigned int i1, const unsigned int i2, const unsigned int i3, const unsigned int i4) const{
    return (Dim[0] * Dim[1] * Dim[2] * i1 +
	    Dim[1] * Dim[2] * i2 +
	    Dim[2] * i3 +
	    i4);
  };

  TensorOrder4::TensorOrder4(const TensorOrder4& t1){
    t1.AssertGoodSize(4);
    this->Name = t1.Name;
    this->Order = t1.Order;
    for(unsigned int i = 0; i < t1.Element.size(); i++)
      this->Element.push_back(t1.Element[i]);
    for(unsigned int i = 0; i < t1.Dim.size(); i++)
      this->Dim.push_back(t1.Dim[i]);
  };


  TensorOrder4::TensorOrder4(const char * Name, unsigned int SizeDim1, unsigned int SizeDim2, unsigned int SizeDim3, unsigned int SizeDim4){
    this->Name = Name;
    Dim.clear();
    Dim.reserve(4);
    Dim.push_back(SizeDim1);
    Dim.push_back(SizeDim2);
    Dim.push_back(SizeDim3);
    Dim.push_back(SizeDim4);
    Order = Dim.size();
    Element.clear();
    Element.reserve(Dim[0]*Dim[1]*Dim[2]*Dim[3]);
    TComplex Zero(0.,0.);
    for(unsigned int i1 = 0; i1 < SizeDim1; i1++){
      for(unsigned int i2 = 0; i2 < SizeDim2; i2++){
	for(unsigned int i3 = 0; i3 < SizeDim3; i3++){
	  for(unsigned int i4 = 0; i4 < SizeDim4; i4++)
	    Element.push_back(Zero);
	}
      }
    }
  };
  


  TensorOrder4::TensorOrder4(const char * Name, unsigned int SizeDim){
    this->Name = Name;
    Dim.clear();
    Dim.reserve(4);
    Dim.push_back(SizeDim);
    Dim.push_back(SizeDim);
    Dim.push_back(SizeDim);
    Dim.push_back(SizeDim);
    Order = Dim.size();
    Element.clear();
    Element.reserve(Dim[0]*Dim[1]*Dim[2]*Dim[3]);
    TComplex Zero(0.,0.);
    for(unsigned int i1 = 0; i1 < SizeDim; i1++){
      for(unsigned int i2 = 0; i2 < SizeDim; i2++){
	for(unsigned int i3 = 0; i3 < SizeDim; i3++){
	  for(unsigned int i4 = 0; i4 < SizeDim; i4++)
	    Element.push_back(Zero);
	}
      }
    }
  };

  TComplex TensorOrder4::At(const unsigned int i,
			  const unsigned int j,
			  const unsigned int k,
			  const unsigned int l) const{
    if(!CheckIndex(i, j, k, l))
      exit(1);

    return Element[GetGlobalIndex(i, j, k, l)];
  };

  void TensorOrder4::Set(unsigned int i, unsigned int j, unsigned int k, unsigned int l, TComplex c){
    if(!CheckIndex(i, j, k, l))
      exit(1);
    
    Element[GetGlobalIndex(i, j, k, l)] = c;
    
  };

  void TensorOrder4::Print() const{
    std::cout << "______________________________________________________________________" << std::endl;
    std::cout << "Tensor Order" << this->Dim.size() << " with name '" << Name << "'" << std::endl;
    for(unsigned int i = 0; i < Order; i++){
      std::cout << "Dimension(" << i << ") size: " << Dim[i] << std::endl;
    }

    std::cout << "Element(<position0>, <position1>, <position2> <position3> = <globalposition>)" << std::endl;
    for(unsigned int i = 0; i < Dim[0]; i++){
      for(unsigned int j = 0; j < Dim[1]; j++){
	for(unsigned int k = 0; k < Dim[2]; k++){
	  for(unsigned int l = 0; l < Dim[3]; l++){
	    std::cout << "Element(" << i << ", " << j << " = " << this->GetGlobalIndex(i,j,k,l) << ") = " << this->At(i,j,k,l) << std::endl;
	  }
	}
      }
    }
    
    std::cout << "______________________________________________________________________" << std::endl;
    
  };


  TComplex& TensorOrder4::operator()(const unsigned int i1,
				   const unsigned int i2,
				   const unsigned int i3,
				   const unsigned int i4){
    if(!CheckIndex(i1, i2, i3, i4))
      exit(1);

    return Element[GetGlobalIndex(i1, i2, i3, i4)];
  };

  TComplex TensorOrder4::operator()(const unsigned int i1,
				  const unsigned int i2,
				  const unsigned int i3,
				  const unsigned int i4) const{
    if(!CheckIndex(i1, i2, i3, i4))
      exit(1);

    return Element[GetGlobalIndex(i1, i2, i3, i4)];
  };


  TensorOrder4 Real(const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Real();
    return t;
  };

  TensorOrder4 Imaginary(const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Imaginary();
    return t;
  };

  TensorOrder4 Abs(const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Abs();
    return t;
  };
  
  TensorOrder4 operator-(const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Minus();
    return t;
  };
  
  TensorOrder4 operator+(const int c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder4 operator-(const int c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder4 operator*(const int c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder4 operator/(const int c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder4 operator+(const TensorOrder4 &t1, const int c){
    TensorOrder4 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder4 operator-(const TensorOrder4 &t1, const int c){
    TensorOrder4 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder4 operator*(const TensorOrder4 &t1, const int c){
    TensorOrder4 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder4 operator/(const TensorOrder4 &t1, const int c){
    TensorOrder4 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder4 operator+(TComplex c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder4 operator-(TComplex c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder4 operator*(TComplex c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder4 operator/(TComplex c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder4 operator+(const TensorOrder4 &t1, TComplex c){
    TensorOrder4 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder4 operator-(const TensorOrder4 &t1, TComplex c){
    TensorOrder4 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder4 operator*(const TensorOrder4 &t1, TComplex c){
    TensorOrder4 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder4 operator/(const TensorOrder4 &t1, TComplex c){
    TensorOrder4 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder4 operator+(double c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder4 operator-(double c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Minus();
    t.Add(c);
    return t;
  };
  TensorOrder4 operator*(double c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder4 operator/(double c, const TensorOrder4 &t1){
    TensorOrder4 t(t1);
    t.OneOverElementWise();
    t.Multiply(c);
    return t;
  };

  TensorOrder4 operator+(const TensorOrder4 &t1, double c){
    TensorOrder4 t(t1);
    t.Add(c);
    return t;
  };
  TensorOrder4 operator-(const TensorOrder4 &t1, double c){
    TensorOrder4 t(t1);
    t.Subtract(c);
    return t;
  };
  TensorOrder4 operator*(const TensorOrder4 &t1, double c){
    TensorOrder4 t(t1);
    t.Multiply(c);
    return t;
  };
  TensorOrder4 operator/(const TensorOrder4 &t1, double c){
    TensorOrder4 t(t1);
    t.Divide(c);
    return t;
  };

  TensorOrder4 operator+(const TensorOrder4 &t1, const TensorOrder4 &t2){
    TensorOrder4 t(t1);
    t.Add(t2);
    return t;
  };

  TensorOrder4 operator-(const TensorOrder4 &t1, const TensorOrder4 &t2){
    TensorOrder4 t(t1);
    t.Subtract(t2);
    return t;
  };

  bool operator==(const TensorOrder4 &t1, const TensorOrder4 &t2){ return t1.IsEqual(t2);          };
  bool operator!=(const TensorOrder4 &t1, const TensorOrder4 &t2){ return t1.IsDifferent(t2);      };
  bool operator>=(const TensorOrder4 &t1, const TensorOrder4 &t2){ return t1.IsBiggerOrEqual(t2);  };
  bool operator<=(const TensorOrder4 &t1, const TensorOrder4 &t2){ return t1.IsSmallerOrEqual(t2); };
  bool operator> (const TensorOrder4 &t1, const TensorOrder4 &t2){ return t1.IsBigger(t2);         };
  bool operator< (const TensorOrder4 &t1, const TensorOrder4 &t2){ return t1.IsSmaller(t2);        };

  bool operator==(double d, const TensorOrder4 &t1){ return t1.IsEqual(d);          };
  bool operator!=(double d, const TensorOrder4 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(double d, const TensorOrder4 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(double d, const TensorOrder4 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (double d, const TensorOrder4 &t1){ return t1.IsSmaller(d);        };
  bool operator< (double d, const TensorOrder4 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder4 &t1, double d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder4 &t1, double d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder4 &t1, double d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder4 &t1, double d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder4 &t1, double d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder4 &t1, double d){ return t1.IsSmaller(d);        };

  bool operator==(TComplex d, const TensorOrder4 &t1){ return t1.IsEqual(d);          };
  bool operator!=(TComplex d, const TensorOrder4 &t1){ return t1.IsDifferent(d);      };
  bool operator>=(TComplex d, const TensorOrder4 &t1){ return t1.IsSmallerOrEqual(d); };
  bool operator<=(TComplex d, const TensorOrder4 &t1){ return t1.IsBiggerOrEqual(d);  };
  bool operator> (TComplex d, const TensorOrder4 &t1){ return t1.IsSmaller(d);        };
  bool operator< (TComplex d, const TensorOrder4 &t1){ return t1.IsBigger(d);         };

  bool operator==(const TensorOrder4 &t1, TComplex d){ return t1.IsEqual(d);          };
  bool operator!=(const TensorOrder4 &t1, TComplex d){ return t1.IsDifferent(d);      };
  bool operator>=(const TensorOrder4 &t1, TComplex d){ return t1.IsBiggerOrEqual(d);  };
  bool operator<=(const TensorOrder4 &t1, TComplex d){ return t1.IsSmallerOrEqual(d); };
  bool operator> (const TensorOrder4 &t1, TComplex d){ return t1.IsBigger(d);         };
  bool operator< (const TensorOrder4 &t1, TComplex d){ return t1.IsSmaller(d);        };

  
}

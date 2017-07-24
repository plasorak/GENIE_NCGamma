#include "NuGamma/TensorInc.h"

#include "TComplex.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdlib>

namespace TensorUtils{

  // Deletion function for mem leaks
  Tensor::~Tensor(){
    Element.clear();
    Dim.clear();
  };


  // Get the order of the tensor, technically this is just the class,
  // but it is always interesting to have it for tensor-wide function
  unsigned int Tensor::GetOrder() const{
    // Check that everything is behaving as expected
    this->AssertGoodSize(this->Order);
    return Order;
  };

  // Get the dimension of every order (which is could be different for every order...)
  unsigned int Tensor::GetOrderDim(const unsigned int i) const{
    // Check that everything is behaving as expected
    this->AssertGoodSize(this->Order);

    if(i >= 0 && i < Order)
      return Dim[i];
    else{ // insult the user if there is a problem and leave
      std::cout << "Trying to access the size of dimension which doesn't exist in the tensor."<<  std::endl;
      std::cout << "Index = " << i << std::endl;
      std::cout << "Dim.size() = " << Dim.size() << std::endl;
      std::cout << "Dim[" << i << "] = " << Dim[i] << std::endl;
      std::cout << "Program will exit" << std::endl;
      exit(1);
    }
  };

  // usefull for dot product, sum over all the element of the tensor
  TComplex Tensor::SumOver() const{
    TComplex Zero(0,0);
    for(unsigned int i = 0; i < Element.size(); i++)
      Zero = Zero + Element[i];
    return Zero;
  }
  
  // "Assign" operator, one of the tensor-wide function
  Tensor& Tensor::operator=(const Tensor &t1){
    this->Name = t1.Name;
    this->Order  = t1.Order;
    for(unsigned int i = 0; i < t1.Element.size(); i++)
      this->Element[i] = t1.Element[i];
    for(unsigned int i = 0; i < t1.Dim.size(); i++)
      this->Dim[i] = t1.Dim[i];
    return *this;
  };

  // Copy constructor
  Tensor::Tensor(const Tensor& t1){
    this->Name = t1.Name;
    this->Order  = t1.Order;
    for(unsigned int i = 0; i < t1.Element.size(); i++)
      this->Element.push_back(t1.Element[i]);
    for(unsigned int i = 0; i < t1.Dim.size(); i++)
     this->Dim.push_back(t1.Dim[i]);
  };

  // Access operator (write)
  TComplex& Tensor::operator()(const unsigned int GlobalIndex){
    return Element[GlobalIndex];
  };

  // More access operator (this time for reading)
  TComplex  Tensor::operator()(const unsigned int GlobalIndex) const{
    return Element[GlobalIndex];
  };

  // Function to check that everything is consistent
  // this checks that the order is the same as the dimension of the vector of dimension
  bool Tensor::GoodSize(unsigned int iDim) const{
    if(this->Order != Dim.size()) { return false; }
    if(this->Order != iDim) { return false; }
    return true;
  };

  // Check that the tensor (this) has the dimension size Dim in the dimension Dim (eugh)
  bool Tensor::GoodSize(unsigned int iDim, unsigned int Dim) const{
    if(this->Order != this->Dim.size()) { return false; }
    if(this->Order <= iDim) { return false; }
    if(this->Dim[iDim] != Dim) { return false; }
    return true;
  };

  // Check that both tensor (this and t1) have same dimension sizes
  bool Tensor::GoodSize(const Tensor &t1) const{
    this->GoodSize(t1.GetOrder());
    t1.GoodSize(this->GetOrder());
    for(unsigned int i = 0; i < this->GetOrder(); i++){
      this->GoodSize(i, t1.GetOrderDim(i));
    }
    return true;
  };

  // Check that the tensor (this) has the dimension Dim
  bool Tensor::AssertGoodSize(unsigned int iDim) const{
    if(this->Order != Dim.size()){
      std::cout << "bool Tensor::AssertGoodSize(int iDim) const" << std::endl;
      std::cout << "Dimension of the tensor has changed!" << std::endl;
      std::cout << "Order            = " << Order << std::endl;
      std::cout << "Dim.size() = " << Dim.size() << std::endl;
      std::cout << "Program will exit" << std::endl;
      exit(1);
    }
    if(this->Order != iDim){
      std::cout << "bool Tensor::AssertGoodSize(int iDim) const" << std::endl;
      std::cout << "The number of dimensions don't match!" << std::endl;
      std::cout << "Order       = " << iDim << std::endl;
      std::cout << "this->Order = " << this->Order << std::endl;
      std::cout << "Program will exit" << std::endl;
      exit(1);
    }
    return true;
    
  };

  // Check that the tensor (this) has the dimension size Order in the dimension Dim (eugh)
  bool Tensor::AssertGoodSize(unsigned int iDim, unsigned int Dim) const{

    if(this->Order != this->Dim.size()){
      std::cout << "bool Tensor::AssertGoodSize(int iDim, int Dim) const" << std::endl;
      std::cout << "Dimension of the tensor has changed!" << std::endl;
      std::cout << "this->Order            = " << this->Order << std::endl;
      std::cout << "this->Dim.size() = " << this->Dim.size() << std::endl;
      std::cout << "Program will exit" << std::endl;
      exit(1);
    }
    if(this->Order <= iDim){
      std::cout << "bool Tensor::AssertGoodSize(int iDim, int Dim) const" << std::endl;
      std::cout << "The number of dimensions don't match!" << std::endl;
      std::cout << "Order       = " << iDim << std::endl;
      std::cout << "this->Order = " << this->Order << std::endl;
      std::cout << "Program will exit" << std::endl;
      exit(1);
    }
    if(this->Dim[iDim] != Dim){
      std::cout << "bool Tensor::AssertGoodSize(int iDim, int Dim) const" << std::endl;
      std::cout << "The dimension size of dimension " << iDim  << " don't match!" << std::endl;
      std::cout << "Dim   = " << iDim << std::endl;
      std::cout << "this->Dim[" << iDim << "] = " << this->Dim[iDim] << std::endl;
      std::cout << "Program will exit" << std::endl;
      exit(1);
    }
      
    return true;

  };

  // Check that both tensor (this and t1) have same dimension sizes
  bool Tensor::AssertGoodSize(const Tensor &t1) const{
    this->AssertGoodSize(t1.GetOrder());
    t1.AssertGoodSize(this->GetOrder());

    for(unsigned int i = 0; i < this->GetOrder(); i++){
      this->AssertGoodSize(i, t1.GetOrderDim(i));
    }
    return true;
  };

  void Tensor::Conjugate(){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = TComplex::Conjugate(this->Element[i]);
  };

  void Tensor::Minus(){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = - this->Element[i];
  };

  void Tensor::OneOverElementWise(){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = 1. / this->Element[i];
  };

  void Tensor::Real(){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i].Re();
  };

  void Tensor::Imaginary(){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i].Im();
  };

  void Tensor::Abs(){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = TComplex::Abs(this->Element[i]);
  };


  // Scalar operations
  void Tensor::Add(const TComplex c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] + c;
  };
  
  void Tensor::Subtract(const TComplex c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] - c;
  };

  void Tensor::Multiply(const TComplex c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] * c;
  };

  void Tensor::Divide(const TComplex c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] / c;
  };


  void Tensor::Add(const double d){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] + d;
  };

  void Tensor::Subtract(const double d){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] - d;
  };

  void Tensor::Multiply(const double d){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] * d;
  };

  void Tensor::Divide(const double d){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] / d;
  };

  // Entry by entry operations
  void Tensor::Add(const Tensor &t1){

    // First check the dimension/sizes of the tensor
    this->AssertGoodSize(t1);

    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] + t1.Element[i];

  };

  void Tensor::Subtract(const Tensor &t1){

    // First check the dimension/sizes of the tensor
    this->AssertGoodSize(t1);

    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] - t1.Element[i];
  };

  void Tensor::MultiplyElementWise(const Tensor &t1){

    // First check the dimension/sizes of the tensor
    this->AssertGoodSize(t1);
  
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] * t1.Element[i];
  };

  void Tensor::DivideElementWise(const Tensor &t1){

    // First check the dimension/sizes of the tensor
    this->AssertGoodSize(t1);

    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] / t1.Element[i];
  };

  void Tensor::Add(const int c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] + (double)c;
  };

  void Tensor::Subtract(const int c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] - (double)c;
  };

  void Tensor::Multiply(const int c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] * (double)c;
  };

  void Tensor::Divide(const int c){
    for(unsigned int i = 0; i < this->Element.size(); i++)
      this->Element[i] = this->Element[i] / (double)c;
  };

  bool Tensor::IsEqual(double d) const{
    for(unsigned int i = 0; i < this->Element.size(); i++)
      if(TComplex::Abs(this->Element[i]) != d)
	return false;
    return true;
  };

  // Boolean operations
  bool Tensor::IsEqual(const Tensor &t1) const{

    // First check the dimension/sizes of the tensor
    this->GoodSize(t1);

    // Then check the element by element
    for(unsigned int i = 0; i < t1.Element.size(); i++){
      if(this->Element[i].Re() != t1.Element[i].Re() || this->Element[i].Im() != t1.Element[i].Im())
	return false;
    }
    return true;
  };

  bool Tensor::IsDifferent(const Tensor &t1) const{
    return !IsEqual(t1);
  };

  bool Tensor::IsBiggerOrEqual(const Tensor &t1) const{

    // First dimensions the dimension/sizes of the tensor
    this->GoodSize(t1);

    bool flag = true;
    unsigned int i = 0;

    // Then check the element by element
    while(flag && i < t1.Element.size()){
      if(TComplex::Abs(this->Element[i]) <= TComplex::Abs(t1.Element[i]))
	flag = false;
    }
    return flag;
  };

  bool Tensor::IsSmaller(const Tensor &t1) const{
    return !IsBiggerOrEqual(t1);
  };
  
  bool Tensor::IsSmallerOrEqual(const Tensor &t1) const{
    return IsSmaller(t1) || IsEqual(t1);
  };
  
  bool Tensor::IsBigger(const Tensor &t1) const{
    return !IsSmallerOrEqual(t1);
  };
  
  bool Tensor::IsDifferent(double d) const{
    return !IsEqual(d);
  };

  bool Tensor::IsBiggerOrEqual(double d) const{
    for(unsigned int i = 0; i < this->Element.size(); i++)
      if(TMath::Abs(d) <= TComplex::Abs(this->Element[i]))
	return false;
    return true;
  };

  bool Tensor::IsSmallerOrEqual(double d) const{
    for(unsigned int i = 0; i < this->Element.size(); i++)
      if(TMath::Abs(d) >= TComplex::Abs(this->Element[i]))
	return false;
    return true;
  };

  bool Tensor::IsBigger (double d) const{
    for(unsigned int i = 0; i < this->Element.size(); i++)
      if(TMath::Abs(d) < TComplex::Abs(this->Element[i]))
	return false;
    return true;
  };

  bool Tensor::IsSmaller (double d) const{
    for(unsigned int i = 0; i < this->Element.size(); i++)
      if(TMath::Abs(d) > TComplex::Abs(this->Element[i]))
	return false;
    return true;
  };

  bool Tensor::IsEqual(TComplex d) const{
    for(unsigned int i = 0; i < this->Element.size(); i++){
      if(this->Element[i].Re() == d.Re() && this->Element[i].Im() == d.Im())
	return false;
    }
    return true;
  };

  bool Tensor::IsDifferent(TComplex d) const{
    return !IsEqual(d);
  };

  bool Tensor::IsBiggerOrEqual(TComplex d) const{
    return IsBiggerOrEqual(TComplex::Abs(d));
  };

  bool Tensor::IsSmallerOrEqual(TComplex d) const{
    return IsSmallerOrEqual(TComplex::Abs(d));
  };

  bool Tensor::IsBigger (TComplex d) const{
    return IsBigger(TComplex::Abs(d));
  };

  bool Tensor::IsSmaller (TComplex d) const{
    return IsSmaller(TComplex::Abs(d));
  };

}

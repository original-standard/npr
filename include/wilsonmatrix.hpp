
#ifndef _WILSONMATRIX_H
#define _WILSONMATRIX_H

#include <complex>
#include "gamma_matrix.hpp"
#include <iostream>
using std::complex;

class wilsonmatrix{
private:
	std::complex<double> * rawdata;
public:
  void destruction() // Do not call this function !
  {
	rawdata = nullptr;
  }
  
  inline std::complex<double> * ret_data(){
    return rawdata;
  }

  inline  std::complex<double>  ret_data_const(int i) const{
    return rawdata[i];
  }
  
  wilsonmatrix(wilsonmatrix && wmat)
  {
    rawdata = wmat.ret_data();
    wmat.destruction();
  }

  wilsonmatrix(wilsonmatrix & wmat)
  {
    rawdata = new complex<double>[144];
    auto dat = wmat.ret_data();
    for(int j = 0;j < 144;j++)
      rawdata[j] = *(dat+j);
  }


  wilsonmatrix(const wilsonmatrix & wmat)
  {
    rawdata = new complex<double>[144];
    for(int j = 0;j < 144;j++)
      rawdata[j] = ret_data_const(j);
  }

  wilsonmatrix operator=(wilsonmatrix wmat)
  {
    rawdata = new complex<double>[144];
    auto dat = wmat.ret_data();
    for(int j = 0;j < 144;j++)
      rawdata[j] = *(dat+j);
    return wmat;
  }
  wilsonmatrix(){
    rawdata = new complex<double>[144];
  }

  wilsonmatrix(double d){
    rawdata = new complex<double>[144];
    for(int j = 0;j < 144;j++)
      rawdata[j] = complex<double>(d,d);
  }

  wilsonmatrix(double *d){
    rawdata = new complex<double>[144];
    for(int j = 0;j < 144;j++)
      rawdata[j] = complex<double>(*(d+j),0);
  }

  wilsonmatrix(complex<double> *c){
    rawdata = new complex<double>[144];
    for(int j = 0;j < 144;j++)
      rawdata[j] = *(c+j);
  }
  
  virtual ~wilsonmatrix(){
    delete[] rawdata;
  }
  void element(int d,int c,int d2, int c2, complex<double> C)
  {
    rawdata[d + 4 * c + 12 * d2 + 48 * c2] = C;
  }

  inline complex<double> operator()(int d, int c, int d2, int c2) const
  {
    return  rawdata[d + 4 * c + 12 * d2 + 48 * c2];
  }
  
  void conj(void)
  {
    for(int i = 0;i < 144;i++)
      rawdata[i] = std::conj(rawdata[i]);
  }



  //operator overloads

  
  wilsonmatrix operator* ( wilsonmatrix &rhs) const
  {
    wilsonmatrix mat(0.);

    auto copy(*this);
    auto lhs = mat.ret_data();
    auto dat = copy.ret_data();
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	for(int a = 0;a < 12;a++) // ???
	  lhs[i + k * 12] += rawdata[i + a * 12] * dat[a + k * 12];
       
      
    return mat;
  }


  wilsonmatrix operator- (const wilsonmatrix &obj)
  {

    wilsonmatrix mat(*this);
    wilsonmatrix copy(*this);
    auto lhs = mat.ret_data();
    auto dat = copy.ret_data();
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	  lhs[i + k * 12]  -= dat[i + k * 12];
       
      
    return mat;
  }

    wilsonmatrix operator+ (const wilsonmatrix &obj) const
  {

    wilsonmatrix mat(*this);
    wilsonmatrix copy(obj);

    auto lhs = mat.ret_data();
    auto dat = copy.ret_data();
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	  lhs[i + k * 12] += dat[i + k * 12];
       
      
    return mat;
  }


  wilsonmatrix operator+= (wilsonmatrix &obj)
  {



    auto lhs = rawdata;
    auto dat = obj.ret_data();
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	  lhs[i + k * 12] += dat[i + k * 12];
       
      
    return *this;
  }


  wilsonmatrix operator* (const double scalar) const
  {

    wilsonmatrix mat(*this);
    auto lhs = mat.ret_data();
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	lhs[i + k * 12] *= scalar;
       
      
    return mat;
  }



  // multiply

  void gr(int dir) // B = A * GAMMA_dir
  {
    
    auto A = rawdata;
    complex<double> * matrix = new complex<double>[16];
    complex<double> * B = new complex<double>[144];
    for(int i = 0;i < 16;i++)
      matrix[i] = complex<double>(gamma_matrix[dir][0][i], gamma_matrix[dir][0][i]);
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	for(int a = 0;a < 12;a++) // ???
	  B[i + k * 12] += A[i + a * 12] * matrix[a % 3 + (k % 3) * 4];
    
  }

  void gl(int dir) // B = GAMMA_dir * A
  {
    
    auto A = rawdata;
    complex<double> * matrix = new complex<double>[16];
    complex<double> * B = new complex<double>[144];
    for(int i = 0;i < 16;i++)
      matrix[i] = complex<double>(gamma_matrix[dir][0][i], gamma_matrix[dir][0][i]);
    for(int i = 0;i < 12;i++) // thread safe ?
      for(int k = 0;k < 12; k++) // thread safe ?
	for(int a = 0;a < 12;a++) // ???
	  B[i + k * 12] += matrix[i % 3 + (a % 3) * 4] * A[a + k * 12];
    
  }

  complex<double> Trace(void)
  {
    complex<double> ret(0.,0.);
    for(int i = 0;i < 12;i++)
      ret += rawdata[i + i * 12];
    return ret;
  }
};










#endif

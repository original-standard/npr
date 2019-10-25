
#ifndef _WILSONMATRIX_H
#define _WILSONMATRIX_H

#include <complex>

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

    wilsonmatrix operator+ (wilsonmatrix &obj)
  {

    wilsonmatrix mat(*this);

    auto lhs = mat.ret_data();
    auto dat = obj.ret_data();
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

  
};
#endif

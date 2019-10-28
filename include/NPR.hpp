#ifndef _NPR_H
#define _NPR_H

#include <stdio.h>
#include <stdlib.h>

#include <complex>
#include <fstream>
#include <sstream>
#include "wilsonmatrix.hpp"
#include <eigen3/Eigen/Dense>
#include "templates.hpp"
#include "operator_overload.hpp"
#include <map>
#include <array>

#ifndef _X_SITES
#define _X_SITES 24
#endif

#ifndef _Y_SITES
#define _Y_SITES 24
#endif

#ifndef _Z_SITES
#define _Z_SITES 24
#endif

#ifndef _T_SITES
#define _T_SITES 64
#endif


#define BIN 1
class NPR{
private:
  int volume_inv[4];
  int volume[4] = {_X_SITES,_Y_SITES,_Z_SITES,_T_SITES};
  int V;
  int loop_N;
  int n_src;
  int mom;
  int count;
  std::map<int,std::vector<wilsonmatrix>> wmat;
  std::map<int,std::vector<wilsonmatrix>> wmato;
  std::map<int,std::vector<wilsonmatrix>> wmati;
public:
  NPR(int,char **);
  ~NPR(){};
  int inline ret_V(){return V;};
  std::vector<std::vector<int>> vectors_sub(std::vector<std::vector<std::complex<double>>> &);
  std::vector<std::vector<int>>  axials_sub(std::vector<std::vector<std::complex<double>>> &);
  // std::vector<std::vector<int>> stensor_sub(std::vector<std::vector<std::complex>> &);
  std::vector<std::vector<int>> sscalar_sub(std::vector<std::vector<std::complex<double>>> &);
};

#endif


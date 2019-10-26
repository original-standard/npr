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





#define BIN 1
class NPR{
private:
  int volume_inv[4];
        int volume[4] = {24,24,24,64};
  //    int volume[4] = {32,32,32,64};
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
  std::vector<std::vector<int>> vectors_sub(std::vector<std::vector<std::complex<double>>> &);
  std::vector<std::vector<int>> axials_sub(std::vector<std::vector<std::complex<double>>> &);
  // std::vector<std::vector<int>> stensor_sub(std::vector<std::vector<std::complex>> &);
  // std::vector<std::vector<int>> sscalar_sub(std::vector<std::vector<std::complex>> &);
};

#endif


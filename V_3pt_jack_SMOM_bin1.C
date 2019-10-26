
#include "NPR.hpp"

using namespace std;

main(int argc,char *argv[])
{
  NPR npr(argc,argv);
  int V = 12;//npr.ret_V();

  vector<vector<complex<double>>> vec;
  auto mom = npr.vectors_sub(vec);
   vector<double> Dev;
   for(int i = 0;i < vec.size();i++)
    {
      

      auto & P =  vec[i];
      complex<double> ev(0.,0.),dev(0.,0.);
      complex<double> W(1.,0.);
      for(auto & j : P){
	ev = ev + W * j;
          
      }

      ev = ev / P.size();
      for(auto & j : P)
	dev = dev + (ev - W * j).real() * (ev - W * j).real();

      dev = dev * (P.size() - 1.) / P.size();
      Dev.push_back(sqrt(dev.real()));
    }

 
  for(int i = 0; i < vec.begin()->size();i++)
    for(int j = 0; j < Dev.size();j++)
      {
	cout << (double)mom[j][0] / V / V * (2. * M_PI) * (2. * M_PI);
	cout << " " << mom[j][1];
	cout << " " << mom[j][2];
	cout << " " << vec[j][i].real() << " " << Dev[j] << endl;
	
      }
  //wmat.resize(0);
}


#include "NPR.hpp"

using std::cout;
using std::endl;
using std::vector;
using std::complex;

int main(int argc,char *argv[])
{
  NPR npr(argc,argv);
  int V = npr.ret_V();

  vector<vector<complex<double>>> sc;
  auto mom = npr.sscalar_sub(sc);
   vector<double> Dev;
   for(int i = 0;i < sc.size();i++)
    {
      

      auto & P =  sc[i];
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

 
  for(int i = 0; i < sc.begin()->size();i++)
    for(int j = 0; j < Dev.size();j++)
      {
	cout << (double)mom[j][0] / V / V * (2. * M_PI) * (2. * M_PI);
	cout << " " << mom[j][1];
	cout << " " << mom[j][2];
	cout << " " << sc[j][i].real() << " " << Dev[j] << endl;
	
      }

  return 0;
}

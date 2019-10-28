#include "NPR.hpp"
#include "wilsonmatrix.hpp"

using std::vector;
using std::map;
vector<vector<int>> NPR::sscalar_sub(vector<vector<complex<double>>> &prop){
  vector<vector<int>> prop_n;

  for(auto & mom1 : wmat)
   for(auto & mom2 : wmat)
    {
      int k = mom1.first;
      int j = mom2.first;

      int p1[4];
      int l = k;
      for(int i = 0;i < 4;i++)
	{
	  p1[i] = l % (2 * mom + 1) - mom;
	  l /= (2 * mom + 1);
	}
      
      int p2[4];
      l = j;
      for(int i = 0;i < 4;i++)
        {
          p2[i] = l % (2 * mom + 1) - mom;
          l /= (2 * mom + 1);
        }

      int p_1 = 0;
      for(int i = 0;i < 4;i++) p_1 += p1[i] * p1[i] * volume_inv[i] * volume_inv[i];
      int p_2 = 0;
      for(int i = 0;i < 4;i++) p_2 += p2[i] * p2[i] * volume_inv[i] * volume_inv[i];
      int p_3 = 0;
      for(int i = 0;i < 4;i++) p_3 += (p1[i]-p2[i]) * (p1[i] -p2[i]) * volume_inv[i] * volume_inv[i];


      
if(p_1 == 0)
continue;
      if(p_2 != p_3)
continue;
        if(p_1 != p_2)
continue;

	//      if(hyper(p1,p2))
	//	continue;

	vector<complex<double>> P(loop_N,complex<double>(0.,0.));


	  {
	    vector<wilsonmatrix> Gamma_insert(loop_N,wilsonmatrix(0.0));
	    wilsonmatrix median(0.0);
	    for(int i = 0;i < loop_N;i++)
	      {
		Gamma_insert[i] = 0.;
		for(int m = 0;m < n_src;m++)
		  {
		    //		auto & gamma = Gamma_insert[i];
		    auto gamma = wmato[j][i];
		    gamma.conj();
		    gamma.gr(5);
		    gamma.gl(5);
		    

		
		    Gamma_insert[i] =  Gamma_insert[i] +  wmato[k][i] * gamma;
		  }
		median = median + Gamma_insert[i];
	      }
	    median = median * (1. /  loop_N);
	    std::transform(Gamma_insert.begin(),Gamma_insert.end(),Gamma_insert.begin(),
			   Jackknife<wilsonmatrix>(loop_N,median));
	    for(int i = 0;i < loop_N;i++)
	      {
		wilsonmatrix Si1,Si2;
		Si1 = wmati[k][i];
		Si2 = wmati[j][i];
		Si2.conj();
		Si2.gr(5);
		Si2.gl(5);
		
		  
		complex<double> C = (Si1 * Gamma_insert[i] * Si2).Trace() / 12. / n_src;
		
		P[i] = P[i] + C;
				
	      }
	  }
	prop.push_back(P);
	vector<int> tmp;
	tmp.push_back(p_1);
	tmp.push_back(mom1.first);
	tmp.push_back(mom2.first);
	prop_n.push_back(tmp);


      }
  return prop_n;
}

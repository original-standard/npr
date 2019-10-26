#include "NPR.hpp"

using std::vector;
using std::map;


vector<vector<int>> NPR::vectors_sub(vector<vector<complex<double>>> &prop){
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



	vector<complex<double>> P(loop_N,complex<double>(0.,0.));
 	vector<wilsonmatrix> Gamma_insert(loop_N,wilsonmatrix(0.0));
	for(int dir1 = 0;dir1 < 4;dir1++)
	  {
	    vector<wilsonmatrix> Gamma_insert_dir(loop_N,wilsonmatrix(0.0));
	    wilsonmatrix median(0.0);
	    for(int i = 0;i < loop_N;i++)
	      {
		//for(int m = 0;m < n_src;m++)
		  {
		    // auto gamma = wmato[m][j][i];
		    auto gamma = wmato[j][i];

		    gamma.conj();

		    gamma.gr(5);
		    gamma.gl(5);
		    //		    gamma.gl(-5);
		    gamma.gl(dir1); //operator

		    Gamma_insert_dir[i] = wmato[k][i] * gamma * (p1[dir1] -p2[dir1]) * volume_inv[dir1];
		    //Gamma_insert_dir[i] += wmato[k][i] * gamma;



		  }
		median = median + Gamma_insert_dir[i];
	      }

	    

	    median = median * (1. /  loop_N);
	    std::transform(Gamma_insert_dir.begin(),Gamma_insert_dir.end(),Gamma_insert_dir.begin(),
			   Jackknife<wilsonmatrix>(loop_N,median));
	    Gamma_insert = Gamma_insert + Gamma_insert_dir;



	  }
	for(int i = 0;i < loop_N;i++){
	  for(int dir2 = 0;dir2 < 4;dir2++){
	    wilsonmatrix Si1,Si2;
	    Si1 = wmati[k][i];
	    Si2 = wmati[j][i];
	    Si2.conj();
	    Si2.gr(5);
	    Si2.gl(5);
	    
	    //	    Si2.gr(-5);
	    Si2.gr(dir2); //projector
	    
	    complex<double> C = (Si1 * Gamma_insert[i] * Si2).Trace()  * (p1[dir2] -p2[dir2]) * volume_inv[dir2] / 12. / n_src / p_3;


	    //	    Float R = 0.,I = 0.;
	    

	    P[i] = P[i] + C;

	  }
	}
	//	prop[p_1] = P + prop[p_1];
	//	prop_n[p_1]++;
	prop.push_back(P);
	vector<int> tmp;
	tmp.push_back(p_1);
	tmp.push_back(mom1.first);
	tmp.push_back(mom2.first);
	prop_n.push_back(tmp);
    }

 return prop_n;
//  for(auto &i : prop)
//	i.second = i.second / prop_n[i.first];
}

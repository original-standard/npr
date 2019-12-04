#include "NPR.hpp"

using std::vector,std::map;
using Eigen::Matrix4cd;
typedef std::complex<double> cd;
vector<vector<int>> NPR::threequark_sub(vector<vector<complex<double>>> &prop,int op1, int op2){
  
  //void NPR::threequark_sub(map<unsigned int,vector<cd>> &prop){


  
  Matrix4cd gamma_x;
  gamma_x << cd(0,0) , cd(0,0) , cd(0,0) , cd(1,0)
    ,cd(0,0) , cd(0,0) , cd(-1,0) , cd(0,0)
    ,cd(0,0) , cd(-1,0) , cd(0,0) , cd(0,0)
    ,cd(1,0) , cd(0,0) , cd(0,0) , cd(0,0);

  Matrix4cd gamma_y;
  gamma_y << cd(0,0) , cd(0,0) , cd(0,-1) , cd(0,0)
    ,cd(0,0) , cd(0,0) , cd(0,0) , cd(0,1)
    ,cd(0,1) , cd(0,0) , cd(0,0) , cd(0,0)
    ,cd(0,0) , cd(0,-1) , cd(0,0) , cd(0,0);

  Matrix4cd gamma_z;
  gamma_z << cd(1,0) , cd(0,0) , cd(0,0) , cd(0,0)
    ,cd(0,0) , cd(1,0) , cd(0,0) , cd(0,0)
    ,cd(0,0) , cd(0,0) , cd(-1,0) , cd(0,0)
    ,cd(0,0) , cd(0,0) , cd(0,0) , cd(-1,0);

  Matrix4cd gamma_t;
  gamma_t << cd(0,0) , cd(0,0) , cd(0,0) , cd(0,-1)
    ,cd(0,0) , cd(0,0) , cd(0,-1) , cd(0,0)
    ,cd(0,0) , cd(0,1) , cd(0,0) , cd(0,0)
    ,cd(0,1) , cd(0,0) , cd(0,0) , cd(0,0);

  Matrix4cd gamma_5;
  gamma_5 << cd(0,0) , cd(0,0) , cd(1,0) , cd(0,0)
    ,cd(0,0) , cd(0,0) , cd(0,0) , cd(1,0)
    ,cd(1,0) , cd(0,0) , cd(0,0) , cd(0,0)
    ,cd(0,0) , cd(1,0) , cd(0,0) , cd(0,0);
  vector<vector<int>> prop_n;
  for(auto &mom1 : wmat)
    {
      int k = mom1.first;
      int j = k;
      int p1[4];
      int p2[4];
      int p3[4];
      int l = k;
      for(int i = 0;i < 4;i++)
	{
	  p1[i] = l % (2 * mom + 1) - mom;
	  l /= (2 * mom + 1);
	}
      
      l = k;
      for(int i = 0;i < 4;i++)
        {
          p2[i] = l % (2 * mom + 1) - mom;
          l /= (2 * mom + 1);
        }

      l = k;
      for(int i = 0;i < 4;i++)
        {
          p3[i] = l % (2 * mom + 1) - mom;
          l /= (2 * mom + 1);
        }
      
      int color_in [][3] = {
	{0,1,2},
	{0,2,1},
	{1,0,2},
	{1,2,0},
	{2,1,0},
	{2,0,1}
      }; 
      int p_1 = 0;
      for(int i = 0;i < 4;i++) p_1 += p1[i] * p1[i] * volume_inv[i] * volume_inv[i];
      int p_2 = 0;
      for(int i = 0;i < 4;i++) p_2 += p2[i] * p2[i] * volume_inv[i] * volume_inv[i];
      int p_3 = 0;
      for(int i = 0;i < 4;i++) p_3 += p3[i] * p3[i] * volume_inv[i] * volume_inv[i];
      if(p_1 != p_2)
	continue;
      if(p_2 != p_3)
	continue;
      {
	vector<cd> X(loop_N,cd(0.,0.));
	vector<vector<cd>> P(loop_N,vector<cd>(12 * 12 * 12 * 4,(0.,0.)));
	vector<cd> median(12 * 12 * 12 * 4,cd(0.,0.));

	Matrix4cd G1 = gamma_5 * gamma_5;
	Matrix4cd G2 = gamma_5 * gamma_5;

	Matrix4cd G1p = gamma_5 * gamma_5;
	Matrix4cd G2p = gamma_5 * gamma_5;
	switch (op1 % 10){
	case 0: // SS
	  G1 = gamma_t * gamma_y;
	  G2 = gamma_5 * gamma_5;
	  break;
	case 1: // PP
	  G1 = gamma_t * gamma_y  * gamma_5;
	  G2 = gamma_5;
	  break;
	case 2: // SP
	  G1 = gamma_t * gamma_y;
	  G2 = gamma_5;
	  break;
	case 3: // PS
	  G1 = gamma_t * gamma_y  * gamma_5;
	  G2 = gamma_5 * gamma_5;
	  break;
	case 4: // left left
	  G1 = gamma_t * gamma_y * (gamma_5 * gamma_5 - gamma_5);
	  G2 = (gamma_5 * gamma_5 - gamma_5);
	  break;
	case 5: // right left
	  G1 = gamma_t * gamma_y * (gamma_5 * gamma_5 + gamma_5);
	  G2 = (gamma_5 * gamma_5 -  gamma_5);
	  break;
	  //to be continued...
	}

	switch (op2 % 10){
	case 0: // SS
	  G1p = gamma_t * gamma_y;
	  G2p = gamma_5 * gamma_5;
	  break;
	case 1: // PP
	  G1p = gamma_t * gamma_y  * gamma_5;
	  G2p = gamma_5;
	  break;
	case 2: // SP
	  G1p = gamma_t * gamma_y;
	  G2p = gamma_5;
	  break;
	case 3: // PS
	  G1p = gamma_t * gamma_y  * gamma_5;
	  G2p = gamma_5 * gamma_5;
	  break;
	case 4: // left left
	  G1p = gamma_t * gamma_y * (gamma_5 * gamma_5 - gamma_5) / 2;
	  G2p = (gamma_5 * gamma_5 - gamma_5) / 2;
	  break;
	case 5: // right left
	  G1p = gamma_t * gamma_y * (gamma_5 * gamma_5 + gamma_5) / 2;
	  G2p = (gamma_5 * gamma_5 - gamma_5) / 2;
	  break;
	  //to be continued...
	}

	//	auto G1 = gamma_t * gamma_y * gamma_5;
	//	auto G2 = gamma_5 * gamma_5;
	for (int i = 0;i < loop_N;i++)
	  {
	    auto &R = P[i];
	    
	   
	    //	  for(int c1 = 0;c1 < 3;c1++)
	    //	  for(int c2 = 0;c2 < 3;c2++)
	    //	  for(int c3 = 0;c3 < 3;c3++)
	    



#pragma omp parallel for
	    for(int color_spin = 0;color_spin < 3 * 3 * 3 * 4 * 4 * 4 ;color_spin++)
	      {
		int color = color_spin / (4 * 4 * 4);
		int c1 = color % 3;
		int c2 = (color / 3) % 3;
		int c3 = (color / (3 * 3)) % 3;
		int spin = color_spin % (4 * 4 * 4);
		
		int a = spin % 4;
		int b = (spin / 4) % 4;
		int c = (spin / (4 * 4)) % 4;				   
		auto prop = wmato[j][i];
		//for(int a = 0;a < 4;a++)
		//for(int b = 0;b < 4;b++)
		//	  for(int c = 0;c < 4;c++)
		//		      for(int spin1 = 0; spin1 < 4 * 4 * 4;spin1++)
		//		cd G1lhs = cd(G1(a,b).real(),G1(a,b).imag());
		//if(abs(G1lhs) == 0)
		//  continue;
		//		if(abs(G1rhs) == 0)
		//  continue;
		for(int m = 0;m < 4;m++)
		  for(int n = 0;n < 4;n++)
		    {
		      cd G1rhs = cd(G1(m,n).real(),G1(m,n).imag());
		      if(abs(G1rhs) == 0)
		        continue;
		      
		      for(int x = 0;x < 4;x++)
			for(int d = 0;d < 4;d++)
			  for(int c4 = 0;c4 < 6;c4++)
			    {
			      cd G2rhs = cd(G2(c,d).real(),G2(c,d).imag());
			      if(abs(G2rhs) == 0)
				continue;
			      cd p1 = prop(a,c1,m,color_in[c4][0]);
			      cd p2 = prop(b,c2,n,color_in[c4][1]);
			      cd p3 = prop(d,c3,x,color_in[c4][2]);
			      
			      //R = R + p1 * p2 * G1lhs * G1rhs / 3.;
			      R[c1 * 4 + a + (c2 * 4 + b) * 12 + (c3 * 4 + c) * 12 * 12 + x * 12 * 12 * 12 ] =
				R[c1 * 4 + a + (c2 * 4 + b) * 12 + (c3 * 4 + c) * 12 * 12 + x * 12 * 12 * 12] + p1 * p2 * p3 * G1rhs * G2rhs;
			    }
		    }
	      }
	    median = median + R;
	  }
	
	
		   
	median = median * (1. /  loop_N);
	std::transform(P.begin(),P.end(),P.begin(),
		       Jackknife<vector<cd>>(loop_N,median));
	


#pragma omp parallel for
	for (int i = 0;i < loop_N;i++)
	  {
	    auto &R = P[i];
	    auto prop = wmati[j][i];
	    
	    for(int a = 0;a < 4;a++)
	      for(int b = 0;b < 4;b++)
		for(int m = 0;m < 4;m++)
		  for(int n = 0;n < 4;n++)
		    {
		      
		      cd G1rhs = cd(G1p(m,n).real(),G1p(m,n).imag());
		      //cd G1rhs = cd(G1p(n,m).real(),G1p(n,m).imag());
		      //if(abs(G1lhs) == 0)
		      //	continue;
		      if(abs(G1rhs) == 0)
			continue;
		      for(int c4 = 0;c4 < 6;c4++)
			for(int c1 = 0;c1 < 3;c1++)
			  for(int c2 = 0;c2 < 3;c2++)
			    {
			      cd p1 = prop(m,color_in[c4][0],a,c1);
			      cd p2 = prop(n,color_in[c4][1],b,c2);
			      for(int c3 = 0;c3 < 3;c3++)
				for(int x = 0;x < 4;x++)
				  for(int c = 0;c < 4;c++)
				    for(int d = 0;d < 4;d++)
				      {
					cd G2rhs = cd(G2p(c,d).real(),G2p(c,d).imag());
					//					cd G2rhs = cd(G2p(d,x).real(),G2p(d,x).imag());
					cd p3 = prop(x,color_in[c4][2],d,c3);
					X[i] = X[i] + R[x * 12 * 12 * 12  + c1 * 4 + a + (c2 * 4 + b) * 12 + (c3 * 4 + c) * 12 * 12] * p1 * p2 * p3 *  G1rhs * G2rhs / 96.;
					//X[i] = X[i] + prop2(a,c1,n,c3) * wmato[j][i](a,c1,n,c3);
				      }
			    }
		    }
	  }
	prop.push_back(X);
	vector<int> tmp;
	tmp.push_back(p_1);
	tmp.push_back(mom1.first);
	tmp.push_back(mom1.first);
	//	tmp.push_back(mom2.first);
	prop_n.push_back(tmp);
	
      }
    }
  
  return prop_n;
}

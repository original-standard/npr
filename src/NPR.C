
#include "NPR.hpp"
#include <iostream>

using std::cout,std::endl;

using std::vector,std::map,Eigen::MatrixXcd;

using std::ifstream,std::stringstream,std::string;
NPR :: NPR(int argc, char *argv[])
  : mom(16) {
  n_src = 1;


  char * conf_file = argv[1];  
  int loop_s = atoi(argv[2]);
  int loop_e = atoi(argv[3]);
  int loop_skp = atoi(argv[4]);

  loop_N = (loop_e - loop_s) / loop_skp + 1;



  int loop_i = 0;

  map<int,wilsonmatrix> WSi;
  map<int,wilsonmatrix>  WS; 


  int V_inv = 1;
  for(int i = 1;i < volume[3];i++)
    {
    for(int j = 0;j < 4;j++)
      if(volume[j] % i)
	goto escape;
        V_inv = i;
    escape:
    ;
    }
  V = 1;
  for(int i = 0;i < 4; i++)
    V *= volume[i] / V_inv;
  for(int i = 0;i < 4; i++){
    volume_inv[i] = V / (volume[i] / V_inv);
  }
  V *= V_inv;
  //  V *= V_inv;

  for(int loop=loop_s; loop<=loop_e; loop+=loop_skp){


    char lat_fname[200];

    sprintf(lat_fname,"%s.%d",conf_file,loop) ;

    ifstream in(lat_fname);
    printf("reading %s ...\n",lat_fname);
    std::string buf;
    char status = 0; // 0 , 1 standby, 2 reading                                                                                    
    int x[4];
    int p[4];


    map<int,wilsonmatrix> WS_l;

    int src = 0;
      while(getline(in,buf))
	{
	  stringstream line(buf);
	  line >> buf;
          if(!buf.compare("Cg")){
            status = 0;
          }
          if(!buf.compare("Fourier"))
            {
              status = 2;
              for(int i = 0;i < 4;i++){
                line >> buf >> p[i];
		//		std::cout << p[i] << std::endl;
	      }
                int index = 0;
                for(int i = 0;i < 4;i++)
                  index = index * (2 * mom + 1) + p[3 - i] + mom;
		auto dot = [](auto *p1,auto *p2) -> double
		  {
		    double prod = 0.;
		    for(int i = 0; i < 4;i++){
		      prod += p1[i] * p2[i];
		    }
		    return prod;
		  };
		double elements[8][4] = {
		  {1./2.,1./2.,1./2.,1./2.},
		  {-1./2.,1./2.,1./2.,1./2.},
		  {1./2.,-1./2.,1./2.,1./2.},
		  {1./2.,1./2.,-1./2.,1./2.},
		  {1./2.,1./2.,1./2.,-1./2.},
		  {-1./2.,-1./2.,1./2.,1./2.},
		  {1./2.,-1./2.,-1./2.,1./2.},		  
		  {1./2.,-1./2.,1./2.,-1./2.}
		}; 
		
		char flag = 0;
		flag = 1; // always pass
		for(int j = 0;j < 8;j++){
		  auto prod = dot(p,elements[j]);
		  double q[4];
		  for(int i = 0;i < 4;i++)
		    q[i] = p[i] - prod * elements[j][i];
		  if(dot(q,q) / dot(p,p)   < .5){
		    flag = 1;
		    //		    for(int i = 0;i < 4;i++)
		      //		      std::cout << p[i] << " ";
		    
		      //    std::cout << prod << " " << j <<  " " << std::endl;
		    break;
		  }
		}

                for(int s1=0; s1<4; ++s1) 
                  for(int c1=0; c1<3; ++c1)
                    for(int s2=0; s2<4; ++s2) 
                      for(int c2=0; c2<3; ++c2){
                        getline(in,buf);
                        stringstream line(buf);
			double r,i;
			int a,b,c,d;
			line >> a >> b >> c >> d >> r >> i;
			wilsonmatrix temp(0.);
			if(flag){
			  temp.element(a,b,c,d,complex<double>(r  * 2 / (5.-1.8) ,i  * 2 / (5.-1.8)));
			  WS_l[index] = WS_l[index] + temp;
			}
		      }

	    }
	}
      for(auto &j : WS_l)
	{
	  auto i = j.first;
	    WS[i] += WS_l[i];
	    wmat[i].push_back(WS_l[i]);
	    wmato[i].push_back(WS_l[i]);
	    wmati[i].push_back(WS_l[i]);
	}
      loop_i++;    
  }
  //    auto wmato(wmat);
  for(auto &i : WS)
        {
          auto j = i.first;
    WS[j] = WS[j] * (1. / loop_N);
    auto &wm = wmat[j];
    auto &wmi = wmati[j];
    std::transform(wm.begin(),wm.end(),wm.begin(),Jackknife<wilsonmatrix>(wm.size(),WS[j]));
    for(int i = 0;i < wm.size();i++){
      MatrixXcd temp(12,12);
      auto &W = wm[i];
      auto &Wi = wmi[i];

      for(int color_sink = 0;color_sink < 3;color_sink++)
	for(int dirac_sink = 0;dirac_sink < 4;dirac_sink++)
	  for(int color_source = 0;color_source < 3;color_source++)
	    for(int dirac_source = 0;dirac_source < 4;dirac_source++)
	      {
		int inner = color_sink * 4 + dirac_sink;
		int outer = color_source * 4 + dirac_source;
		temp(inner,outer) = std::complex<double>(W(dirac_sink,color_sink,dirac_source,color_source).real(),
							 W(dirac_sink,color_sink,dirac_source,color_source).imag());
	      }
      
      //      MatrixXd real_tmp = real;
      //MatrixXd imag_tmp = imag;
      //      real = real_tmp.inverse() + imag_tmp.inverse() * real_tmp * imag_tmp.inverse();
      //imag = -real_tmp.inverse() * (real_tmp * imag_tmp.inverse() + imag_tmp * real_tmp.inverse()).inverse();
      //real = -imag * real_tmp * imag_tmp.inverse();
      MatrixXcd tmp_inv = temp.inverse();
      for(int color_sink = 0;color_sink < 3;color_sink++)
	for(int dirac_sink = 0;dirac_sink < 4;dirac_sink++)
	  for(int color_source = 0;color_source < 3;color_source++)
	    for(int dirac_source = 0;dirac_source < 4;dirac_source++)
	      {
		int inner = color_sink * 4 + dirac_sink;
		int outer = color_source * 4 + dirac_source;
		Wi.element(dirac_sink,color_sink,dirac_source,color_source,complex<double>(tmp_inv(inner,outer).real(),
										    tmp_inv(inner,outer).imag()));
		//std::cout << i << " " << tmp_inv(inner,outer).real() << std::endl;
	      }
    }
  }


}


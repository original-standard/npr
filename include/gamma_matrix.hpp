
#ifndef _GAMMA_MATRIX_H
#define _GAMMA_MATRIX_H

static double gamma_matrix[6][2][16] =
  {
    {
      { // G0 real
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.
      },
      { // G0 imag
	0.,0.,0.,-1.,
	0.,0.,-1.,0.,
	0.,1.,0.,0.,
	1.,0.,0.,0.
      }
    },
    {
      { // G1 real
	0.,0.,0.,1.,
	0.,0.,-1.,0.,
	0.,-1.,0.,0.,
	1.,0.,0.,0.
      },
      { // G1 imag
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.
      }
    },
    {
      { // G2 real
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.	
      },
      { // G2 imag
	0.,0.,-1.,0.,
	0.,0.,0.,1.,
	1.,0.,0.,0.,
	0.,-1.,0.,0.
      }
    },
    {
      { //G3 real
	1.,0.,0.,0.,
	0.,1.,0.,0.,
	0.,0.,-1.,0.,
	0.,0.,0.,-1.
      },
      { //G3 imag
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.
      }
    },
    {
      { //G4???? real
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.
      },
      { //G4????
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.
      }
    },
    {
      { //G5 real
	0.,0.,1.,0.,
	0.,0.,0.,1.,
	1.,0.,0.,0.,
	0.,1.,0.,0.
      },
      { //G5 imag
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.,
	0.,0.,0.,0.
      }
    }
  };


#endif

///////////////////////////////////////////////////////////////////////
// Linear Algebra for B-Spline (LABS)
//
// This external function calculate the matrices and vectors 
// to use B-Spline
//
// It needs to be written in CPP file because B-spline module in GSL
//
// Developer : M.K. Lee
// Last updated : 2013.01.03
//
// Bug report : mk@ices.utexas.edu
///////////////////////////////////////////////////////////////////////


#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <string>
#include <fstream>
#define PI 3.141592653589793238462643

using namespace std;

void bspline_integration_coefficients(const size_t nderiv, double *coeffs, size_t inc,  gsl_matrix *dB, gsl_bspline_workspace *w, gsl_bspline_deriv_workspace *dw);

extern "C" {
  void LABS_fill_full(double *M0, double *M1, double *M2, double *M3, double *M4, double *yc, double *yknot, double *int_coeff, double *center_coeff, int *K, int *N, double *STRETCH, double *Ly);
}


void LABS_fill_full (double *M0, double *M1, double *M2, double *M3, double *M4, double *yc, double *yknot, double *int_coeff, double *center_coeff, int *K, int *N, double *STRETCH, double *Ly){

  int i, j;

  // nknots = number of B-splines + (k-1)
  //        = nbreak - 1 + (k-1)
  //        = nbreak + k - 2
  // nbreak = nknots -k + 2
  const double pi = PI;
  const size_t k = *K;                           // the order of B-Spline, eg. k=4 cubic
  const size_t nknots = *N;                      // Number of knots
  const size_t nbreak = nknots - k + 2;
  const double stretch = *STRETCH;
  const double ly = *Ly;
  // allocation GSL variables
  gsl_bspline_workspace *BW;                    //The information of B-Splines  
  gsl_bspline_deriv_workspace *DBW;             //The workspace of B-Splines Derivatives
  gsl_matrix *BM;                               //Matrix of Bsplines at collocation points
  gsl_vector *y;                                //breakpoint
  gsl_vector *B1;                               //Bsplines at a given collocation point
  gsl_vector *ycol;                             //Collocation Points
  gsl_matrix *DB;                               //Derivative of Bspline at collocation points
  gsl_matrix *DB_i;                             //Derivative of Bspline at collocation points_integration
  gsl_matrix *DB1;                              //Matrix of 1st order derivative at collocation points
  gsl_matrix *DB2;                              //Matrix of 2nd order derivative at collocation points
  gsl_matrix *DB3;                              //Matrix of 3rd order derivative at collocation points
  gsl_matrix *DB4;                              //Matrix of 4th order derivative at collocation points
  

  BW   = gsl_bspline_alloc(k,nbreak);
  DBW  = gsl_bspline_deriv_alloc(k);
  B1   = gsl_vector_alloc(nknots);
  y    = gsl_vector_alloc(nbreak);
  BM   = gsl_matrix_alloc(nknots, nknots);
  ycol = gsl_vector_alloc(nknots);
  DB   = gsl_matrix_alloc(nknots,5);
  DB1  = gsl_matrix_alloc(nknots,nknots);
  DB2  = gsl_matrix_alloc(nknots,nknots);
  DB3  = gsl_matrix_alloc(nknots,nknots);
  DB4  = gsl_matrix_alloc(nknots,nknots);
  DB_i   = gsl_matrix_alloc(k,3);


  // calculate the location of breakpoints
  // cout << "The location of break points" << endl;
  if (abs(stretch) <= 0.0000000000000001) 
    for (i=0; i<nbreak; ++i){
      double yi = -1.0 + 2.0/(nbreak-1)*i;
      gsl_vector_set(y,i,yi);
      yknot[i] = yi;
    }
  else{
    for (i=0; i<nbreak;++i){
      // double yi = pi*(1+cos(-pi*(nbreak-i-1)/(nbreak-1)));
      //  double yi = 0.5*ly*tanh(stretch*(-1.0+(2.0*i)/(nbreak-1)))/tanh(stretch);

      // Jimenez Grid with Stretch Factor
      double yi = sin(stretch*PI/2.0*(-1.0+(2.0*i)/(nbreak-1)))/sin(stretch*PI/2.0); 

      // Chebyshev Polynomial Grid Spacing
      // double yi = -cos(0.5*PI*((2.0*i)/(nbreak-1))); 

      // printf("%24.16e \n",yi);
      gsl_vector_set(y,i,yi);
      yknot[i]=yi;
    }
  }
  // Calculate the location of knots (internally)
  gsl_bspline_knots(y,BW);

  //  cout << "The location of collocation points" << endl;
  for (i=0; i<nknots; ++i){
    double yci = gsl_bspline_greville_abscissa (i,BW);
    // cout << yci << endl;
    gsl_vector_set(ycol,i,yci);
    gsl_bspline_deriv_eval(yci,4,DB,BW,DBW);
    for (j=0; j<nknots; ++j){
      gsl_matrix_set(DB1,i,j,gsl_matrix_get(DB,j,1));
      gsl_matrix_set(DB2,i,j,gsl_matrix_get(DB,j,2));
      gsl_matrix_set(DB3,i,j,gsl_matrix_get(DB,j,3));
      gsl_matrix_set(DB4,i,j,gsl_matrix_get(DB,j,4));
    }
  }
  

  for(i=0;i<nknots;++i){
    double yci = gsl_vector_get(ycol,i);
    // Calculate the B-Spline Basis Function at given location
    // in this case, the location is breakpoint
    gsl_bspline_eval(yci,B1,BW);
   
    // Store the 1 column of B-spline Basis Functions into a Matrix
    for (j=0;j<nknots;++j){
      double Bj = gsl_vector_get(B1,j);
      gsl_matrix_set(BM,i,j,Bj);
    }
  }
  bspline_integration_coefficients(0,int_coeff,1,DB_i,BW,DBW);

  for(i=0;i<nknots;++i){
    for(j=0;j<nknots;++j){
      M0[j+i*nknots] = gsl_matrix_get(BM,i,j);
      M1[j+i*nknots] = gsl_matrix_get(DB1,i,j);
      M2[j+i*nknots] = gsl_matrix_get(DB2,i,j);
      M3[j+i*nknots] = gsl_matrix_get(DB3,i,j);
      M4[j+i*nknots] = gsl_matrix_get(DB4,i,j);
    }
    yc[i] = gsl_vector_get(ycol,i);
  }


  gsl_bspline_eval(0.0,ycol,BW);
  for(i=0;i<nknots;++i){
    center_coeff[i] = gsl_vector_get(ycol,i);
  }



  /*
  // print out the matrix
  printf("\nThe matrix of value of B-Spline basis functions at given collocation points\n"); 
  for(i=0;i<nknots;++i){
    for(j=0;j<nknots;++j){
      //      printf ("%f ", M0[j+i*nknots]);
        printf ("%14.4e ", gsl_matrix_get(BM,i,j));
    }
    printf ("\n");
  }

  printf("\nThe matrix of value of 1st order derivative of B-Spline basis functions at given collocation points\n"); 
  for(i=0;i<nknots;++i){
    for(j=0;j<nknots;++j){
      // printf ("%f ", M1[j+i*nknots]);
      printf ("%14.4e ", gsl_matrix_get(DB1,i,j));
    }
    printf ("\n");
  }
    
  printf("\nThe matrix of value of 2nd order derivative of B-Spline basis functions at given collocation points\n"); 
  for(i=0;i<nknots;++i){
    for(j=0;j<nknots;++j){
      //printf ("%f ", M2[j+i*nknots]);
      printf ("%14.4e ", gsl_matrix_get(DB2,i,j));
    }
    printf ("\n");
  }
  */
  gsl_bspline_free(BW);
  gsl_bspline_deriv_free(DBW);	 
  gsl_vector_free(B1);
  gsl_vector_free(y);
  gsl_matrix_free(BM);
  gsl_vector_free(ycol);
  gsl_matrix_free(DB);
  gsl_matrix_free(DB_i);
  gsl_matrix_free(DB1);
  gsl_matrix_free(DB2);
  gsl_matrix_free(DB3);
  gsl_matrix_free(DB4);

} // end of function




// This routine is written by Rhys Ulerich
// This routine calculate the coefficient vector to calculate the integration over given collocation points
// Example
// U = Bu      <----- U: Original vector / B: Bspline Mass Matrix / u: Bspline coeff for U
// int(U) = (coeff,u) <--- inner product

void bspline_integration_coefficients(
					  const size_t nderiv,
					  double * coeffs,
					  size_t inc,
					  gsl_matrix *dB,
					  gsl_bspline_workspace *w,
					  gsl_bspline_deriv_workspace *dw)
{
  /* Obtain an appropriate order Gauss-Legendre integration rule */
    gsl_integration_glfixed_table * const tbl
      = gsl_integration_glfixed_table_alloc((w->k - nderiv + 1)/2);
    //    if (tbl == NULL) {
    //SUZERAIN_ERROR("failed to obtain Gauss-Legendre rule from GSL",
    //		     SUZERAIN_ESANITY);
    //}

    /* Zero integration coefficient values */
    for (size_t i = 0; i < w->n; ++i) {
      coeffs[i * inc] = 0.0;
    }

    /* Accumulate the breakpoint-by-breakpoint contributions to coeffs */
    double xj = 0, wj = 0;
    for (size_t i = 0; i < (w->nbreak - 1); ++i) {

      /* Determine i-th breakpoint interval */
      const double a = gsl_bspline_breakpoint(i,   w);
      const double b = gsl_bspline_breakpoint(i+1, w);

      for (size_t j = 0; j < tbl->n; ++j) {

	/* Get j-th Gauss point xj and weight wj */
	gsl_integration_glfixed_point(a, b, j, &xj, &wj, tbl);

	/* Evaluate basis functions at point xj */
	size_t kstart, kend;
	gsl_bspline_deriv_eval_nonzero(xj, nderiv,
				       dB, &kstart, &kend, w, dw);

	/* Accumulate weighted basis evaluations into coeffs */
	for (size_t k = kstart; k <= kend; ++k) {
	  coeffs[k * inc] += wj * gsl_matrix_get(dB,
						 k - kstart, nderiv);
	}
      }
    }

    /* Free integration rule resources */
    gsl_integration_glfixed_table_free(tbl);

}

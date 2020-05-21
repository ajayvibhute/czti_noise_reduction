/* 
 * File:   cztmatrix.h
 * Author: tanul
 *
 * Created on July, 5, 2013, 09:33 AM
 */

#ifndef CZTMATRIX_H
#define CZTMATRIX_H


//**********************Include files***********************
#include<linalg.h>
#include<ap.h>
#include "glog/logging.h"


using namespace std;
using namespace alglib;
/******************************************************************************************************
 MATRIX FUNCTIONS
 ******************************************************************************************************/
/**
 * To multiply two matrices A and B with sizes mxn and nxp to give output matrix C with size mxp
 * @param A
 * @param B
 * @param C
 * @param m
 * @param n
 * @param p
 * @return 
 */
int matrixproduct(real_2d_array &A,real_2d_array &B,real_2d_array &C,
        ae_int_t m,ae_int_t n,ae_int_t p);

/**
 * Multiplies inverse of A with B to get the result in X array  
 * @param Aarr
 * @param Barr
 * @param m
 * @param n
 * @param p
 * @param X
 * @return 
 */
int multiply_AinvB(float *Aarr,float *Barr,int m,int n,int p,float *X);

/******************************************************************************************************
 MATRIX FUNCTIONS OVER
 ******************************************************************************************************/
#endif	/* CZTMATRIX_H */
/* 
 * File:   ExpMap.h
 * Author: preeti
 *
 * Created on January 30, 2012, 4:08 PM
 */

#ifndef EXPMAP_H
#define	EXPMAP_H

#include<fitsio.h>
#include<cstring>
#include<iostream>
#include<cmath>
#include<fstream>

#include "utils.h"
#include "glog/logging.h"

#define DET_MASKSIZE 1950
#define QUAD_MASKSIZE 8100
#define NPIXELS_GAP 100
#define PIX_SIZ_MASK  20     //20 micrometer
#define CLOSE 0

#define PIXELS_PER_COL_DET 16
#define PIXELS_PER_ROW_DET 16

#define PIXELS_PER_ROW_QUAD 64
#define PIXELS_PER_COL_QUAD 64

#define NUM_DET_PER_QUADRANT 16     //Number of detectors per quadrant
#define NUM_DET_PER_QUAD_PER_ROW 4

#define MASKHEIGHT 481         //in mm
#define DETECTORWIDTH 39       //in mm
#define TOTAL_DET_MODULES 64

/**
 * Class to retrieve required data from mask file
 * Stores data from one quadrant at a time in the object
 */
/*****************************************************************************
//Is is assumed that fits file for Mask will have 4 extensions for 4 quadrants
 * The masks for 4 quadrants would be in hdu 2, 3, 4 and 5
 *****************************************************************************/

//int writeImg(char *file,float img[XSIZE][YSIZE],int m,int n);  //Just to check

class Mask{
    char maskfile[FLEN_FILENAME];                //mask file name
    float **maskquad;   //size is QUAD_MASKSIZE
public:
    Mask();                   //default constructor   
    Mask(char *maskfile);     //parametrized constructor
    Mask(const Mask &m);
    ~Mask();   
  //int getMaskDet(int detid,float mask[DET_MASKSIZE][DET_MASKSIZE]);
    int getMaskDet(int detid,float **mask);
    int getMaskQuad(int quadid);
    
};

/**
 * Function to shift the mask array by dx and dy (in mm) in a detector
 * @param mask
 * @param delta_x  (in mm)
 * @param delta_y  (in mm) 
 * @return Returns the shifted mask array in the input array itself
 */
//int shiftMaskDet(float mask[DET_MASKSIZE][DET_MASKSIZE],double delta_x,double delta_y);
int shiftMaskDet(float **mask,double delta_x,double delta_y);
/**
 * Function to get open fraction of mask for each pixel in a detector
 * @param shiftedmask
 * @param openFraction
 * @return 
 */
//int getOpenFraction(float shiftedmask[DET_MASKSIZE][DET_MASKSIZE],
//        float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_ROW_DET]);
  
int getOpenFraction(float **shiftedmask,float **openFraction);
/**
 * Function to get weight for pixels in one detector
 * @param openFraction
 * @param weight
 * @return 
 */
//int getWeightDet(float openFraction[PIXELS_PER_COL_DET][PIXELS_PER_COL_DET],
//                float weight[PIXELS_PER_COL_DET][PIXELS_PER_COL_DET]);
//  
int getWeightDet(float **openFraction,float **weight);

/**
 * Function to get area of each pixel in 128 x 128 array and total area
 * @param area
 * @return 
 */
void getArea(float **area,double *totarea);
    

/**
 * Function to assign weight to each pixel depending on shift of mask. 
 * It also return the open fraction for each pixel
 * @param maskweight : Mask weight is returned in this buffer
 * @param openf : OPen Fraction is returned in this buffer
 * @param maskfile : maskfile path to be given
 * @param dx : shift in x direction
 * @param dy : shift in y direction
 * @return 
 */
//int getWeight(float maskweight[YSIZE][XSIZE],float openf[YSIZE][XSIZE],
 //       char *maskfile,float dx,float dy);

/**
 * Function to generate shadow of known sources given its camera coordinates and flux.
 * Uses mask file for shadow generation
 * @param thetax
 * @param thetay
 * @param sourceflux
 * @param maskfile
 * @param shadow
 * @return 
 */
int cztgenshadow(double thetax,double thetay,double sourceflux,char *maskfile,
        float shadow[YSIZE][XSIZE]);

/**
 * Function to generate exposure map for a source direction by creating array of transmission fraction
 * @param theta_x
 * @param theta_y
 * @param maskfile
 * @param emap
 * @return 
 */
int getExpMap(double theta_x,double theta_y,char *maskfile,float emap[XSIZE*YSIZE]);


int getWeight(float **maskweight,float **openf,char *maskfile,float dx,float dy);


#endif	/* EXPMAP_H */


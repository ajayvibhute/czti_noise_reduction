/* 
 * File:   ExpMap.h
 * Author: preeti, Tanul Gupta
 *
 * Created on January 30, 2012, 4:08 PM
 * Modified on November 14, 2014
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
#include "coordinateTransformation.h"

/*******
#define DET_MASKSIZE 1950

#define COMP_FIELD_SIZE 25

#define COMPCOLS 327
#define COMPROWS 8175

#define NPIXELS_GAP 125
#define PIX_SIZ_MASK  20     //20 micrometer
#define CLOSE 0

#define PIXELS_PER_COL_DET 16
#define PIXELS_PER_ROW_DET 16

#define PIXELS_PER_ROW_QUAD 64
#define PIXELS_PER_COL_QUAD 64

#define NUM_DET_PER_QUAD 16     //Number of detectors per quadrant
#define NUM_DET_PER_QUAD_PER_ROW 4

#define MASKHEIGHT 477         //in mm
#define DETECTORWIDTH 39       //in mm
#define TOTAL_DET_MODULES 64

#define BIG_PIXEL_WIDTH 2.46 // in mm
#define SMALL_PIXEL_WIDTH 2.31 // in mm

#define AVG_PIXEL_WIDTH 2.4375  //in mm
 * *************/


/**
 * Class to retrieve required data from mask file
 * Stores data from one quadrant at a time in the object
 */
/*****************************************************************************
Is is assumed that fits file for Mask will have 4 extensions for 4 quadrants
 * The masks for 4 quadrants would be in hdu 2, 3, 4 and 5
 *****************************************************************************/

/******************************************************************************
 maskfile - variable to store mask size
 
 uncompress - function which takes maskpattern(storing all 8100x270 in a single row) 
 * along with the pixel number i as the input. each pixel is of integer class representing 
 * 4 bytes each i.e. 32 bits. uncompress now take one by one first thirty bits and store them
 * in a character array uncompressed field of size 1x30. Since uncompress is selecting only 1 bit from 32 bits
 * from the maskpattern[i] so the data that gets stored in uncompress_field is either 0 or 1.
 
 mask - 2D array to store mask for 4 quadrants. its size is NUMQUAD x (QUADMASKSIZExQUADMASKSIZE) i.e. 4 x (8100x8100)
 * 
 Mask - function to copy NULL in maskfile during initialization
 * 
 read - function to allocate memory and read mask file in mask array
 * 
 getDetMask - function to extract the detector mask from uncompressed mask genrated using read().
 * It takes quadrant id (0-3) and detector id (0-15) as the input and then store mask of that particular quadrant in 2D array detmask
 * which is of size 1950x1950
 * 
 ~Mask -  destructor
 * 
 write - fits file write. use to write uncompressed mask files mainly for debugging purposes.
 ******************************************************************************/
class Mask{
protected:
    char maskfile[FLEN_FILENAME];
    int uncompress(int *maskpattern,unsigned char *uncompressed_field, int i); 
public:
    unsigned char **mask;
    Mask();
    int read(char *maskfile); 
    int getDetMask(int quadid,int detid,float **detmask);
    ~Mask();
    int write(char *file);
};


/**
* Class containing functions to perform mask weighting
* Computes mask weight for each of the detector
* Computes shadow for the given direction of source
* Computes Exposure map for given direction of source
*/

/******************************************************************************************************
 thetax, thetay - variable to store angular shift of mask pattern
 dx, dy - variable to store shift of mask pattern in mm.     dx=MASKHEIGHT*tan(thetax); dy=MASKHEIGHT*tan(thetay);
 * where MASKHEIGHT=481mm.
 sourceflux -  float variable to store flux of source for shadow computation
 openf - 1D array of size 128x128 to store weight for each detector and its open fraction
 area - array of size 128x128 to store area for each detector
 
 MaskWeigting - function to initialize MaskWeighting class.
 computeArea() - function to compute area of each 128x128 pixels in sq. mm and store it in array area. Along with individual area this function
 * also computes total area and store it in variable totalarea.
 run()- it computes the area, then reads the mask file and then compute open fraction for each pixel in a quadrant and then store it in an array of size
 * [4][64x64] and also make full open fraction of size 128x128 and store it in variable openf.
 computeOpenFraction()- function to compute the open fraction of the all quadrants and store it in an array openfraction[4][64x64]
 shiftMaskDet() - function that shifts the detector mask according to dx and dy and returns the mask array in its shifted form.
 getOpenFractionDet() - function that provides open fraction for each pixel in a detector in a 2D array of size 16x16.
 getMaskWeight()- function to get maskweight of each detector and store it in an array of size 128x128. It calculates the mask weight
 * using the formula maskweight=2*(openf-1)
 * getExposureMap()-  its simply copying open fraction in the variable expmap.
 
 ******************************************************************************************************/
class MaskWeighting:public Mask{
private:
    double thetax,thetay;        //angular shift of mask patterm
    double dx,dy;            //shift of mask pattern in mm
    float sourceflux;                //flux of source for shadow computation
    float *openf;        //size is XSIZE*YSIZE , arrays to store weight for each detector and its open fraction
    float **openf_quad; //size is [4][64*64], array to store weight of each detector pixel quadrantwise.
    void computeArea();
    int computeOpenFraction();
    int shiftMask();    // shifts mask for all 4 quadrants by dx and dy
    int shiftMaskDet(float **detmask);    //to shift mask of single detector  
    int getOpenFractionDet(float **shifteddetmask,float **openfractiondet);
    
public:
    //MaskWeighting();                   //default constructor   
    MaskWeighting(char *maskfile,double thetax,double thetay,float flux=1.0);     //parametrized constructor
    ~MaskWeighting();   
 
    int run();
    int getShadow(float *shadow,int size=XSIZE*YSIZE);                 //size is XSIZE*YSIZE for shadow
    int getExposureMap(float *expmap,int size=XSIZE*YSIZE);        //size is XSIZE*YSIZE for expmap   Open Fraction and Exposure map are same
    int getMaskWeight(double *maskweight,int size=XSIZE*YSIZE);
    int getMaskWeightQuadWise(double **mskwt_quad, int size=PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD);
    int getExposureMapQuadWise(float **openfraction_quad);
     
    float *area;                              //size is XSIZE*YSIZE, to store area of each 64 detector
    double totalarea;                 //total area of all 64 detectors
    
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


int allToQuad(float *area, float **area_quad);

/*****************************************************************************************************************
 * Description: This function finds out the RA, DEC, FLUX, thetax and thetay (camera coordinates) of known sources in the field of view.
 * Input Required:
 * 1. Catalog filename
 * 2. Aspect filename
 * 3. Extension number of catalog file to be read.
 * 4. RA column name
 * 5. DEC column name
 * 6. FLUX column name
 * 
 * Output generated;
 * 1. RA vector: vector containing RA of known sources in the FOV. (in Radians)
 * 2. DEC vector: vector containing DEC of known sources in the FOV. (in Radians)
 * 3. FLUX vector: vector containing FLUX of known sources in the FOV (obtained from Catalog)
 * 4. thetax vector: vector containing camera coordinate x of known sources in the field of view. (in Radians)
 * 5. thetay vector: vector containing camera coordinate y of known sources in the field of view. (in RAdians)
 * 
 * Description:
 * This program reads catalogfile and aspect file. It then finds shift in X and y direction of all catalog sources.
 * Sources for which shift is less than detector width are then selected and included in list of known sources.
 */
int findKnownSources(char *catalog_filename, Mkf mkf, Teldef teldef, int extnumCatalog, 
        char* RA_colname, char* DEC_colname, char* FLUX_colname,
        vector<float> &RA, vector<float> &DEC, vector <float> &FLUX, 
        vector <float> &thetax, vector <float> &thetay);


int findKnownSources(char *catalog_filename, string aspectFileName, int extnumCatalog,
        char* RA_colname, char* DEC_colname, char* FLUX_colname,
        vector<float> &RA, vector<float> &DEC, vector <float> &FLUX,
        vector <float> &thetax, vector <float> &thetay);
#endif	/* EXPMAP_H */


/* 
 * File:   cztimage.h
 * Author: preeti
 *
 * Created on October 16, 2012, 12:15 PM
 */

#ifndef CZTIMAGE_H
#define	CZTIMAGE_H

#include<pil.h>
#include<linalg.h>
#include<ap.h>
#include"utils.h"
#include"ExpMap.h"
#include "glog/logging.h"
#include "cztstring.h"

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define FFT 0
#define IFFT 1

#define NUMQUAD_ROW 2
#define NUMQUAD_COL 2

#define PEAKLIMIT 20

#define DPI 1
#define MASK 0

using namespace std;
using namespace alglib;

class cztimage{
private:
    //Input parameters
    char modulename[NAMESIZE];
    char dpifile[PIL_LINESIZE];   //Input DPI file for image creation
    char maskfile64x64[PIL_LINESIZE];
    char maskfile[PIL_LINESIZE];
    char aspectfile[PIL_LINESIZE];
    char outfile[PIL_LINESIZE];
    char quad_no[25];
    int oversamplingfactor;
    int debugmode;             //if yes/true, creates intermediate file 
    int clobber;
    int history;
    
    //data parameters
    float *dpi,*mask,*rawskyimage,*oversampledmask,*oversampleddpi;
    int size;            //size after oversampling
    int skyimgsize;    //size for skyimage
    fitsfile *fout;
     
    /**
     * read RA, DEC & TWIST from Aspect file 
     * @return 
     */
     int addWCS();
     
     /**
    * Function to read the (128 x 128) mask array from mask file and 
    * store it in the single dimnsion array mask of size 128x128.
    * @return - Returns 0 on success
    */
    int getMaskArray();
    
     /**
    * Function to read dpi for each quadrant from dpi file and 
    * store it in the single dimnsion array dpi of size 128x128.
    * @return - Returns 0 on success
    */
    int getDPI();
     
      
   //Function to oversample mask and DPI array by oversamplingfactor and stores them in oversampleddpi and oversampledmask.
    void oversample();
    
    int writeSkyImage();
    
    /**
    * Generates required number of background shadows maximum is six
    * Shadows should have size [nBkgdCoef][XSIZE*YSIZE]
    * Shadows are stored in the array 'shadows'
    * @return 
    */
 


    int getHistory(vector<string> &vhistory);
    
    /**
     * Function to subtract known source shadows from the DPI
     * @return known source subracted DPI
     */    
    
public:  
    cztimage();
    ~cztimage();
    int read(int argc,char **argv);
    int read(char *dpifile,char *maskfile64x64,char *maskfile,char *aspectfile,
        int oversamplingfactor, char *outfile, char* quad_no, int debugmode=NO,
        int clobber=YES, int history=YES);
    void display();
  
    /**
    * Function to create image from DPI
    * @param par
    * @return 
    */
    int cztimageProcess();
};



void doFft(float data[], unsigned long nn, int isign);

/**
 * Function that retains  the central (XSIZE/4)*(YSIZE/4) of image which is given in correlation matrix
 * factor gives the oversampling factor
 * The new size of image is returned in centralsize;
 * @param correlationMatrix
 * @param TOTALROWS
 * @param TOTALCOLS
 * @param factor
 * @param centralsize
 * @return 
 */
int keepCentral(float *correlationMatrix,int TOTALROWS,int TOTALCOLS,int factor,
        int *centralsize);
/**
 * Function to swap alternate wuadrants i.e. data of quadrant 0 with 3rd and data of quadrant 1st with second in any matrix. It divides a matrix of size n*n into 4 quadrants of 
 * size n/2 * n/2 and then swap them.
 * @param data
 * @param TOTALROWS
 * @param TOTALCOLS
 * @param factor
 */
void swapComplex(float *data,int TOTALROWS,int TOTALCOLS,int factor);


void generateCorrelationMatrix_2d(float *dph,float *mask,float *correlationMatrix,
        int TOTALROWS,int TOTALCOLS,int factor);

int applySigma(vector<float> &intensity,vector<int> &row,vector<int> &col,float threshold);

/**
 * Function to remove contribution of known sources from DPI
 * @param knownSourcesShadowfile
 * @param dpi
 * @param rows
 * @param XSIZE
 * @return 
 */
//int removeKnownSources(char *knownSourcesShadowfile,float *dpi,int rows,int XSIZE);

//svd related



#endif	/* CZTIMAGE_H */


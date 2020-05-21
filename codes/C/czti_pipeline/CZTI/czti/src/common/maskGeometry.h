/* 
 * @file  maskGeometry.h
 * @author Tanul Gupta
 * @date Created on July 21, 2015, 3:44 PM
 * @brief This class can read and store both compressed mask and uncompressed 64x64 mask.
 * @details This class can read and store both compressed mask and uncompressed 64x64 mask.
 * @version 1.0
 */

#ifndef MASKGEOMETRY_H
#define MASKGEOMETRY_H

#include <iostream>
#include <vector>
#include <string>
#include "glog/logging.h"
#include "fitsio.h"
#include "utils.h"
#include "detectorGeometry.h"
#include "Exposure.h"




class TanMaskGeometry {
private:
    string mask64x64Filename;
    string compMaskFilename; //name of compressed mask file
    
    vector <float> xmask0; //x coordinates of 18197 [8175+8175+1847] sub-pixels
    vector <float> ymask0; //y coordinates of 18197 [8175+8175+1847] sub-pixels
    /* To store 64x64 mask pattern */
    vector < vector <unsigned char> > q0Mask;
    vector < vector <unsigned char> > q1Mask;
    vector < vector <unsigned char> > q2Mask;
    vector < vector <unsigned char> > q3Mask;
    vector < vector <unsigned char> > fullMask;
    
    /* To store 8175x8175 uncompressed mask */
    vector < vector <unsigned char> > q0UnMask;    
    vector < vector <unsigned char> > q1UnMask;    
    vector < vector <unsigned char> > q2UnMask;    
    vector < vector <unsigned char> > q3UnMask;    
    vector < vector <unsigned char> > fullUnMask;
    
    int uncompress(unsigned int *maskpattern, unsigned char *uncompressed_field, int index);

public:
    TanMaskGeometry();
    /**
     * Reads and stores 64x64 mask pattern.
     * @param maskFilename
     * @return 
     */
    
    /**
     * Defines locations (lower left corners) of the mask sub-pixels (sampled at 0.002 cm).
     * All lengths are in cm.
     * @return 
     */
    int init_full_finer_mask();
    
    int read_64x64_mask(string maskFilename);
    /**
     * Reads finer mask from compressed mask file in CALDB. It uncompresses it and stores in
     * class variables.
     * @param compMaskFilename
     * @return 
     */
    int read_compressed_mask(string compMaskFilename);
    
    /**
     * Generates full 64x64 i.e. 128x128 mask image by combining individual quadrant masks.
     * @return 
     */
    int get_full64x64_mask();
    /**
     * Generates full finer mask image by combining individual quadrant masks.
     * @return 
     */
    int get_full_finer_mask();
    /**
     * Gets finer detector mask for a particular quadrant and specific detector Number.
     * This is how detector is visualized in this particular function. 
     * ---------------------
     * | 12 | 13 | 14 | 15 |
     * ---------------------
     * | 8  | 9  | 10 | 11 |
     * ---------------------
     * | 4  | 5  | 6  | 7  |
     * ---------------------
     * | 0  | 1  | 2  | 3  |
     * ---------------------
     * @param quadid
     * @param tanMask
     * @param detNo
     * @param detMask
     * @return 
     */
     friend vector < vector <unsigned char> > get_detector_mask(int quadid, int detNo, TanMaskGeometry tanMask, int &status);
    
    /**
     * Gets shifted finer detector mask based on shifts in X and Y direction.
     * @param quadid
     * @param detNo
     * @param TanMaskGeometry
     * @param dX : shift in X direction in cm
     * @param dY : shift in Y direction in cm
     * @param status 
     * @return 
     */
    friend vector < vector <unsigned char> > get_shifted_detector_mask(int quadid, int detNo,
                                    TanMaskGeometry tanMask, float dX, float dY, int &status);
    //GETTERS
    vector <float> get_xmask0(){return xmask0;}
    vector <float> get_ymask0(){return ymask0;}
    vector <vector <unsigned char> > get_fullMask(){ return fullMask;} 
    vector <vector <unsigned char> > get_q0Mask(){ return q0Mask;} 
    vector <vector <unsigned char> > get_q1Mask(){ return q1Mask;} 
    vector <vector <unsigned char> > get_q2Mask(){ return q2Mask;} 
    vector <vector <unsigned char> > get_q3Mask(){ return q3Mask;} 
    vector <vector <unsigned char> > get_fullUnMask(){ return fullUnMask;}
    vector <vector <unsigned char> > get_q0UnMask(){ return q0UnMask;}
    vector <vector <unsigned char> > get_q1UnMask(){ return q1UnMask;}
    vector <vector <unsigned char> > get_q2UnMask(){ return q2UnMask;}
    vector <vector <unsigned char> > get_q3UnMask(){ return q3UnMask;}
    
    //SETTERS
    void reset_full_uncomressed_mask(){fullUnMask.clear(); q0UnMask.clear();
                                            q1UnMask.clear(); q2UnMask.clear();
                                            q3UnMask.clear();}
    
    
};

//Mask weighting by using the principle of shifting. 
//NOTE: used in earlier codes.
class TanMaskWeighting {
private:
    double thetaX, thetaY;
    double dX, dY;
    float sourceflux;
    vector < vector <float> > openFractionDet; 
    vector < vector <float> > openFractionQ0; 
    vector < vector <float> > openFractionQ1; 
    vector < vector <float> > openFractionQ2; 
    vector < vector <float> > openFractionQ3; 
    vector < vector <float> > openFractionFull; 
public:
    TanMaskWeighting(double thetaX=0.0, double thetaY=0.0, float flux=1.0);
    
    
    /**
     * Computes open fraction for all quadrants
     * @param tanMask
     * @return 
     */
    int compute_open_fraction(TanMaskGeometry tanMask);
    /**
     * Computes open fraction for entire quadrant
     * @return 
     */
    int compute_open_fraction_quad(int quadNo, TanMaskGeometry tanMask);
    /**
     * Calculates open fraction of all 16x16 pixels of a detector 
     *          open fraction = 1 - close_count/total_pixels;
     * @param detMask: detector Mask for which open fraction needs to be calculated
     * @return 
     */
    int calculate_open_fraction_det(vector < vector <unsigned char> > detMask );
    
    //GETTERS
    vector < vector <float> > get_openfraction_det(){return openFractionDet;}
    vector < vector <float> > get_openfractionQ0(){return openFractionQ0;}
};

struct MaskOpenFractionTable{
    vector <int> detX; //detX coordinate (in cm)
    vector <int> detY; //detY coordinate (in cm)
    vector <int> pixX; //pixX coordinate (in cm)
    vector <int> pixY; //pixY coordinate (in cm)
    vector <float> openfrac;    //open fraction assuming blocking means no transmission at all
    vector <float> oftrans;     //open fraction with transmission considered
    
    long get_nrows(){return detX.size();}
    int pushback_openfraction(int detx, int dety, int pixx, int pixy, 
            float openfraction, float oftransmission);
    int reset();
};

class TransMask {
private:
    double thetaX, thetaY; //in degrees
    double dx, dy;
    vector <vector <unsigned char> > fullUnMask;
    vector < vector <float> > openFractionFull; //128x128 open fraction as a result of mask transmission
    vector < vector <float> > oftransFull; // 128x128 open fraction taking into account transmission
                                           // from closed mask pixels.
    MaskOpenFractionTable oftable;
    
public:
    TransMask(double thetaX=0.0, double thetaY=0.0);
    /**
     * Computes mask openfraction by dividing oversampling each detector pixel.
     * @param detgeom
     * @param energy: in keV
     * @param oversampling
     * 
     * @return 
     */
    int compute_mask_openfraction(CztDetectorGeometry detgeom,int oversampling, float energy=10.0);
    /**
     * Computes whether mass pixel at x,y is open or not
     * @param x
     * @param y
     * @param pixstate: CLOSE:0 OPEN:1
     * @return 
     */
    int transmask(float x, float y, int &pixstate);
    
    int generate_full_openfraction_image();
    int generate_full_oftrans_image();
    
    //SETTERS
    void set_fullUnMask(vector <vector <unsigned char> > fullUnMask){this->fullUnMask = fullUnMask;}
    void set_thetaXY(double tx, double ty){this->thetaX=tx; this->thetaY=ty;}
    //GETTERS
    vector < vector <float> > get_full_openfraction(){return openFractionFull;} 
    vector < vector <float> > get_full_oftrans(){return oftransFull;} 
    
};


//Independent functions

/**
 * calculates shift of mask pattern based on source direction.
 * @param thetaX: in degrees
 * @param thetaY: in degrees
 * @param dx: in cm
 * @param dy: in cm
 * @return 
 */
template <class T> int compute_shift(T thetaX, T thetaY, T &dx, T &dy){
    int status=0;
    float radThetaX = thetaX * M_PI/180;
    float radThetaY = thetaY * M_PI/180;
    
    dx = MASKHEIGHT*tan(radThetaX);
    dy = MASKHEIGHT*tan(radThetaY);
    return EXIT_SUCCESS;
}

#endif /* MASKGEOMETRY_H */
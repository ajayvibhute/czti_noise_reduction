/* 
 * @file  Exposure.h
 * @author Tanul Gupta
 * @date Created on August 23, 2015, 3:50 PM
 * @brief Exposure and Shadow Calculation
 * @details This class contains function to calculate exposure and shadow.
 * @version 0.2.0
 */
#ifndef EXPOSURE_H
#define	EXPOSURE_H

//Include section
#include "glog/logging.h"
#include "utils.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include <fitsio.h>
#include <exception>
#include "level2validation.h"
#include "detectorGeometry.h"
#include "maskGeometry.h"
#include "caldbHandler.h"
#include "badpixCALDB.h"

//Macro section


//Constants
#define CALDBVERSION "v01"

using namespace std;



class Shadow{
private:
    ShadowTable shadowTab;
    float thetaXr; //source angle x in radians
    float thetaYr; //source angle y in radians
    string cameraGeomfile; //camera geometry filename
    float sourceFlux; //flux of source in counts/area
    int includeMask; //set it to 1 for including mask, 0 mask is not included
                     //in exposure calculation
    int maskOversampling; //mask oversampling factor
    float energyStart; //in keV
    float energyEnd; //in keV
    int nBins; //number of energy bins required for energy calculation
    
public:
    Shadow(string cameraGeomFile="", float sourceFlux=1.0, float thetaXr=0.0, float thetaYr=0.0,
        int nBins=10, float energyStart=10.0, float energyEnd=100.0, 
        int includeMask=1, int maskOversampling=10);

    int compute_shadow(vector <vector <float> > &fullnormalizedEffarea,
            vector <vector <unsigned char> > &fullfinermask);
    
    //Parameter validation
    int check_shadow_par();
    //Display
    int display();
    
    //Getters
    ShadowTable get_shadow_table(){return shadowTab;}
    
    //Setters
    void set_thetaXYr(float thetaXr, float thetaYr){this->thetaXr=thetaXr; this->thetaYr=thetaYr;};
    void set_camerageomfile(string cameraGeomfile){this->cameraGeomfile=cameraGeomfile;}
    void set_sourceFlux(float sourceFlux){this->sourceFlux=sourceFlux;}
    void set_maskParameter(int includeMask, int maskOversampling=10){this->includeMask=includeMask; this->maskOversampling=maskOversampling;}
    void set_energyrange(int nBins, float energyStart, float energyEnd){this->nBins=nBins; this->energyStart=energyStart; this->energyEnd=energyEnd;}
    void set_energyvalue(int nBins, float energyStart){this->nBins=nBins, this->energyStart=energyStart;}
    void reset_shadowtab(){shadowTab.reset();}
};

//INDEPENDENT FUNCTIONS

/**
 * Returns absorption coefficient of Aluminium in units of cm^(-1)
 * Polynomial fits in logarithmic space is used to interpolate opacity
 * including coherent and incoherent scattering.
 * 
 * Range of Validity: 1-1000 keV
 * @param ekev: Energy in keV
 * @return 
 */
float absco_Al(float ekev);
/**
 * Returns absorption coefficient of Tantalum in units of cm^(-1)
 * Polynomial fits in logarithmic space is used to interpolate opacity
 * including coherent and incoherent scattering.
 * 
 * Range of Validity: 2-1000 keV
 * @param ekev: Energy in keV
 * @return 
 */
float absco_Ta(float ekev);
/**
 * Returns absorption coefficient of Cadmium Zinc Telluride in units of cm^(-1)
 * Polynomial fits in logarithmic space is used to interpolate opacity
 * including coherent and incoherent scattering.
 * 
 * Range of Validity: 5-1000 keV
 * @param ekev: Energy in keV
 * @return 
 */
float absco_CZT(float ekev);
/**
 * Returns absorption coefficient of Printec Circuit Board material on the CZTI
 * detectors used in Astrosat CZT imager in units of cm^(-1)
 * Polynomial fits in logarithmic space is used to interpolate opacity
 * including coherent and incoherent scattering.
 * 
 * Range of Validity: 1-100 MeV
 * @param ekev: Energy in keV
 * @return 
 */
float absco_pcb(float ekev);

/**
 * Returns effective absorption coefficient of material by taking average of absorption
 * coefficients calculated at interval of dekev = (maxEnergy-minEnergy/nIter).
 * For ex. if minEnergy=5, maxEnergy=95 and nIter=10, then absco are evaluated at 
 * 5,15,25,...95. 
 * @param material: "Al" | "Tan" | "Czt" | "Pcb"
 * @param minEnergy: 
 * @param maxEnergy
 * @param nIter: number of distinct coefficient for which average has to be calculated
 * @param status
 * @return 
 */
float effAbsco(string material, float minEnergy, float maxEnergy, int nIter, int &status);

/**
 * Creates exposure array in the energy range min energy to max energy
 * @param thetax (in degrees) camera coordinate x
 * @param thetay (in degrees) camera coordinate y
 * @param cameraGeomFilename 
 * @param fullfinermask
 * @param expTable
 * @param includeMask
 * @param maskOversampling
 * @param eboundsFilename: if present, then minEnergy and maxEnergy provided by user will be ignored.
 * @param minEnergy: minimum energy (keV) NOTE: cannot be less than 1 keV
 * @param maxEnergy: maximum energy (keV) NOTE: cannot be greater than 300 keV
 * @return 
 */
int create_exposure_array(float thetaxd, float thetayd, string cameraGeomFilename, 
        vector <vector <unsigned char> > &fullfinermask, ExposureTable &expTable, int includeMask, int maskOversampling, 
        string eboundsFilename, float minEnergy=1.0, float maxEnergy=301.0, int nchannels=300);
/**
 * @param ekev: Energy in keV
 * @param thetax: source x direction (in degrees)
 * @param thetay: source y direction (in degrees)
 * @param cameraGeomFileName: Name of CALDB camera geometry file
 * @param expTable: exposure Table
 * @return 
 */
int calculate_pixexposure(float ekev, float thetax, float thetay, 
        string cameraGeomFileName, vector <vector <float> > &ofmaskFull,
        int includeMask, ExposureTable &expTable);

/**
 * Calculates average exposure of all pixels in a given energy range
 * @param thetax: in degrees
 * @param thetay: in degrees
 * @param cameraGeomFilename
 * @param compMaskFile
 * @param expTable
 * @param includeMask
 * @param maskOversampling
 * @param nbins
 * @param minEnergy
 * @param maxEnergy
 * @return 
 */
int calculate_average_pixexposure(float thetax, float thetay, string cameraGeomFilename, 
        vector <vector <unsigned char> > &fullfinermask, ExposureTable &expTable, int includeMask, int maskOversampling, 
        int nbins=1,float minEnergy=10.0, float maxEnergy=100.0);


int calculate_renormalized_weights(ExposureTable &expTable, string badpixFilename, int badpixThreshold,
        vector <int> quadsToProcess, string effareaFilename, bool includeCztPcbAbsorption,string evtfilename);

#endif	/* EXPOSURE_H */

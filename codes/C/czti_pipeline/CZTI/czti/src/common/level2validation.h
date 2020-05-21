/* 
 * @file  level2validation.h
 * @author Tanul Gupta
 * @date Created on May 8, 2015, 1:03 PM
 * @brief Verifies formatting of level-2 data products
 * @details This module contains classes and functions to verify level2data.
 * @version 0.2.0
*/


#ifndef LEVEL2VALIDATION_H
#define	LEVEL2VALIDATION_H

#include "glog/logging.h"
#include "caldbHandler.h"
#include "badpixCALDB.h"
#include "l1evtdecode.h"
#include "utils.h"
#include "cztHeaderParam.h"
#include "level1handler.h"
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <numeric>
#include "errorHandler.h"

struct MKFrecord; //declared in level1handler

struct FitsInfo{
    string hduName;
    int hduType;
    long rows;
    int cols;
};
class FitsFileInfo {
private:
    int numHDU;
    map<string, FitsInfo> fitsFileMap;
public:
    /**
     * This function performs basic validation of fits file. It checks the number
     * of HDUs, HDUnames, nrows and cols in the fits file are as specified by the 
     * user.
     * If not, it will generate a warning message.
     * @param numHDU: number of HDUs expected in input fits file.
     * @param hduNames: name of HDUs expected in input fits file.
     * @param nrows: [OPTIONAL] number of rows expected in fits file.
     * @param ncols: [OPTIONAL] number of columns expected in fits file.
     * @return 
     */
    int validate_fitsfile(int numHDUs, vector<string>hduNames, vector<int> hduTypes, int* status, vector<int> ncols, vector<long> nrows);
    int get_fitsinfo(fitsfile *fptr);
    
    /**
     * Function to get information stored in fitsFileMap
     * 
     * @return 
     */
    

    int display_fitsinfo();
    
    // GETTER FUNCTIONS
    /**
     * Function to get information stored in fitsFileMap
     * @param key: [INPUT] key(hduname) value for which information is to be extracted
     * @param hdutype: [OUTPUT] stores the hdutype of the corresponding fits extension
     * @param ncols: [OUTPUT] number of columns in the hdu
     * @param nrows: [OUTPUT] number of rows in the hdu
     * @return 
     */
    int get_hdu_info(string key, int &hdutype, int &ncols, long &nrows);
};

struct eventFileQuad{
    double time;
    double localTime;
    unsigned char detID;
    unsigned char pixID;
    unsigned short veto;
    unsigned char alpha;
    unsigned char detX;
    unsigned char detY;
    unsigned short cztNtick;
    float energy;
    unsigned short PHA;
    unsigned short PI;
};

class EventFileHandler{
private:
    long nrowsQuad;
    string eventFilename;
    vector<double> vecTimeTemp;
    vector<double> vecLocalTimeTemp;
    vector<double> vecTimeQuad;
    vector<double> vecLocalTimeQuad;
    vector<unsigned short> vecPHAQuad;
    vector<unsigned short> vecPIQuad;
    vector<unsigned char> vecDetIDQuad;
    vector<unsigned char> vecPixIDQuad;
    vector<unsigned short> vecVetoQuad;
    vector<unsigned char> vecAlphaQuad;
    vector<unsigned char> vecDetXQuad;
    vector<unsigned char> vecDetYQuad;
    vector<unsigned short> vecCztntickQuad;
    vector<float> vecEnergyQuad;
    vector<vector <float> > vecTempQuadrant;

	long temp_start_index;
	     
    //Keywords
    long tstarti;
    long tstopi;
    double tstartf;
    double tstopf;
    
public:
    /**
     * Creates empty level-2 event file using event template
     * @param l2eventFilename
     * @param eventTemplate
     * @param l1eventFilename
     * @return 
     */
    int create_l2event_file(string l2eventFilename, string eventTemplate, string l1eventFilename);
    
    /**
     * Extracts and writes level-1 data 
     * @param l2eventFilename
     * @param l2bunchFilename
     * @param l1eventFilename
     * @return 
     */
    int write_l2event_file(string l2eventFilename,string bunchfile, L1evtData l1data);
    
    /**
     * Creates empty level-2 header file
     * @param l2hdrFilename
     * @param hdrTemplate
     * @param l1eventFilename
     * @return 
     */
    int create_l2hdr_file(string l2hdrFilename, string hdrTemplate, string l1eventFilename);
    
    /**
     * Extracts and writes level-1 header data.
     * @param l2hdrFilename
     * @param l1data
     * @return 
     */
    int write_l2_hdr_file(string l2hdrFilename, L1evtData l1data);
    /**
     * Function to read Temp extension of event file (any kind). This function reads
     * temperature in 2D temperature vector which contains vectors of size 16 representing
     * temperature values for 16 detectors in that quadrant.
     * @param fptr: file pointer to event file
     * @param quad_no: quadrant for which temperature values have to be read.
     * @return 
     */
    int create_l2_gti_file(string l2gtiFilename, string gtiTemplate, string l1eventFilename);
    
    int write_l2_gti_file(string l2gtiFilename, vector <double> tstart, vector <double> tstop);
    
    int create_l2_bti_file(string l2btiFilename, string btiTemplate, string l1eventFilename);

    int write_l2_bti_file(string l2btiFilename, vector <double> tstart, vector <double> tstop);
  
    /**
     * New function to copy 5 GTI extensions to the 
     * event file
     */
    int copy_gti_extensions(string l2eventFilename, string l1scienceFilename);

    int read_event_temp(fitsfile *fptr, int quad_no);
    /**
     * Function to read Q_ extension of event file (any kind). This function reads values of
     * all columns in event file and store them in class variables. 
     * @param fptr: file pointer to event file
     * @param quad_no: quadrant number to decide which extension to read.
     * @return 
     */
    int read_quad_extension(fitsfile *fptr, int quad_no);
    
    int read_quad_extension(fitsfile *fptr, int quad_no, string colname, long startRowno, long endRowno);

    /**
     * returns tstart & tstop of a particular quadrant
     * @param quadNo
     * @return 
     */
    int find_tstart_tstop(fitsfile *fptr, int quadNo, double *tstart, 
                        double *tstop, string l2evtFilename="");

    /**
     * returns minimum start and maximum stop time for the quadrants specified by quadsToProcess
     * @param fptr
     * @param quadsToProcess
     * @param tmin
     * @param tmax
     * @param l2evtFilename
     * @return 
     */
    int get_min_max_time(fitsfile *fptr, vector <int> quadsToProcess, 
                        double *tmin, double *tmax, string l2evtFilename="");
    
    /**
     * Returns common tstart and tstop time for all quadrants specified by quadsToProcess
     * @param fptr
     * @param quadsToProcess
     * @param tmin
     * @param tmax
     * @param l2evtFilename
     * @return 
     */
    int get_common_time(fitsfile *fptr, vector <int> quadsToProcess, double *tmin, double *tmax,
                        string l2evtFilename="");
    /**
     * returns true if event lies in time range specified by tstart and tstop & in
     * energy range specified by estart and estop. This event data is then stored in
     * evtquad.
     * @param recordNo
     * @param tstart
     * @param tstop
     * @param estart
     * @param estop
     * @param evtquad
     * @return 
     */
    bool event_satisfy_criteria(int recordNo, vector <double> tstart, 
                                    vector <double> tstop, vector <float> estart,
                                    vector <float> estop, eventFileQuad* evtquad,int tfilter,int efilter);
    
    /**
     * returns 1: column is empty
     *         0: column is not empty
     *        -1: column does not exist
     * @param colname
     * @return 
     */
    int is_column_empty(string colname);
    //Getters
    long get_nrowsQuad(){return nrowsQuad;}
    vector<double> get_vecTimeTemp(){return vecTimeTemp;}
    vector<double> get_vecLocalTimeTemp(){return vecLocalTimeTemp;}
    vector<double> get_vecTimeQuad(){return vecTimeQuad;}
    vector<double> get_vecLocalTimeQuad(){return vecLocalTimeQuad;}
    vector<unsigned short> get_vecPHAQuad(){return vecPHAQuad;}
    vector<unsigned char> get_vecDetIDQuad(){return vecDetIDQuad;}
    vector<unsigned char> get_vecPixIDQuad(){return vecPixIDQuad;}
    vector<unsigned short> get_vecVetoQuad(){return vecVetoQuad;}
    vector<unsigned char> get_vecDetXQuad(){return vecDetXQuad;}
    vector<unsigned char> get_vecDetYQuad(){return vecDetYQuad;}
    vector<unsigned short> get_vecCztntickQuad(){return vecCztntickQuad;}
    vector<vector <float> > get_vecTempQuadrant(){return vecTempQuadrant;}
    vector <float> get_vecEnergyQuad(){return vecEnergyQuad;}
    vector<float> get_quadrant_temperature(double eventTime, int &status);
    
};

/*********************EXPOSURE MAP*****************************************/
struct ExposureTable{
    float thetaxd; //thetax in degrees
    float thetayd; //thetay in degrees
    string energyRange; //energy range (in keV)
    int nchannels; //number of channels
    vector <unsigned char> detX; //detX coordinate 
    vector <unsigned char> detY; //detY coordinate
    vector <unsigned char> pixX; //pixX coordinate
    vector <unsigned char> pixY; //pixY coordinate
    vector <unsigned char> locX; //actual location X on 128x128 array
    vector <unsigned char> locY; //actual location Y on 128x128 array
    vector <float> openfrac;    //open fraction for single energy or average open fraction for an energy range
    vector <float> weights;     //weights for all pixels at a single energy value
    vector < vector <float> > fullExposure;
    vector < vector <float> > openfracArray; // open fraction for all pixels in multiple channels
                                             // It will be used to store open fraction array of both response 
                                             // and spectrum channels.
    vector < vector <float> > weightsArray;  // weights array for all pixels in multiple channels
    vector <unsigned short> PIs;             // PI bins (0-511 if exposure table is generated for Ebounds data)
    vector <float> energies;                 // energies for which openfracArray is generated.
    vector <unsigned char> badpixFlag;        // Badpixel flag of pixel as read from badpixel file
    vector <float> renormOffsetArray;        // renormalization offset values at multiple channels
    vector <float> areaArray;                // Area of detector at different energy values.
    vector < vector <long> > indexTable;     // Index table that stores (locX, locY) => recordNo 
                                             // of corresponding pixel in openfracArray, weightsArray etc.
    bool eboundFlag; //set to true if generated for ebounds energy range.
    
    ExposureTable();
    int generate_full_exposure_image();
    /**
     * Generates index table to locate pixels faster in the table.
     * It is generated only if earlier not available.
     * For optimization.
     * @regenerateFlag if true, then table is generated again.
     */
    int generate_index_table(bool regenerateFlag=false);
    
    long get_nrows(){return detX.size();}
    int set_exposure(unsigned char detx, unsigned char dety, unsigned char pixx, unsigned char pixy, float openfraction, float weight);

    void set_openfrac(vector<float> openfraction){this->openfrac= openfraction;}
    void set_openFracArray(vector <vector <float> > openfracArray){this->openfracArray = openfracArray;}
    void set_weightsArray(vector <vector <float> > weightsArray){this->weightsArray = weightsArray;}
    
    //Getters
    float get_weight(unsigned char locx, unsigned char locy, unsigned short PIchannel);
    //Reset vector sizes to 0 for all exposure table vectors
    int reset();
};

class ExposuremapFileHandler{
private:
    long nrows; //number of rows in exposure table
    ExposureTable expTab;
public:
    ExposuremapFileHandler();
    int read_exposure_map_file(string expMapFilename);
    /**
     * Writes data in exposure map file.
     * @param expMapFilename
     * @param expMapTemplate
     * @return 
     */
    
    int write_exposure_map_file(string expMapFilename, string expMapTemplate);

    int read_exposure_array_file(string expArrFilename);
    
    int write_exposure_array_file(string expArrFilename, string expArrTemplate);
    //SETTERS
    int set_exposuretable(ExposureTable expTable);
    //GETTERS
    vector <unsigned char> get_vecDetX(){return expTab.detX;}
    vector <unsigned char> get_vecDetY(){return expTab.detY;}
    vector <unsigned char> get_vecPixX(){return expTab.pixX;}
    vector <unsigned char> get_vecPixY(){return expTab.pixY;}
    vector <float> get_vecOpenfrac(){return expTab.openfrac;}
    ExposureTable get_exposure_table(){return expTab;}
    

};


/*******************END EXPOSURE MAP****************************************/

/******************ASPECT FILE HANDLER *************************************/
class AspectFileHandler{
private:
    long nrows;
    float RA;
    float DEC;
    float TWIST;
    string teldefFilename;
    vector<double> vecTime;
    vector<float> vecNx;
    vector<float> vecNy;
    vector<float> vecNz;
    vector<float> vecNxt;
    vector<float> vecNyt;
    vector<float> vecNzt;
public:
    AspectFileHandler();
    /**
     * Read a single aspect file and stores following information
     * 1. RA, DEC, TWIST
     * 2. vecTime, vecNx, vecNy, vecNz, vecNxt, vecNyt, vecNzt
     * @param aspectFilename
     * @return 
     */
    int read_aspect_file(string aspectFilename);
    
    /**
     * Reads multiple aspect file and store following information
     * 1. vecTime, vecNx, vecNy, vecNz, vecNxt, vecNyt, vecNzt
     * @param aspectFilenames
     * @return 
     */
    int read_multiple_aspect_files(vector <string> aspectFilenames);
    
    /**
     * Take average of both pointing and twist vectors, and provides the mean value of nx,
     * ny, nz, nxt, nyt, nzt to the user.
     * @param avg_nX
     * @param avg_nY
     * @param avg_nZ
     * @param avg_nXt
     * @param avg_nYt
     * @param avg_nZt
     * @param status
     * @return 
     */
    int get_avg_pointing_and_twist_vector(float &avg_nX, float &avg_nY, float &avg_nZ, float &avg_nXt, 
                                            float &avg_nYt, float &avg_nZt, int &status);
    long find_nrows(){ return vecTime.size();}
    
    //Getters
    float get_RA(){ return RA;}
    float get_DEC(){ return DEC;}
    float get_TWIST(){ return TWIST;}
    string get_teldefFilename(){return teldefFilename;}

};
/*****************END ASPECT FILE HANDLER***************************************/

/******************DPI FILE HANDLER******************************************/
class DPI{
private:
    vector <double> vecTime; //UTC time
    vector <unsigned char> vecDetID;
    vector <unsigned char> vecPixID;
    vector <unsigned char> vecDetX;
    vector <unsigned char> vecDetY;
    vector <float> vecEnergy; //in keV
    vector < vector <long> > dph2D;
    vector < vector <long> > fullDPH2D;
    vector < vector <float> > dpi2D;
    vector < vector <float> > fullDPI2D;
    float eMin;
    float eMax;
    double tStart;
    double tStop;
    long totalDPHcount;
    float totalDPIcount;
    string badpixFile;
    int badpixThreshold;
    bool badpixFlag;

public:
    
    cztHeaderParam headerParam;	
    DPI();
    int read_dpi_file(string dpiFilename, string extname);
    /**
     * Creates empty DPI file with header values copied from event file.
     * @param outDPIfilename
     * @param DPItemplate
     * @param inEvtFilename
     * @return 
     */
    int create_dpi_file(string outDPIfilename, string DPItemplate, string inEvtFilename);
    
    int write_dpih_quad(string DPIfilename, string filetype,string extname);
        
    int write_dpih_full(string DPIfilename, string filetype, string extname);

    /**
     * Creates DPH for a quadrant based on energy range and time range supplied by user.
     * @param eventFilename
     * @param quad_no
     * @param emin
     * @param emax
     * @param tstart
     * @param tstop
     * @return 
     */
    int create_dph2d(string eventFilename, int quad_no, float emin, float emax, 
        double tstart, double tstop,int timefilter,int energyfilter);
    
    /**
     * Places quadDPI/H into full DPI/H based on quadrant number
     * @param quadNo
     * @return 
     */
    int add_quad_to_fulldpih(int quadNo, string fileType);

    /**
     * Creates DPI from DPH based on effective area file and eMin & eMax
     * @param effAreaFilename
     * @param quad_no
     * @param eMin
     * @param eMax
     * @return 
     */
    int create_dpi2d(string effAreaFilename, int quad_no);
    //GETTERS
    vector <unsigned char> get_detID(){return vecDetID;}
    vector <unsigned char> get_pixID(){return vecPixID;}
    vector <unsigned char> get_detX(){return vecDetX;}
    vector <unsigned char> get_detY(){return vecDetY;}
    vector <float> get_energy(){return vecEnergy;}
    //Reads DPI image based on extension name provided by user
    vector <vector <float> > get_dpi_image(string dpiFilename, string extname);
    //Reads DPH image based on extension name provided by user
    vector <vector <long> > get_dph_image(string dphFilename, string extname);
    
    //SETTERS
    void set_badpixfile(string badpixelFile){this->badpixFile = badpixelFile;}
    void set_badpixThreshold(int badpixThreshold){this->badpixThreshold=badpixThreshold;}
    
};
/******************SHADOW FILE HANDLER***********************/
struct ShadowTable {
    vector <int> detX; //detX coordinate
    vector <int> detY; //detY coordinate
    vector <int> pixX; //pixX coordinate on detector
    vector <int> pixY; //pixY coordinate on detector
    vector <int> locx; //actual pixel x coordinate
    vector <int> locy; //actual pixel y coordinate
    vector <float> shadow; //shadow
    vector < vector <float> > shadow2D;
    //Information parameters
    float RA;
    float DEC;
    float thetaXr; //camera coordinate in radians
    float thetaYr;
    float energyMin; //in keV
    float energyMax;
    float sourceFlux;
    int nBins; //number of energy bins

    long get_nrows() {
        return detX.size();
    }
    
    int is_empty();
    //SETTERS
    int set_shadow(int detx, int dety, int pixx, int pixy,
            int x, int y, float shadow);
    int reset();
};

class ShadowFileHandler{
private:
    ShadowTable shadowtab;
    string shadowFilename; //shadow file for which data is stored
    vector <float> shadowImg;
    string extname; //extension of shadow file for which data is stored in shadowtab
    bool headerReadFlag;
public:
    //ShadowFileHandler();
    /**
     * Reads a particular shadow file and store data read in shadowtab 
     * @param shadowFilename
     * @param extName
     * @return 
     */
    int read_shadow_file(string shadowFilename, string extName);
    /**
     * Creates empty shadow file based on template provided by user
     * @param shadowFilename
     * @param extName
     * @return 
     */
    int create_shadow_file(string shadowFilename, string shadowTemplate, vector<string> extNames);
    /**
     * writes shadow file
     * @return 
     */
    int write_shadow_file(string shadowFilename, string extName);
    
    //GETTERS
    vector < vector <float> > get_shadow(string quadname="QALL");
    //SETTERS
    void set_shadowtab(ShadowTable shadowtab){this->shadowtab = shadowtab;}
    int reset();
};

/**********SPECTRUM CLASS********************************/
class SpectrumLc{
public:
    vector <float> flux;   //in counts/sec
    vector <float> energy; //in keV
    vector <unsigned int> PI; //PI channel(for spectrum binning)
	vector <unsigned int> channels;
	vector <double> UT; //Universal time (for light curve binning)
    vector <float> fracexp;
    double timedel; //binsize for lc
    vector <float> error;
    float thetaxd; //thetax in degrees
    float thetayd; //thetay in degrees
	float theta;
    string energyRange; //energy range (in keV)
	string gtitype;

  	string expmapfile;
	string badpixfile;
	string eventfile;
	int badpixthresh;
	int maskweighted;
	    
    int quadrantid;
    int allquad;
	 
    float exposure_time; //Exposure time

    int nchannels; //number of channels
    
    int bin_spectrum(string eventFileName, vector<double> tstart, vector<double> tstop, 
                     vector <float> estart, vector <float> estop, vector <int> quadsToProcess,
                        ExposureTable &expTab,int tfilter,int efilter,int maskWeight);

    int bin_lc(string eventFileName, vector<double> tstart, vector<double> tstop,
            vector <float> estart, vector <float> estop, vector <int> quadsToProcess,
            ExposureTable &expTab, double binsize,int tfilter,int efilter, char *livetimefile,int maskWeight);
    
    int verify_spectrum();
    int verify_lc();
    int clear();
};
/**********SPECTRUM CLASS END***************************/
/**********Spectrum File Handler*************************/
class SpectrumFileHandler {
private:
    string spectrumFilename;
    SpectrumLc spec;

public:

    int read_spectrum_file();

    int write_spectrum_file(string specFilename, string specTemplate,string infile);

    
    //Setters
    void set_spectrum_filename(string specFilename){ this->spectrumFilename = specFilename;} 
    void set_spectrum(SpectrumLc spec){ this->spec = spec;}

};
/************Spectrum file handler end********************/

class LightCurveFileHandler {
    string lcFilename;
    SpectrumLc lc;
public:
    int write_lc_file(string lcFilename, string lcTemplate,string infile);
    //Setters
    void set_lc_filename(string lcFilename){ this->lcFilename = lcFilename;}   
    void set_lc(SpectrumLc lc){ this->lc = lc;}
};
//Independent Functions
int verify_num_hdus(fitsfile *fptr, int expectedNumHDUs, int status);

bool is_in_range();

/***************Image File Handler*******************/
class Image {
public:
    vector < vector <float> > crossCorrelationImage;
    long xsize;
    long ysize;
    string teldefFile;
    string aspectFile;
    int oversamplingfactor;
    //Keywords in DPI
    float eMin;
    float eMax;
    double tStart;
    double tStop;
    string badpixFile;
    int badpixThreshold;
    //WCS KEYWORDS
    //RA DEC
    float crpix1;
    float crpix2;
    float crval1;
    float crval2;
    float cdelt1;
    float cdelt2;
    float crota1;
    float crota2;
    float cd1_1;
    float cd1_2;
    float cd2_1;
    float cd2_2;
    //Camera coordinates
    float crpix1b;
    float crpix2b;
    float crval1b;
    float crval2b;
    float cdelt1b;
    float cdelt2b;
    float crota1b;
    float crota2b;
    float cd1_1b;
    float cd1_2b;
    float cd2_1b;
    float cd2_2b;
    

    Image();
    int is_empty();
    void evaluate_wcs_coordinates(int oversamplingFactor=1, string aspectFile="");
    void generate_cross_correlation_image(
                    vector <vector <float> > *dpi, 
                    vector <vector <float> > *mask,
                    int oversampfactor,float xshift,float yshift);
    //Getters
    long get_imgxsize(){return xsize;}
    long get_imgysize(){return ysize;}
    
};

class ImageFileHandler{
private:
    string imageFilename;
    Image imgq0;
    Image imgq1;
    Image imgq2;
    Image imgq3;
    Image imgqall;
  
    
public:
    cztHeaderParam headerParam;	
    ImageFileHandler();
    int read_image_file();
    void create_image_file(string outImageFilename, string imageTemplate);
    int write_image_file(string imgFilename, string extname, Image *img);
    
    //setters
    //void set_image(Image img){this->img=img;}
    
};

// New functions (Mithun)
int getAvgExposureTime(fitsfile *fevt,vector <int>quadsToProcess,float *exp_time);
//int create_livetime_file(char *livetimefile);Commented by Mayuri,16th Dec 2017
int create_livetime_file(char *livetimefile,double livetime_binsize);//Added by Mayuri,16th Dec 2017
int writelivetime(double *timearray,double*livetime,long num_timebins,int qid,char*livetimefile);

#endif	/* LEVEL2VALIDATION_H */


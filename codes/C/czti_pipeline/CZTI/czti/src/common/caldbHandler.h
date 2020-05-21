/* 
 * @file  caldbHandler.h
 * @author Tanul Gupta
 * @date Created on April 23, 2015, 3:50 PM
 * @brief CALDB files handler
 * @details This class contains function to read data from CALDB files.
 * @version 0.2.0
 */

#ifndef CALDBHANDLER_H
#define	CALDBHANDLER_H

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
#include <map>
#include "level2validation.h"
#include "errorHandler.h"
//Macro section


//Constants
#define CALDBVERSION "v01"

using namespace std;

class CaldbHandler{
private:
    string basepath;
    map<string,string> CaldbFiles;
public:
    void generate_caldb_filepaths(string basepath);
    //GETTERS
    string get_caldb_file_path(string filetype);
    map<string, string> get_Caldb_files(){return CaldbFiles;}
};
class Teldef {
private:
    string fileType; // Type of calibration data
    string coord0; // Name of First coordinate system
    string coord1; // Name of Second coordinate system
    string coord2; // Name of Third coordinate system

    string rawXcol; // Name of raw X column in event files
    float rawXsize; // RAW address space X size (pixels)
    float rawXpix1; // RAW address space X first pixel number (pixel)
    float rawXscl; // RAW X scale
    string rawYcol; // Name of raw Y column in event files
    float rawYsize; // RAW address space Y size (pixels)
    float rawYpix1; // RAW address space Y first pixel number (pixel)
    float rawYscl; // RAW Y scale
    string rawUnit; // Physical unit of RAW coordinates
    
    string detXcol; // Name of detector X column in event files
    float detXsize; // DET address space X size (pixels)
    float detXpix1; // DET address space X first pixel number (pixel)
    float detXscl; // DET X scale
    string detYcol; // Name of detector Y column in event files
    float detYsize; // DET address space Y size (pixels)
    float detYpix1; // DET address space Y first pixel number (pixel)
    float detYscl; // DET Y scale
    string detUnit; // Physical unit of detector coordinates

    float coeXA; // Coefficients to translate RAW into intermediate
    float coeXB; // coordinate systems; The new values so obtained
    float coeXC; // are used to transform them into DET coordiante
    float coeYA; // systems.
    float coeYB;
    float coeYC;

    float detXoff; // X offset between intermediate and DET coords
    float detYoff; // Y offset between intermediate and DET coords
    float detXflip; // flip x-axis in RAW -> DET 0:No 1:Yes
    float detYflip;
    float detScal; // detector scaling
    float detRotd; // detector rotation


    string skyXcol; // Name of SKY X column in event files
    float skyXsize; // SKY address space X size (pixels)
    float skyXpix1; // SKY address space X first pixel number (pixel)
    string skyYcol; // Name of SKY Y column in event files
    float skyYsize; // SKY address space Y size (pixels)
    float skyYpix1; // SKY address space Y first pixel number (pixel)
    string skyUnit; // physical unit of SKY coordinates
    string skyFrom; // coordinates from which sky coordinates are calculated
    
    float alignM11; // CZTI -> SAT coordinates alignment matrix Mij
    float alignM12;
    float alignM13;
    float alignM21;
    float alignM22;
    float alignM23;
    float alignM31;
    float alignM32;
    float alignM33;
    
    float focalLength;  // Telescope focal length (mm)
    float optAxisY;     // optical axis Y in detector coordinates
    
public:
    Teldef();
    /**
     * Function to read TELDEF CALDB FILE and store the information extracted from header
     * in class variables.
     * @param teldefFilename
     * @return 
     */
    int read_teldef_file(string teldefFilename);
    int display();
    
    //Getter functions
    float get_alignM11(){ return alignM11;}
    float get_alignM12(){ return alignM12;}
    float get_alignM13(){ return alignM13;}
    float get_alignM21(){ return alignM21;}
    float get_alignM22(){ return alignM22;}
    float get_alignM23(){ return alignM23;}
    float get_alignM31(){ return alignM31;}
    float get_alignM32(){ return alignM32;}
    float get_alignM33(){ return alignM33;}
};

class Ebounds {
private:
    vector <unsigned short> channel;   //Stores channel numbers from ebounds file
    vector <float> eMin;    //Stores minimum energy values for corresponding channel number
    vector <float> eMax;    //Stores maximum energy values for corresponding channel number
public:
    Ebounds();
    /**
     * Function to read EBOUNDS CALDB FILE and store the data extracted in vectors of the object.
     * @param eboundsFilename
     * @return status: 0:Successful else Failed
     */
    int read_ebounds_file(string eboundsFilename);
    unsigned short find_channel_from_energy(float energy, int &status);
    
    /**
     * Returns mean of emin and emax for the user specified channel number.
     * @param channelNo
     * @return 
     */
    float find_energy_from_channel(unsigned short channelNo);
    int display();
    
    //Getters
    long get_nrows(){return channel.size();}
    float get_energy(int rowno);
    unsigned int get_channel(int rowno);
};

class GainOffset {
private:
    vector <unsigned char> detID;
    vector <unsigned char> pixID;
    vector <float> gain00;
    vector <float> gain05;
    vector <float> gain10;
    vector <float> gain15;
    vector <float> gain20;
    vector <float> offset00;
    vector <float> offset05;
    vector <float> offset10;
    vector <float> offset15;
    vector <float> offset20;
public:
    GainOffset();
    /**
     * This function reads GAIN file provided by user. It then extracts and stores data
     * in class variables.
     * @param gainsFilename
     * @param extname: HDU extension name of gains file to be read.
     * @return 
     */
    int read_gainoffset_file(string gainsFilename, string extname);
    int display();
    
    // Getter Functions
    vector<unsigned char> get_detID(){return detID;}
    vector<unsigned char> get_PixID(){return pixID;}
    
    /**
     * Function that returns all gain values for a particular temperature
     * @param temperature: "T1", "T2", "T3", "T4", "T5"
     * @param status: 0: success else failure
     * @return Gain vector for that temperature
     */
    vector<float> get_gain(int temperature, int &status);
    vector<float> get_offset(int temperature, int &status);

};


class ResponsePar {
private:
    vector <unsigned char> detID;
    vector <unsigned char> pixID;
    vector <unsigned char> detX;
    vector <unsigned char> detY;
    vector <float> mRes;
    vector <float> cRes;
public:
    ResponsePar();
    int read_response_par_file(string responseParsFilename, string extname);
    int display();
    
};

class EffArea {
private:
    vector <unsigned char> detID;
    vector <unsigned char> pixID;
    vector <unsigned char> detX;
    vector <unsigned char> detY;
    vector <vector <float> > area;
    vector <vector <float> > normalizedEffArea;
    vector <vector <float> > fullNormalizedEffArea;
    float maxEffArea;
    
    /**
     * Adds effective area covering energy range eMin-eMax for the given 
     * area one-dimensional vector.
     * @param area1d
     * @param eMin
     * @param eMax
     * @return 
     */
    float add_effarea(vector <float> area1d, float eMin, float eMax);
public:
    EffArea();
    
    /**
     * Reads a particular extension of effective area file
     * @param effareaFilename
     * @param extname
     * @return 
     */
    int read_effarea_file(string effareaFilename, string extname);
    
    
    int display();
    
    /**
     * Function to calculate normalized effective area for a single quadrant
     * @param eMin
     * @param eMax
     * @return 
     */
    int calculate_normalized_effarea(float eMin, float eMax);
    
    /**
     * Function to calculate full normalized effective area i.e. for all quadrants.
     * @param eMin
     * @param eMax
     * @return 
     */
    int calculate_full_normalized_effarea(string effareaFilename,
                                        float eMin, float eMax);
    /**
     * Interpolates effective area of a particular pixel for a user specified energy value.
     * @param pixx
     * @param pixy
     * @param ekeV
     * @return 
     */
    vector<float> get_effective_area(unsigned char pixx, unsigned char pixy, vector<float> ekeV);
    
    //GETTERS
    vector < vector <float> > get_normalized_effarea(){ return normalizedEffArea;}
    vector < vector <float> > get_fullnormalized_effarea(){ return fullNormalizedEffArea;}
    float get_max_effarea(){return maxEffArea;}
    long get_nrows(){return detID.size();}
    
    //SETTERS
    void reset();
    
};

/**
 * Defines the CZTI geometry as a collection of plane surfaces. Each surface
 * has following attributes:
 * zmin, zmax : Vertical extent
 * xmin, xmax, ymin, ymax : horizontal extents in x and y direction, may depend
 *                          on z, described by a straight line equation.
 * coefs xminz, xmaxz, yminz, ymaxz
 * intercepts xminc, xmacx, yminc, ymaxc
 * xnp, ynp, znp: direction cosines of the normal to the surface
 * dnp: constant to define surface location
 * 
 * xnp*x+ ynp*y + znp*z = dnp is the equation of the surface
 * 
 * Material: Each surface is composed of either Aluminium or Tantalum or both.
 * thickAL, thickTa: effective normal thickness of the Al and Ta layers in cm.
 * 
 * The Ta layers may have holes cut into it for creating coded mask. If so,
 * then openfract quantifies the fraction of the Ta surface area that has been
 * cut out and is hence open.
 * 
 * The top coded mask plate is not included in the definitions provided here.
 * Function transmask below should be used to model transmission through the mask
 * plate.
 * 
 * Coordinate system conventions:
 *      z direction is parallel to the roll axis.
 *      x direction is normal to the radiator plate, directed from the radiator
 *      plate towards the camera.
 *      y direction is parallel to the radiator plate, directed towards the SXT.
 * 
 * A flag, "disable", attached to every surface controls whether or not it is used
 * in the computation. 
 * Disable=0: retain the surface
 * Disable=1: ignore the surface
 */
class CztiCameraGeometry {
private:
    vector <unsigned short> indexSurf;
    vector <bool> disable;
    vector <float> zmin, zmax;
    vector <float> xminz, xminc, xmaxz, xmaxc;
    vector <float> yminz, yminc, ymaxz, ymaxc;
    vector <float> xnp, ynp, znp, dnp;
    float chsm, chsc, chcos, chsin, chdcom;
    vector <float> thickAl, thickTa, openfrac;
    
public:
    CztiCameraGeometry(); 
    
    /**
     * Function to define all 62 surfaces of CZTI in terms of variables define above.
     * @return 0: success 1: failure
     */
    int define_camera_geometry();
    
    /**
     * Function to create camera geometry file from corresponding template. 
     * This file stores the camera geometry evaluated by define_camera_geometry in a FITS file.
     * @param cameraGeomFileName: Name of camera geometry file.
     * @param geomTemplate: Name of camera geometry template file.
     * @return 
     */
    static int create_camera_geometry_file(string cameraGeomFileName, string geomTemplate);
    
    /**
     * Function that writes camera geometry evaluated by define_camera_geometry in a empty camera geometry
     * fits file created by by create_camera_geometry_file.
     * @param cameraGeomFileName: Name of camera geometry file.
     * @return 
     */
    int write_camera_geometry_file(string cameraGeomFileName, string geomTemplate);
    
    /**
     * Function to read camera geometry file and store data in object variables.
     * @return 
     */
    int read_camera_geometry_file(string cameraGeomFileName);

    //GETTERS
    vector <unsigned short> get_indexSurf(){return indexSurf;}
    vector <bool> get_disable(){return disable;}
    vector <float> get_zmin(){return zmin;}
    vector <float> get_zmax(){return zmax;}
    vector <float> get_xmaxz(){return xmaxz;}
    vector <float> get_xmaxc(){return xmaxc;}
    vector <float> get_ymaxz(){return ymaxz;}
    vector <float> get_ymaxc(){return ymaxc;}
    vector <float> get_xminz(){return xminz;}
    vector <float> get_xminc(){return xminc;}
    vector <float> get_yminz(){return yminz;}
    vector <float> get_yminc(){return yminc;}
    vector <float> get_xnp(){return xnp;}
    vector <float> get_ynp(){return ynp;}
    vector <float> get_znp(){return znp;}
    vector <float> get_dnp(){return dnp;}
    vector <float> get_thickAl(){return thickAl;}
    vector <float> get_thickTa(){return thickTa;}
    vector <float> get_openfrac(){return thickAl;}
};

struct Catalog{
    vector<string> sourceNameBAT;
    vector<string> counterpartName;
    vector<float> RA; //in degrees (read from catalog)
    vector<float> DEC; //in degrees (read from catalog)
    vector<float> thetaxr; //in radians
    vector<float> thetayr; //in radians
    vector<float> SNR;
    vector<float> ctptRA;
    vector<float> ctptDEC;
    vector<float> Flux;
    vector<float> FluxLO;
    vector<float> FluxHI;
    
    int get_nsources(){return sourceNameBAT.size();}
    //Setters
    int set_catalog_values(string srcname, string cprtname, float ra, float dec,
                    float thetaxr, float thetayr, float snr, float ctptra, float ctptdec, float flux,
                    float fluxlo, float fluxhi);
    
    
};

class StarCatalog{
private:
    string catalogFilename;
    Catalog catalog;
public:
    int read_catalog_file(string catalogFilename, int catalogExtnum);
    int find_sources_in_FOV(Catalog &catalog, double RA, double DEC, double TWIST, 
                double thetaxr, double thetayr);
    //Getters
    Catalog get_catalog(){return this->catalog;};
};
/**
 * Function to generate full path of CALDB file. It takes the value of CALDB 
 * directory from environment variable "CALDBnew" and then appends "/bcf/filename"
 * to it.
 * @param filename: Name of CALDB file
 * @return 
 */
string caldb_full_path_generator(string filename="");

/**
 * Function to generate CALDB filename based on user input. This definition is 
 * in accordance with Astrosat_CALDB_convention_v0.1
 * @param content: content of caldb file (descriptor with no spaces)
 * @param version: version number
 * @param quadNo: quadrant number for which CALDB file is valid. If "" then CALDB
 *                  file is valid for all quadrants.
 * @param creationTime: time of creation of file 
 * @param extname: extension name of file
 * @return 
 */
string caldb_filename_generator(string content, string version, string creationTime, string extname="fits", string quadNo="");

/**
 * Writes checksum and datasum of each HDU of fits file & updates following CALDB header keywords
 * 1. Origin: Origin of CALDB file | ex. CZTI-POC
 * 2. Creator: Creator of CALDB file | ex. cztiCALDBcreate_v0.2
 * 3. Filename: name of file
 * 4. CVSD0001: Validity Start Date in UTC (yyyy-mm-dd)
 * 5. CVST0001: Validity Start Time in UTC (hh:mm:ss)
 * 
 * 
 * @param filename
 * @param creator
 * @return 
 */
int update_CALDB_keywords(char* filename, char *creator);



//template <class Class_T,class Array_T> 
//int create_2D_array(Class_T o, vector <vector <Array_T> > *rearranged_array){
//    int i,j =0; //counter variables
//    int status=0; //status variable
//    for (i=0; i<= o.detID.size(); i++){
//            (*rearranged_array)[(int) o.detID[i]][(int) o.pixID[i]] = o.offsetT2[i];
//    }
//    
//    return EXIT_SUCCESS;
//}


//Functions to query caldb files

int QueryCaldb(string telescope,string instrument,string detname,string codename,double tstart,double tstop,string &caldb_filename,int &extnum);

int getcif(string telescope,string instrument,string cif);

#endif	/* CALDBHANDLER_H */


/* 
 * @file  mkfRegeneration.h
 * @author Tanul Gupta
 * @date Created on October 13, 2015, 2:51 PM
 * @brief Regenerates mkf file for level-2 analysis
 * @details
 * @version 2.0
 */

#ifndef MKFREGENERATION_H
#define	MKFREGENERATION_H

#include <string>
#include <fitsio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "glog/logging.h"
#include "utils.h"
#include "errorHandler.h"
#include "level1handler.h"
#include "macrodef.h"
#include "Mvector.h"
//#include <my_global.h>
//#include <mysql.h>
#include<stdlib.h>

using namespace std;



/**************MKF_thresholds HANDLER***********************/
struct MKFthrStruct {
    string paramaeterName;
    short int parFlag;
    string minValue;
    string maxValue;
};

class MKFthreshold {
public:
    vector <MKFthrStruct> mkfthr;
    int read_mkf_thresholds(string MKFthrFilename);

    bool verify_mkf_record(MKFrecord mkfr);

    /**
     * Displays values read from mkf threshold file
     */
    void display();

    /**
     * Checks whether a particular value is in range or not.
     * @param val
     * @param minValue
     * @param maxValue
     * @return 
     */
    template <class T>
    bool is_in_range(T val, string minValue, string maxValue);
    //Getter

    unsigned int get_npar() {
        return mkfthr.size();
    }

};

template <class T>
bool MKFthreshold::is_in_range(T val, string minValue, string maxValue) {
    bool x = true;
    double minVal;
    double maxVal;
    if (minValue == "-" && maxValue != "-") {
        maxVal = atof((char*) maxValue.c_str());
        x = (val <= (T) maxVal) ? true : false;
    } else if (minValue != "-" && maxValue == "-") {
        minVal = atof((char*) minValue.c_str());
        x = (val >= (T) minVal) ? true : false;
    } else if (minValue != "-" && maxValue != "-") {
        maxVal = atof((char*) maxValue.c_str());
        minVal = atof((char*) minValue.c_str());
        x = (val >= (T) minVal && val <= (T) maxVal) ? true : false;
    } else {
        x = true;
    }
    return x;
}

/************* END MKF THRESHOLD ***********************/

struct MKFrecord {
    double time;
    vector<double>  Qsat;
    float rollRA;
    float rollDEC;
    float pitchRA;
    float pitchDEC;
    float yawDEC;
    float yawRA;
    float rollRot;
    float posX;
    float posY;
    float posZ;
    float velX;
    float velY;
    float velZ;
    float earthLAT;
    float earthLON;
    float altitude;
    float sunAngle;
    float moonAngle;
    float elv;
    float sunelev;
    float angoffset;
    float timeSinceSAA;
    unsigned long cpmRate;
    unsigned char q1modeID;
    unsigned char q2modeID;
    unsigned char q3modeID;
    unsigned char q4modeID;
    float q1Pos5vMonitor;
    float q2Pos5vMonitor;
    float q3Pos5vMonitor;
    float q4Pos5vMonitor;
    unsigned long q1CZTcounter;
    unsigned long q2CZTcounter;
    unsigned long q3CZTcounter;
    unsigned long q4CZTcounter;
    float q1Pos2_5vMonitor;
    float q2Pos2_5vMonitor;
    float q3Pos2_5vMonitor;
    float q4Pos2_5vMonitor;
    float q1CZTHVMonitor;
    float q2CZTHVMonitor;
    float q3CZTHVMonitor;
    float q4CZTHVMonitor;
    float q1VetoHVMonitor;
    float q2VetoHVMonitor;
    float q3VetoHVMonitor;
    float q4VetoHVMonitor;
    float q1DVDD;
    float q2DVDD;
    float q3DVDD;
    float q4DVDD;
    float q1VetoLLD;
    float q2VetoLLD;
    float q3VetoLLD;
    float q4VetoLLD;
    unsigned long q1VetoCounter;
    unsigned long q2VetoCounter;
    unsigned long q3VetoCounter;
    unsigned long q4VetoCounter;
    unsigned long q1AlphaCounter;
    unsigned long q2AlphaCounter;
    unsigned long q3AlphaCounter;
    unsigned long q4AlphaCounter;
    float q1temp;
    float q2temp;
    float q3temp;
    float q4temp;    
};

class Mkf{
public:
    vector <MKFrecord> mkf;
    vector <double> vec_time;
    vector < vector <double> > vec_Qsat;
    vector <float> vec_rollRA;
    vector <float> vec_rollDEC;
    vector <float> vec_pitchRA;
    vector <float> vec_pitchDEC;
    vector <float> vec_yawDEC;
    vector <float> vec_yawRA;
    vector <float> vec_rollRot;
    vector <float> vec_posX;
    vector <float> vec_posY;
    vector <float> vec_posZ;
    vector <float> vec_velX;
    vector <float> vec_velY;
    vector <float> vec_velZ;
    vector <float> vec_earthLAT;
    vector <float> vec_earthLON;
    vector <float> vec_altitude;
    vector <float> vec_sunAngle;
    vector <float> vec_moonAngle;
    vector <float> vec_elv;
    vector <float> vec_sunelev;
    vector <float> vec_angoffset;
    vector <float> vec_timeSinceSAA;
    vector <unsigned long> vec_cpmRate;

    vector <float> vec_q1Pos5vMonitor;
    vector <float> vec_q2Pos5vMonitor;
    vector <float> vec_q3Pos5vMonitor;
    vector <float> vec_q4Pos5vMonitor;
    vector <unsigned long> vec_q1CZTcounter;
    vector <unsigned long> vec_q2CZTcounter;
    vector <unsigned long> vec_q3CZTcounter;
    vector <unsigned long> vec_q4CZTcounter;
    vector <float> vec_q1Pos2_5vMonitor;
    vector <float> vec_q2Pos2_5vMonitor;
    vector <float> vec_q3Pos2_5vMonitor;
    vector <float> vec_q4Pos2_5vMonitor;
    vector <float> vec_q1CZTHVMonitor;
    vector <float> vec_q2CZTHVMonitor;
    vector <float> vec_q3CZTHVMonitor;
    vector <float> vec_q4CZTHVMonitor;
    vector <float> vec_q1VetoHVMonitor;
    vector <float> vec_q2VetoHVMonitor;
    vector <float> vec_q3VetoHVMonitor;
    vector <float> vec_q4VetoHVMonitor;
    vector <float> vec_q1DVDD;
    vector <float> vec_q2DVDD;
    vector <float> vec_q3DVDD;
    vector <float> vec_q4DVDD;
    vector <float> vec_q1VetoLLD;
    vector <float> vec_q2VetoLLD;
    vector <float> vec_q3VetoLLD;
    vector <float> vec_q4VetoLLD;
    vector <unsigned long> vec_q1VetoCounter;
    vector <unsigned long> vec_q2VetoCounter;
    vector <unsigned long> vec_q3VetoCounter;
    vector <unsigned long> vec_q4VetoCounter;
    vector <unsigned long> vec_q1AlphaCounter;
    vector <unsigned long> vec_q2AlphaCounter;
    vector <unsigned long> vec_q3AlphaCounter;
    vector <unsigned long> vec_q4AlphaCounter;
    vector <float> vec_q1temp;
    vector <float> vec_q2temp;
    vector <float> vec_q3temp;
    vector <float> vec_q4temp;
    vector <unsigned char> vec_q1modeid;
    vector <unsigned char> vec_q2modeid;
    vector <unsigned char> vec_q3modeid;
    vector <unsigned char> vec_q4modeid;
    //hdr parameters
    //added by shrikant
	float hkparam[8];
	int error_count,boot_page_no,modeid;
	long alphacounter,vetocounter,start_indexQ0,start_indexQ1,start_indexQ2,start_indexQ3;
	double RAsun,DecSun;
	


    //Mkf parameters
    long nrowsMkf; //number of records in MKF file.
    long mkfLastIndex; //index of last record read.
    //Keywords
    double tstart,tstop;
    long tstarti;
    double tstartf;
    long tstopi;
    double tstopf;
    double RAPnt; //RA pointing (degrees)
    double DECPnt; //DEC pointing (degrees)
    string obsID;
    long lastRecord;
    bool fileSorted; //if true then file is sorted with respect to time
    //in ascending order.
    bool headerReadFlag; //set to true if all header parameters have 
    //been read.
    string mkfFilename;
    bool optionalParRead;

    //Files read to regenerate mkf
    string attFilename;
    string lbtFilename;
    string orbFilename;

    Mkf();
    int write_header(string,string);
    /**
     * Reads mkf file.
     * @param mkffilename
     * @return 
     */
    int read_mkf_file(string mkffilename);
    
    /**
     * Regenerates mkf values from attitude, lbt and orbit file.
     * @param attFilename
     * @param lbtFilename
     * @param orbFilename
     * @return 
     */
    int regenerate_mkf_values(string attFilename, string lbtFilename, string orbFilename,string hdrFilename);
    
    /**
     * Writes regenerated mkf file.
     * @param mkffilename
     * @param rMkfTemplate
     * @return 
     */
    int write_regenerated_mkffile(string mkffilename, string rMkfTemplate);
    //convert class data into vectors
    int get_vectors();
    

    int read_mkf_time(fitsfile *fptr);
    
    /**
     * Reads rollRA, rollDEC, rollROT, pitchRA, pitchDEC, yawDEC & yawRA columns of MKF file.
     * @param fptr: file pointer to mkf file.
     * @return 
     */
    int read_mkf_rpy(fitsfile *fptr); //reads roll, pitch & yaw RA and Dec information.
    
    /**
     * From mkf file it calculates start and stop time of gti satisfying mkf threshold criteria.
     * @param mkfThrFilename
     * @param tstart
     * @param tstop
     * @return 
     */
    int get_filtered_gti(string mkfThrFilename, vector <double> &tstart, vector <double> &tstop);
    
    int display();
   
    
    //Getters
    vector<double> get_vecTime(){ return vec_time;}
    vector<float> get_vecRollRA(){return vec_rollRA;}
    vector<float> get_vecRollDEC(){return vec_rollDEC;}
    vector<float> get_vecPitchRA(){return vec_pitchRA;}
    vector<float> get_vecPitchDEC(){return vec_pitchDEC;}
    vector<float> get_vecYawRA(){return vec_yawRA;}
    vector<float> get_vecYawDEC(){return vec_yawDEC;}
    int get_last_index(){return mkfLastIndex;}
    long get_nrows(){return vec_time.size();}
    
    //Setters
    int set_last_index(long lastIndex);
    
    //Get Interpolated data
    
    /**
     * Function to get interpolated RA,DEC values for roll pitch and yaw axis at instantaneous time.
     * @param time: Instantaneous time provided by user
     * @param rollRA
     * @param rollDEC
     * @param pitchRA
     * @param pitchDEC
     * @param yawRA
     * @param yawDEC
     * @return 
     */
    int get_inst_rpy(double time, float &rollRA, float &rollDEC, float &pitchRA, 
                        float &pitchDEC, float &yawRA, float &yawDEC, long startIndex=0);
    /**
     * Function to get average of rollRA, rollDEC, pitchRA, pitchDEC, yawRA, yawDEC from mkf file.
     * @param rollRA
     * @param rollDEC
     * @param pitchRA
     * @param pitchDEC
     * @param yawRA
     * @param yawDEC
     * @param status
     * @return 
     */
    int get_avg_rpy(float &avg_rollRA, float &avg_rollDEC, float &avg_pitchRA, float &avg_pitchDEC, float &avg_yawRA, 
                float &avg_yawDEC);
    

    /**
     * gets mkf record for a particular row number.
     * @param record_no
     * @return 
     */
    int get_mkf_record(long record_no, MKFrecord &mkfr);
    //VERIFICATION
    
    /**
     * Function to check repetition of time values. If repetition is greater than 4 then corresponding
     * time values along with their repetition rate is printed on terminal.
     * If user provides outputFilename then these values are also written in this file.
     * @param outputFilename
     * @return 
     */
    int check_time_repetition(string outputFilename="");
};


/**
 * //Angular offset calculation
 * @param RApnt: degrees
 * @param DECpnt: degrees
 * @param rollRA
 * @param rollDEC
 * @return 
 */
double calculate_angular_offset(double RApnt, double DECpnt, double rollRA, double rollDEC);

/**
 * //Elevation calculation
 * @param RApnt: degrees
 * @param DECpnt: degrees
 * @param x
 * @param y
 * @param z
 * @return 
 */
double calculate_elevation(double RApnt, double DECpnt, double x, double y, double z);
//Sun angle & moon angle calculation

/**
 * Calculates 
 * @param utTime
 * @param RAsun: in degrees
 * @param DecSun: in degrees
 * @param RAmoon: in degrees
 * @param DecMoon: in degrees
 * @return 
 */
int sunmoon(double utTime, double *RAsun, double* DecSun, double* RAmoon, double *DecMoon);

int calculate_sun_moon_angle(double utTime,
         double rollRA, double rollDEC, double *sunAngle, double* moonAngle,double *OutRAsun,double *OutDecSun);
//added by shrikant
int gethk(double Time[],int DataID[],int HKChannelNo[],int ModeID[],int ErrorCount[],int BootPageNo[],long ADCOutput[],long AlphaCount[],long VetoCount[], int quadrant,double search_time,float hkparam[],int *error_count,int *boot_page_no,int *modeid,long *alphacounter,long *vetocounter,long nrowshdr,long *start_indexQ0,long *start_indexQ1,long *start_indexQ2,long *start_indexQ3);
long getNumrows(const char *infile, int hdunum);
#endif	/* MKFREGENERATION_H */


/* 
 * @file level1handler.h
 * @author Tanul Gupta
 * @date Created on May 20, 2015, 9:20 AM
 * @brief This file contains classes that will act as interface between level1 data and the program.
 * @details Getting, Setting and handling of level-1 data.
 */
#ifndef LEVEL1HANDLER_H
#define	LEVEL1HANDLER_H

#include <iostream>
#include <cstring>
#include <vector>
#include <algorithm>
#include "level2validation.h"
#include "glog/logging.h"
#include "utils.h"
#include "macrodef.h"
#include <numeric>


using namespace std;

struct TctSlopeOffsets {
    double tstart;
    double tstop;
    double slope;
    double offset;
};

struct tctStruct{
    double spsTime;
    double spsObt;
    double instrumentTime;
};

class Tct {
private:
    vector <double> vec_spsTime; 
    vector <double> vec_spsObt; 
    vector <double> vec_instrumentTime;
    vector <tctStruct> vec_tctdata;
    vector <TctSlopeOffsets> vec_tctSlopeOffsets;
    bool resetFlag;
    vector <long> resetIndex;
    long nrecOutTCT;
	long last_index;

   
public:

    Tct();
    /**
     * Reads level-1 tct file and stores data in object variables.
     * It also evaluates slope and offset values after linear fitting of spsTime &
     * instrumentTime at intervals of 100 seconds.
     * @param tctFilename
     * @return 
     */
    int read_tct_file(string tctFilename);
    int display();
    int is_reset();
    /**
     * Interpolates instrument time using TCT file and returns UTC time.
     * @param instrumentTime
     * @param intTime: interpolated time
     * @param isReset: whether this time value is after reset or not (1: reset; 0: not reset)
     * NOTE: isReset needs to be used to handle resetting of TCT time.
     * @return 
     */
    int interpolate_time(double instrumentTime,  double &intTime, bool isReset=0);
    
    int get_clockCorrectionFactor(double &corFac);

    int clean_tct_data();
    //Getters
    long get_nrecOutTCT(){ return nrecOutTCT;}
    vector <double> getSpsTime(){return vec_spsTime;}
    vector <double> getSpsObtTime(){return vec_spsObt;}
    vector <double> getInstrumentTime(){return vec_instrumentTime;}
    
};

//ATTITUDE FILE
struct attStruct{
   double time;
   vector <double> qSat;
   double rollRA;
   double rollDEC;
   double rollRot;
   
   attStruct();
    //Natural ordering is ascending in time.
    friend bool operator<(const attStruct &a, const attStruct &b) {
        return a.time < b.time;
    }
};

//Comparison function to check sorting
//These template function will be used to check whether attitude data is sorted or not.

bool whether_ascending_att(attStruct a, attStruct b);
bool whether_descending_att(attStruct a, attStruct b);

class Attitude{
private:
    string attFilename;
    vector <attStruct> att;
    //Keywords
    long tstarti;
    double tstartf;
    long tstopi;
    double tstopf;
    double RAPnt; //RA pointing (degrees)
    double DECPnt; //DEC pointing (degrees)
    string obsID;
    
    bool attFileRead;
    long lastRecord;
    long nrows;
    bool fileSorted; //if true then file is sorted with respect to time
    //in ascending order.
    bool headerReadFlag; //set to true if all header parameters have 
    //been read.
    double dataMinTime;
    double dataMaxTime;
    double minSamplingInterval; 
    double maxSamplingInterval;
    
public:
    Attitude();
    int read_attitude_file(string attFilename);
    attStruct get_interpolated_attitude(double time, bool extrapolateFlag = false);
    //GETTERS
    vector <attStruct> get_att(){return att;}
    double get_RAPnt(){return RAPnt;}
    double get_DECPnt(){return DECPnt;}
    void get_min_max_sampling_interval();
    void get_min_max_dataTime();
    double get_minDataTime(){return dataMinTime;}
    double get_maxDataTime(){return dataMaxTime;}
    //DISPLAY
    void display_header_keywords();
};
//ATTITUDE FILE END

//ORBIT FILE
struct orbStruct{
    double time; //seconds
    float x; //in km
    float y; // in km
    float z; // in km
    float vX; //in km/s
    float vY; // in km/s
    float vZ; // in km/s
    float lat; // in degrees
    float lon; // in degrees
    float altitude; //in km
    float semimajorAxis; //in km
    float ecentricity;  //spelled incorrectly in orbit file
    float inclination; //in degrees
    float aop; //in degrees
    float raan; //in degrees
    float trueAnomaly; //in degrees
    
    orbStruct();

    //Natural ordering is ascending.
    friend bool operator<(const orbStruct &a, const orbStruct &b){
        return a.time < b.time;
    }
};

//Comparison function to check sorting
//These template function will be used to check whether orbit data is sorted or not.

bool whether_ascending_orb(orbStruct a, orbStruct b);
bool whether_descending_orb(orbStruct a, orbStruct b);

class Orbit{
private:
    string orbFilename;
    vector <orbStruct> orb;
    //Keywords
    long tstarti;
    double tstartf;
    long tstopi;
    double tstopf;
    double RAPnt; //RA pointing (degrees)
    double DECPnt; //DEC pointing (degrees)
    string obsID;
    
    
    bool orbFileRead;
    long lastRecord;
    long nrows;
    bool fileSorted; //if true then file is sorted with respect to time
                     //in ascending order.
    bool headerReadFlag; //set to true if all header parameters have 
                         //been read.
    double dataMinTime;
    double dataMaxTime;
    double minSamplingInterval;
    double maxSamplingInterval;
public:
    Orbit();
    int read_orbit_file(string orbFilename);
    void sort_orbit_data();
    
    //GETTERS
    vector <orbStruct> get_orb(){return orb;}
    void get_min_max_dataTime();
    void get_min_max_sampling_interval();
    orbStruct get_interpolated_orbit(double time, bool extrapolateFlag=false);
    double get_minDataTime(){return dataMinTime;}
    double get_maxDataTime(){return dataMaxTime;}
    //DISPLAY
    void display_header_keywords();
};


//LBT file
struct lbtStruct{
    double time;
    //Quadrant 1
    unsigned char q1modeID;
    float q1Pos5VMonitor;
    unsigned long q1CZTCounter;
    float q1Temperature1;
    float q1VetoHVMonitor;
    float q1Pos2dot5VMonitor;
    float q1CZTHVMonitor;
    unsigned long q1VetoCounter;
    float q1vetoLLD;
    float q1DVDD;
    unsigned long q1AlphaCounter;
    //Quadrant 2
    unsigned char q2modeID;
    float q2Pos5VMonitor;
    unsigned long q2CZTCounter;
    float q2Temperature1;
    float q2VetoHVMonitor;
    float q2Pos2dot5VMonitor;
    float q2CZTHVMonitor;
    unsigned long q2VetoCounter;
    float q2vetoLLD;
    float q2DVDD;
    unsigned long q2AlphaCounter;
    //Quadrant 3
    unsigned char q3modeID;
    float q3Pos5VMonitor;
    unsigned long q3CZTCounter;
    float q3Temperature1;
    float q3VetoHVMonitor;
    float q3Pos2dot5VMonitor;
    float q3CZTHVMonitor;
    unsigned long q3VetoCounter;
    float q3vetoLLD;
    float q3DVDD;
    unsigned long q3AlphaCounter;
    //Quadrant 4
    unsigned char q4modeID;
    float q4Pos5VMonitor;
    unsigned long q4CZTCounter;
    float q4Temperature1;
    float q4VetoHVMonitor;
    float q4Pos2dot5VMonitor;
    float q4CZTHVMonitor;
    unsigned long q4VetoCounter;
    float q4vetoLLD;
    float q4DVDD;
    unsigned long q4AlphaCounter;
    
    unsigned char cztMemLvl;
    unsigned char cztBootPage;
    unsigned char cztPEErrCnt;
    unsigned long cztLastCmd;
    unsigned long cpmRate;

    lbtStruct();
    //Natural ordering is ascending.
    friend bool operator<(const lbtStruct &a, const lbtStruct &b) {
        return a.time < b.time;
    }
};

//Comparison function to check sorting
//These template function will be used to check whether lbt data is sorted or not.

bool whether_ascending_lbt(lbtStruct a, lbtStruct b);
bool whether_descending_lbt(lbtStruct a, lbtStruct b);
bool is_event_modeID(unsigned char modeID);
class LBT{
private:
    string lbtFilename;
    vector <lbtStruct> lbt;
    //Keywords
    long tstarti;
    double tstartf;
    long tstopi;
    double tstopf;
    double RAPnt; //RA pointing (degrees)
    double DECPnt; //DEC pointing (degrees)
    string obsID;
    
    bool lbtFileRead;
    long lastRecord;
    long nrows;
    bool fileSorted; //if true then file is sorted with respect to time
    //in ascending order.
    bool headerReadFlag; //set to true if all header parameters have 
    //been read.
    double dataMinTime;
    double dataMaxTime;
    double minSamplingInterval;
    double maxSamplingInterval;
public:
    LBT();
    int read_lbt_file(string lbtFilename);
    
    //GETTERS
    vector <lbtStruct> get_lbt(){return lbt;}
    lbtStruct get_lbt(double time, bool interpolateFlag=false);
    void get_min_max_sampling_interval();
    void get_min_max_dataTime();
    double get_minDataTime(){return dataMinTime;}
    double get_maxDataTime(){return dataMaxTime;}
    //DISPLAY
    void display_header_keywords();
};


// SORT CHECKING FUNCTIONS END

#endif	/* LEVEL1HANDLER_H */


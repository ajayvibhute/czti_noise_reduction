
#ifndef CZTHEADERPARAM_H
#define CZTHEADERPARAM_H

#include<utils.h>

// For reading and writing header parameters
class cztHeaderParam{

private:
    double tstart;
    double tstop;
    long tstarti;
    long tstopi;
    double tstartf;
    double tstopf;
    string date_obs;
    string time_obs;
    string date_end;
    string time_end;
    string object;
    float ra_pnt;
    float dec_pnt;
    float ra_obj;
    float dec_obj;
    string obs_id;
    string obs_mode;
    string origin;
    string mission;
    string telescop;
    string instrume;
    string filename;
    double timedel;
    double exp_time;
	double telapse;
    long mjdrefi;
    long mjdreff;
    string timesys;
    string timeunit;
    int equinox;

    string radecsys;

public:
    
    int quadid;
    int set_default_values();
    int readFromHeader(fitsfile *fptr);
    int writeToHeader(fitsfile *fptr);
    int writeTimekey(double tstart, double tstop,fitsfile *fptr); 
    int getTelapse(double &telapse);
};


#endif 

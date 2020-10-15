/* 
 * File:   cztgtigen.h
 * Author: preeti
 *
 * Created on September 11, 2012, 10:16 AM
 */

#ifndef CZTGTIGENV2_H
#define	CZTGTIGENV2_H

#include<pil.h>
#include "utils.h"
#include "glog/logging.h"
#include "validations.h"
#include "level1handler.h"
#include "caldbHandler.h"
#include "level2validation.h"
#include "gtiHandler.h"
#include "macrodef.h"
#include "mkfRegeneration.h"
#include "Mvector.h"
#include <fstream>
#include "cztHeaderParam.h"

/*#define MKF 1
#define CUSTOM 2
#define MKF_REFINED 3
#define CUSTOM_REFINED 4 
#define COPY 5
*/

using namespace std;

class cztgtigen{
private:
    char modulename[NAMESIZE];
    char eventfile[PIL_LINESIZE];
    char mkffile[PIL_LINESIZE];
    char thresholdfile[PIL_LINESIZE];  //ASCII file containing valid ranges for each parameter 
    char outfile[PIL_LINESIZE];
    char usergtifile[PIL_LINESIZE];
	int clobber;
    int history;
public:
    cztgtigen();
    void display();
    int read(int argc,char **argv);
    int read(char *eventfile,char *mkffile,char *thresholdfile,char *outfile,char *usergtifile,int clobber=NO,int history=NO);
    
    int cztgtigenProcess();
    void getHistory(vector<string> &vhistory);
    
};


class MKF_Threshold{
 public:
    string parameter;
    string range;
    float LV,UV;              //lower value and upper value; 
    void display();
    
};

#endif	/* CZTGTIGEN_H */


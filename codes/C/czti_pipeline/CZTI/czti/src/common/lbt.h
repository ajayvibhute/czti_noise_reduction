/* 
 * File:   lbt.h
 * Author: preeti
 *
 * Created on November 28, 2011, 5:59 PM
 */

#ifndef LBT_H
#define	LBT_H

#include<string>
#include<fitsio.h>
#include<vector>
#include "glog/logging.h"

using namespace std;

class LBT_Threshold{
 public:
    string parameter;
    string range;
    int LV,UV;              //lower value and upper value; 
    void display();
    
};

int readLBT_Threshold(char *file,vector<LBT_Threshold> &lbt_thr);

class LBTparam{
public:
    double time;
    unsigned char q1_modeid;
    unsigned char q1_5v_monitor;
    unsigned char q1_msb_czt_counter;
    unsigned char q1_temperature;
    unsigned char q1_lsb_czt_counter;
    unsigned char q1_2p5v_monitor;
    unsigned char q1_vetoHVmonitor;   
    unsigned char q1_cztHVmonitor;
    unsigned char q1_msb_vetocounter;
    unsigned char q1_veto_lld;
    unsigned char q1_dvdd;
    unsigned char q1_msb_alphacounter;
    unsigned char q1_lsb_alphacounter;
    unsigned char q1_io9_15_8;
    unsigned char q1_io9_7_0;
    unsigned char q1_lsb_vetocounter;
    unsigned char q2_modeid;
    unsigned char q2_5v_monitor;
    unsigned char q2_msb_czt_counter;
    unsigned char q2_temperature;
    unsigned char q2_lsb_czt_counter;
    unsigned char q2_2p5v_monitor;
    unsigned char q2_vetoHVmonitor;
    unsigned char q2_cztHVmonitor;
    unsigned char q2_msb_vetocounter;
    unsigned char q2_veto_lld;
    unsigned char q2_dvdd;
    unsigned char q2_msb_alphacounter;
    unsigned char q2_lsb_alphacounter;
    unsigned char q2_io9_15_8;
    unsigned char q2_io9_7_0;
    unsigned char q2_lsb_vetocounter;
    unsigned char q3_modeid;
    unsigned char q3_5v_monitor;
    unsigned char q3_msb_czt_counter;
    unsigned char q3_temperature;
    unsigned char q3_lsb_czt_counter;
    unsigned char q3_2p5v_monitor;
    unsigned char q3_vetoHVmonitor;
    unsigned char q3_cztHVmonitor;
    unsigned char q3_msb_vetocounter;
    unsigned char q3_veto_lld;
    unsigned char q3_dvdd;
    unsigned char q3_msb_alphacounter;
    unsigned char q3_lsb_alphacounter;
    unsigned char q3_io9_15_8;
    unsigned char q3_io9_7_0;
    unsigned char q3_lsb_vetocounter;
    unsigned char q4_modeid;
    unsigned char q4_5v_monitor;
    unsigned char q4_msb_czt_counter;
    unsigned char q4_temperature;
    unsigned char q4_lsb_czt_counter;
    unsigned char q4_2p5v_monitor;
    unsigned char q4_vetoHVmonitor;
    unsigned char q4_cztHVmonitor;
    unsigned char q4_msb_vetocounter;
    unsigned char q4_veto_lld;
    unsigned char q4_dvdd;
    unsigned char q4_msb_alphacounter;
    unsigned char q4_lsb_alphacounter;
    unsigned char q4_io9_15_8;
    unsigned char q4_io9_7_0;
    unsigned char q4_lsb_vetocounter;
    unsigned char syncword;
        
    int read(fitsfile *fptr,int row);
    void display();

};

int readLBT(char *file,char *lbtextname,vector<LBTparam> &lbtparam);

#endif	/* LBT_H */


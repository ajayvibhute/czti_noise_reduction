
/* 
 * @file cztgaas_v2.h
 * @author Tanul Gupta
 * @date Created on May 17, 2015, 9:28 AM
 * @brief This module computes the average aspect of satellite corresponding to a specific event file.
 * @details This module computes the average aspect of satellite corresponding to a specific event file.
 */

#ifndef CZTGAAS_V2_H
#define	CZTGAAS_V2_H

#include "caldbHandler.h"
#include <iomanip>
#include <pil.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>
#include <utils.h>
#include "glog/logging.h"
#include "cztstring.h"
#include "level1handler.h"
#include "level2validation.h"
#include "fitsio.h"
#include "coordinateTransformation.h"
#include "mkfRegeneration.h"
#include "Mvector.h"
#include <numeric>

using namespace std;

class Cztgaas{
private:
    char modulename[NAMESIZE];
    char eventfile[PIL_LINESIZE];
    /*char teldeffile0[PIL_LINESIZE];
    char teldeffile1[PIL_LINESIZE];
    char teldeffile2[PIL_LINESIZE];
    char teldeffile3[PIL_LINESIZE];*/
    char mkffile[PIL_LINESIZE];
    char outAspectFile[PIL_LINESIZE];
    int clobber;
    int history;
public:    
    Cztgaas();
    int read(int argc,char **argv);
    int read(char *eventfile,char *outAspectFile, char* mkffile, int clobber=NO,int history=NO);
    void display();
    
    int cztgaas_process();
    
    /**
    * Creates history for the module to be written to the output files
    * @param par : Object of cztgaas_param class
    * @param vhistory : vector ho hold history strings
    * @return 
    */
    int get_history(vector<string> &vhistory);

};

int create_aspect_file(string outFileName);



#endif	/* CZTGAAS_V2_H */


#ifndef CZTDPIGEN_V2_H
#define	CZTDPIGEN_V2_H

#include <iostream>
#include <cstdlib>
#include "utils.h"
#include "glog/logging.h"
#include "level2validation.h"
#include <vector>
#include <string>
#include "caldbHandler.h"
#include "cztstring.h"
#include "validations.h"
#include "fitsio.h"

class Cztdpigen{
private:
    char modulename[NAMESIZE];
    char infile[PIL_LINESIZE];    //input event file
    char outDPHfile[PIL_LINESIZE];   //output filename, will serve as prefix if more than one DPH is to be produced
    char outDPIfile[PIL_LINESIZE];   //output filename, will serve as prefix if more than one DPI is to be produced
    string effAreafile;	
    char badpixFile[PIL_LINESIZE]; //Badpixel file
    char ebins[PIL_LINESIZE];
    char timerange[PIL_LINESIZE];
    char quadsToProcess[25]; //new addition to incorporate quadrant wise functionality
    int badpixThreshold;
    int ntbins;
    int nebins;  //number of energy bins, decides number of output files to be produced
    int ndpi;
    double *emin,*emax;
    double *tstart,*tstop;
    int history;
    int clobber;


public:    
    

    Cztdpigen();
    void display();
    int read(int argc,char **argv);
    int read(char* infile,  char* badpixelFile, char* outDPHfile,
            char* outDPIfile, char *quadsToProcess, int badpixThreshold,
            char* ebins, char *timerange,
            int clobber=YES, int history=YES);
    
    /**
    * Function to generate DPI from event file
    * @param par
    * @return 
    */
    int cztdpigenProcess();   
    
    /**
    * Function to generate history for the module to be written to the output file
    * @param vhistory
    * @return 
    */
    int getHistory(vector<string> &vhistory);
    
};

/**
 * Function to generate the names of output DPI files based on number of DPIs to
 * be generated.
 * @param ndpi: number of DPIs to be generated
 * @param outfile: name of output file provided by user
 * @param outputFileNames: string vector containing names of output DPI files.
 * @return 
 */
int getNoutputfiles(int ndpi,char *outDPIfile, char* outDPHfile,
        vector<string> &outputDPIFileNames, vector<string> &outputDPHFileNames);
/**
 * Function to get tstart and tstop values from timerange string
 * Reads time from input event data file if timerange is '-' and provide the 
 * maximum tstart and minimum tstop i.e. tstart and tstop which is
 * which is common to all quadrants.
 * @param fptr
 * @param timerange
 * @param tstart
 * @param tstop
 * @return 
 */
int getTime(fitsfile *fptr,char *timerange,int ntbins,double *tstart,double *tstop, vector <int> quadsToProcess);

/**
 * Function to get emin and emax values for each ebin entered. If ebins="-", then the program selects
 * maximum and minimum energy based on ecol. If ecol="PI"or"PHA" then emin=0 and emax=511, otherwise if 
 * ecol="ENERGY" then emin=10 and emax=150;
 * @param ebins
 * @param ecol
 * @param nebins
 * @param emin
 * @param emax
 * @return 
 */
int getEnergy(char *ebins,int nebins,double *emin,double *emax, string eboundsFile);

#endif	/* CZTDPIGEN_V2_H */


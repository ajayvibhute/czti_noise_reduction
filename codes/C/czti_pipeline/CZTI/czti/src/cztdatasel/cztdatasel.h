

#ifndef CZTDATASEL_H
#define	CZTDATASEL_H

#include<pil.h>
#include<cstdlib>
#include "glog/logging.h"

#include"utils.h"
#include "cztHeaderParam.h"

using namespace std;

class cztdatasel{
private:
    char modulename[NAMESIZE];
    char infile[PIL_LINESIZE];        //Input Event Data file
    char gtifile[PIL_LINESIZE];        // gti file for filtering
    char gtitype[PIL_LINESIZE]; //edited by Mithun to ask which gti type
    char outfile[PIL_LINESIZE];
    int clobber;                      //Overwrite Existing file
    int history;
public:
    cztdatasel();
    int read(int argc, char **argv);
    int read(char *infile,char *gtitype, char *outfile,int clobber=YES,int history=NO);
    void display();
    /**
    * Function to select event data based on GTIs. The GTI file should contain four 
    * extensions each for four quadrants or one extension to be used for all four quadrants.
    * The function will update the input file if output filename is either '-' or same as input filename
    * The function will add GTI extensions to the output file
    * @param par
    * @return 
    */
    int cztdataselProcess();
    /**
    * Function to create history for the module to be added to input file
    * @param par
    * @param vhistory
    * @return 
    */
    int getHistory(vector<string> &vhistory);
    
};


//Function to compute exposure time  and tstart,tstop from GTI
int computeExptime(fitsfile *fgti,double *exp_time, double *tstart, double *tstop);

/**
 * Function to remove any existing GTI extensions from the input file
 * @param fptr
 * @return 
 */
int removeGTIext(fitsfile *fptr);


#endif	/* CZTDATASEL_H */


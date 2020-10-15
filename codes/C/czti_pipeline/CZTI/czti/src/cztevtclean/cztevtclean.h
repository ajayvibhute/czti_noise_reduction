

#ifndef CZTEVTCLEAN_H
#define	CZTEVTCLEAN_H


#include <iostream>
#include <pil.h>
#include <fitsio.h>
#include <vector>
#include <cstring>
#include <string>
#include "utils.h"
#include "glog/logging.h"
#include "errorHandler.h"

class cztevtclean{
private:
     char modulename[NAMESIZE];
     char infile[PIL_LINESIZE];          //Input Event Data file
     char outfile[PIL_LINESIZE];         //Output selected file
     char alphaval[10];                   //alpha values to accept, it can be '0','1','0,1'
     char vetorange[512];                 //veto range/ranges to accept 
     int isdoubleEvent;
     int clobber;                        //Overwrite Existing file
     int history;

     //if outfile and infile are same or if outfile is " " then infile is updated 
public:
     cztevtclean();
     int read(int argc, char **argv);
     int read(char *infile,char *outfile,char *alphaval,char *vetorange,int isdoubleEvent=NO,
        int clobber=NO,int history=NO);
     void display();
     
     int cztevtcleanProcess();
     
     /**
    * Creates history for the module to be written to the output files
    * @param par : Object of cztevtclean_param class
    * @param vhistory : vector ho hold history strings
    * @return 
    */
    int getHistory(vector<string> &vhistory);
    
};



/**
 * Function to check alpha satisfies given criteria in alphaval
 * @param alpha
 * @param alphaval
 * @return - returns 0 if alpha satisfies the criteria
 */
int checkAlpha(unsigned char alpha,char *alphaval);

/**
 * Function to check whether veto is within vetorange
 * @param veto
 * @param vetorange
 * @return : returns 0 if alpha satisfies the criteria
 */
int checkVeto(unsigned char veto,char *vetorange);

#endif	/* CZTEVTCLEAN_H */



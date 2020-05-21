/* 
 * File:   validations.h
 * Author: Tanul Gupta
 *
 * Created on November 11, 2014,10:26 AM
 * 
 * Description:
 * It contains functions to validate different type of input arguments.
 */

#ifndef VALIDATIONS_H
#define	VALIDATIONS_H

#include <iostream>
#include "glog/logging.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <string>
#include "cztstring.h"
#include "regex.h"
#include <fitsio.h>

/*macro-variable definitions*/
#define EVENTDATA_HDU_FIRST 2
#define EVENTDATA_HDU_LAST 5
#define EVENTDATA_MIN_COLS 10
/****definition ends here****/
using namespace std;

int isStringValid(string inString, vector<string> sampleStrings);

int isOversamplingfactorValid(int oversamplingfactor);


/******************************************************************************************************
 *Function to check the validity of quadrant string passed. It checks that the user has passed no other value 
 * than 0,1,2,3 or '-'.
 * @param quadString - variable in which user defined character array of quadrants to be processed is passed.
 ******************************************************************************************************/
int isQuadrantNumbersValid(char* quadString);

/******************************************************************************************************
 *Function to check the validity of number of background coefficients passed. It checks that the user has passed no other value 
 * than 0,1,2,3,4,5,6.
 * @param nBkgdCoef - number of Background coefficients.
 ******************************************************************************************************/
int isNBkgdCoefValid(int nBkgdCoef);

/******************************************************************************************************
 *Function to check the validity of resolution limit. It checks that the value in the range 2-12.
 * @param resolution_limit - resolution limit in arcminutes..
 ******************************************************************************************************/
int isResolutionLimitValid(int resolution_limit);

/******************************************************************************************************
 *Function to check the validity of bins range entered by the user.
 * @param ebins- bin range, which could be time bins or energy bins.
 *****************************************************************************************************/
int isBinsValid(char* bins);

int isOuttypeValid(char* outtype);

/**
 * To check that the parameter lies in the specified range
 * @param minval: minimum value that parameter can attain
 * @param maxval: maximum value that parameter can attain
 * @param val: value of parameter
 * @return 
 */
template<class T> int checkRange(T minval, T maxval, T val);
template<class T> int checkRange(T minval, T maxval, T val){
    if(val<minval || val>maxval){
        LOG(ERROR)<<"***Allowed range is [" << minval << "," << maxval << "] ***";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool isEventFileValid(fitsfile *fptr);


#endif	/* VALIDATIONS_H */
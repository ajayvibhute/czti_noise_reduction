#include <vector>

#include "validations.h"

int isStringValid(string inString, vector<string> sampleStrings){
    int istring=0;
    int badFlag=TRUE;
    
    for(istring=0; istring<sampleStrings.size(); istring++){
        if(inString == sampleStrings[istring]){
            badFlag=FALSE;
            break;
        }
    }
    return badFlag;
}
int isOversamplingfactorValid(int oversamplingfactor){
    float intpart=0, power=0, fraction=0;
    power = (float) log(oversamplingfactor) / (float) (log(2));
    fraction = modf(power, &intpart);
    if(oversamplingfactor<=0){
        LOG(ERROR) << "*** Oversampling factor cannot be 0 or negative value***";
        return (EXIT_FAILURE);
    }
    
    if (fraction > 0) {
        LOG(ERROR) << "***Over sampling must be provided in powers of 2***";
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}

int isQuadrantNumbersValid(char* quadString) {    
    if (strcmp(quadString, "-") == 0) {
        strcpy(quadString, "0,1,2,3");
    }
  // LOG(INFO)<<quadString<<"\n"; 
    // REGULAR EXPRESSION VALIDATION
    regex_t reg;
    string pattern = "^((([0-3]){1})(,[0-3]){1,3})*";//added by mayuri for fedora compatibility purpose
//"([0-3]){1}(,[0-3])*";//commented by mayuri 
    size_t nmatch=1;
    regmatch_t matches[2];
    int result=regcomp(&reg, pattern.c_str(),REG_EXTENDED|REG_ICASE);
//LOG(INFO)<<"Regex result : "<<result<<"\n";
//LOG(INFO)<<"Regexec result : "<<regexec(&reg, quadString,nmatch,matches ,0)<<"\n";
    if (regexec(&reg, quadString,nmatch,matches ,0)==0){
        if(strlen(quadString)!= matches[0].rm_eo-matches[0].rm_so){
            LOG(INFO) << "Valid Quadrant IDs entered by the user: " <<((string)quadString).substr(matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so);
            LOG(ERROR) << "*** Although the quadrant IDs are partially correct, but the entire string entered by the user does not conform to the standards. ***";
            LOG(ERROR) << "*** Quadrant IDs can only be 0 or 1 or 2 or 3***";
            LOG(ERROR) << "Example Input: 0,2,3"; 
            return (EXIT_FAILURE);
        }
    } else {
            LOG(ERROR) << "*** Quadrant IDs can only be 0 or 1 or 2 or 3***";
            LOG(ERROR) << "Example Input: 0,2,3"; 
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);

}

int isBinsValid(char* bins){

    regex_t reg;
    string pattern = "(([0-9]+)-([0-9]+)){1}(,([0-9]+)-([0-9]+))*";
    regmatch_t matches[1];
    regcomp(&reg, pattern.c_str(), REG_EXTENDED | REG_ICASE);
    if (!(strcmp(bins, "-") == 0)) {
        if (regexec(&reg, bins, 2, matches, 0) == 0) {
            if (strlen(bins) != matches[0].rm_eo - matches[0].rm_so) {
                LOG(INFO) << "Valid bins entered by the user: " << ((string) bins).substr(matches[0].rm_so, matches[0].rm_eo - matches[0].rm_so);
                LOG(ERROR) << "*** Although the bins entered are partially correct, but the entire string entered by the user does not conform to the standards. ***";
                LOG(ERROR) << "*** Bins should be entered in the following manner***";
                LOG(ERROR) << "Example Input: 0-2,2-5,3-9";
                return (EXIT_FAILURE);
            }
        } else {
            LOG(ERROR) << "*** Bins should be entered in the following manner***";
            LOG(ERROR) << "Example Input: 0-2,2-5,3-9";
            return (EXIT_FAILURE);
        }
    }
    
//    regex_t reg_extract;
//    string extract_pattern = "([0-9]+)-([0-9]+)+";
//    regmatch_t matches_extract[512];
//    regcomp(&reg_extract, extract_pattern.c_str(), REG_EXTENDED | REG_ICASE);
//    int length_bins = strlen(bins);
//    string bins_left;
//   regexec(&reg_extract, bins, 2, matches_extract, 0);
//        LOG(INFO) << ">>DBG " << ((string) bins).substr(matches_extract[0].rm_so, matches_extract[0].rm_eo - matches_extract[0].rm_so); 
//        bins_left = ((string) bins).substr(matches_extract[0].rm_eo, length_bins-matches_extract[0].rm_eo);
//        string binss(bins_left) ;
//        LOG(INFO) << binss;
//        regexec(&reg_extract, binss.c_str(), 2, matches_extract, 0);
//        LOG(INFO) << ">>DBG " << (binss).substr(matches_extract[0].rm_so, matches_extract[0].rm_eo - matches_extract[0].rm_so); 
    
    
    return EXIT_SUCCESS;
}

int isNBkgdCoefValid(int nBkgdCoef){
    if (nBkgdCoef > 6 || nBkgdCoef < 1) {
            LOG(ERROR) << "***Error reading number of background coefficients --Allowed range is 1-6***";
            return EXIT_FAILURE;
        }
    return EXIT_SUCCESS;
}

int isResolutionLimitValid(int resolution_limit){
    if(resolution_limit<2 || resolution_limit>16){
        LOG(ERROR)<< "*** Error reading resolution limit. Allowed range is <2-16>***";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}


bool isEventFileValid(fitsfile *fptr) {
    int status = 0;
    int hdutype = 0;
    int numhdu = 0;
    int ncols = 0;
    fits_get_num_hdus(fptr, &numhdu, &status);
    if (status) {
        LOG(ERROR) << "***Error in getting number of HDUs - isEventFileValid()***";
        return false;
    }
    LOG(INFO) << "Number of HDUs in event file are " << numhdu;
    if (numhdu < 5) {
        LOG(ERROR) << "Least number of HDUs in event file must be 5 -isEventFileValid()******";
        return false;
    }
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    if (status) {
        LOG(ERROR) << "***Error moving to HDU 1 in event file - isEventFileValid()***";
        return false;
    }
    if (hdutype != IMAGE_HDU) {
        LOG(ERROR) << "***Primary Extension is not Image HDU - is EventFileValid()***";
        return false;
    }
    for (int i = EVENTDATA_HDU_FIRST; i <= EVENTDATA_HDU_LAST; i++) {
        fits_movabs_hdu(fptr, i, &hdutype, &status);
        if (status) {
            LOG(ERROR) << "***Error moving to HDU " << i << " in event file***";
            return false;
        }

        if (hdutype != BINARY_TBL) {
            LOG(ERROR) << "***Extension " << i << " must be Binary Table - is EventFileValid()***";
            return false;
        }

        fits_get_num_cols(fptr, &ncols, &status);
        if (status) {
            LOG(ERROR) << "***Error reading number of columns - is EventFileValid()***";
            return false;
        }

        if (!(ncols >= EVENTDATA_MIN_COLS)) {
            LOG(ERROR) << "***NUmber of columns in Binary HDUs cannot be less than "<< EVENTDATA_MIN_COLS << "***";
            return false;
        }
    }
    return true;
}

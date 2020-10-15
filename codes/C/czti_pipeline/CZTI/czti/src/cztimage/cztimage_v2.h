
#ifndef CZTIMAGEV2_H
#define CZTIMAGEV2_H

#include <pil.h>
#include "utils.h"
#include "validations.h"
#include "glog/logging.h"
#include "caldbHandler.h"
#include "coordinateTransformation.h"
#include "level2validation.h"
#include <vector>
#include "caldbHandler.h"
#include "maskGeometry.h"

#include "fft.h"

#define DPIFILE 0
#define DPHFILE 1
#define SHADOWFILE 2
using namespace std;
class Cztimage{
    private:
    //Input parameters
    char modulename[NAMESIZE];
    char inputType[NAMESIZE];
    char infile[PIL_LINESIZE]; //Input DPI/DPH/SHADOW file for image creation
    //char maskfile64x64[PIL_LINESIZE];
    //char maskfile8100[PIL_LINESIZE];
    char shadowfile[PIL_LINESIZE];
    char aspectfileQ0[PIL_LINESIZE];
    char aspectfileQ1[PIL_LINESIZE];
    char aspectfileQ2[PIL_LINESIZE];
    char aspectfileQ3[PIL_LINESIZE];
    char catalogfile[PIL_LINESIZE];
    int catalogExtnum;
    char RA_colname[25];
    char DEC_colname[25];
    char FLUX_colname[25];
    float resolutionLimit;
    float threshold;
    char sourcelist[PIL_LINESIZE];
    char outImgFile[PIL_LINESIZE];
    char quadsToProcess[25];
    int oversamplingfactor;
    float sourceStrengthThr;
    int nBkgdCoef; //number of background coefficients
    int debugmode; //if yes/true, creates intermediate file 
    int clobber;
    int history;
    
    int inFlag; //0:DPI 1:DPH 2: SHADOW
    
public:
    Cztimage();
    int read(int argc, char** argv);
    void display();
    int cztimageProcess();
};

#endif /*CZTIMAGEV2_H*/

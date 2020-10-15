
#ifndef CZTRSPGEN_H
#define CZTRSPGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include<fitsio.h>
#include<pil.h>
#include<pthread.h>
#include<iomanip>
#include <string>
#include <utils.h>
#include "level2validation.h"
#include "errorHandler.h"
#include "Exposure.h"

#define COMP_COLS 327
#define TOTALROWS 8175
#define DET_NUMPIXELS 64
#define COLS 1950
#define ROWS 1950
#define MATHPI 3.14159265


int checkbit(int *,int ,int);
void calculateTrans(float,float ,float ,int *,float *);
void getModule(int ,int *);
void printerror( int status);
void getMaskPattern(char*,int *,int);
int getShadow(float tx,float ty,int qid,float *shadow_pixels,int *maskElements, ExposureTable &exptable);

extern"C" {
    void fcompute_tx_ty_(double *RAroll,double *Decroll,double *roll_rot,double *RAsrc,double *Decsrc,double *thetax,double *thetay);
}

int calculate_renormalized_weights(ExposureTable &expTable,Badpix &badpix,int badpixThreshold,int quadID,string effareaFilename,string evtfilename);

using namespace std;

class cztrspgen{

    private:
        char modulename[NAMESIZE];
        char phafile[PIL_LINESIZE];    
        char rspfile[PIL_LINESIZE];
		char evtfile[PIL_LINESIZE];
        char badpixfile[PIL_LINESIZE];

	    char compmaskfile[PIL_LINESIZE];
    	char cameraGeomFile[PIL_LINESIZE];
	    char effectiveAreafile[PIL_LINESIZE];
	    char eboundsfile[PIL_LINESIZE];
    	char LLDfile[PIL_LINESIZE];
		char pixrespFile[PIL_LINESIZE];
		char respparFile[PIL_LINESIZE];

        int clobber;
        int history;

    public:
        
        cztrspgen();
        ~cztrspgen();
        int cztrspgenProcess();
        int read(int argc,char **argv);
		void display();
		int read(char *phafile,char *rspfile,char *evtfile, char *compmaskfile,char *eboundsfile, char *cameraGeomFile,
        char *effectiveAreafile,char *LLDfile, char *respparFile,char*pixrespFile,char *badpixfile,int clobber,int history);

        
        int write_rsp(float *Elo,float *Ehi,int *ngrp,int *fchan,int *nchan,float **rspmatrix,int n_Ebin,int n_PIbin,char *febounds_name);
        int getHistory(vector<string> &vhistory);
};

int compute_exposure_frac(double a[16][5],int quad);

#endif

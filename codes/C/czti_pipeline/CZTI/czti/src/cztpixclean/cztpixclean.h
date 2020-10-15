/*cztpixclean.h

Mithun NPS
(30/11/15)
*/

#ifndef CZTPIXCLEAN_H
#define CZTPIXCLEAN_H

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
#include"utils.h"
#include "errorHandler.h"
#include "cztstring.h"
#include "caldbHandler.h"
#include "level2validation.h"

class cztpixclean{
private:
	char modulename[100];
	char infile[PIL_LINESIZE];
	char inlivetimefile[PIL_LINESIZE];
    char outlivetimefile[PIL_LINESIZE];
    char outfile1[PIL_LINESIZE];
    char outfile2[PIL_LINESIZE];
	int clobber;                      //Overwrite Existing file
    	int history;	
	//char caldb_badpix[PIL_LINESIZE];
	char badpixfile[PIL_LINESIZE];
	int nsigma;
	float det_tbinsize,pix_tbinsize;
	int det_count_thresh,pix_count_thresh;
	int writedblevt;

public:
	cztpixclean();
	~cztpixclean();	
	int cztpixcleanProcess();
	void display();
	int read(int argc,char **argv);
    int read(char *infile,char *inlivetimefile,int writedblevt,char *outfile1, char *outfile2,char *outlivetimefile, char *badpixfile, int nsigma,float det_tbinsize,float pix_tbinsize,int det_count_thresh,int pix_count_thresh);

};

//Independant functions
int remove_events(char* infile,int writedblevt,char* outfile1, char *outfile2,char *evt_flag,int qid,double **pixel_exposure);
double compute_mean(long *dph,int *pixflag);
double compute_stdev(long *dph,int *pixflag,double avg);
int findbadpix(long *dph,int *pixflag,int nsig);
int read_lldfile(string caldb_lld,int qid,int *lld);
int read_badpixfile(string caldb_badpix,int qid, int *pixflag);
int write_badpix(char *fname,int **flg);
int create_primaryfits(fitsfile *fout);

#endif

/*cztbunchclean.h


Mithun NPS
(30/11/15)
*/

#ifndef CZTBUNCHCLEAN_H
#define CZTBUNCHCLEAN_H

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
#include "level2validation.h"

class cztbunchclean{
private:
    char modulename[NAMESIZE];
    char infile[PIL_LINESIZE];
	char bunchfile[PIL_LINESIZE];
    char outfile[PIL_LINESIZE];
    char livetimefile[PIL_LINESIZE];
	int bunchdeftime;
	float skipT1,skipT2,skipT3;
	int bunch_length_thresh;
    int clobber;                      
    int history;
	double livetime_binsize;

public:
	cztbunchclean();
	~cztbunchclean();	
	int cztbunchcleanProcess();
	int read(int argc,char **argv);
    int read(char *infile,char *bunchfile,char *outfile,char * livetimefile,int bunchdeftime,int bunch_length_thresh,float skipT1 ,float skipT2,float skipT3,double livetime_binsize, int clobber, int history);
    int getHistory(vector<string> &vhistory);
	void display();

 //	int calculateLivetime(char*infile,int qid);
//	int writeOutfile(double *timearray,float*livetime,int qid,char*filename);
    

};

//Independant function
int remove_events(char* infile,char* outfile,char *evt_flag,int qid,double **pixel_exposure);


#endif

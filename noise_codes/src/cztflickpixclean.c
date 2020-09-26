/*
cztflickpixclean.c
* Author : Ajay Ratheesh, TIFR Mumbai
* 		   Ajay Vibhute, IUCAA Pune
* 		   Mayuri Shinde, IUCAA Pune

Code is developed for removing flickering pixels

inputs:- 1.Event file name 2.Badpix file name
output:- 1.Flickering pix cleaned level2 event file 2.modified bad pixel file 3.BTI file
* 
* Updated by Ajay Ratheesh to include GTI extentions on 20 Dec, 2017
* 
* 
*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "pil.h"
/*#define flick_energy_thresh 50000
#define sigma_thresh_1 3
#define flick_count_thresh_1 200000
#define sigma_thresh_2 5
#define flick_count_thresh_2 1*/
#define BUFFSIZE 1024
int dimen=64;
void processFlickPixReduction();
void printerror( int status);
void createEventFile(char *outputfile,char *infile);
void createBadpixFile(char *outputfile);
void createBTIFile(char *outputfile);
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum, double exposure);
//void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int bufsize,int hdunum);
void writeBadpix(unsigned char *detid,unsigned char *pixid,unsigned char *pixx,unsigned char *pixy,unsigned char *pix_flag,char *outputfile,int hdunum);
void writeBTI(unsigned char *detid,unsigned char *pixid,unsigned char *detx,unsigned char *dety,double *tstart,double *tstop,int writesize,char *outputfile,int hdunum);
void modifyEventHeaderParams(char *outputfile);
void writeBadpixExtension(char *outputfile);
void modifyExposure(char *outputfile, double **pix_exposure,int qid);
int pixThresholdFix(int countrate, int num_pix);
void modifycommonGTI(char *outputfile);

int read_input_parameters(int r, char *infile, char *inbadpixfile,char *thresholdfile,char *outfile,char *outbadpixfile,char *outbtifile,int *clobber,int *history);
int display_input_parameters(char *infile,char *inbadpixfile, char *thresholdfile,char *outfile,char *outbadpixfile,char *outbtifile,int clobber,int history);

char *CZTNOISECLEAN,infile[BUFFSIZE],inbadpixfile[BUFFSIZE],parfilename[BUFFSIZE], outbadpixfile[BUFFSIZE],outbtifile[BUFFSIZE],outfile[BUFFSIZE],thresholdfile[BUFFSIZE];
int clobber=1, history=1;
int flick_energy_thresh,sigma_thresh_1,flick_count_thresh_1,sigma_thresh_2,flick_count_thresh_2;

int main(int argc, char **argv)
{
	/*if(argc==3)
        {
		processFlickPixReduction(argv[1],argv[2]);
	}
	if(argc!=3)
	{
		printf("Enter all command line arguments\n1:Event file name\n2.Badpix file name\n");
		exit(-1);
	}*/

	CZTNOISECLEAN = getenv ("CZTNOISECLEAN");
	if(CZTNOISECLEAN==NULL)
	{
		printf("CZTNOISECLEAN Variable is not set\n");
		exit(0);
	}
	strcpy(parfilename,CZTNOISECLEAN);
	strcat(parfilename,"/paramfiles/cztflickpixclean.par");
	PILSetModuleName(parfilename);
	//Read all inputs through Parameter Interface Library
	int r=PILInit(argc,argv);
	
	read_input_parameters(r,infile,inbadpixfile,thresholdfile,outfile,outbadpixfile,outbtifile,&clobber,&history);
	display_input_parameters(infile,inbadpixfile,thresholdfile,outfile,outbadpixfile,outbtifile,clobber,history);	
	processFlickPixReduction();
	
	
	return 0;
}

void processFlickPixReduction()
{
	fitsfile *caldbfptr,*evtfptr,*badpixfptr; 
	FILE *logfile;      
	int ii,flag,status, hdunum, hdutype,   anynull,i,intnull,badpixdetidcolnum,badpixpixidcolnum,pixflagcolnum,size=1;
	long frow, felem, nelem,nrows,evtnrows, longnull,badpixnrows;
	double *evttime,*cztseccnt,*finalevttime,*finalcztseccnt,doublenull,*btitstart,*btitstop;
	int j,*pi,*finalpi;
	unsigned short *cztntick,*veto,*finalcztntick,*finalveto,*pha,*finalpha;
	unsigned char *detid,*pixid,*detx,*dety,*alpha,*finaldetid,*finalpixid,*finaldetx,*finaldety,*finalalpha,*bpdetx,*bpdety,*bpdetid,*bppixid,*bppix_flag,*btidetx,*btidety,*btidetid,*btipixid;
	float *energy,floatnull,*finalenergy;
	unsigned char *caldb_detid,*caldb_pixid,*pix_flag,*final_pix_flag,bytenull,*caldb_detx,*caldb_dety;
	int quadstart=0,quadend=4;
	char outlogfile[1000];
	int *detpixid,qid,bppix_flag_count=0;
	char extname[20]; 

		evttime  = (double*)malloc(size * sizeof(double));
		cztseccnt  = (double*)malloc(size * sizeof(double));
		pha  = (unsigned short*)malloc(size * sizeof(unsigned short));
		cztntick  = (unsigned short*)malloc(size * sizeof(unsigned short));
		veto  = (unsigned short*)malloc(size * sizeof(unsigned short));
		detid  = (unsigned char*)malloc(size * sizeof(unsigned char));
		pixid  =(unsigned char*) malloc(size * sizeof(unsigned char));
		detx  = (unsigned char*)malloc(size * sizeof(unsigned char));
		dety  = (unsigned char*)malloc(size * sizeof(unsigned char));
		alpha = (unsigned char*)malloc(size * sizeof(unsigned char));
		energy = (float*)malloc(sizeof(float)*size);
		pi = (int*)malloc(sizeof(int)*size);

		finalevttime  = (double*)malloc(size * sizeof(double));
		finalcztseccnt  = (double*)malloc(size * sizeof(double));
		finalpha  = (unsigned short*)malloc(size * sizeof(unsigned short));
		finalcztntick  = (unsigned short*)malloc(size * sizeof(unsigned short));
		finalveto  = (unsigned short*)malloc(size * sizeof(unsigned short));
		finaldetid  = (unsigned char*)malloc(size * sizeof(unsigned char));
		finalpixid  =(unsigned char*) malloc(size * sizeof(unsigned char));
		finaldetx  = (unsigned char*)malloc(size * sizeof(unsigned char));
		finaldety  = (unsigned char*)malloc(size * sizeof(unsigned char));
		finalalpha = (unsigned char*)malloc(size * sizeof(unsigned char));
		finalenergy = (float*)malloc(sizeof(float)*size);
		finalpi = (int*)malloc(sizeof(int)*size);

		bpdetid  = (unsigned char*)malloc(size * sizeof(unsigned char));
		bppixid  =(unsigned char*) malloc(size * sizeof(unsigned char));
		bpdetx  = (unsigned char*)malloc(size * sizeof(unsigned char));
		bpdety  = (unsigned char*)malloc(size * sizeof(unsigned char));
		bppix_flag  = (unsigned char*)malloc(size * sizeof(unsigned char));
		detpixid  = (int*)malloc(size * sizeof(int));

		btidetid  = (unsigned char*)malloc(size * sizeof(unsigned char));
		btipixid  =(unsigned char*) malloc(size * sizeof(unsigned char));
		btidetx  = (unsigned char*)malloc(size * sizeof(unsigned char));
		btidety  = (unsigned char*)malloc(size * sizeof(unsigned char));
		btitstart  = (double*)malloc(size * sizeof(double));
		btitstop  = (double*)malloc(size * sizeof(double));

	status = 0;
	hdunum = 2;
	FILE *thrfile;
	thrfile = fopen(thresholdfile, "r");

        if (thrfile == NULL) {
                printf("Couldn't open the threshold file.");
                return ;
        }

	const size_t label_size = 300;
	char* label = malloc(label_size);

	while (fgets(label, label_size, thrfile) != NULL)  
	{
		char *val;
		strtok_r (label, " ", &val);
		if (strcmp(label,"flick_energy_thresh") == 0)
			flick_energy_thresh=atoi(val);
		else if (strcmp(label,"sigma_thresh_1") == 0)
			sigma_thresh_1=atoi(val);
		else if (strcmp(label,"flick_count_thresh_1") == 0)
			flick_count_thresh_1=atoi(val);
		else if (strcmp(label,"sigma_thresh_2") == 0)
			sigma_thresh_2=atoi(val);
		else if (strcmp(label,"flick_count_thresh_2") == 0)
			flick_count_thresh_2=atoi(val);	                
	}
	free(label);
	char* tempevt = calloc(strlen(infile)+1, sizeof(char));
	strcpy(tempevt, infile);

	status=0;
	if ( fits_open_file(&evtfptr, infile, READONLY, &status) )
	{
	    	printf("Error in opening a file : %d",status);
	   	printerror( status );
	}

	status=0;
	if ( fits_open_file(&badpixfptr, inbadpixfile, READONLY, &status) )
	{
	    	printf("Error in opening a file : %d",status);
	   	printerror( status );
	}
	
	//create event and badpix file
	char *file = strtok(infile, ".");
	if(clobber==1)
	{ 
		remove(outfile);
		remove(outbadpixfile);
		remove(outbtifile);
        }

	createEventFile(outfile,tempevt);
	createBadpixFile(outbadpixfile);
	createBTIFile(outbtifile);
		
	sprintf(outlogfile, "%s_fc.log",file);
	remove(outlogfile);
	logfile=fopen(outlogfile,"a"); 

	fprintf(logfile,"Inputs for flickering pixel clean\n1.Event file : %s\n2.Bad pixel file : %s\n\n",tempevt,inbadpixfile);
	fprintf(logfile,"Important thresholds ....\n1.flick_energy_thresh = %d\n2.sigma_thresh_1 = %d\n3.flick_count_thresh_1 = %d\n4.sigma_thresh_2 : %d\n5.flick_count_thresh_2 : %d\n\n",flick_energy_thresh,sigma_thresh_1,flick_count_thresh_1,sigma_thresh_2,flick_count_thresh_2);

	frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;
	floatnull=0.0;
	
	
	double **pixel_exposure, **pixel_badtime;
	
	pixel_exposure=(double**)malloc(sizeof(double*)*4);
	pixel_badtime=(double**)malloc(sizeof(double*)*4);
	for(i=0;i<4;i++) 
	{
		pixel_exposure[i]=(double*)malloc(sizeof(double)*4096);
		pixel_badtime [i]=(double*)malloc(sizeof(double)*4096);
	}
    	for(qid=0;qid<4;qid++)
    	{
		for(j=0;j<4096;j++) {
			//~ pixel_exposure[qid][j]=1.0;
			pixel_badtime[qid][j]=0.0;}
		
	}
	
	for(qid=quadstart; qid<quadend; qid++)
	{
		
		printf("Q%d started\n",qid);
		
		status=0;
		fprintf(logfile,"Quad %d is processing...\n",qid);
		fits_movabs_hdu(badpixfptr, qid+2, NULL, &status);
		fits_get_num_rows(badpixfptr, &badpixnrows, &status);

		bpdetid  	= (unsigned char*)realloc(bpdetid,badpixnrows * sizeof(unsigned char));
		bppixid  	=(unsigned char*) realloc(bppixid,badpixnrows * sizeof(unsigned char));
		bpdetx 		= (unsigned char*)realloc(bpdetx,badpixnrows * sizeof(unsigned char));
		bpdety 		= (unsigned char*)realloc(bpdety,badpixnrows * sizeof(unsigned char));
		bppix_flag  = (unsigned char*)realloc(bppix_flag,badpixnrows * sizeof(unsigned char));

		fits_read_col(badpixfptr, TBYTE, 1, frow, felem,badpixnrows, &bytenull, bpdetid, &anynull, &status);  
		fits_read_col(badpixfptr, TBYTE, 2, frow, felem,badpixnrows, &bytenull, bppixid, &anynull, &status);        
		fits_read_col(badpixfptr, TBYTE, 3, frow, felem, badpixnrows, &bytenull, bpdetx, &anynull, &status);   
		fits_read_col(badpixfptr, TBYTE, 4, frow, felem,badpixnrows, &bytenull, bpdety,&anynull, &status);
		fits_read_col(badpixfptr, TBYTE, 5, frow, felem,badpixnrows, &bytenull, bppix_flag,&anynull, &status);
		
			
		status=0;
		fits_movabs_hdu(evtfptr, qid+2, NULL, &status);
		fits_get_num_rows(evtfptr, &evtnrows, &status);
		
		evttime  = (double*)realloc(evttime,evtnrows * sizeof(double));
		cztseccnt  = (double*)realloc(cztseccnt,evtnrows * sizeof(double));
		pha  = (unsigned short*)realloc(pha,evtnrows * sizeof(unsigned short));
		cztntick  = (unsigned short*)realloc(cztntick,evtnrows * sizeof(unsigned short));
		veto  = (unsigned short*)realloc(veto,evtnrows * sizeof(unsigned short));
		detid  = (unsigned char*)realloc(detid,evtnrows * sizeof(unsigned char));
		pixid  =(unsigned char*) realloc(pixid,evtnrows * sizeof(unsigned char));
		detx  = (unsigned char*)realloc(detx,evtnrows * sizeof(unsigned char));
		dety  = (unsigned char*)realloc(dety,evtnrows * sizeof(unsigned char));
		alpha = (unsigned char*)realloc(alpha,evtnrows * sizeof(unsigned char));
		energy = (float*)realloc(energy,sizeof(float)*evtnrows);
		pi = (int*)realloc(pi,sizeof(int)*evtnrows);
		detpixid  = (int*)realloc(detpixid,evtnrows * sizeof(int));


		finalevttime  = (double*)realloc(finalevttime,evtnrows * sizeof(double));
		finalcztseccnt  = (double*)realloc(finalcztseccnt,evtnrows * sizeof(double));
		finalpha  = (unsigned short*)realloc(finalpha,evtnrows * sizeof(unsigned short));
		finalcztntick  = (unsigned short*)realloc(finalcztntick,evtnrows * sizeof(unsigned short));
		finalveto  = (unsigned short*)realloc(finalveto,evtnrows * sizeof(unsigned short));
		finaldetid  = (unsigned char*)realloc(finaldetid,evtnrows * sizeof(unsigned char));
		finalpixid  =(unsigned char*) realloc(finalpixid,evtnrows * sizeof(unsigned char));
		finaldetx  = (unsigned char*)realloc(finaldetx,evtnrows * sizeof(unsigned char));
		finaldety  = (unsigned char*)realloc(finaldety,evtnrows * sizeof(unsigned char));
		finalalpha = (unsigned char*)realloc(finalalpha,evtnrows * sizeof(unsigned char));
		finalenergy = (float*)realloc(finalenergy,sizeof(float)*evtnrows);
		finalpi = (int*)realloc(finalpi,sizeof(int)*evtnrows);


		fits_read_col(evtfptr, TDOUBLE, 1, frow, felem, evtnrows, &doublenull, evttime,&anynull, &status);
		fits_read_col(evtfptr, TDOUBLE, 2, frow, felem, evtnrows, &doublenull, cztseccnt, &anynull, &status);        
		fits_read_col(evtfptr, TUSHORT, 3, frow, felem, evtnrows, &intnull, cztntick, &anynull, &status); 
		fits_read_col(evtfptr, TUSHORT, 4, frow, felem, evtnrows, &intnull, pha,&anynull, &status);
		fits_read_col(evtfptr, TBYTE, 5, frow, felem,evtnrows, &bytenull, detid, &anynull, &status);  
		fits_read_col(evtfptr, TBYTE, 6, frow, felem,evtnrows, &bytenull, pixid, &anynull, &status);        
		fits_read_col(evtfptr, TBYTE, 7, frow, felem, evtnrows, &bytenull, detx, &anynull, &status);   
		fits_read_col(evtfptr, TBYTE, 8, frow, felem,evtnrows, &bytenull, dety,&anynull, &status);
		fits_read_col(evtfptr, TUSHORT, 9, frow, felem,evtnrows, &intnull, veto,&anynull, &status);
		fits_read_col(evtfptr, TBYTE, 10, frow, felem, evtnrows, &bytenull, alpha,&anynull, &status); 
		fits_read_col(evtfptr, TINT, 11, frow, felem, evtnrows, &intnull, pi,&anynull, &status);
		fits_read_col(evtfptr, TFLOAT, 12, frow, felem, evtnrows, &floatnull, energy,&anynull, &status);


	
		

		for(i=0;i<evtnrows;i++)
		{
			detpixid[i]=detid[i]*256+pixid[i];
		}
		
		int numbin_100s=(int)(evttime[evtnrows-1]-evttime[0])/100;
		//printf("%lf\t%lf\t%d\n",evttime[evtnrows-1],evttime[0],evtnrows);
		int numbin_10s=(int)(evttime[evtnrows-1]-evttime[0])/10;
		int numbin_1s=(int)(evttime[evtnrows-1]-evttime[0]);
		//printf("%d\t%d\t%d\n",numbin_100s,numbin_10s,numbin_1s);
			
		double **pix_lc_lt, **pix_lc_gt;
		pix_lc_lt=(double**)malloc(sizeof(double*)*4096);
		pix_lc_gt=(double**)malloc(sizeof(double*)*4096);
		for(i=0;i<4096;i++) 
		{
			pix_lc_lt[i]=(double*)malloc(sizeof(double)*(numbin_100s+1));
			pix_lc_gt[i]=(double*)malloc(sizeof(double)*(numbin_100s+1));
		}
			
		
		double *lc_lt,*lc_gt,*flivetime;
		
		lc_lt=(double*)malloc(sizeof(double)*(numbin_100s+1));
		lc_gt=(double*)malloc(sizeof(double)*(numbin_100s+1));
		flivetime=(double*)malloc(sizeof(double)*(numbin_100s+1));
		
		//double pix_lc_lt[4096][numbin_100s+1],pix_lc_gt[4096][numbin_100s+1];
		//printf("TTTT\n");
		//double lc_lt[numbin_100s+1],lc_gt[numbin_100s+1],flivetime[numbin_100s+1];
		
		
		//float pix_mean_lt,pix_sd_lt,pix_mean_gt,pix_sd_gt;
		double pix_mean_lt,pix_sd_lt,pix_mean_gt,pix_sd_gt;
		int ignore_pix_lt,ignore_pix_lt_2,ignore_pix_gt,ignore_pix_gt_2,*pix_flag,*evt_flag,flick_all_1s=0,flick_all_10s=0,flick_thresh_10s,flick_thresh_1s;
		double *flick_pix_time_1s,*flick_pix_time_10s;
		
		int *flick_pix_1s,*flick_pix_10s,*flick_detpixid_1s,*flick_detpixid_10s;
		flick_pix_1s=(int*)malloc(sizeof(int)*(4096));
		flick_pix_10s=(int*)malloc(sizeof(int)*(4096));
		flick_detpixid_1s=(int*)malloc(sizeof(int)*(4096));
		flick_detpixid_10s=(int*)malloc(sizeof(int)*(4096));
		
		
		pix_flag=(int*)calloc(sizeof(int),4096);
		evt_flag=(int*)calloc(sizeof(int),evtnrows);

		
		//flick_pix_time_1s=(double*)calloc(sizeof(double),numbin_1s+1);
		//flick_pix_time_10s=(double*)calloc(sizeof(double),numbin_10s+1);
		
		flick_pix_time_1s=(double*)calloc(sizeof(double),4096);
		flick_pix_time_10s=(double*)calloc(sizeof(double),4096);
		
		
		//flick_pix_1s=(int*)calloc(sizeof(int),numbin_1s+1);
		//flick_pix_10s=(int*)calloc(sizeof(int),numbin_10s+1);

		btidetid  = (unsigned char*)realloc(btidetid,4096 * sizeof(unsigned char));
		btipixid  =(unsigned char*) realloc(btipixid,4096 * sizeof(unsigned char));
		btidetx  = (unsigned char*)realloc(btidetx,4096 * sizeof(unsigned char));
		btidety  = (unsigned char*)realloc(btidety,4096 * sizeof(unsigned char));
		btitstart  = (double*)realloc(btitstart,4096 * sizeof(double));
		btitstop  = (double*)realloc(btitstop,4096 * sizeof(double));



		for(i=0;i<numbin_100s+1;i++)
		{
			lc_lt[i]=0.0;
			lc_gt[i]=0.0;
			flivetime[i] = 0.0;
			for(j=0;j<4096;j++)
			{
				pix_lc_lt[j][i]=0.0;
				pix_lc_gt[j][i]=0.0;
			}
		}
		
		long gtinrows=0;
		double tstart,tstop;
		double *gtitstart,*gtitstop;
		fits_read_key(evtfptr,TDOUBLE,"TSTART",&tstart,NULL, &status);
		fits_read_key(evtfptr,TDOUBLE,"TSTOP",&tstop,NULL, &status);
		
		
		//~ fits_movabs_hdu(evtfptr, qid+10, &hdutype, &status);
		sprintf(extname, "Q%d_GTI",qid);
		fits_movnam_hdu(evtfptr, BINARY_TBL, extname, 0, &status);
		
		fits_get_num_rows(evtfptr, &gtinrows, &status);
		
		gtitstart=(double*)malloc(sizeof(double)*gtinrows);
		gtitstop=(double*)malloc(sizeof(double)*gtinrows);
		
		frow      = 1;
		felem     = 1;
		doublenull = 0.;
		intnull = 0;
		bytenull = 0;
		floatnull=0.0;

		fits_read_col(evtfptr, TDOUBLE, 1, frow, felem, gtinrows, &doublenull, gtitstart,&anynull, &status);
		fits_read_col(evtfptr, TDOUBLE, 2, frow, felem, gtinrows, &doublenull, gtitstop, &anynull, &status);
		
		
		int time_gti_1,time_gti_2;
		//printf("True---2\n");
		for(i=0;i<gtinrows;i++)
		{
			time_gti_1 = (int)(gtitstart[i]-evttime[0])/100;
			time_gti_2 = (int)(gtitstop[i]-evttime[0])/100;
			
			//printf("%d %d \n",time_gti_1,time_gti_2);
			
			if(time_gti_1==time_gti_2)
			{
				flivetime[time_gti_1]=flivetime[time_gti_1]+(gtitstop[i]-gtitstart[i]);
				//printf("%lf \n",flivetime[time_gti_1]);
			}
			else if(time_gti_2>time_gti_1)
			{
				flivetime[time_gti_1] = flivetime[time_gti_1] + (evttime[0]+100.0*(double)(time_gti_1+1)-gtitstart[i]);
				flivetime[time_gti_2] = flivetime[time_gti_2] + (gtitstop[i]-(evttime[0]+100.0*(double)(time_gti_2)));
				
				//printf("%lf %lf\n",flivetime[time_gti_1],flivetime[time_gti_2]);
				
				for(j=time_gti_1+1;j<time_gti_2;j++)
				{
					flivetime[j]=flivetime[j]+100.0;
					//printf("%lf \n",flivetime[j]);
				}
			}
			
		}
		
	
		for(i=0;i<evtnrows;i++)
		{
			int tempnumbin_100s=(int)(evttime[i]-evttime[0])/100;
			//printf("%d\n",tempnumbin_100s);
					
			if(energy[i]<flick_energy_thresh)
			{	
				lc_lt[tempnumbin_100s]++;
				pix_lc_lt[detpixid[i]][tempnumbin_100s]++;
			}
			if(energy[i]>=flick_energy_thresh)
			{
				lc_gt[tempnumbin_100s]++;
				pix_lc_gt[detpixid[i]][tempnumbin_100s]++;
			}

		}
		
		int flick_pix_counter_1s=0,flick_pix_counter_10s=0,k=0;
		double tot_exposure;
		
		fits_read_key(evtfptr,TDOUBLE,"EXPOSURE",&tot_exposure,NULL, &status);
		//~ fits_movabs_hdu(evtfptr,13, NULL, &status);
		fits_movnam_hdu(evtfptr, BINARY_TBL, "EXPOSURE", 0, &status);
		fits_read_col(evtfptr, TDOUBLE, qid+1, frow, felem, 4096, &doublenull, pixel_exposure[qid],&anynull, &status);
		

		
	
		for(j=0;j<numbin_100s;j++)
		{
			//printf("%lf %lf\n",lc_lt[j],lc_gt[j]);
			//printf("%lf \n",flivetime[j]);
			if(flivetime[j]<20.0){flivetime[j]=0.0;}
			
			//printf("%lf\n",flivetime[j]);
			
			flivetime[j]=flivetime[j]/100.0;
			//printf("%lf \n",flivetime[j]);
		}
		
		
		for(i=0;i<4096;i++)
		{	
			//printf("%f\n",pixel_exposure[qid][i]);
			pix_mean_lt=0.0;
			pix_sd_lt=0.0;
			ignore_pix_lt=0;
			ignore_pix_lt_2=0;
			ignore_pix_gt_2=0;
			pix_mean_gt=0.0;
			pix_sd_gt=0.0;
			ignore_pix_gt=0;
			bppix_flag_count=0;
			
			int ig_count=0;
			for(j=0;j<numbin_100s;j++)
			{
				if(flivetime[j]==0.0)
				{
					ig_count++;
					continue;
				}
				//pix_lc_lt[i][j] = pix_lc_lt[i][j]/flivetime[j];
				//pix_lc_gt[i][j] = pix_lc_gt[i][j]/flivetime[j];
				//lc_lt[j] = lc_lt[j]/flivetime[j];
				//lc_gt[j] = lc_gt[j]/flivetime[j];
				/*if(pix_lc_lt[i][j]==0.0)
				{
					printf("TRUE   %d %d\n",i,j);
				}*/
				
				//printf("%lf %lf %lf %lf\n",pix_lc_lt[i][j],pix_lc_gt[i][j],lc_lt[j],lc_gt[j]);
				
				pix_lc_lt[i][j]=pix_lc_lt[i][j]/lc_lt[j];
				pix_lc_gt[i][j]=pix_lc_gt[i][j]/lc_gt[j];
				pix_mean_lt=pix_mean_lt+pix_lc_lt[i][j];
				pix_mean_gt=pix_mean_gt+pix_lc_gt[i][j];
				
				
			}
		

			pix_mean_lt=pix_mean_lt/(double)numbin_100s;
			pix_mean_gt=pix_mean_gt/(double)numbin_100s;

			//printf("mean = %lf %lf\n",pix_mean_lt,pix_mean_gt);
			
			
			for(j=0;j<numbin_100s;j++)
			{
				if(flivetime[j]==0.0)
				{
					continue;
				}
				
				pix_sd_lt=pix_sd_lt+pow((pix_lc_lt[i][j]-pix_mean_lt),2);
				pix_sd_gt=pix_sd_gt+pow((pix_lc_gt[i][j]-pix_mean_gt),2);
			}
			pix_sd_lt = sqrt(pix_sd_lt/((double)numbin_100s-(double)ig_count));
			pix_sd_gt = sqrt(pix_sd_gt/((double)numbin_100s-(double)ig_count));
			
			//printf("sd = %lf %lf\n",pix_sd_lt,pix_sd_gt);
		
			
			
			
			for(j=0;j<numbin_100s;j++)
			{
				if(flivetime[j]==0.0)
				{
					continue;
				}
				
	
				if(pix_lc_lt[i][j] > pix_mean_lt+sigma_thresh_1*pix_sd_lt)
				{
					//printf("LOW 1 %d %d %lf %lf %lf\n",i,j,pix_lc_lt[i][j],pix_mean_lt,pix_sd_lt);
					ignore_pix_lt++;
				}
				if(pix_lc_gt[i][j]> pix_mean_gt+sigma_thresh_1*pix_sd_gt)
				{
					//printf("HIGH 1 %d %d %lf %lf %lf\n",i,j,pix_lc_gt[i][j],pix_mean_lt,pix_sd_lt);
					ignore_pix_gt++;
				}
				if(pix_lc_lt[i][j] > pix_mean_lt+sigma_thresh_2*pix_sd_lt)
				{
					//printf("LOW 2 %d %d %lf %lf %lf\n",i,j,pix_lc_lt[i][j],pix_mean_lt,pix_sd_lt);
					ignore_pix_lt_2++;
				}
				if(pix_lc_gt[i][j]> pix_mean_gt+sigma_thresh_2*pix_sd_gt)
				{
					//printf("HIGH 2 %d %d %lf %lf %lf\n",i,j,pix_lc_gt[i][j],pix_mean_lt,pix_sd_lt);
					ignore_pix_gt_2++;
				}
			}
			
			
			//printf("True---3\n");
			/*
			if(ignore_pix_lt >= flick_count_thresh_1 || ignore_pix_gt >= flick_count_thresh_1)
			{
				//printf("1 FLAG\n");
				pix_flag[i]=1;
				bppix_flag[i]=2;
				bppix_flag_count++;	
				pixel_exposure[qid][i]=0.0;
				
			}
			 */
			if(ignore_pix_lt_2 >= flick_count_thresh_2 || ignore_pix_gt_2 >= flick_count_thresh_2)
			{
				//printf("2 FLAG-----  %d %d\n",ignore_pix_lt_2,i);
				pix_flag[i]=1;
				bppix_flag[i]=2;
				bppix_flag_count++;	
				pixel_exposure[qid][i]=0.0;
			}

		 }
		//printf("Count of flickering pixels : %d\n",bppix_flag_count);
		 int flag;
		 for(i=0;i<evtnrows;i++)
		 {
			if(pix_flag[detid[i]*256+pixid[i]]==1)
			{
				//printf("row\n");
				evt_flag[i]=1;	
	
			}
		}
		
		int valid_pix_counter=0, valid_pix_counter_temp=0;
		
		for(i=0;i<badpixnrows;i++)
		{
			if(bppix_flag[i]==0 || bppix_flag[i]==1)
			{
				valid_pix_counter++;
			}
			if(bppix_flag[i]==2)
			{
				valid_pix_counter_temp++;
			}
			
		}
		
		printf("No of valid pixels = %d   No of flickering pixels excluded = %d\n",valid_pix_counter,valid_pix_counter_temp);
			


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		double bintime;
		int tempind1=0,tempind2,jj=0;
		
		for(bintime=evttime[0];bintime< evttime[evtnrows-1];bintime+=1.0)
		{
			//printf("%lf\n",bintime);

			flick_all_1s = 0;
			for(j=0;j<4096;j++){flick_detpixid_1s[j]=0;}
			for(j=tempind1;j<evtnrows;j++)
			{
				if(evt_flag[j]==1)
				{
					continue;
				}
				
				
				if(bintime<=evttime[j] && bintime+1.0 > evttime[j])
				{	
					flick_detpixid_1s[detpixid[j]]++;
					flick_all_1s++;
					tempind1=j;
					
				}
				if(evttime[j]-bintime >= 3.0){break;}
			}
			//if(qid==1){printf("%d\t%d\n",flick_detpixid_1s[230], flick_all_1s++);}
			flick_pix_counter_1s=0;
			flick_thresh_1s = pixThresholdFix(flick_all_1s,valid_pix_counter);
			//if(bintime>211223974.15 && bintime<211223990.15)
			//printf("flick_thresh_1s : %d\n",flick_thresh_1s);

			for(jj=0;jj<4096;jj++)
			{

				//if(qid==1 && j==230){ printf("%d\t%d\n",flick_detpixid_1s[230], flick_thresh_1s);}
				if(flick_detpixid_1s[jj]>flick_thresh_1s)
				{
					
					pixel_badtime[qid][jj] = pixel_badtime[qid][jj]+1.0; 
					//flick_pix_time_1s[flick_pix_counter_1s]=bintime;
					flick_pix_1s[flick_pix_counter_1s]=jj;
					//if(qid==1 && j==230){printf("1    %d\t%lf\n",j,flick_pix_time_1s[flick_pix_counter_1s]);}
					flick_pix_counter_1s++;
					
				}


			}
			//printf("%lf\t%d\n",bintime,flick_pix_counter_1s);
			//~ if(flick_pix_counter_1s>0) printf("%d\t%d\t%d\n",flick_pix_counter_1s,flick_thresh_1s,qid);
			for(k=0;k<flick_pix_counter_1s;k++)
			{
				for(j=tempind2;j<evtnrows;j++)
				{
					//if(qid==1 && detpixid[j]==230 && flick_pix_1s[k]==230){ printf("True %d\t%lf\n",flick_pix_1s[k],flick_pix_time_1s[k]);}
					
					if((evttime[j] >= bintime && evttime[j] <= bintime+1.0) && detpixid[j] == flick_pix_1s[k])
					{
						//if(qid==1 && detpixid[j]==230 && flick_pix_1s[k]==230){ printf("True %d\t%lf\t%lf\t%d\n",flick_pix_1s[k],bintime,evttime[j],j);}
						//printf("1s row\n");
						evt_flag[j]=1;
						tempind2=j;
						//break;

					}
					if(evttime[j]-bintime >= 3.0){break;}
				}
		
				btipixid[k]=flick_pix_1s[k]/256;
				btidetid[k]=flick_pix_1s[k]%256;
				btidetx[k]=((btidetid[k]%4)*16)+(btipixid[k]%16);
				btidety[k]=((btidetid[k]/4)*16)+(btipixid[k]/16);
			
				if(qid==0 || qid==3)
					btidety[k]=63-btidety[k];
				else
					btidetx[k]=63-btidetx[k];

				btitstart[k]=bintime;
				btitstop[k]=bintime+1.0;
			}
			
			writeBTI(btidetid,btipixid,btidetx,btidety,btitstart,btitstop,flick_pix_counter_1s,outbtifile,qid+2);


			
		}
		
		//printf("%lf\n",bintime);
		
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		
		
		/*
		tempind1=0;tempind2=0,jj=0;
		for(bintime=evttime[0];bintime<evttime[evtnrows-1];bintime+=10.0)
		{
			

			flick_all_10s = 0;
			for(j=0;j<4096;j++){flick_detpixid_10s[j]=0;}
			for(j=tempind1;j<evtnrows;j++)
			{
				
				if(evt_flag[j]==1)
				{
					continue;
				}
				
				if(bintime<=evttime[j] && bintime+10.0 > evttime[j])
				{	
					flick_detpixid_10s[detpixid[j]]++;
					
					
					flick_all_10s++;
					
					tempind1=j;
					
				}
				if(evttime[j]-bintime >= 15.0){break;}
			}
			flick_pix_counter_10s=0;
			flick_thresh_10s = pixThresholdFix(flick_all_10s,valid_pix_counter);
			//if(bintime>211223974.15 && bintime<211223990.15)
			//printf("flick_thresh_10s : %d\n",flick_thresh_10s);

			for(jj=0;jj<4096;jj++)
			{
				if(flick_detpixid_10s[jj]>flick_thresh_10s)
				{
					
					pixel_badtime[qid][jj] = pixel_badtime[qid][jj]+10.0;
					flick_pix_time_10s[flick_pix_counter_10s]=bintime;
					flick_pix_10s[flick_pix_counter_10s]=jj;
					//if(qid==1 && i/256==0 && i%256==230){printf("10    %d\t%lf\n",j,flick_pix_time_10s[flick_pix_counter_10s]);}
					flick_pix_counter_10s++;
					
				}
			}
			
			//if(flick_pix_counter_10s>0) printf("%d\t%d\t10\t%d\n",flick_pix_counter_10s,flick_thresh_10s,qid);
			for(k=0;k<flick_pix_counter_10s;k++)
			{
				for(j=tempind2;j<evtnrows;j++)
				{
					
					if((evttime[j]>=flick_pix_time_10s[k] && evttime[j]<flick_pix_time_10s[k]+10.0) && detpixid[j]== flick_pix_10s[k])
					{
						//printf("10s row\n");
						evt_flag[j]=1;
						tempind2=j;
						//break;

					}
					if(evttime[j]-flick_pix_time_10s[k] >= 15.0){break;}
				}
				btipixid[k]=flick_pix_10s[k]/256;
				btidetid[k]=flick_pix_10s[k]%256;
				btidetx[k]=((btidetid[k]%4)*16)+(btipixid[k]%16);
				btidety[k]=((btidetid[k]/4)*16)+(btipixid[k]/16);
			
				if(qid==0 || qid==3)
					btidety[k]=63-btidety[k];
				else
					btidetx[k]=63-btidetx[k];

				btitstart[k]=flick_pix_time_10s[k];
				btitstop[k]=flick_pix_time_10s[k]+10.0;
			}
			
			writeBTI(btidetid,btipixid,btidetx,btidety,btitstart,btitstop,flick_pix_counter_10s,outbtifile,qid+2);
				
		} */
		
		
		
		int l=0;
		fprintf(logfile,"Events removed after flickering pixels detection\nIndex\tTime\tDETX\tDETY\n");
		for(i=0;i<evtnrows;i++)
			{
				
				if(evt_flag[i]==0)
				{
					finalevttime[l]=evttime[i];
					finalcztseccnt[l]  = cztseccnt[i];
					finalpha[l]  = pha[i];
					finalcztntick[l]  = cztntick[i];
					finalveto[l]  = veto[i];
					finaldetid[l]  = detid[i];
					finalpixid[l]  = pixid[i];
					finaldetx[l]  = detx[i];
					finaldety[l]  = dety[i];
					finalalpha[l] = alpha[i];
					finalenergy[l] = energy[i];
					finalpi[l] = pi[i];
					l++;
				}
				else if(evt_flag[i]==1)
					fprintf(logfile,"%d\t%f\t%u\t%u\n",i,evttime[i],detx[i],dety[i]);
				
		
			}
				
			
		for(i=0;i<4096;i++)
		{
			//~ printf("%f\n",pixel_exposure[qid][i]);
			if(pixel_exposure[qid][i]==0.0){continue;}
			
			
			else
			{
				pixel_exposure[qid][i] = (tot_exposure*pixel_exposure[qid][i]-pixel_badtime[qid][i])/tot_exposure;
			}
			
				
		}	
			
		
		modifyExposure(outfile,pixel_exposure,qid);
	
		
		writeBadpix(bpdetid,bppixid,bpdetx,bpdety,bppix_flag,outbadpixfile,qid+2);
		writeEvent(finalevttime,finalcztseccnt,finalcztntick,finalpha,finaldetid,finalpixid,finaldetx,finaldety,finalveto,finalalpha,finalpi,finalenergy,outfile,l,qid+2,tot_exposure);
		printf("Q%d end\n",qid);
		for(i=0;i<4096;i++) 
		{
			free(pix_lc_lt[i]);free(pix_lc_gt[i]);
		}
		free(pix_lc_lt);free(pix_lc_gt);	
			
		free(pix_flag);free(evt_flag);
		free(flick_pix_time_1s);free(flick_pix_time_10s);
		free(flick_pix_1s);free(flick_pix_10s);free(flick_detpixid_1s);free(flick_detpixid_10s);
		free(gtitstart);free(gtitstop);
		free(lc_lt);free(lc_gt);free(flivetime);
		
	}
	//~ modifycommonGTI(outfile);
	
	modifyEventHeaderParams(outfile);
	
	writeBadpixExtension(outbadpixfile);
	
	
	
	free(evttime);free(cztseccnt);free(finalevttime);free(finalcztseccnt);
	free(cztntick);free(veto);free(finalcztntick);free(finalveto);free(pha);free(finalpha);
	free(detid);free(pixid);free(detx);free(dety);free(alpha);free(finaldetid);free(finalpixid);
	free(finaldetx);free(finaldety);free(finalalpha);
	free(energy);free(finalenergy);
	free(pi);free(finalpi);
	free(bpdetx);free(bpdety);free(bpdetid);free(bppixid);free(bppix_flag);
	free(btidetx);free(btidety);free(btidetid);free(btipixid);free(btitstart);free(btitstop);
	free(detpixid);

	for(i=0;i<4;i++) 
	{
		free(pixel_exposure[i]);
		free(pixel_badtime[i]);
	}
	free(pixel_exposure);
	free(pixel_badtime);
	free(tempevt);
	fclose(logfile);
	fclose(thrfile);

	//if ( fits_close_file(caldbfptr, &status) )       
	//        printerror( status );
	if ( fits_close_file(evtfptr, &status) )       
	        printerror( status );
	if ( fits_close_file(badpixfptr, &status) )       
	        printerror( status );
	printf("\nFLICKERING PIXEL REDUCTION COMPLETED SUCCESSFULLY.\n");
	
	return;
}
void writeBadpix(unsigned char *detid,unsigned char *pixid,unsigned char *pixx,unsigned char *pixy,unsigned char *pix_flag,char *outputfile,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	status=0;
	int writesize=4096;
	printf("Moving to hdu :%d\n",hdunum);
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	frow      = 1;//firstindex;
	felem     = 1;
	
		
        fits_write_col(fptrOut, TBYTE, 1, frow, felem, writesize,detid,&status);
        fits_write_col(fptrOut, TBYTE, 2,frow, felem, writesize, pixid,&status);
	fits_write_col(fptrOut, TBYTE, 3,frow, felem, writesize, pixx,&status);
        fits_write_col(fptrOut, TBYTE, 4,frow, felem, writesize, pixy,&status);
        fits_write_col(fptrOut, TBYTE, 5,frow, felem, writesize,pix_flag,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

int pixThresholdFix(int countrate, int num_pix)
{
	double mean=(double)countrate/(4096.0-(double)num_pix),chance;//4096-bad+dead pix count
	//printf("Mean = %lf\n",mean);
	int x=2,i;
	int factx;
	while(1)
	{
		factx=1;
		for(i=1;i<=x;i++){factx=factx*i;}
		chance = pow(2.718,-mean)*pow(mean,(double)x)/(double)factx*(4096.0-(double)num_pix);//4096-bad+dead pix count 
		if(chance<1.0){break;}
		x++;
	}

	return x;

}

//creating event file 
void createEventFile(char *outputfile,char *eventfile)
{
	
	fitsfile *fptrOut,*fptrevt;      
	int status, hdutype,tfields=12,i,hdunum=2;
	long frow, felem;
	int mjdrefi=55197,mjdreff=0,equinox=2000;
	float ra_pnt,dec_pnt;
	double timedel,telapse,doublenull;
	char object[20],obs_id[20],obs_mode[20],date_obs[20],time_obs[20],date_end[20],time_end[20],date[20],creator[20],filename[70],checksum[20],datasum[20],chksumcomm[50],datasumcomm[50];
      
	char extname[20]; 
	          
	
	char *ttype[] = { "TIME", "CZTSECCNT","CZTNTICK","PHA","DetID","pixID","DETX","DETY","veto","alpha","PI","ENERGY"};
	char *tform[] = { "D","D","I","I","B","B","B","B","I","B","I","E"};
	char *tunit[] = {"s","s","micro-sec","counts","","","","","counts","counts","",""};
       
	status=0;
        if (fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

	//if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	//        printerror( status );

	status=0;
	
	if ( fits_open_file(&fptrevt,eventfile, READONLY, &status) ) 
	         printerror( status );

	// Copy Primary
        if(fits_movabs_hdu(fptrevt, 1, 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	
    	for(i=0;i<4;i++)
    	{      
		status=0;
		sprintf(extname, "Q%d",i);
		strcpy(chksumcomm,"HDU checksum updated ");
		strcpy(datasumcomm,"Data unit checksum updated ");
		
		if(fits_movabs_hdu(fptrevt,i+hdunum, 0, &status)) 
			printerror( status );  

		
		fits_read_key(fptrevt,TSTRING,"OBJECT",&object,NULL, &status);
		fits_read_key(fptrevt,TFLOAT,"RA_PNT",&ra_pnt,NULL, &status);
		fits_read_key(fptrevt,TFLOAT,"DEC_PNT",&dec_pnt,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"OBS_ID",&obs_id,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"OBS_MODE",&obs_mode,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"DATE-OBS",&date_obs,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"TIME-OBS",&time_obs,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"DATE-END",&date_end,NULL, &status);		
		fits_read_key(fptrevt,TSTRING,"TIME-END",&time_end,NULL, &status);
		fits_read_key(fptrevt,TDOUBLE,"TIMEDEL",&timedel,NULL, &status);
		fits_read_key(fptrevt,TDOUBLE,"TELAPSE",&telapse,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"DATE",&date,NULL , &status);	
		fits_read_key(fptrevt,TSTRING,"FILENAME",&filename,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"CHECKSUM",&checksum,NULL, &status);
		fits_read_key(fptrevt,TSTRING,"DATASUM",&datasum,NULL, &status);
		//printf("DATE : %s\n",)

		status=0;

		if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        	printerror( status );
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields, ttype, tform,tunit, extname, &status) )
			printerror( status );
		if ( fits_movabs_hdu(fptrOut, i+hdunum, &hdutype, &status) ) 
	         	printerror( status );

		fits_write_key(fptrOut,TSTRING,"OBJECT",object,"Target name", &status);
		fits_write_key(fptrOut,TFLOAT,"RA_PNT",&ra_pnt,"Nominal pointing RA", &status);
		fits_write_key(fptrOut,TFLOAT,"DEC_PNT",&dec_pnt,"Nominal pointing DEC", &status);
		fits_write_key(fptrOut,TSTRING,"OBS_ID",obs_id,"Observation ID", &status);
		fits_write_key(fptrOut,TSTRING,"OBS_MODE",obs_mode,NULL, &status);
		fits_write_key(fptrOut,TSTRING,"DATE-OBS",date_obs,"Start date of observation", &status);
		fits_write_key(fptrOut,TSTRING,"TIME-OBS",time_obs,"Start time of observation", &status);
		fits_write_key(fptrOut,TSTRING,"DATE-END",date_end,"End date of observation", &status);		
		fits_write_key(fptrOut,TSTRING,"TIME-END",time_end,"End time of observation", &status);
		fits_write_key(fptrOut,TDOUBLE,"TIMEDEL",&timedel,"Time resolution", &status);
		fits_write_key(fptrOut,TDOUBLE,"TELAPSE",&telapse,"Elapsed time", &status);	
		fits_write_key(fptrOut,TSTRING,"DATE",date,"File creation date(YY-MM-DDThh:mm:ss UT)" , &status);	
		fits_write_key(fptrOut,TSTRING,"CREATOR","dataselection","Module that created this file", &status);
		fits_write_key(fptrOut,TSTRING,"FILENAME",filename,NULL, &status);
		fits_write_key(fptrOut,TSTRING,"CHECKSUM",checksum,strcat(chksumcomm,date), &status);
		fits_write_key(fptrOut,TSTRING,"DATASUM",datasum,strcat(datasumcomm,date), &status);

		fits_write_key(fptrOut, TSTRING,"EXTNAME", extname,"Name of this binary table extension", &status);
		fits_write_key(fptrOut, TINT,"QUADID",&i,"Quadrant Number", &status);
		fits_write_key(fptrOut, TSTRING,"MISSION", "ASTROSAT","Name of the mission/satellite", &status);
		fits_write_key(fptrOut, TSTRING,"TELESCOP","ASTROSAT","Name of the mission/satellite", &status);
		fits_write_key(fptrOut, TSTRING,"INSTRUME","CZTI","Name of the instrument/detector", &status);
		fits_write_key(fptrOut, TSTRING,"ORIGIN","CZTI POC","Source of FITS FILE", &status);
		fits_write_key(fptrOut, TINT,"MJDREFI", &mjdrefi,"MJDREF Integer part", &status);
		fits_write_key(fptrOut, TINT,"MJDREFF", &mjdreff,"MJDREF Fractional part", &status);
		fits_write_key(fptrOut, TSTRING,"TIMESYS","UTC","Time in UTC", &status);
		fits_write_key(fptrOut, TINT,"EQUINOX", &equinox,"J2000", &status);
		fits_write_key(fptrOut, TSTRING,"RADECSYS", "ICRS","Reference Frame", &status);
		fits_write_key(fptrOut, TSTRING,"TIMEUNIT","s","Time is in seconds", &status);
		
		if ( fits_close_file(fptrOut, &status) )       
	        	printerror( status );
	}
	
	status=0;
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        	printerror( status );
	// Copy VETOSPECTRUM
        if(fits_movnam_hdu(fptrevt, BINARY_TBL, "VETOSPECTRUM", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );
	
	// Copy SSM data
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "SSM Data", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	// Copy Temp
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "TEMP", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	// Copy GTI
	/*
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "GTI", 0, &status)) 
		printerror( status );         
	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );
	*/
	
	char *ttype_gti[] = { "START", "STOP"};
	char *tform_gti[] = { "1D","1D"};
	char *tunit_gti[] = {"sec","sec"};
	status=0;
	sprintf(extname,"GTI");
	int tfields_gti=2;
	

	
	//~ if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields_gti, ttype_gti, tform_gti,tunit_gti, extname, &status) )
			//~ printerror( status );
	//~ if(fits_movabs_hdu(fptrOut,9, 0, &status)) 
		//~ printerror( status );
		
	
	// Copy Q0_GTI
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q0_GTI", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	// Copy Q1_GTI
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q1_GTI", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	// Copy Q2_GTI
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q2_GTI", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	// Copy Q3_GTI
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q3_GTI", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	// Copy Exposure
	/*
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "EXPOSURE", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );
	
	*/
	
	
	int tfields_exp=4;
	
	char *ttype_exp[] = {"EXPOSURE_Q0","EXPOSURE_Q1","EXPOSURE_Q2","EXPOSURE_Q3"};
	char *tform_exp[] = { "D","D","D","D"};
	char *tunit_exp[] = {"","","",""};
	
	 
	sprintf(extname, "EXPOSURE");

	
	if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields_exp, ttype_exp, tform_exp,tunit_exp, extname, &status) )
	printerror( status );
		
	//~ if ( fits_movabs_hdu(fptrOut, 14, &hdutype, &status) ) 
	         	//~ printerror( status );
	if ( fits_movnam_hdu(fptrOut, BINARY_TBL, "EXPOSURE", 0, &status) ) 
	      	printerror( status );

	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	if ( fits_close_file(fptrevt, &status) )       
	        printerror( status );
	return;
	
	

}

void createBadpixFile(char *outputfile)
{
	fitsfile *fptrOut;      
	int status, hdutype;
	long frow, felem;
     	int tfields,i;
	char extname[20];
	
	tfields=5;
	char *ttype[] = {"DetID","PixID","DetX","DetY","PIX_FLAG"};
	char *tform[] = { "B","B","B","B","B"};
	char *tunit[] = {"","","","",""};
	        
	int hdunum=2;
	
	status=0;
	
        if (fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

	//if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	//         printerror( status );

    	for(i=0;i<4;i++)
    	{       
		sprintf(extname, "Q%d",i);
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields, ttype, tform,tunit, extname, &status) )
			printerror( status );
		if ( fits_movabs_hdu(fptrOut, i+hdunum, &hdutype, &status) ) 
	         	printerror( status );
	}
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

void createBTIFile(char *outputfile)
{
	fitsfile *fptrOut;      
	int status, hdutype;
	long frow, felem;
      	int tfields,i;
	char extname[20];
	
	tfields=6;
	char *ttype[] = {"DetID","PixID","DetX","DetY","TStart","TStop"};
	char *tform[] = { "B","B","B","B","D","D"};
	char *tunit[] = {"","","","","s","s"};
	        
	int hdunum=2;
	
	status=0;
	
        if (fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

//	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
//	         printerror( status );

    	for(i=0;i<4;i++)
    	{       
		sprintf(extname, "Q%d_BTI",i);
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields, ttype, tform,tunit, extname, &status) )
			printerror( status );
		if ( fits_movabs_hdu(fptrOut, i+hdunum, &hdutype, &status) ) 
	         	printerror( status );
	}
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum, double exposure)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	status=0;
	int tstarti,tstopi;
	double tstart,tstop,tstartf,tstopf;//,exposure;
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	tstart=evttime[0];
	tstop=evttime[writesize-1];

	tstarti=(int)tstart;
	tstopi=(int)tstop;

	tstartf=tstart-tstarti;
	tstopf=tstop-tstopi;
	//exposure=tstop-tstart;

	frow      = 1;
	felem     = 1;
	
	fits_write_key(fptrOut, TDOUBLE,"TSTART", &tstart,"Start time of observation", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP", &tstop,"Stop time of observation", &status);
	fits_write_key(fptrOut, TINT,"TSTARTI", &tstarti,"Start time of observation Integer part", &status);
	fits_write_key(fptrOut, TINT,"TSTOPI", &tstopi,"Stop time of observation Integer part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTARTF", &tstartf,"Start time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOPF", &tstopf,"Stop time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"EXPOSURE", &exposure,"Exposure time", &status);


    fits_write_col(fptrOut, TDOUBLE, 1, frow, felem, writesize, evttime,&status);
	fits_write_col(fptrOut, TDOUBLE, 2, frow, felem, writesize, cztseccnt,&status);
    fits_write_col(fptrOut, TUSHORT, 3,frow, felem, writesize, cztntick,&status);  
    fits_write_col(fptrOut, TUSHORT, 4,frow, felem, writesize, pha,&status);
	fits_write_col(fptrOut, TBYTE, 5, frow, felem, writesize,detid,&status);
	fits_write_col(fptrOut, TBYTE, 6,frow, felem, writesize, pixid,&status);
	fits_write_col(fptrOut, TBYTE, 7,frow, felem, writesize, detx,&status);
    fits_write_col(fptrOut, TBYTE, 8,frow, felem, writesize, dety,&status);
	fits_write_col(fptrOut, TUSHORT, 9, frow, felem, writesize, veto,&status);   
    fits_write_col(fptrOut, TBYTE, 10,frow, felem, writesize,alpha,&status);
	fits_write_col(fptrOut, TINT, 11,frow, felem, writesize, pi, &status);
	fits_write_col(fptrOut, TFLOAT, 12,frow, felem, writesize,energy, &status);

	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}



void writeBadpixExtension(char *outputfile)
{
	fitsfile *fptrOut;      
	int status, ii, jj,hdunum=2,i=0,j=0,xstart,xend,ystart,yend;
	long  fpixel, nelements,nrows;
	short *badpix;
	unsigned char *pix_flag,*detx,*dety,**pix_flag_2d;

	badpix=(short*)malloc((128*128)*sizeof(short));
	pix_flag_2d=(unsigned char**)malloc(sizeof(unsigned char*)*64);
	for(i=0;i<64;i++) 
	{
		pix_flag_2d[i]=(unsigned char*)malloc(sizeof(unsigned char)*64);
	}


	int bitpix   =  SHORT_IMG; 
	long naxis    =   2;                           
	long naxes[2] = { 128, 128 };

	status = 0;         /* initialize status before calling fitsio routines */

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
		printerror( status );

	for(ii=0;ii<128;ii++)	
		for(jj=0;jj<128;jj++)
			badpix[ii*128+jj]=0;
	for(i=0;i<4;i++)
	{   
	 for (ii = 0 ; ii <64 ; ii++)
	{  
		 for (jj =0 ; jj <64; jj++)
		 {
		    pix_flag_2d[ii][jj] = 0;
		 }
	} 
	if ( fits_movabs_hdu(fptrOut, i+hdunum, NULL, &status) ) 
		printerror( status );
	fits_get_num_rows(fptrOut, &nrows, &status);
	pix_flag  = (unsigned char*)malloc(nrows * sizeof(unsigned char));
	detx  = (unsigned char*)malloc(nrows * sizeof(unsigned char));
	dety  = (unsigned char*)malloc(nrows * sizeof(unsigned char));

	fits_read_col(fptrOut, TBYTE, 3, 1, 1, nrows, NULL, detx,NULL, &status);   
	fits_read_col(fptrOut, TBYTE, 4, 1, 1,nrows, NULL, dety,NULL, &status);
	fits_read_col(fptrOut, TBYTE, 5, 1, 1,nrows, NULL, pix_flag,NULL, &status);

	if(i==0)
	{
		xstart=0;	xend=63;
		ystart=64;	yend=127;
	}
	if(i==1)
	{
		xstart=64;	xend=127;
		ystart=64;	yend=127;
	}
	if(i==2)
	{
		xstart=64;	xend=127;
		ystart=0;	yend=63;
	}
	if(i==3)
	{
		xstart=0;	xend=63;
		ystart=0;	yend=63;
	}

	for (j = 0; j < nrows; j++)
	{
		pix_flag_2d[detx[j]][dety[j]]=pix_flag[j];
	}
	//printf("xstart %d xend %d ystart %d yend %d\n",xstart,xend,ystart,yend);
	int k=0,l=0;
	for(ii = ystart,k=0 ; ii <=yend ; ii++,k++)
	{  
		 for (jj =xstart,l=0 ; jj <=xend; jj++,l++)
		 {
		    badpix[ii*128+jj] = pix_flag_2d[l][k];
		 }
	}
	}

	if ( fits_create_img(fptrOut,  bitpix, naxis, naxes, &status) )
	 printerror( status );  
	status=0;
	fpixel = 1;                              
	nelements = naxes[0] * naxes[1];
	fits_update_key(fptrOut, TSTRING, "EXTNAME","BADPIX", NULL, &status);

	if ( fits_write_img(fptrOut, TSHORT, fpixel, nelements, &badpix[0], &status) )
	printerror( status );
	
	free(pix_flag);free(detx);free(dety);          
	free(badpix);
	for(i=0;i<64;i++) 
	{
	free(pix_flag_2d[i]);
	}
	free(pix_flag_2d);

	if ( fits_close_file(fptrOut, &status) )
	printerror( status ); 
	
	

	return;
}

void writeBTI(unsigned char *detid,unsigned char *pixid,unsigned char *detx,unsigned char *dety,double *tstart,double *tstop,int writesize,char *outputfile,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem,nrows;
	status=0;
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );
	fits_get_num_rows(fptrOut, &nrows, &status);

	frow      = nrows+1;//firstindex;
	felem     = 1;

    	fits_write_col(fptrOut, TBYTE, 1, frow, felem, writesize,detid,&status);
    	fits_write_col(fptrOut, TBYTE, 2,frow, felem, writesize, pixid,&status);
	fits_write_col(fptrOut, TBYTE, 3,frow, felem, writesize, detx,&status);
    	fits_write_col(fptrOut, TBYTE, 4,frow, felem, writesize, dety,&status);
    	fits_write_col(fptrOut, TDOUBLE, 5,frow, felem, writesize,tstart,&status);
	fits_write_col(fptrOut, TDOUBLE, 6,frow, felem, writesize,tstop,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

/*
void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i;
	double tstart,tstop,tstartf,tstopf,start[4],stop[4],largest,smallest,exposure;
	int tstarti,tstopi;
	status=0;

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         printerror( status );

	
    	for(i=0;i<=3;i++)
    	{       
		fits_movabs_hdu(fptrOut, i+2, NULL, &status);
		fits_read_key(fptrOut,TDOUBLE,"TSTART",&tstart,NULL, &status);
		fits_read_key(fptrOut,TDOUBLE,"TSTOP",&tstop,NULL, &status);
		
		start[i]=tstart;
		stop[i]=tstop;
	}

	smallest = start[0];
	for (i = 1; i < 4; i++)
		if (smallest > start[i])
			smallest = start[i];

	largest = stop[0];
	for (i = 1; i < 4; i++)
		if (largest < stop[i])
			largest = stop[i];

	tstarti=(int)smallest;
	tstopi=(int)largest;

	tstartf=smallest-tstarti;
	tstopf=largest-tstopi;
	exposure=largest-smallest;


	if ( fits_movabs_hdu(fptrOut, 1, NULL, &status) ) 
	      	printerror( status );

	fits_update_key(fptrOut, TDOUBLE,"TSTART",&smallest,"Start time of observation",&status);
	fits_update_key(fptrOut, TDOUBLE,"TSTOP",&largest,"Stop time of observation",&status);
	fits_update_key(fptrOut, TINT,"TSTARTI", &tstarti,"Start time of observation Integer part", &status);
	fits_update_key(fptrOut, TINT,"TSTOPI", &tstopi,"Stop time of observation Integer part", &status);
	fits_update_key(fptrOut, TDOUBLE,"TSTARTF", &tstartf,"Start time of observation Fractional part", &status);
	fits_update_key(fptrOut, TDOUBLE,"TSTOPF", &tstopf,"Stop time of observation Fractional part", &status);
	fits_update_key(fptrOut, TDOUBLE,"EXPOSURE", &exposure,"Exposure time", &status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	return;
}
 */

void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i,tstarti,tstopi;
	double tstart,tstop,start[4],stop[4],largest,smallest,texposure, exposure[4],tstartf,tstopf;
	status=0;
	double smallest_exposure;
	char *comment; 

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         printerror( status );
    	for(i=0;i<4;i++)
    	{       
		fits_movabs_hdu(fptrOut, i+2, NULL, &status);
		fits_read_key(fptrOut,TDOUBLE,"TSTART",&tstart,NULL, &status);
		fits_read_key(fptrOut,TDOUBLE,"TSTOP",&tstop,NULL, &status);
		fits_read_key(fptrOut,TDOUBLE,"EXPOSURE",&texposure,NULL, &status);
		
		start[i]=tstart;
		stop[i]=tstop;
		exposure[i] = texposure;
	}
	smallest_exposure = exposure[0];
	for (i = 1; i < 4; i++)
		if (smallest_exposure > exposure[i])
			smallest_exposure = exposure[i];


	smallest = start[0];
	for (i = 1; i < 4; i++)
		if (smallest > start[i])
			smallest = start[i];

	largest = stop[0];
	for (i = 1; i < 4; i++)
		if (largest < stop[i])
			largest = stop[i];
			

	if ( fits_movabs_hdu(fptrOut, 1, NULL, &status) ) 
	      	printerror( status );
	tstarti=(int)smallest;
	tstopi=(int)largest;

	tstartf=smallest-tstarti;
	tstopf=largest-tstopi;
	//exposure=largest-smallest;
	printf("Time : %f\t%f\t%d\t%d\t%10f\t%10f\n",smallest,largest,tstarti,tstopi,tstartf,tstopf);
	
	fits_update_key(fptrOut, TDOUBLE,"TSTART",&smallest,"Start time of observation",&status);
	fits_update_key(fptrOut, TDOUBLE,"TSTOP",&largest,"Stop time of observation",&status);
	fits_update_key(fptrOut, TINT,"TSTARTI", &tstarti,"Start time of observation Integer part", &status);
	fits_update_key(fptrOut, TINT,"TSTOPI", &tstopi,"Stop time of observation Integer part", &status);
	fits_update_key(fptrOut, TDOUBLE,"TSTARTF", &tstartf,"Start time of observation Fractional part", &status);
	fits_update_key(fptrOut, TDOUBLE,"TSTOPF", &tstopf,"Stop time of observation Fractional part", &status);
	fits_update_key(fptrOut, TDOUBLE,"EXPOSURE", &smallest_exposure,"Exposure time", &status);
	
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	return;
}


void modifyExposure(char *outputfile, double **pix_exposure,int qid)
{
	fitsfile *fptrOut; 
	int status, hdutype;
	long frow, felem;
	status=0;
	//double tstartf,tstopf;//,exposure=0.0;
	//int tstarti, tstopi;
	frow      = 1;
	felem     = 1;
	//tstarti = (int)tstart;
	//tstopi = (int)tstop;
	//tstartf = tstart-(double)tstarti;
	//tstopf = tstop-(double)tstopi;
	//printf("%d\t%d\n",tstarti,tstopi);
	

	status=0;
	
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	        
	
	//~ if ( fits_movabs_hdu(fptrOut, 13, &hdutype, &status) ) 
	      	//~ printerror( status );
	      	
	if ( fits_movnam_hdu(fptrOut, BINARY_TBL, "EXPOSURE", 0, &status) ) 
	      	printerror( status );
	      	
	      	
	      	
		   	
	fits_write_col(fptrOut, TDOUBLE, qid+1, frow, felem,4096,pix_exposure[qid],&status);
	
	
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	
	return;
		
}


void modifycommonGTI(char *outputfile)
{
	fitsfile *fptrOut;
	double **gtitstart,**gtitstop;
	int anynull;
	long gtinrowsf=0,gtinrows_temp;
	long gtinrows[4];
	double tstart,tstop,tstart_temp,tstop_temp;
	int i,j,k;
	int status, hdutype,intnull;
	long frow, felem,nrows;
	double doublenull;
	unsigned char bytenull;
	status=0;
	frow      = 1;
	felem     = 1;
	doublenull=0.0;
	anynull=0;
	
	gtitstart =	(double**)malloc(sizeof(double*)*4);
	gtitstop  =	(double**)malloc(sizeof(double*)*4);	
	
	if (fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
			printerror( status );
	
	int initial_flag=0;
	
	for(i=0;i<4;i++)
    	{
		fits_movabs_hdu(fptrOut,i+10, &hdutype, &status);
		fits_read_key(fptrOut,TDOUBLE,"TSTART",&tstart_temp,NULL, &status);
		fits_read_key(fptrOut,TDOUBLE,"TSTOP",&tstop_temp,NULL, &status);
		
		if(initial_flag < 1)
		{
			tstart = tstop_temp;
			tstop  = tstop_temp;
		}
		
		if(tstart_temp < tstart) {tstart=tstart_temp;}
		
		if(tstop_temp > tstop) {tstop=tstop_temp;}


		fits_get_num_rows(fptrOut, &gtinrows[i] , &status);
		gtitstart[i]= (double*)malloc(sizeof(double)*gtinrows[i]);
		gtitstop[i]= (double*)malloc(sizeof(double)*gtinrows[i]);
		fits_read_col(fptrOut, TDOUBLE, 1, frow, felem, gtinrows[i], &doublenull, gtitstart[i],&anynull, &status);
		fits_read_col(fptrOut, TDOUBLE, 2, frow, felem, gtinrows[i], &doublenull, gtitstop[i], &anynull, &status);
		gtinrowsf= gtinrowsf+gtinrows[i];
	}
	
	double *gtitstart_f,*gtitstop_f,*gtitstart_temp,*gtitstop_temp;
	gtitstart_f = (double*)malloc(sizeof(double)*gtinrowsf);
	gtitstop_f  = (double*)malloc(sizeof(double)*gtinrowsf);
	gtitstart_temp = (double*)malloc(sizeof(double)*gtinrowsf);
	gtitstop_temp  = (double*)malloc(sizeof(double)*gtinrowsf);
	
	gtinrows_temp = gtinrows[0];
	for(i=0;i<gtinrows[0];i++)
	{
		if( gtitstart[0][i] >= gtitstop[0][i])
		{
			printf("GTI ERROR PRESENT in Quadrant 0\n");
		}
		
		
		gtitstart_temp[i] = gtitstart[0][i];
		gtitstop_temp[i]  = gtitstop[0][i];
	}
	long gti_counter=0;
	
	for(i=1;i<4;i++)
    	{
		gti_counter=0;
		for(j=0;j<gtinrows_temp;j++)
		{
			for(k=0;k<gtinrows[i];k++)
			{
				if( gtitstart[i][k] >= gtitstop[i][k])
				{
					printf("GTI ERROR PRESENT in Quadrant %d\n",i);
				}
				
				if( (gtitstart[i][k]>=gtitstart_temp[j] &&  gtitstart[i][k]<gtitstop_temp[j])  &&  (gtitstop[i][k]>gtitstart_temp[j] &&  gtitstop[i][k]<=gtitstop_temp[j]) )
				{
					gtitstart_f[gti_counter] = gtitstart[i][k];
					gtitstop_f[gti_counter]  = gtitstop[i][k];
					gti_counter++;
				}
				else if( !(gtitstart[i][k]>=gtitstart_temp[j] &&  gtitstart[i][k]<gtitstop_temp[j])  &&  (gtitstop[i][k]>gtitstart_temp[j] &&  gtitstop[i][k]<=gtitstop_temp[j]) )
				{
					gtitstart_f[gti_counter] = gtitstart_temp[j];
					gtitstop_f[gti_counter]  = gtitstop[i][k];
					gti_counter++;				
				}
				else if( (gtitstart[i][k]>=gtitstart_temp[j] &&  gtitstart[i][k]<gtitstop_temp[j])  &&  !(gtitstop[i][k]>gtitstart_temp[j] &&  gtitstop[i][k]<=gtitstop_temp[j]) )
				{
					gtitstart_f[gti_counter] = gtitstart[i][k];
					gtitstop_f[gti_counter]  = gtitstop_temp[j];
					gti_counter++;
				}
				else if(gtitstart[i][k]<=gtitstart_temp[j] &&   gtitstop[i][k]>gtitstop_temp[j]  )
				{
					gtitstart_f[gti_counter] = gtitstart_temp[j];
					gtitstop_f[gti_counter]  = gtitstop_temp[j];
					gti_counter++;					
				}
				
			}
		}
		
		for(j=0;j<gti_counter;j++)
		{
			gtitstart_temp[j] = gtitstart_f[j];
			gtitstop_temp[j]  = gtitstop_f[j];
		}
		gtinrows_temp=gti_counter;
		
		
	}
	
	double exposure=0.0;
	for(j=0;j<gti_counter;j++)
	{
		exposure = exposure + (gtitstop_f[j]-gtitstart_f[j]);
		//printf("%lf\t%lf\n",gtitstart_f[j],gtitstop_f[j]);
    
	}
	
	printf("EXPOSURE = %lf\n",exposure);
	
	double tstartf,tstopf;
	int tstarti, tstopi;
	frow      = 1;
	felem     = 1;
	tstarti = (int)tstart;
	tstopi = (int)tstop;
	tstartf = tstart-(double)tstarti;
	tstopf = tstop-(double)tstopi;
	
	if ( fits_movabs_hdu(fptrOut,9, &hdutype, &status) ) 
	     	printerror( status );
	     	
	
	fits_write_key(fptrOut, TDOUBLE,"TSTART", &tstart,"Start time of observation", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP", &tstop,"Stop time of observation", &status);
	fits_write_key(fptrOut, TINT,"TSTARTI", &tstarti,"Start time of observation Integer part", &status);
	fits_write_key(fptrOut, TINT,"TSTOPI", &tstopi,"Stop time of observation Integer part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTARTF", &tstartf,"Start time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOPF", &tstopf,"Stop time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"EXPOSURE", &exposure,"Exposure time", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTART",&tstart,"Start time of observation",&status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP",&tstop,"Stop time of observation",&status);
	
	//printf("True\n");
	fits_write_col(fptrOut, TDOUBLE, 1, frow, felem,gti_counter,gtitstart_f,&status);
    	fits_write_col(fptrOut, TDOUBLE, 2,frow, felem,gti_counter,gtitstop_f,&status);
	
	free(gtitstart_f);free(gtitstop_f);free(gtitstart_temp);free(gtitstop_temp);
	
	for(i=0;i<4;i++) {free(gtitstart[i]);free(gtitstop[i]);}
	free(gtitstart);free(gtitstop);

	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
	
	

	/*
	for(i=0;i<gti_counter;i++)
	{
		for(j=0;j<ii;j++)
		{
			if(((GTITSTART_merge[j]>=gtitstart_sel[i])&&(GTITSTART_merge[j]<gtitstop_sel[i]))&&((GTITSTOP_merge[j]>=gtitstart_sel[i])&&(GTITSTOP_merge[j]<gtitstop_sel[i])))
			{
				
				gtitstart_f[gtinrows_f] =  GTITSTART_merge[j];
				gtitstop_f[gtinrows_f] =  GTITSTOP_merge[j];
				gtinrows_f++;
			}
			
			else if((!(GTITSTART_merge[j]>=gtitstart_sel[i])&&(GTITSTART_merge[j]<gtitstop_sel[i]))&&((GTITSTOP_merge[j]>=gtitstart_sel[i])&&(GTITSTOP_merge[j]<gtitstop_sel[i])))
			{
			
				gtitstart_f[gtinrows_f] =  gtitstart_sel[i];
				gtitstop_f[gtinrows_f] =  GTITSTOP_merge[j];
				gtinrows_f++;
			}
			
			else if(((GTITSTART_merge[j]>=gtitstart_sel[i])&&(GTITSTART_merge[j]<gtitstop_sel[i]))&& !((GTITSTOP_merge[j]>=gtitstart_sel[i])&&(GTITSTOP_merge[j]<gtitstop_sel[i])))
			{
			
				gtitstart_f[gtinrows_f] =  GTITSTART_merge[j];
				gtitstop_f[gtinrows_f] =  gtitstop_sel[i];
				gtinrows_f++;
			}
			else if((GTITSTART_merge[j]<=gtitstart_sel[i])&& (GTITSTOP_merge[j]>=gtitstop_sel[i]))
			{
			
				gtitstart_f[gtinrows_f] =  gtitstart_sel[i];
				gtitstop_f[gtinrows_f] = gtitstop_sel[i];
				gtinrows_f++;
			}
			
			
		}
		
	}*/

	
}
	
	



void printerror( int status)
{
	if (status)
	{
		fits_report_error(stderr, status); 	
		exit( status );    	
    }
  	return;
}

int read_input_parameters(int r, char *infile,char *inbadpixfile, char *thresholdfile,char *outfile,char *outbadpixfile,char *outbtifile,int *clobber,int *history)
{

	int status=0,hdutype=0;
	FILE *fp;
	fitsfile *fptr;
	if(r<0)
	{
		printf("Error(%s:%d) : Error while loading par file\n",__FILE__,__LINE__);
		exit(0);
	}
	
	r=PILGetFname("infile",infile);
	if (fits_open_file(&fptr, infile, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,infile);
		printerror( status );
	}
	fits_close_file(fptr,&status);

	r=PILGetFname("inbadpixfile",inbadpixfile);
	if (fits_open_file(&fptr, inbadpixfile, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,inbadpixfile);
		printerror( status );
	}
	fits_close_file(fptr,&status);

	r=PILGetFname("thresholdfile",thresholdfile);
	if((fp = fopen(thresholdfile, "r")) == NULL) {
		 printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,thresholdfile);
	}
	fclose(fp);
	r=PILGetFname("outfile",outfile);
	r=PILGetFname("outbadpixfile",outbadpixfile);
	r=PILGetFname("outbtifile",outbtifile);

	//---------------------------------------------------------Overwrite files and remember history------------------------------------------------------------
	r=PILGetBool("clobber", clobber);
	r=PILGetBool("history", history);
	//-------------------------------------------------------------End of input paramters----------------------------------------------------------------------	
	PILClose(r);	//Closing PIL file
	
	return(status);
}
int display_input_parameters(char *infile, char *inbadpixfile,char *thresholdfile,char *outfile,char *outbadpixfile,char *outbtifile,int clobber,int history)
{
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf("                                                    CZTFLICKPIXCLEAN PARAMETERS \n");
    	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Input file              		: %s\n",infile);
	printf(" Input badpix file              	: %s\n",inbadpixfile);
    	printf(" Threshold file	              		: %s\n",thresholdfile);
	printf(" Output badpix file 			: %s\n",outbadpixfile);
	printf(" Output BTI file             		: %s\n",outbtifile);
	printf(" Output file 				: %s\n",outfile);
	printf(" Clobber				: %d\n",clobber);
    	printf(" History				: %d\n",history);
    	printf("----------------------------------------------------------------------------------------------------------------------------\n\n");
    	return 0;
}

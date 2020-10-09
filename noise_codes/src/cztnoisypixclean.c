/*
cztnoisypixclean.c
Ajay Vibhute,Mayuri Shinde,Ajay Ratheesh.

code to filter the noisy pixels and generate noise cleaned event file. 

inputs:- 1:Event file name 2.CALDB badpix file 3.Bunch file name 4.Threshold file 
output:- 1.Noise cleaned level2 event file 2.level2 livetime file 3.bad pixel file
* 
* Updated by Ajay Ratheesh on 20th Dec 2017
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "pil.h"
#define BUFFSIZE 1024
int dimen=64;
long size=1;
void processNoiseReduction();
void printerror( int status);
void createEventFile(char *outputfile,char *infile);
void createBadpixFile(char *outputfile);
void createLivetimeFile(char *outputfile);
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char *detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,long bufsize,int hdunum, double exposure);
void writeBadpix(unsigned char *detid,unsigned char *pixid,unsigned char *pixx,unsigned char *pixy,unsigned char *pix_flag,char *outputfile,int hdunum);
void writeLivetime(double *time,float *livetime,int livetimesize,char *outputfile,int hdunum);
void livetimeGeneration(double *evttime,double*,long evtnrows,double *bunchtime,long bunchnrows,char *outputfile,int hdunum);
void modifyEventHeaderParams(char *outputfile);
void writegtiextension(double *gtistart, double *gtistop, double gtinrows,char *outputfile, int hdunum, double tstart, double tstop, double exposure);
void modifyExposure(char *outputfile, double **pix_exposure,int qid);
float * calculateSD(double **dph,float *mean,int pix_flag_2d[dimen][dimen]);
float * calculateMean(double **dph,int pix_flag_2d[64][64] );
int read_input_parameters(int r, char *infile, char *caldb_badpix_file, char *thresholdfile,char *outbadpixfile,char *outfile,int *clobber,int *history);
int display_input_parameters(char *infile, char *caldb_badpix_file, char *thresholdfile,char *outbadpixfile,char *outfile,int clobber,int history);

char *CZTNOISECLEAN,infile[BUFFSIZE],parfilename[BUFFSIZE], outfile[BUFFSIZE],outbadpixfile[BUFFSIZE],caldb_badpix_file[BUFFSIZE],thresholdfile[BUFFSIZE];
int clobber=1, history=1;

int main(int argc, char *argv[])
{

	/*if(argc!=4)
	{
		printf("Enter all command line arguments\n1:Event file name\n2.CALDB badpix file\n3.Threshold file \n");
		exit(-1);
	}
	
	if(argc==4)
    	{
		processNoiseReduction(argv[1],argv[2],argv[3]);
	}*/

	CZTNOISECLEAN = getenv ("CZTNOISECLEAN");
	if(CZTNOISECLEAN==NULL)
	{
		printf("CZTNOISECLEAN Variable is not set\n");
		exit(0);
	}
	strcpy(parfilename,CZTNOISECLEAN);
	strcat(parfilename,"/paramfiles/cztnoisypixclean.par");
	PILSetModuleName(parfilename);
	//Read all inputs through Parameter Interface Library
	int r=PILInit(argc,argv);
	//Read all input files & all the selection parameters
	
	read_input_parameters(r, infile, caldb_badpix_file, thresholdfile,outbadpixfile,outfile,&clobber,&history);
	display_input_parameters(infile, caldb_badpix_file, thresholdfile,outbadpixfile,outfile,  clobber,  history);

	processNoiseReduction();

	//return 0;
}

void processNoiseReduction()
{
	fitsfile *caldbfptr,*evtfptr,*bunchfptr;  
	FILE *logfile;
    long ii;     
	int flag,status, hdunum, hdutype,   anynull,i,intnull,badpixdetidcolnum,badpixpixidcolnum,pixflagcolnum;
	long frow, felem, nelem,nrows,evtnrows,bunchnrows, longnull;
	double *evttime,*cztseccnt,*finalevttime,*finalcztseccnt,doublenull,*bunchtime;
	int j;
    	unsigned short *cztntick,*veto,*finalcztntick,*finalveto,*pha,*finalpha;
	unsigned char *detid,*pixid,*detx,*dety,*alpha,*finaldetid,*finalpixid,*finaldetx,*finaldety,*finalalpha;
	float *energy,floatnull,*finalenergy;
    	int *pi,*finalpi;
	unsigned char *caldb_detid,*caldb_pixid,*pix_flag,*final_pix_flag,bytenull,*caldb_detx,*caldb_dety;
    	int pix_flag_2d[64][64],ignore_count[2];
	//int energycolnum,picolnum,timecolnum,cztseccntcolnum,cztntickcolnum,phacolnum,detidcolnum,pixidcolnum,detxcolnum,detycolnum,vetocolnum,alphacolnum;
	int quadstart=0,quadend=4;
	//float **dph;
	double **dph;
	//double dph[64][64];
	float *sd,*mean;
	char outlivetimefile[1000],outlogfile[1000];
	int DEAD_PIX_THRESHOLD,DEAD,NOISY;
	float THRESHOLD;

	FILE *fp;
	
	//create event and badpix file
	//char *file = strtok(infile, ".");

	//sprintf(outevtfile,"%s_nc.evt",file);

	//sprintf(outbadpixfile,"%s_badpix_nc.fits",file);
	//sprintf(outlivetimefile,"%s_nc_livetime.fits",file);
	
	if(clobber==1)
	{
		remove(outfile);
		remove(outbadpixfile);
		//remove(outlivetimefile);

	}

    	createEventFile(outfile,infile);
	createBadpixFile(outbadpixfile);
	//createLivetimeFile(outlivetimefile);

	sprintf(outlogfile, "%s_nc.log",infile);
	logfile=fopen(outlogfile,"a"); 


	caldb_detid  = (unsigned char*)malloc(size * sizeof(unsigned char));
	caldb_pixid  = (unsigned char*)malloc(size * sizeof(unsigned char));
	pix_flag  = (unsigned char*)malloc(size * sizeof(unsigned char));
	caldb_detx  = (unsigned char*)malloc(size * sizeof(unsigned char));
	caldb_dety  = (unsigned char*)malloc(size * sizeof(unsigned char));
	final_pix_flag  = (unsigned char*)malloc(size * sizeof(unsigned char));

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
	
	bunchtime  = (double*)malloc(size * sizeof(double));
	
	status = 0;
	hdunum = 2;


	fprintf(logfile,"Inputs for Noisy Pixel Clean \n1.Event file : %s\n2.Caldb file : %s\n3.Threshold file : %s\n\n",infile,caldb_badpix_file,thresholdfile);

	if((fp = fopen(thresholdfile, "r")) == NULL) {
		 printf("Failed to open Noise Reduction Threshold file\n");
	}
	
	const size_t label_size = 300;
	char* label = malloc(label_size);
	
	while (fgets(label, label_size, fp) != NULL)  
	{
		char *value;
		strtok_r (label, " ", &value);
		if (strcmp(label,"DEAD_PIX_THRESHOLD") == 0)
			DEAD_PIX_THRESHOLD=atoi(value);
		else if (strcmp(label,"DEAD") == 0)
			DEAD=atoi(value);
		else if (strcmp(label,"NOISY") == 0)
			NOISY=atoi(value);
		else if (strcmp(label,"THRESHOLD") == 0)
			THRESHOLD=atof(value);	
	}
	free(label); 
	if ( fits_open_file(&caldbfptr, caldb_badpix_file, READONLY, &status) )
	{
	    printf("Error in opening a file : %d",status);
	    printerror( status );
	}


	status=0;
	if ( fits_open_file(&evtfptr, infile, READONLY, &status) )
	{
	    	printf("Error in opening a file : %d",status);
	   	printerror( status );
	} 


	/*
	status=0;
	if ( fits_open_file(&bunchfptr, bunchfile, READONLY, &status) )
	{
	    printf("Error in opening a file : %d",status);
	    printerror( status );
	}*/
	frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;
	floatnull=0.0;
	
	double **pixel_exposure;
	pixel_exposure=(double**)malloc(sizeof(double*)*4);
	for(i=0;i<4;i++) pixel_exposure[i]=(double*)malloc(sizeof(double)*4096);
	for(i=0;i<4;i++){for(j=0;j<4096;j++) pixel_exposure[i][j]=0.0;}

	for(i=quadstart; i<quadend; i++)
	{

		fprintf(logfile,"************************************Quad %d is processing**********************************\n",i);
		status=0;
		fits_movabs_hdu(caldbfptr, i+2, NULL, &status);
	    fits_get_num_rows(caldbfptr, &nrows, &status);
		
		caldb_detid  = (unsigned char*)realloc(caldb_detid,nrows * sizeof(unsigned char));
		caldb_pixid  = (unsigned char*)realloc(caldb_pixid,nrows * sizeof(unsigned char));
		pix_flag  = (unsigned char*)realloc(pix_flag,nrows * sizeof(unsigned char));
		caldb_detx  = (unsigned char*)realloc(caldb_detx,nrows * sizeof(unsigned char));
		caldb_dety  = (unsigned char*)realloc(caldb_dety,nrows * sizeof(unsigned char));
		final_pix_flag  = (unsigned char*)realloc(final_pix_flag,nrows * sizeof(unsigned char));
	


		if( fits_get_colnum(caldbfptr, CASEINSEN, "DETID",&badpixdetidcolnum, &status) )
		{
			printf("Error(%s:%d): Unable to get column no for DETID in file %s\nExiting...",__FILE__,__LINE__,caldb_badpix_file);
			exit(-1);
		}
		if(fits_get_colnum(caldbfptr, CASEINSEN, "PIXID", &badpixpixidcolnum, &status))
		{
	
			printf("Error(%s:%d): Unable to get column no for PIXID\nExiting...",__FILE__,__LINE__);
			exit(-1);
		}
		if(fits_get_colnum(caldbfptr, CASEINSEN, "PIX_FLAG",&pixflagcolnum, &status))
		{
	
			printf("Error(%s:%d): Unable to get column no for PIX_FLAG\nExiting...",__FILE__,__LINE__);
			exit(-1);
		}


/*

	  	fits_read_col(caldbfptr, TBYTE, badpixdetidcolnum, frow, felem, nrows, &bytenull, caldb_detid,&anynull, &status);    
		fits_read_col(caldbfptr, TBYTE, badpixpixidcolnum, frow, felem, nrows, &bytenull, caldb_pixid,&anynull, &status);
		fits_read_col(caldbfptr, TBYTE, pixflagcolnum, frow, felem, nrows, &bytenull, pix_flag,&anynull, &status); 

*/

	  	fits_read_col(caldbfptr, TBYTE, 1, frow, felem, nrows, &bytenull, caldb_detid,&anynull, &status);    
		fits_read_col(caldbfptr, TBYTE, 2, frow, felem, nrows, &bytenull, caldb_pixid,&anynull, &status);
		fits_read_col(caldbfptr, TBYTE, 3, frow, felem, nrows, &bytenull, pix_flag,&anynull, &status); 

/*
		fits_get_colnum(evtfptr, CASEINSEN, "Time", &timecolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "CZTSECCNT", &cztseccntcolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "CZTNTICK",&cztntickcolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "PHA",&phacolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "DetID", &detidcolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "pixID",&pixidcolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "DETX",&detxcolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "DETY", &detycolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "veto",&vetocolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "alpha",&alphacolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "pi",&picolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "energy",&energycolnum, &status);

*/
		int gcount=0,bcount=0;
		//detid,pixid to detx,dety conversion
		for (ii = 0; ii < nrows; ii++)
		{
			
			caldb_detx[ii]=((caldb_detid[ii]%4)*16)+(caldb_pixid[ii]%16);
			caldb_dety[ii]=((caldb_detid[ii]/4)*16)+(caldb_pixid[ii]/16);
			
			if(i==0 || i==3)
				caldb_dety[ii]=63-caldb_dety[ii];
			else
				caldb_detx[ii]=63-caldb_detx[ii];
			pix_flag_2d[caldb_detx[ii]][caldb_dety[ii]]=pix_flag[ii];
			if(pix_flag_2d[caldb_detx[ii]][caldb_dety[ii]]==0)
					gcount++;
			if(pix_flag_2d[caldb_detx[ii]][caldb_dety[ii]]==1)
					bcount++;
		}
		printf("Good Pix : %d\tSpectroscopically Bad Pix : %d\n",gcount,bcount);
		
		fits_movabs_hdu(evtfptr, i+2, NULL, &status);
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

		//printf("%ld\n",evtnrows);
		

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
		
		
		//fits_movabs_hdu(bunchfptr, i+2, NULL, &status);
		//fits_get_num_rows(bunchfptr, &bunchnrows, &status);
		
		//bunchtime  = (double*)realloc(bunchtime,bunchnrows * sizeof(double));
		//fits_read_col(bunchfptr, TDOUBLE, 1, frow, felem, bunchnrows, &doublenull, bunchtime,&anynull, &status);

		

		dph=(double**)calloc(64,sizeof(double*));//for valgrind it converted to calloc from malloc

		for(ii=0;ii<64;ii++)
		{
			dph[ii]=(double*)calloc(64,sizeof(double));
		}
		//dph = malloc(dimen * sizeof(float *));
		//for(ii = 0; ii < dimen; ii++)
		//{
		//	dph[ii] = malloc(dimen * sizeof(int));
		//}
						
		

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
		fits_read_col(evtfptr, TINT, 11, frow, felem, evtnrows, &floatnull, pi,&anynull, &status);
		fits_read_col(evtfptr, TFLOAT, 12, frow, felem, evtnrows, &floatnull, energy,&anynull, &status);
	    
	   
		//livetimeGeneration(evttime,cztseccnt,evtnrows,bunchtime,bunchnrows,outlivetimefile,i+2);
		//printf("RUN TILL HERE\n");   
		
		
		
		
		
		//generate dph i.e.count events in each pixel
		for (ii = 0; ii < evtnrows; ii++)
		{
			dph[detx[ii]][dety[ii]]++;		
		}
		
		//find dead pixels
		fprintf(logfile,"****************************************Dead pixels****************************************\nDETX\tDETY\n");
		for(ii=0;ii<dimen;ii++)
		{
			for(j=0;j<dimen;j++)
			{
				
				//printf("%lf\n",dph[ii][j]);
				
				
				if(dph[ii][j]<DEAD_PIX_THRESHOLD)
				{
					pix_flag_2d[ii][j]=DEAD;				
					fprintf(logfile,"%d\t%d\n",ii,j);
				}			
			}
		}
		//return ;
		//printf("%ld\n",nrows);
		//correct for area
		
		
		
		float width=0.0,height=0.0,area=0.0;
		for (ii = 0; ii < nrows; ii++)
		{
			if((caldb_detx[ii]+1)%16==0 || (caldb_detx[ii])%16==0 || (caldb_dety[ii]+1)%16==0 || (caldb_dety[ii])%16==0)
			{
				height=2.28;
				width=2.28;
			}
			else
			{
				height=2.46;
				width=2.46;
			}
			
			area=width*height;
			//printf("%u\t%u AREA === %f DPH ==== %f Division product=====%f\n",caldb_detx[ii],caldb_dety[ii],area,dph[caldb_detx[ii]][caldb_dety[ii]],dph[caldb_detx[ii]][caldb_dety[ii]]/area);
			dph[caldb_detx[ii]][caldb_dety[ii]] = dph[caldb_detx[ii]][caldb_dety[ii]]/area;
			//printf("%f\n",dph[caldb_detx[ii]][caldb_dety[ii]]);
			
		}
		
		
		
		
		
		
		//for(ii=0;ii<64;ii++)
		//{
		//	for(j=0;j<64;j++)
		//	{
		//		printf("%f\n",dph[ii][j]);
		//	}
		//}
		mean=calculateMean(dph,pix_flag_2d);
        sd = calculateSD(dph,mean,pix_flag_2d);	
        
        //printf("%f\t%f\t%f\t%f\t%f\n",mean[0],sd[0],mean[1],sd[1]);
        
		fprintf(logfile,"****************************************NOISY pixels****************************************\nDETX\tDETY\n");
		//divide pixels in two categories good and banana and find noisy pixels
		
			
		while(1)
		{
			ignore_count[0]=0;
			ignore_count[1]=0;
			for(ii=0; ii<dimen; ++ii)
	    		{
				for(j=0; j<dimen; ++j)
	    			{
						
					
					if(pix_flag_2d[ii][j]==0)
					{
						if(dph[ii][j]>mean[0]+THRESHOLD*sd[0] || dph[ii][j]<mean[0]-THRESHOLD*sd[0])
						{
							dph[ii][j]=-1;
							pix_flag_2d[ii][j]=NOISY;
							fprintf(logfile,"%d\t%d\n",ii,j);
							ignore_count[0]++;
						}
					}	
					if(pix_flag_2d[ii][j]==1)
					{
						if(dph[ii][j]>mean[1]+THRESHOLD*sd[1] || dph[ii][j]<mean[1]-THRESHOLD*sd[1])
						{
							dph[ii][j]=-1;
							pix_flag_2d[ii][j]=NOISY;
							fprintf(logfile,"%d\t%d\n",ii,j);
							ignore_count[1]++;
						}
					}
				}
			}
			if(ignore_count[0]==0 && ignore_count[1]==0)
				break;
			
			//recompute mean and sd of good and banana pixels
			mean=calculateMean(dph,pix_flag_2d);
                	sd = calculateSD(dph,mean,pix_flag_2d);			
		}
		
        
        int k;
        long l=0;


	
		//write noise cleaned event file
		for (ii = 0; ii < evtnrows; ii++)
		{

			if(pix_flag_2d[detx[ii]][dety[ii]]!=NOISY && pix_flag_2d[detx[ii]][dety[ii]]!=DEAD)
			{
				finalevttime[l]=evttime[ii];
				finalcztseccnt[l]  = cztseccnt[ii];
				finalpha[l]  = pha[ii];
				finalcztntick[l]  = cztntick[ii];
				finalveto[l]  = veto[ii];
				finaldetid[l]  = detid[ii];
				finalpixid[l]  = pixid[ii];
				finaldetx[l]  = detx[ii];
				finaldety[l]  = dety[ii];
				finalalpha[l] = alpha[ii];
				finalenergy[l] = energy[ii];
				finalpi[l] = pi[ii];
				l++;
			}
		
		}
		
		int jj,detx_t,dety_t;
		
		
		
		for(k=0;k<16;k++)
		{
			for(jj=0;jj<256;jj++)
			{
				detx_t = ((k%4)*16)+(jj%16);
				dety_t = ((k/4)*16)+(jj/16);
				
				if(i ==0 ||i ==3)
					dety_t=63-dety_t;
				else
					detx_t=63-detx_t;
				
				//pix_exposure_det_pix[i][j] = pix_exposure[detx_t][dety_t];
				
				if(pix_flag_2d[detx_t][dety_t]<2)
				{
					pixel_exposure[i][k*256+jj] = 1.0;
					
				}
				
				//printf("%lf\t%lf\n",pix_exposure_det_pix[i][j],pixel_exposure[qid][i*256+j]);
				
			}
			
		}

		int agcount=0,abcount=0;

		//write badpix file		
		for (ii = 0; ii < nrows; ii++)
		{
			final_pix_flag[ii]=pix_flag_2d[caldb_detx[ii]][caldb_dety[ii]];
			if(pix_flag_2d[caldb_detx[ii]][caldb_dety[ii]]==0)
					agcount++;
			if(pix_flag_2d[caldb_detx[ii]][caldb_dety[ii]]==1)
					abcount++;
		}
		fprintf(logfile,"\nFrom CALDB...\n");
		fprintf(logfile,"Good Pix : %d\tSpectroscopically Bad Pix : %d\n",gcount,bcount);
		fprintf(logfile,"After processing...\n");
		fprintf(logfile,"Good Pix : %d\tSpectroscopically Bad Pix : %d\n",agcount,abcount);
		int gdiff=gcount-agcount;
		int bdiff=bcount-abcount;
		printf("good_diff = %d  Banana diff = %d \n",gdiff,bdiff);
		fprintf(logfile,"Difference is...\n");
		fprintf(logfile,"Good Pix : %d\tSpectroscopically Bad Pix : %d\n\n",gdiff,bdiff);
		writeBadpix(caldb_detid,caldb_pixid,caldb_detx,caldb_dety,final_pix_flag,outbadpixfile,i+2);
		
		
			
		double *gtitstart,*gtitstop,*new_gtitstart,*new_gtitstop;
		double tstart,tstop;
		long gtinrows,ngtinrows=0;
	
		fits_read_key(evtfptr,TDOUBLE,"TSTART",&tstart,NULL, &status);
		fits_read_key(evtfptr,TDOUBLE,"TSTOP",&tstop,NULL, &status);
		
		//printf("STATUS %d\n",status);
		
		
		fits_movabs_hdu(evtfptr, i+10, &hdutype, &status);
		fits_get_num_rows(evtfptr, &gtinrows, &status);
		
		gtitstart=(double*)malloc(sizeof(double)*gtinrows);
		gtitstop= (double*)malloc(sizeof(double)*gtinrows);
		
		new_gtitstart	= (double*)malloc(sizeof(double)*gtinrows);
		new_gtitstop	= (double*)malloc(sizeof(double)*gtinrows);
		
		
		
		frow      = 1;
		felem     = 1;
		doublenull = 0.;
		intnull = 0;
		bytenull = 0;
		floatnull=0.0;
	
		fits_read_col(evtfptr, TDOUBLE, 1, frow, felem, gtinrows, &doublenull, gtitstart,&anynull, &status);
		fits_read_col(evtfptr, TDOUBLE, 2, frow, felem, gtinrows, &doublenull, gtitstop, &anynull, &status); 
		
		//long gti_counter=0;
		int j;
		for(j=0;j<gtinrows;j++)
		{
			
			if(((gtitstart[j]>=tstart)&&(gtitstart[j]<tstop))&&((gtitstop[j]>tstart)&&(gtitstop[j]<=tstop)))
			{
			
				new_gtitstart[ngtinrows]= gtitstart[j];
				new_gtitstop[ngtinrows]= gtitstop[j];
				ngtinrows++;
			}
			else if(((gtitstart[j]>=tstart)&&(gtitstart[j]<tstop)) && !((gtitstop[j]>tstart)&&(gtitstop[j]<=tstop)))
			{
			
				new_gtitstart[ngtinrows]= gtitstart[j];
				new_gtitstop[ngtinrows]= tstop;
				ngtinrows++;
			}
			else if(!((gtitstart[j]>=tstart)&&(gtitstart[j]<tstop)) && ((gtitstop[j]>tstart)&&(gtitstop[j]<=tstop)))
			{
				
				new_gtitstart[ngtinrows]= tstart;
				new_gtitstop[ngtinrows]= gtitstop[j];
				ngtinrows++;
			}
			else if(((gtitstart[j]<=tstart) && (gtitstop[j]>=tstop)))
			{
				new_gtitstart[ngtinrows]= tstart;
				new_gtitstop[ngtinrows]= tstop;
				ngtinrows++;
				
			}
		
		}
		
		double exposure=0.0;
			
		for(j=0;j<ngtinrows;j++)
		{
			if( new_gtitstart[j] >= new_gtitstop[j])
			{
				printf("ERROR-------GTI  Tstart >= Tstop---------\n");
			}
			//printf("%lf\t%lf\n",new_gtitstart[j],new_gtitstop[j]);
			exposure = exposure+(new_gtitstop[j]-new_gtitstart[j]);
		}
			
			
		//printf("EXPOSURE ---%lf\n",exposure);
			
		modifyExposure(outfile,pixel_exposure,i);	
		writegtiextension(new_gtitstart, new_gtitstop,ngtinrows,outfile,i+10,tstart,tstop,exposure);
		writeEvent(finalevttime,finalcztseccnt,finalcztntick,finalpha,finaldetid,finalpixid,finaldetx,finaldety,finalveto,finalalpha,finalpi,finalenergy,outfile,l,i+2,exposure);
	
		free(gtitstart);free(gtitstop);free(new_gtitstart);free(new_gtitstop);
		for(ii=0;ii<64;ii++){free(dph[ii]);} free(dph);

	}
	modifyEventHeaderParams(outfile);
	//modifyEventHeaderParams(outlivetimefile);
	//free memory
	free(evttime);free(cztseccnt);free(finalevttime);free(finalcztseccnt);free(bunchtime);
	free(cztntick);free(veto);free(finalcztntick);free(finalveto);free(pha);free(finalpha);
	free(detid);free(pixid);free(detx);free(dety);free(alpha);free(finaldetid);free(finalpixid);
	free(finaldetx);free(finaldety);free(finalalpha);
	free(energy);free(finalenergy);
	free(pi);free(finalpi);free(caldb_detid);free(caldb_pixid);
	free(pix_flag);free(final_pix_flag);free(caldb_detx);free(caldb_dety);
	for(i=0;i<4;i++) free(pixel_exposure[i]);
	free(pixel_exposure);
		free(sd);free(mean);
	
	if ( fits_close_file(caldbfptr, &status) )       
	        printerror( status );
	if ( fits_close_file(evtfptr, &status) )       
	        printerror( status );
	//if ( fits_close_file(bunchfptr, &status) )       
	//        printerror( status );

	fclose(fp);fclose(logfile);

	printf("NOISE REDUCTION COMPLETED SUCCESSFULLY.\n");

	return;
}


void livetimeGeneration(double *evttime,double *cztseccnt,long evtnrows,double *bunchtime,long bunchnrows,char *outputfile,int hdunum)
{
	double *time;
        float *livetime;
	long bintime,k=0,temp2=0;
	int j,temp1=0;
	int *counter_evt,*counter_bunch,bunchcount;
	int livetimesize;
	double lastevt,lastbunch,ssmlc=0.0;
	long lasttime=0,tmod=0,lasttmod=0;
	livetimesize=(int)(evttime[evtnrows-1])-(int)(evttime[0])+1;//+1;
	time=(double*)malloc(sizeof(double)*livetimesize);
	livetime=(float*)malloc(sizeof(float)*livetimesize);
	double *bunchtime_cztsec;
	bunchtime_cztsec = (double*)malloc(sizeof(double)*bunchnrows);
	double time_thresh;
	counter_evt=(int*)calloc(sizeof(int),livetimesize);
	counter_bunch=(int*)calloc(sizeof(int),livetimesize);
	double last_event_czt, last_bunch_czt, last_event_ut, last_bunch_ut, corr=0.0, last_evt;
	int temp3=0;
	float ssmlc1=0;
	int isfracremain=-1;
	double frac=0;
	long i=0;
	int bunch_counter=0,event_counter=0;
	lasttime=(long)cztseccnt[0];
	time_thresh = evttime[0] - cztseccnt[0];
	for(bintime=(long)evttime[0];bintime<(long)evttime[evtnrows-1]+1;bintime++)//read tstart and tend of that quadrant instead off (int)evttime[0]
	{
		
		lastbunch=bintime;
		lastevt=bintime;
		bunchcount=0;
		ssmlc=0;

		int mm=0;
		for(j=temp1;j<evtnrows;j++)
		{

			if(evttime[j]>(double)bintime && evttime[j]<(double)bintime+1)
			{
			
				
				tmod=(long)cztseccnt[j]%100;
				/*if(evttime[j]>221154093.0 && evttime[j]<221154094.5)
				{
					printf("%f\t%f\t%d\t%ld\t%ld\n",evttime[j],cztseccnt[j],tmod,(long)evttime[j],(long) (evttime[j]-0.4));
				}*/


				if(lasttmod!=tmod && tmod==0 &&j!=0 )
				{
					


					//if 100 ms drop comes on the second boundry
					if((long)evttime[j]==(long)(evttime[j]-0.3))
					{
						
						ssmlc=0.3;
						ssmlc1=0.0;

					}
					else
					{

						mm=1;
					//	mm=0;
						/*ssmlc=0.4-((evttime[j]+0.4)-(long)(evttime[j]+0.4));	
						ssmlc1=0.4-ssmlc;
						isfracremain=0;
						//printf("%f\t%f\t%f\t%ld\n",ssmlc,ssmlc1,evttime[j]+0.4, (long)(evttime[j]+0.4));
						//scanf("%d");
						break;*/
					}

				}
				else
				{
					/*if(isfracremain==0)
					{
						ssmlc=ssmlc1;
						isfracremain=-1;
						ssmlc1=0.0;
					}
					else	*/
					{
					//	ssmlc=0.0;
					//	ssmlc1=0.0;
					}
				}
				lasttmod=tmod;

				
	
				counter_evt[i]++;
				temp1=j;
				lastevt=evttime[j];
				time_thresh = evttime[j] - cztseccnt[j];
			}
			if(evttime[j]-(double)bintime >= 5.0){break;}
	
			lasttime=(long)cztseccnt[j];

		}
		for(k=temp2;k<bunchnrows;k++)
		{
			bunchtime_cztsec[k] = bunchtime[k]-time_thresh;
			if(bunchtime[k]>(double)bintime && bunchtime[k]<(double)bintime+1)
			{	
					
				counter_bunch[i]++;
				temp2=k;
				lastbunch=bunchtime[k];
			
			}
			if(bunchtime[k]-(double)bintime >= 10){break;}
		}
		time[i]=bintime+0.5;
		isfracremain--;
		if(mm==1)
		{
			livetime[i-1]-=0.3;
			livetime[i]=1.0-((float)counter_bunch[i]*20.0/1000000.0);
		}
		else
			livetime[i]=1.0-((float)counter_bunch[i]*20.0/1000000.0)-ssmlc;
		
			i++;
	}
	temp1=0;
	temp2=0;
	k=0;	

	for(bintime=(long)cztseccnt[0];bintime<(long)cztseccnt[evtnrows-1]+1;bintime++)
	{
		event_counter=0;bunch_counter=0;
		for(j=temp1;j<evtnrows;j++)
		{
			if(cztseccnt[j]>(double)bintime && cztseccnt[j]<(double)bintime+1.0)

			{
				event_counter++;
				last_event_ut = evttime[j];
				last_event_czt = cztseccnt[j];
				temp1=j;
			}
			if(cztseccnt[j]-(double)bintime >= 10.0){break;}
		}

		for(k=temp2;k<bunchnrows;k++)
		{
			if(bunchtime_cztsec[k]>(double)bintime && bunchtime_cztsec[k]<(double)bintime+1.0)
			{
				bunch_counter++;
				last_bunch_ut = bunchtime[k];
				last_bunch_czt = bunchtime_cztsec[k];
				temp2=k;
			}
			if(bunchtime_cztsec[k]-(double)bintime >= 10.0){break;}

		}
		corr =0.0;
	//	printf("no of evt = %d\n",bunch_counter+event_counter);
		if(((double)bintime+1.0-last_event_czt)>0.05 || ((double)bintime+1.0-last_bunch_czt)>0.05)
		{
			
			if(last_bunch_czt<last_event_czt) {corr = (double)bintime+1-last_event_czt;last_evt=last_event_czt;}
			else{corr = (double)bintime+1.0-last_bunch_czt;last_evt=last_bunch_czt;}
			//printf("no of events = %d %lf\n",bunch_counter+event_counter, corr);
			for(j=temp3; j < i ; j++)
			{

				//printf("%lf %lf\n",time[j],last_evt);
				if((time[j]+0.5-time_thresh) >= bintime && (time[j]+0.5-time_thresh) <= bintime+1.0)
				{	
					if((time[j]+0.5-time_thresh-last_evt)>0)
					{
						
						livetime[j]-=(time[j]+0.5-time_thresh-last_evt);
						livetime[j+1]-=bintime+1.0-(time[j]+0.5-time_thresh);
						//printf("@@@@@@%lf %lf %lf\n",time[j]+0.5-time_thresh-last_evt,livetime[j],livetime[j+1]);
						if(livetime[j]<0.1){livetime[j]=1.0;}
						if(livetime[j+1]<0.1){livetime[j+1]=1.0;}

					}
					else
					{
						
						livetime[j+1]-=corr;
						//printf("$$$$$$$%lf %lf %lf\n",time[j]+0.5-time_thresh-last_evt,livetime[j],livetime[j+1]);
						if(livetime[j+1]<0.1){livetime[j+1]=1.0;}
					}
					temp3 = j;
				}
			}

		}
		/*if( (bunch_counter+event_counter) >= 3072)
		{

			printf("no of events = %d\n",bunch_counter+event_counter );
			if(last_bunch_czt<last_event_czt) {corr = bintime+1-last_event_czt;last_evt=last_event_czt;}
			else{corr = bintime+1.0-last_bunch_czt;last_evt=last_bunch_czt;}
			for(j=temp3; j < i ; j++)
			{
				if((time[j]+0.5-time_thresh) >= bintime && (time[j]+0.5-time_thresh) <= bintime+1.0)
				{	
					if((time[j]+0.5-time_thresh-last_evt)>0)
					{
						livetime[j]-=(time[j]+0.5-time_thresh-last_evt);
						livetime[j+1]-=bintime+1.0-(time[j]+0.5);
						if(livetime[j]<0.05){livetime[j]=1.0;}
						if(livetime[j+1]<0.05){livetime[j+1]=1.0;}

					}
					else
					{
						livetime[j+1]-=corr;
						if(livetime[j+1]<0.05){livetime[j+1]=1.0;}
					}
					temp3 = j;
				}
			}
		}*/
	}
	
	writeLivetime(time,livetime,i,outputfile,hdunum);
	free(time);
        free(livetime);
	free(counter_evt);free(counter_bunch);
	free(bunchtime_cztsec);

}

//function to calculate mean
float * calculateMean(double **dph,int pix_flag_2d[64][64] )
{
	double goodsum = 0.0,bananasum=0.0;
	double g_count=0.0,b_count=0.0;
    	int i,j;
	float *mean;
	mean = calloc(2,sizeof(float));
	for(i=0; i<dimen; i++)
	{
		for(j=0; j<dimen; j++)
	    	{
			//good pixel
			//printf("%f\n",dph[i][j]);				
			if( pix_flag_2d[i][j]==0 && dph[i][j]>0)
			{
				goodsum += dph[i][j];
				g_count++;
				//printf("%f\t%f\t%f\n",g_count,goodsum,dph[i][j]);
			}
			//banana pixel
			if( pix_flag_2d[i][j]==1 && dph[i][j]>0)
			{
				bananasum += dph[i][j];
				b_count++;
				//printf("%f\t%f\t%f\n",b_count,bananasum,dph[i][j]);
			}
	    	}
	}
	//printf("%f\t%f\t%lf\t%lf\n",g_count,b_count,goodsum,bananasum);
	
	mean[0] = (float)(goodsum/g_count);
	mean[1] = (float)(bananasum/b_count);
	//printf("Mean : %f\t%f\n",mean[0],mean[1]);
	return mean;
}

//function to calculate standard deviation
float * calculateSD(double **dph,float *mean,int pix_flag_2d[dimen][dimen])
{
	double g_count=0.0,b_count=0.0;
    	int i,j;
	double *variance;
	float *sd;
	
		
	variance = calloc(2,sizeof(double));
	sd = calloc(2,sizeof(float));
	for(i=0; i<dimen; ++i)
	{
		
		for(j=0; j<dimen; ++j)
		{
			//good pixel
			if( pix_flag_2d[i][j]==0 && dph[i][j]>0)
			{
				variance[0] += pow(dph[i][j] - (double)mean[0],2);
				g_count++;
			}
			//banana pixel
			if( pix_flag_2d[i][j]==1 && dph[i][j]>0)
			{
				variance[1] += pow(dph[i][j] - (double)mean[1],2);
				b_count++;
			}
	    	}
	}
	sd[0] = (float)sqrt(variance[0]/g_count);
	sd[1] = (float)sqrt(variance[1]/b_count); 
	//printf("SD : %f\t%f\n",sd[0],sd[1]);
	free(variance);
	return sd;
}

//creating event file 


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

    	for(i=0;i<=3;i++)
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

void createLivetimeFile(char *outputfile)
{
	fitsfile *fptrOut;      
	int status, hdutype;
	long frow, felem;
      	int tfields,i;
	char extname[20];
	
	tfields=2;
	char *ttype[] = {"TIME","FRACEXP"};
	char *tform[] = { "D","E"};
	char *tunit[] = {"",""};
	        
	int hdunum=2;
	
	status=0;
	
        if (fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

		//if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         //printerror( status );

    	for(i=0;i<=3;i++)
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

void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,long writesize,int hdunum,double exposure)
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
        fits_write_col(fptrOut, TUSHORT, 3,frow, felem, writesize, cztntick,&status);  //numerical data overflow exception 
        fits_write_col(fptrOut, TUSHORT, 4,frow, felem, writesize, pha,&status); //numerical data overflow exception 
	fits_write_col(fptrOut, TBYTE, 5, frow, felem, writesize,detid,&status);
        fits_write_col(fptrOut, TBYTE, 6,frow, felem, writesize, pixid,&status);
	fits_write_col(fptrOut, TBYTE, 7,frow, felem, writesize, detx,&status);
        fits_write_col(fptrOut, TBYTE, 8,frow, felem, writesize, dety,&status);
	fits_write_col(fptrOut, TUSHORT, 9, frow, felem, writesize, veto,&status); //numerical data overflow exception  
        fits_write_col(fptrOut, TBYTE, 10,frow, felem, writesize,alpha,&status);
	fits_write_col(fptrOut, TINT, 11,frow, felem, writesize, pi, &status);
	fits_write_col(fptrOut, TFLOAT, 12,frow, felem, writesize,energy, &status);	
		
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

void writeBadpix(unsigned char *detid,unsigned char *pixid,unsigned char *pixx,unsigned char *pixy,unsigned char *pix_flag,char *outputfile,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	status=0;
	int writesize=4096;
	
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

void writeLivetime(double *time,float *livetime,int livetimesize,char *outputfile,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	double tstart,tstop;
	status=0;
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	frow      = 1;
	felem     = 1;

	tstart=time[0];
	tstop=time[livetimesize-1];

	fits_write_key(fptrOut, TDOUBLE,"TSTART",&tstart,"Start time of observation",&status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP",&tstop,"Stop time of observation",&status);

        fits_write_col(fptrOut, TDOUBLE, 1, frow, felem,livetimesize ,time,&status);
        fits_write_col(fptrOut, TFLOAT, 2,frow, felem, livetimesize, livetime,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

/*
void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i,tstarti,tstopi;
	double tstart=0,tstop=0,start[4],stop[4],largest,smallest,exposure,tstartf,tstopf;
	status=0;
	double t=0,t1=0;
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         printerror( status );

    	for(i=0;i<=3;i++)
    	{       
		status=0;
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

	if ( fits_movabs_hdu(fptrOut, 1, NULL, &status) ) 
	      	printerror( status );
	tstarti=(int)smallest;
	tstopi=(int)largest;

	tstartf=smallest-tstarti;
	tstopf=largest-tstopi;
	exposure=largest-smallest;
	//printf("Time : %f\t%f\t%d\t%d\t%10f\t%10f\n",smallest,largest,tstarti,tstopi,tstartf,tstopf);
	
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
} */

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


void writegtiextension(double *gtistart, double *gtistop, double gtinrows,char *outputfile, int hdunum, double tstart, double tstop, double exposure)
{
	fitsfile *fptrOut; 
	int status, hdutype;
	long frow, felem;
	status=0;
	double tstartf,tstopf;//,exposure=0.0;
	int tstarti, tstopi;
	frow      = 1;
	felem     = 1;
	tstarti = (int)tstart;
	tstopi = (int)tstop;
	tstartf = tstart-(double)tstarti;
	tstopf = tstop-(double)tstopi;
	//printf("%d\t%d\n",tstarti,tstopi);
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );
	 

	status=0;
	//printf("Final Exposure = %lf\n",exposure);
	//printf("%lf\t%lf\n",exposure,tstop-tstart);
	fits_write_key(fptrOut, TDOUBLE,"TSTART", &tstart,"Start time of observation", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP", &tstop,"Stop time of observation", &status);
	fits_write_key(fptrOut, TINT,"TSTARTI", &tstarti,"Start time of observation Integer part", &status);
	fits_write_key(fptrOut, TINT,"TSTOPI", &tstopi,"Stop time of observation Integer part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTARTF", &tstartf,"Start time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOPF", &tstopf,"Stop time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"EXPOSURE", &exposure,"Exposure time", &status);
	
	fits_write_key(fptrOut, TDOUBLE,"TSTART",&tstart,"Start time of observation",&status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP",&tstop,"Stop time of observation",&status);
	
	
	//printf("FINISHED WRITING HEADER");
	    fits_write_col(fptrOut, TDOUBLE, 1, frow, felem,gtinrows,gtistart,&status);
	    fits_write_col(fptrOut, TDOUBLE, 2,frow, felem,gtinrows,gtistop,&status);
	
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
	
}


void createEventFile(char *outputfile,char *eventfile)
{
	
	fitsfile *fptrOut,*fptrevt;      
	int status, hdutype,anynull,tfields=12,i,hdunum=2;
	long frow, felem;
	int mjdrefi=55197,mjdreff=0,equinox=2000;
	float ra_pnt,dec_pnt;
	double timedel,telapse;
	char object[20],obs_id[20],obs_mode[20],date_obs[20],time_obs[20],date_end[20],time_end[20],date[20],creator[20],filename[70],checksum[20],datasum[20],chksumcomm[50],datasumcomm[50];
      
	char extname[20]; 
	          
	
	char *ttype[] = { "TIME", "CZTSECCNT","CZTNTICK","PHA","DetID","pixID","DETX","DETY","veto","alpha","PI","ENERGY"};
	char *tform[] = { "D","D","U","U","B","B","B","B","U","B","U","E"};
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
	
    	for(i=0;i<=3;i++)
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
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "GTI", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );
	/*
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
	*/
	
	   	
	        	
	
	
	char *ttype_gti[] = { "START", "STOP"};
	char *tform_gti[] = { "1D","1D"};
	char *tunit_gti[] = {"sec","sec"};
	
	
	
	int tfields_gti=2;
	for(i=0;i<=3;i++)
	{
		
		status=0;
		sprintf(extname, "Q%d_GTI",i);
		if(fits_movabs_hdu(fptrevt,i+8+hdunum, 0, &status)) 
			printerror( status );
		
		
		
		
		//printf("%d %d \n",i,hdunum);
		//if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        	//printerror( status );
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields_gti, ttype_gti, tform_gti,tunit_gti, extname, &status) )
			printerror( status );
			
		if ( fits_movabs_hdu(fptrOut, i+10, &hdutype, &status) ) 
	         	printerror( status );
		//if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			//printerror( status );

	}
	
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
		
	if ( fits_movabs_hdu(fptrOut, 14, &hdutype, &status) ) 
	         	printerror( status );
	

	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	if ( fits_close_file(fptrevt, &status) )       
	        printerror( status );
	return;
}










/*
void createEventFile(char *outputfile,char *eventfile)
{
	
	fitsfile *fptrOut,*fptrevt;      
	int status, hdutype,anynull,tfields=12,i,hdunum=2;
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

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );

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
	
    	for(i=0;i<=3;i++)
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
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "GTI", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );
	
	

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

		 
		 
	char *ttype_gti[] = { "START", "STOP"};
	char *tform_gti[] = { "1D","1D"};
	char *tunit_gti[] = {"sec","sec"};
	
	int tfields_gti=2;
	for(i=0;i<=3;i++)
	{
		
		status=0;
		sprintf(extname, "Q%d_GTI",i);
		//if(fits_movabs_hdu(fptrevt,i+8+hdunum, 0, &status)) 
			//printerror( status );
		
		//if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        	//printerror( status );
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields_gti, ttype_gti, tform_gti,tunit_gti, extname, &status) )
			printerror( status );
		//if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			//printerror( status );

	}

}
*/

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
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, 14, &hdutype, &status) ) 
	      	printerror( status );
	      	
	fits_write_col(fptrOut, TDOUBLE, qid+1, frow, felem,4096,pix_exposure[qid],&status);
	
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
		
}

int read_input_parameters(int r, char *infile, char *caldb_badpix_file, char *thresholdfile,char *outbadpixfile,char *outfile,int *clobber,int *history)
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
	
	r=PILGetFname("caldb_badpix_file",caldb_badpix_file);
	if (fits_open_file(&fptr, caldb_badpix_file, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,caldb_badpix_file);
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

	//---------------------------------------------------------Overwrite files and remember history------------------------------------------------------------
	r=PILGetBool("clobber", clobber);
	r=PILGetBool("history", history);
	//-------------------------------------------------------------End of input paramters----------------------------------------------------------------------	
	PILClose(r);	//Closing PIL file
	
	return(status);
}
int display_input_parameters(char *infile, char *caldb_badpix_file, char *thresholdfile,char *outbadpixfile,char *outfile,int clobber,int history)
{
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf("                                                    CZTNOISYPIXCLEAN PARAMETERS \n");
    	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Input file              		: %s\n",infile);
	printf(" Caldb badpix file             		: %s\n",caldb_badpix_file);
    	printf(" Threshold file	              		: %s\n",thresholdfile);
	printf(" Output file 				: %s\n",outfile);
	printf(" Output badpixel file 			: %s\n",outbadpixfile);
    	printf(" Clobber				: %d\n",clobber);
    	printf(" History				: %d\n",history);
    	printf("----------------------------------------------------------------------------------------------------------------------------\n\n");
    	return 0;
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







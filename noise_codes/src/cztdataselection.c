/*
cztdataselection.c
Ajay Vibhute,Mayuri Shinde.

code to generate saa removed event files

inputs:- 1.SAA Config file 2.Event file 3.MKF file
output:- SAA removed event file/files
* 
* 
* Adjusted by Ajay Ratheesh on 22/12/2017 to include GRBs near SAA.. The normal SAA cut doesnot work in this module
*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "pil.h"
#define BUFFSIZE 1024
void processDataSelection();
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outfile,long firstindex,long lastindex,long bufsize,int hdunum);
void printerror( int status);
void createEventFile(char *outfile,char *infile,int flag);
//int AVG_THRESHOLD=30;
void modifyEventHeaderParams(char *outfile);
double pixThresholdFix(int countrate);
int read_input_parameters(int r, char *infile, char *mkffile, char *thresholdfile,char *outfile,int *clobber,int *history);
int display_input_parameters(char *infile, char *mkffile, char *thresholdfile,char *outfile,int clobber,int history);

double SAA_TIME_THRESHOLD,SAA_TIME_THRESHOLD_2;
float SAA_ENERGY_THRESHOLD;
char *CZTNOISECLEAN,infile[BUFFSIZE],parfilename[BUFFSIZE], outfile[BUFFSIZE],mkffile[BUFFSIZE],thresholdfile[BUFFSIZE];
int clobber=1, history=1;

int main(int argc, char **argv)
{
	
	/*if(argc==4)
        {
		processDataSelection(argv[1],argv[2],argv[3]);
	}		
	if(argc!=4)
	{
		printf("Enter all command line arguments\n1:MKF file name\n2.Threshold file \n3.Event file name\n");
		exit(-1);
	}*/
	CZTNOISECLEAN = getenv ("CZTNOISECLEAN");
	if(CZTNOISECLEAN==NULL)
	{
		printf("CZTNOISECLEAN Variable is not set\n");
		exit(0);
	}
	strcpy(parfilename,CZTNOISECLEAN);
	strcat(parfilename,"/paramfiles/cztdataselection.par");
	PILSetModuleName(parfilename);
	//Read all inputs through Parameter Interface Library
	int r=PILInit(argc,argv);
	//Read all input files & all the selection parameters
	
	read_input_parameters(r, infile, mkffile, thresholdfile,outfile,&clobber,&history);
	display_input_parameters(infile, mkffile, thresholdfile,outfile,  clobber,  history);
	processDataSelection();
	
	return 0;
}

void processDataSelection()
{
	fitsfile *mkffptr,*evtfptr; 
	FILE *logfile,*fp;      
	int ii,flag,status, hdunum, hdutype,  nfound, anynull,i,intnull,ncols,cnummkftime,cnumlat,cnumlon,cnumevttime,**dph,size=1;
	long frow, felem, nelem,nrows,evtnrows, longnull;
	float floatnull, *lat,*lon;
	char strnull[10];
	double doublenull,*mkftime,*saastarttime,*saaendtime,exposure=0.0,tstart,tstop,time_diff;
	int SAALatRef,SAALongMin,SAALongMax,nsaa=0,n=-1,k=0,l=0;
	float SAAshift1,SAAshift2,SAASlopeLeft,SAASlopeRight;
	int bufsize,firstindex=0,lastindex=0,count;
	double *evttime,*cztseccnt;
        unsigned short *pha,*cztntick,*veto;
	int j,newindex;
	unsigned char *detid,*pixid,*detx,*dety,*alpha,bytenull;
	float *energy;
	int *pi,*evt_flag,*detpixid,saaflag[100];
	int energycolnum,picolnum,timecolnum,cztseccntcolnum,cztntickcolnum,phacolnum,detidcolnum,pixidcolnum,detxcolnum,detycolnum,vetocolnum,alphacolnum;
	int quadstart=0,quadend=4,total_quad=4;
	int isfilecreated=0,dimen=64;
	char *outtxtfile;
	double thresh_1s;
	float AVG_THRESHOLD;

	int *lcbin_1s,*temp_ind_array,*lcbin_avg_1s,tempindex=0;
	double bintime=0.0,*final_saastarttime,*final_saaendtime,**new_saastarttime,**new_saaendtime;
	int lc_index,boundry_change=0,end_trace=0;
	float avg;
	long *firstindex_array,*lastindex_array;
	
	saastarttime  = calloc(1000,sizeof(double));
	saaendtime  = calloc(1000,sizeof(double));
	outtxtfile=calloc(1000,sizeof(char));

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


  	status = 0;
	hdunum = 2;

	sprintf(outtxtfile, "%s_ds.log",outfile);
	//remove(outtxtfile);
	
	logfile=fopen(outtxtfile,"a"); 

	fprintf(logfile,"Inputs for dataselection\n1.MKF file : %s\n2.Threshold file : %s\n3.Event file : %s\n\n",mkffile,thresholdfile,infile);
	fprintf(logfile,"Important thresholds other than threshold file ....\n1.SAA_TIME_THRESHOLD = %f\n2.SAA_TIME_THRESHOLD_2 = %f\n3.SAA_ENERGY_THRESHOLD = %f\n\n",SAA_TIME_THRESHOLD,SAA_TIME_THRESHOLD_2,SAA_ENERGY_THRESHOLD);

	if((fp = fopen(thresholdfile, "r")) == NULL) {
		 printf("Failed to open SAA Threshold file\n");
	}
	const size_t label_size = 300;
	char* label = malloc(label_size);
	
	while (fgets(label, label_size, fp) != NULL)  
	{
		char *value;
		strtok_r (label, " ", &value);
		if (strcmp(label,"SAALatRef") == 0)
			SAALatRef=atoi(value);
		else if (strcmp(label,"SAALongMin") == 0)
			SAALongMin=atoi(value);
		else if (strcmp(label,"SAALongMax") == 0)
			SAALongMax=atoi(value);
		else if (strcmp(label,"SAASlopeLeft") == 0)
			SAASlopeLeft=atof(value);
		else if (strcmp(label,"SAASlopeRight") == 0)
			SAASlopeRight=atoi(value);
		else if (strcmp(label,"SAA_TIME_THRESHOLD") == 0)
			SAA_TIME_THRESHOLD=atof(value);
		else if (strcmp(label,"SAA_TIME_THRESHOLD_2") == 0)
			SAA_TIME_THRESHOLD_2=atof(value);
		else if (strcmp(label,"SAA_ENERGY_THRESHOLD") == 0)
			SAA_ENERGY_THRESHOLD=atof(value);
	}
	free(label);

	if ( fits_open_file(&mkffptr, mkffile, READONLY, &status) )
	{
	    printf("Error in opening a file : %d",status);
	    printerror( status );
	}
	if ( fits_movabs_hdu(mkffptr, 2, &hdutype, &status) )
	{
	     	printf("Error in moving : %d",status);
	     	printerror( status );	      
	}
	     	
    	fits_get_num_rows(mkffptr, &nrows, &status);
	lat  = malloc(nrows * sizeof(float));
	lon  = malloc(nrows * sizeof(float));
	mkftime  = malloc(nrows * sizeof(double));
	
	fits_get_colnum(mkffptr, CASEINSEN, "TIME",&cnummkftime, &status);
	fits_get_colnum(mkffptr, CASEINSEN, "EARTHLAT", &cnumlat, &status);
	fits_get_colnum(mkffptr, CASEINSEN, "EARTHLON",&cnumlon, &status);
	
	frow      = 1;
    	felem     = 1;
    	longnull  = 0;
    	floatnull = 0.;
	doublenull = 0.0;
	anynull=0;
  	fits_read_col(mkffptr, TDOUBLE, cnummkftime, frow, felem, nrows, &doublenull, mkftime, 
							   &anynull, &status);    
	fits_read_col(mkffptr, TFLOAT, cnumlat, frow, felem, nrows, &floatnull, lat, 
							   &anynull, &status); 
	fits_read_col(mkffptr, TFLOAT, cnumlon, frow, felem, nrows, &floatnull, lon,
							   &anynull, &status); 
	//find non SAA time durations
	for (ii = 0; ii < nrows; ii++){
	//calculate SAAshift1 and SAAshift2
		SAAshift1=SAASlopeLeft*(lat[ii]-SAALatRef);
		SAAshift2=SAASlopeRight*(lat[ii]-SAALatRef);
		//	printf("%d %f %f %f %f %d %d %f %f\n",ii,lat[ii],lon[ii],SAAshift1,SAAshift2,SAALongMin,SAALongMax,SAASlopeLeft,SAASlopeRight);
				
			if((lon[ii]>(SAALongMin+SAAshift1)) &&
				(lon[ii]<(SAALongMax+SAAshift2)) )
			{
				n=0;
				if(k==0)
				{
					saastarttime[l]=mkftime[ii];
					//printf("%f\n",saastarttime[l]);
				//	scanf("%d");
					l++;
				}
				k=1;
				
			}
			else
			{
				k=0;			
				if(n==0 && l>0)
				{
					 saaendtime[nsaa]=mkftime[ii-1];
					 nsaa++;
				}
				n=1;

			}	
	}

	for(ii=0;ii<nsaa;ii++)
		saaflag[ii]=0;


        if(nsaa==0)
	{
		printf("SAA is not present in %s file.\n",infile);
		sprintf(outfile, "%s_0.evt",outfile);
		remove(outfile);
	        createEventFile(outfile,infile,0);	
	}
	else
	{
	printf("SAA STARTTIME \t ENDTIME\n");
	fprintf(logfile,"SAA STARTTIME \t ENDTIME\n");
	for (ii = 0; ii <=nsaa && ii <l; ii++)
	{
		printf("%f\t%f\n",saastarttime[ii],saaendtime[ii]);
		fprintf(logfile,"%f\t%f\n\n",saastarttime[ii],saaendtime[ii]);
	}

	status=0;
	if ( fits_open_file(&evtfptr, infile, READONLY, &status) )
	{
	    	printf("Error in opening a file : %d",status);
	   	printerror( status );
	}   
	
	new_saastarttime = (double **)malloc(total_quad * sizeof(double *));
	new_saaendtime = (double **)malloc(total_quad * sizeof(double *));
        for (i=0; i<total_quad; i++)
	{
        	new_saastarttime[i] = (double *)malloc(nsaa * sizeof(double));
		new_saaendtime[i] = (double *)malloc(nsaa * sizeof(double));
	}
	final_saastarttime  = calloc(nsaa,sizeof(double));
	final_saaendtime  = calloc(nsaa,sizeof(double));
	for(i=quadstart; i<quadend; i++)
	{
	//	printf("QID %d\n",i);		
		fprintf(logfile,"Quad %d is processing...\n",i);
		fits_movabs_hdu(evtfptr, i+2, NULL, &status);
		fits_get_num_rows(evtfptr, &evtnrows, &status);

		evt_flag=(int*)calloc(evtnrows,sizeof(int));

		evttime  = (double*)realloc(evttime,evtnrows * sizeof(double));
		energy = (float*)realloc(energy,sizeof(float)*evtnrows);

		fits_get_colnum(evtfptr, CASEINSEN, "Time", &timecolnum, &status);
		fits_get_colnum(evtfptr, CASEINSEN, "energy",&energycolnum, &status);

		frow      = 1;
		felem     = 1;
		doublenull = 0.;
		intnull = 0;
		bytenull = 0;

		fits_read_col(evtfptr, TDOUBLE, timecolnum, frow, felem, evtnrows, &doublenull, evttime, &anynull, &status);
		fits_read_col(evtfptr, TFLOAT, energycolnum, frow, felem, evtnrows, &floatnull, energy, &anynull, &status);

		lcbin_1s = calloc((int)SAA_TIME_THRESHOLD,sizeof(int));
		temp_ind_array = calloc((int)SAA_TIME_THRESHOLD,sizeof(int));
		lcbin_avg_1s = calloc((int)SAA_TIME_THRESHOLD_2,sizeof(int));
		int nsaa_mem=0;
		nsaa_mem=nsaa+1;
		firstindex_array  = calloc(nsaa_mem,sizeof(long));
		lastindex_array  = calloc(nsaa_mem,sizeof(long));
		
		//Generating 1 sec lightcurve of all events energy>energy_threshold(saastart)
		for(ii=0;ii<nsaa;ii++)
		{

			for(j=0;j<(int)SAA_TIME_THRESHOLD;j++)
				lcbin_1s[j]=0;
		
		   	for(j=0;j<(int)SAA_TIME_THRESHOLD_2;j++)
				lcbin_avg_1s[j]=0;

			lc_index=0;
			avg=0.0;
			tempindex=0;
			boundry_change=0;
			AVG_THRESHOLD=0.0;
		
			//finding average
			for(bintime=(saastarttime[ii]-SAA_TIME_THRESHOLD_2);bintime< (saastarttime[ii]);bintime+=1.0)
			{		
				flag=0;
				
				for(j=tempindex;j<evtnrows;j++)
				{
					if(evttime[j]>=bintime && evttime[j]<(bintime+1.0))
					{					
						flag=1;
						if (energy[j]>SAA_ENERGY_THRESHOLD)
						{
							lcbin_avg_1s[lc_index]++;
						}
					}
					else if(flag==1)
					{
						break;
					}
				}
				tempindex=j;
				avg+=lcbin_avg_1s[lc_index];
				lc_index++;
			}
			
			avg=avg/SAA_TIME_THRESHOLD_2;
			//avg=60;
			//printf("start avg : %f\n",avg);

		//	AVG_THRESHOLD=avg/2.0;
			AVG_THRESHOLD=avg*100;
			lc_index=0;
			tempindex=0;

			//generating 1 sec light curve
			for(bintime=saastarttime[ii];bintime< (saastarttime[ii]+SAA_TIME_THRESHOLD);bintime+=1.0)
			{
				flag=0;
				for(j=tempindex;j<evtnrows;j++)
				{
					if(evttime[j]>=bintime && evttime[j]<(bintime+1.0))
					{
						tempindex=j;
						flag=1;
						if (energy[j]>SAA_ENERGY_THRESHOLD)
						{
							lcbin_1s[lc_index]++;
						}
					}
					else if(flag==1)
					{
						break;
					}
				}
				temp_ind_array[lc_index]=tempindex;
				
				//finding new saastart using average + avg threshold
				if(lcbin_1s[lc_index]>(avg+AVG_THRESHOLD))
				{
					new_saastarttime[i][ii]=evttime[temp_ind_array[lc_index]];//bintime+1;
					boundry_change=1;
					break;
				}
				lc_index++;
			}
			if(boundry_change==0)
			{
				new_saastarttime[i][ii]=saastarttime[ii];
			}
		}
		
		flag=0;
		tempindex=0;
		//Generating 1 sec lightcurve of all events energy>energy_threshold (saaend)
		for(ii=0;ii<nsaa;ii++)
		{
			lc_index=0;
			avg=0.0;
			tempindex=0;
			boundry_change=0;
			end_trace=0;
			for(j=0;j<(int)SAA_TIME_THRESHOLD;j++)
			{
				lcbin_1s[j]=0;
				temp_ind_array[j]=0;
			}
			for(j=0;j<(int)SAA_TIME_THRESHOLD_2;j++)
			{
				lcbin_avg_1s[j]=0;
			}
			
			//finding average
			for(bintime=saaendtime[ii];bintime< (saaendtime[ii]+SAA_TIME_THRESHOLD_2);bintime+=1.0)
			{		
				flag=0;
				for(j=tempindex;j<evtnrows;j++)
				{
					if(evttime[j]>=bintime && evttime[j]<(bintime+1.0))
					{	
						if(tempindex==0)
						{
							end_trace=j;
							tempindex=1;
						}				
						flag=1;
						if (energy[j]>SAA_ENERGY_THRESHOLD)
							lcbin_avg_1s[lc_index]++;
						
					}
					else if(flag==1)
					{
						break;
					}
				}
				tempindex=j;
				avg+=lcbin_avg_1s[lc_index];
				lc_index++;
			}
			
			avg=avg/SAA_TIME_THRESHOLD_2;
			//AVG_THRESHOLD=avg/2.0;
			AVG_THRESHOLD=avg*100;
			//avg=60;
			//printf("end avg : %f\n",avg);
			lc_index=0;
			tempindex=0;
			
			//generating 1 sec light curve
			for(bintime=saaendtime[ii];bintime>(saaendtime[ii]-SAA_TIME_THRESHOLD);bintime-=1.0)
			{
				flag=0;
				for(j=end_trace;j>0;j--)
				{
					if(evttime[j]<bintime && evttime[j]>=(bintime-1.0))
					{
						if(flag==0)
							end_trace=j;
						flag=1;
						if (energy[j]>SAA_ENERGY_THRESHOLD)
							lcbin_1s[lc_index]++;
					}
					else if(flag==1)
					{
						break;
					}
				}

				//finding new saaend using average+avg threshold
				temp_ind_array[lc_index]=end_trace;
				if(lcbin_1s[lc_index]>(avg+AVG_THRESHOLD))
				{
					new_saaendtime[i][ii]=evttime[temp_ind_array[lc_index]];
					boundry_change=1;
					break;
				}
				lc_index++;
			}
			if(boundry_change==0)
			{
				new_saaendtime[i][ii]=saaendtime[ii];
			}
		}
		free(evt_flag);
		free(lcbin_1s);
		free(temp_ind_array);
		free(lcbin_avg_1s);

	}
	for(i=0;i<4;i++)
	{
		printf("Quad %d\n",i);
		for(ii=0;ii<nsaa;ii++)
			printf("%f\t%f\n",new_saastarttime[i][ii],new_saaendtime[i][ii]);
	}
	printf("Final Start ...\nSAA Start\t SAA end\n");
	for (ii = 0; ii <nsaa; ii++)
	{
		//final_saastarttime[ii] = new_saastarttime[0][ii];
		//final_saaendtime[ii] = new_saaendtime[0][ii];
		
		final_saaendtime[ii] = saaendtime[ii]-500.0;
		final_saastarttime[ii] = saastarttime[ii]+500.0;
		
		
		for(i = 1; i<total_quad; i++)
		{
			//if (final_saaendtime[ii] > new_saaendtime[i][ii])
				//final_saaendtime[ii] = new_saaendtime[i][ii];
				final_saaendtime[ii] = saaendtime[ii]-500.0;
				
			//if (final_saastarttime[ii] < new_saastarttime[i][ii])
				//final_saastarttime[ii] = new_saastarttime[i][ii];
				final_saastarttime[ii] = saastarttime[ii]+500.0;
				
		}
		printf("%f\t%f\n",final_saastarttime[ii],final_saaendtime[ii]);

	}

	for(i=quadstart; i<quadend; i++)
	{
		//printf("QID %d\n",i);	
		printf("Quad %d is processing...\n",i);	
		fprintf(logfile,"Quad %d is processing...\n",i);
		fits_movabs_hdu(evtfptr, i+2, NULL, &status);
		fits_get_num_rows(evtfptr, &evtnrows, &status);
		
		evt_flag=(int*)calloc(evtnrows,sizeof(int));

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

		frow      = 1;
		felem     = 1;
		doublenull = 0.;
		intnull = 0;
		bytenull = 0;

		fits_read_col(evtfptr, TDOUBLE, timecolnum, frow, felem, evtnrows, &doublenull, evttime, &anynull, &status);
		fits_read_col(evtfptr, TDOUBLE, cztseccntcolnum, frow, felem, evtnrows, &doublenull, cztseccnt, &anynull, &status);         
		fits_read_col(evtfptr, TUSHORT, cztntickcolnum, frow, felem, evtnrows, &intnull, cztntick,&anynull, &status); 		     fits_read_col(evtfptr, TUSHORT, phacolnum, frow, felem, evtnrows, &intnull, pha, &anynull, &status);
		fits_read_col(evtfptr, TBYTE, detidcolnum, frow, felem,evtnrows, &bytenull, detid, &anynull, &status);  
		fits_read_col(evtfptr, TBYTE, pixidcolnum, frow, felem,evtnrows, &bytenull, pixid, &anynull, &status);        
		fits_read_col(evtfptr, TBYTE, detxcolnum, frow, felem, evtnrows, &bytenull, detx, &anynull, &status);   
		fits_read_col(evtfptr, TBYTE, detycolnum, frow, felem,evtnrows, &bytenull, dety, &anynull, &status);
		fits_read_col(evtfptr, TUSHORT, vetocolnum, frow, felem,evtnrows, &intnull, veto, &anynull, &status);
		fits_read_col(evtfptr, TBYTE, alphacolnum, frow, felem, evtnrows, &bytenull, alpha, &anynull, &status); 
		fits_read_col(evtfptr, TINT, picolnum, frow, felem, evtnrows, &floatnull, pi, &anynull, &status);
		fits_read_col(evtfptr, TFLOAT, energycolnum, frow, felem, evtnrows, &floatnull, energy, &anynull, &status);

                flag=-1;
		double tstart=evttime[0];
		for(j=0;j<nsaa+1;j++)
		{
				firstindex_array[j]=0;
				lastindex_array[j]=0;
		}
		ii=0;
		
		//creating event files of selected data(non saa region)		
		for(k=0;k<evtnrows;k++)
		{
			if(evttime[k]<=final_saastarttime[ii] && evttime[k]>=tstart)
			{		
				if(isfilecreated==0)
				{
					char tempout[BUFFSIZE];
					sprintf(tempout, "%s_%d.evt",outfile,ii);
					if(i==0)
					{
						printf("Creating event file %s\n",tempout);					
						remove(tempout);		
						saaflag[ii]=1;
						createEventFile(tempout,infile,1);													
					}						
					isfilecreated=1;	
					firstindex_array[ii]=k;
				}
				lastindex_array[ii]=k;
			}
			else if(evttime[k]>final_saaendtime[ii])
			{
				tstart=final_saaendtime[ii];
				if(ii==(nsaa-1))
				{
					char tempout[BUFFSIZE];
					sprintf(tempout, "%s_%d.evt",outfile,ii+1);
					if(i==0)
					{						
						remove(tempout);
						saaflag[ii+1]=1;
						createEventFile(tempout,infile,1);													
					}
					isfilecreated=1;
					firstindex_array[ii+1]=k-1;
					lastindex_array[ii+1]=evtnrows-1;
					break;
				}
				ii++;
				isfilecreated=0;
			}		
		}

		for(ii=0;ii<(nsaa+1);ii++)
		{
				//printf("Before writing %d...\n",i);
				char tempout[BUFFSIZE];	
				sprintf(tempout, "%s_%d.evt",outfile,ii);
				fprintf(logfile,"Writing in %s file.....\n",tempout);
				fprintf(logfile,"%d to %d events are selected from input evt file\n",firstindex_array[ii],lastindex_array[ii]);
				printf("Difference %d %d\n",lastindex_array[ii]-firstindex_array[ii],saaflag[ii]);
				if(lastindex_array[ii]-firstindex_array[ii]>0 && saaflag[ii])
					writeEvent(evttime,cztseccnt,cztntick,pha,detid,pixid,detx,dety,veto,alpha,pi,energy,tempout,firstindex_array[ii],lastindex_array[ii],lastindex_array[ii]-firstindex_array[ii],i+2);
				else
					printf("No events found for %d,file will not created\n",ii);
					
		}
		free(evt_flag);		
        }
	for(ii=0;ii<(nsaa+1);ii++)
	{
		if(saaflag[ii])
		{
			char tempout[BUFFSIZE];	
			sprintf(tempout, "%s_%d.evt",outfile,ii);
			modifyEventHeaderParams(tempout);
		}
	}
		//free(evt_flag);
		free(evttime);free(cztseccnt);
		free(pha);free(cztntick);free(veto);
		free(detid);free(pixid);free(detx);
		free(dety);free(alpha);free(energy);free(pi);

	}
	//free(lcbin_1s);free(temp_ind_array);free(lcbin_avg_1s);
	free(firstindex_array);free(lastindex_array);

	free(saastarttime);free(saaendtime);free(lat);free(lon);
        free(mkftime);
	for(i=0;i<4;i++)
	{
		free(new_saastarttime[i]);free(new_saaendtime[i]);	
	}

	free(new_saastarttime);free(new_saaendtime);
	free(final_saastarttime);free(final_saaendtime);

	printf("\nDATA SELECTION COMPLETED SUCCESSFULLY.\n");
	fclose(logfile);fclose(fp);free(outtxtfile);
	if ( fits_close_file(mkffptr, &status) )       
	        printerror( status );
	if ( fits_close_file(evtfptr, &status) )       
	        printerror( status );

	return;
}

//creating event file 
void createEventFile(char *outfile,char *eventfile,int flag)
{
	
	fitsfile *fptrOut,*fptrevt;      
	int status, hdutype,tfields=12,i,hdunum=2;
	long frow, felem;
	int mjdrefi=55197,mjdreff=0,equinox=2000;
	float ra_pnt,dec_pnt;
	double timedel,telapse,exposure[4096];
	char object[20],obs_id[20],obs_mode[20],date_obs[20],time_obs[20],date_end[20],time_end[20],date[20],creator[20],filename[70],checksum[20],datasum[20],chksumcomm[50],datasumcomm[50];
      
	char extname[20]; 
	          
	
	char *ttype[] = { "TIME", "CZTSECCNT","CZTNTICK","PHA","DetID","pixID","DETX","DETY","veto","alpha","PI","ENERGY"};
	char *tform[] = { "D","D","I","I","B","B","B","B","I","B","I","E"};
	char *tunit[] = {"s","s","micro-sec","counts","","","","","counts","counts","",""};

	char *expttype[] = { "EXPOSURE_Q0", "EXPOSURE_Q1","EXPOSURE_Q2","EXPOSURE_Q3"};
	char *exptform[] = { "D","D","D","D"};
	char *exptunit[] = {"","","",""};
       
	status=0;
        if (fits_create_file(&fptrOut, outfile, &status))
	       	 printerror( status );       

//	if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
//	        printerror( status );

	status=0;
	
	if ( fits_open_file(&fptrevt,eventfile, READONLY, &status) ) 
	         printerror( status );

	// Copy Primary
        if(fits_movabs_hdu(fptrevt, 1, 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	if(flag==0)
	{
		// Copy Q0
		if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q0", 0, &status)) 
			printerror( status );         

		if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			printerror( status );

		// Copy Q1
		if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q1", 0, &status)) 
			printerror( status );         

		if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			printerror( status );

		// Copy Q2
		if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q2", 0, &status)) 
			printerror( status );         

		if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			printerror( status );

		// Copy Q3
		if(fits_movnam_hdu(fptrevt, BINARY_TBL, "Q3", 0, &status)) 
			printerror( status );         

		if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			printerror( status );

	}
        if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	if(flag==1)
	{
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

			if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
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
	}
	for(i=0;i<4096;i++)
	{
		exposure[i]=1.0;
	}
	
	
	status=0;
	if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
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

	if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, 4, expttype, exptform,exptunit, "EXPOSURE", &status) )
			printerror( status );

	if(fits_movnam_hdu(fptrOut, BINARY_TBL, "EXPOSURE", 0, &status)) 
	         	printerror( status );

	frow      = 1;
	felem     = 1;

	fits_write_col(fptrOut, TDOUBLE, 1, frow, felem, 4096, exposure,&status);
	fits_write_col(fptrOut, TDOUBLE, 2, frow, felem, 4096, exposure,&status);
	fits_write_col(fptrOut, TDOUBLE, 3, frow, felem, 4096, exposure,&status);
	fits_write_col(fptrOut, TDOUBLE, 4, frow, felem, 4096, exposure,&status);

	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	if ( fits_close_file(fptrevt, &status) )       
	        printerror( status );
	return;
}

//write event file
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outfile,long firstindex,long lastindex,long writesize,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem,nrows;
	status=0;
	int tstarti,tstopi;
	double tstart,tstop,tstartf,tstopf,exposure;

	if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );
	fits_get_num_rows(fptrOut, &nrows, &status);

	tstart=evttime[firstindex];
	tstop=evttime[lastindex];

	tstarti=(int)tstart;
	tstopi=(int)tstop;

	tstartf=tstart-tstarti;
	tstopf=tstop-tstopi;
	exposure=tstop-tstart;
	//printf("Writing %d...\n",hdunum-2);	

	frow      = 1;
	felem     = 1;
	
	fits_write_key(fptrOut, TDOUBLE,"TSTART", &tstart,"Start time of observation", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOP", &tstop,"Stop time of observation", &status);
	fits_write_key(fptrOut, TINT,"TSTARTI", &tstarti,"Start time of observation Integer part", &status);
	fits_write_key(fptrOut, TINT,"TSTOPI", &tstopi,"Stop time of observation Integer part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTARTF", &tstartf,"Start time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"TSTOPF", &tstopf,"Stop time of observation Fractional part", &status);
	fits_write_key(fptrOut, TDOUBLE,"EXPOSURE", &exposure,"Exposure time", &status);

        fits_write_col(fptrOut, TDOUBLE, 1, frow, felem, writesize, &evttime[firstindex],&status);
	fits_write_col(fptrOut, TDOUBLE, 2, frow, felem, writesize, &cztseccnt[firstindex],&status);
        //fits_write_col(fptrOut, TUSHORT, 3,frow, felem, writesize, &cztntick[firstindex],&status);  //numerical data overflow exception 
        fits_write_col(fptrOut, TUSHORT, 4,frow, felem, writesize, &pha[firstindex],&status); 
	fits_write_col(fptrOut, TBYTE, 5, frow, felem, writesize,&detid[firstindex],&status);
        fits_write_col(fptrOut, TBYTE, 6,frow, felem, writesize, &pixid[firstindex],&status);
	fits_write_col(fptrOut, TBYTE, 7,frow, felem, writesize, &detx[firstindex],&status);
        fits_write_col(fptrOut, TBYTE, 8,frow, felem, writesize, &dety[firstindex],&status);
	fits_write_col(fptrOut, TUSHORT, 9, frow, felem, writesize, &veto[firstindex],&status);
        fits_write_col(fptrOut, TBYTE, 10,frow, felem, writesize,&alpha[firstindex],&status);
	fits_write_col(fptrOut, TINT, 11,frow, felem, writesize, &pi[firstindex], &status);
	fits_write_col(fptrOut, TFLOAT, 12,frow, felem, writesize,&energy[firstindex], &status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

void modifyEventHeaderParams(char *outfile)
{
	
	fitsfile *fptrOut;  
	int status,i,tstarti,tstopi;
	double tstart,tstop,start[4],stop[4],largest,smallest,exposure,tstartf,tstopf;
	status=0;

	if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
	         printerror( status );

    	for(i=0;i<4;i++)
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

	if ( fits_movabs_hdu(fptrOut, 1, NULL, &status) ) 
	      	printerror( status );
	tstarti=(int)smallest;
	tstopi=(int)largest;

	tstartf=smallest-tstarti;
	tstopf=largest-tstopi;
	exposure=largest-smallest;
//	printf("Time : %f\t%f\t%d\t%d\t%10f\t%10f\n",smallest,largest,tstarti,tstopi,tstartf,tstopf);
	
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
int read_input_parameters(int r, char *infile, char *mkffile, char *thresholdfile,char *outfile,int *clobber,int *history)
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
	
	r=PILGetFname("mkffile",mkffile);
	if (fits_open_file(&fptr, mkffile, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,mkffile);
		printerror( status );
	}
	fits_close_file(fptr,&status);

	r=PILGetFname("thresholdfile",thresholdfile);
	if((fp = fopen(thresholdfile, "r")) == NULL) {
		 printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,thresholdfile);
	}
	fclose(fp);
	r=PILGetFname("outfile",outfile);
	if (strchr(outfile, '.') != NULL) {
      	//outputfile = strtok(outfile, ".");
      		printf("\nWARNING (%s:%d): Please give the output file name without extension\n",__FILE__,__LINE__);
		exit(0);
  	}
		
	//---------------------------------------------------------Overwrite files and remember history------------------------------------------------------------
	r=PILGetBool("clobber", clobber);
	r=PILGetBool("history", history);
	//-------------------------------------------------------------End of input paramters----------------------------------------------------------------------	
	PILClose(r);	//Closing PIL file
	
	return(status);
}
int display_input_parameters(char *infile, char *mkffile, char *thresholdfile,char *outfile,int clobber,int history)
{
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf("                                                    CZTDATASELECTION PARAMETERS \n");
    	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Input file              		: %s\n",infile);
	printf(" MKF file 	             		: %s\n",mkffile);
    	printf(" Threshold file	              		: %s\n",thresholdfile);
	printf(" Output file 				: %s\n",outfile);
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


/*Author : Ajay Ratheesh, TIFR Mumbai
 * Date	 : 18 Sept 2017
 * This code is used to seperate single and double events
 * 
 * 
*/
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#define double_time 0.00003
#define BUFFSIZE 1024

void createEventFile(char *outputfile,char *eventfile);
void processeventseperation();
void printerror( int status);
void modifyEventHeaderParams(char *outputfile);
//void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum);
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum,double exposure);

int read_input_parameters(int r, char *infile,char *outfile_single,char *outfile_double,int *clobber,int *history);
int display_input_parameters(char *infile, char *outfile_single,char *outfile_double,int clobber,int history);

char *CZTNOISECLEAN,infile[BUFFSIZE],parfilename[BUFFSIZE], outfile_single[BUFFSIZE],outfile_double[BUFFSIZE];
int clobber=1, history=1;

int main(int argc, char **argv)
{
	int r=0;
	/*if(argc==2)
   	{
		//printf("%c\n",argv[2]);
		processeventseperation(argv[1]);
	}
	if(argc!=2)
	{
		printf("Enter all command line arguments\n1:Event file name\n");
		exit(-1);
	}*/
	CZTNOISECLEAN = getenv ("CZTNOISECLEAN");
	if(CZTNOISECLEAN==NULL)
	{
		printf("CZTNOISECLEAN Variable is not set\n");
		exit(0);
	}
	strcpy(parfilename,CZTNOISECLEAN);
	strcat(parfilename,"/paramfiles/czteventsep.par");
	PILSetModuleName(parfilename);
	//Read all inputs through Parameter Interface Library
	r=PILInit(argc,argv);

	read_input_parameters(r, infile, outfile_single, outfile_double, &clobber, &history);
	display_input_parameters(infile, outfile_single, outfile_double, clobber, history);

	
	processeventseperation();
	
	return 0;
}



void processeventseperation()
{
	fitsfile *fevt;
	
	long frow, felem;
	int status=0,hdutype=0;
	int double_flag;
	long evtnrows;
	int qid,i,j,k,l1,l2;
	int intnull,anynull,*evtpi,*finalpi_single,*finalpi_double,*redund_flag;
	unsigned char *evtdetid,*evtpixid,*evtdetx,*evtdety,*evtalpha,*finaldetid_double,*finalpixid_double,*finaldetx_double,*finaldety_double,*finalalpha_double,bytenull;
	unsigned char *finaldetid_single,*finalpixid_single,*finaldetx_single,*finaldety_single,*finalalpha_single;
	unsigned short *evtpha,*finalpha_single,*evtcztntick,*evtveto,*finalcztntick_single,*finalveto_single,*finalcztntick_double,*finalveto_double,*finalpha_double;
	float *evtenergy,floatnull,*finalenergy_single,*finalenergy_double;
	double *evttime,*evtcztseccnt,*finalevttime_single,*finalcztseccnt_single,doublenull,*finalevttime_double,*finalcztseccnt_double;
	//char outfile_single[1000],outfile_double[1000] ;
	
	if ( fits_open_file(&fevt, infile, READONLY, &status) )
	{
	    	printf("Error in opening a file : %d",status);
	   	printerror( status );
	}
	status=0;
	evtnrows=1;
	evttime  = (double*)malloc(evtnrows * sizeof(double));
	evtcztseccnt  = (double*)malloc(evtnrows * sizeof(double));
	evtpha  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	evtcztntick  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	evtveto  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	evtdetid  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	evtpixid  =(unsigned char*) malloc(evtnrows * sizeof(unsigned char));
	evtdetx  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	evtdety  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	evtalpha = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	evtenergy = (float*)malloc(sizeof(float)*evtnrows);
	evtpi = (int*)malloc(sizeof(int)*evtnrows);
	
	finalevttime_single  = (double*)malloc(evtnrows * sizeof(double));
	finalcztseccnt_single  = (double*)malloc(evtnrows * sizeof(double));
	finalpha_single  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	finalcztntick_single  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	finalveto_single  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	finaldetid_single  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finalpixid_single  =(unsigned char*) malloc(evtnrows * sizeof(unsigned char));
	finaldetx_single  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finaldety_single  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finalalpha_single = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finalenergy_single = (float*)malloc(sizeof(float)*evtnrows);
	finalpi_single = (int*)malloc(sizeof(int)*evtnrows);
	
	finalevttime_double  = (double*)malloc(evtnrows * sizeof(double));
	finalcztseccnt_double  = (double*)malloc(evtnrows * sizeof(double));
	finalpha_double  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	finalcztntick_double  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	finalveto_double = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
	finaldetid_double  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finalpixid_double  =(unsigned char*) malloc(evtnrows * sizeof(unsigned char));
	finaldetx_double  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finaldety_double = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finalalpha_double = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
	finalenergy_double = (float*)malloc(sizeof(float)*evtnrows);
	finalpi_double = (int*)malloc(sizeof(int)*evtnrows);

	if(clobber==1){
	remove(outfile_single);
	remove(outfile_double);
	}
	createEventFile(outfile_double,infile);
	createEventFile(outfile_single,infile);
	
	for(qid=0;qid<4;qid++)
	{
		fits_movabs_hdu(fevt, qid+2, &hdutype, &status);
		fits_get_num_rows(fevt, &evtnrows, &status);
	
		redund_flag = (int*)calloc(evtnrows,sizeof(int));
		
		evttime=(double*)realloc(evttime,sizeof(double)*evtnrows);
		evtdetid=(unsigned char*)realloc(evtdetid,sizeof(unsigned char)*evtnrows);
		evtpixid=(unsigned char*)realloc(evtpixid,sizeof(unsigned char)*evtnrows);
		evtdetx=(unsigned char*)realloc(evtdetx,sizeof(unsigned char)*evtnrows);
		evtdety=(unsigned char*)realloc(evtdety,sizeof(unsigned char)*evtnrows);
		evtcztseccnt  = (double*)realloc(evtcztseccnt,evtnrows * sizeof(double));
		evtpha  = (unsigned short*)realloc(evtpha,evtnrows * sizeof(unsigned short));
		evtcztntick  = (unsigned short*)realloc(evtcztntick,evtnrows * sizeof(unsigned short));
		evtveto  = (unsigned short*)realloc(evtveto,evtnrows * sizeof(unsigned short));
		evtalpha = (unsigned char*)realloc(evtalpha,evtnrows * sizeof(unsigned char));
		evtenergy = (float*)realloc(evtenergy,sizeof(float)*evtnrows);
		evtpi = (int*)realloc(evtpi,sizeof(int)*evtnrows);
		
		finalevttime_single  = (double*)realloc(finalevttime_single,evtnrows * sizeof(double));
		finalcztseccnt_single  = (double*)realloc(finalcztseccnt_single,evtnrows * sizeof(double));
		finalpha_single  = (unsigned short*)realloc(finalpha_single,evtnrows * sizeof(unsigned short));
		finalcztntick_single  = (unsigned short*)realloc(finalcztntick_single,evtnrows * sizeof(unsigned short));
		finalveto_single  = (unsigned short*)realloc(finalveto_single,evtnrows * sizeof(unsigned short));
		finaldetid_single  = (unsigned char*)realloc(finaldetid_single,evtnrows * sizeof(unsigned char));
		finalpixid_single  =(unsigned char*) realloc(finalpixid_single,evtnrows * sizeof(unsigned char));
		finaldetx_single  = (unsigned char*)realloc(finaldetx_single,evtnrows * sizeof(unsigned char));
		finaldety_single  = (unsigned char*)realloc(finaldety_single,evtnrows * sizeof(unsigned char));
		finalalpha_single = (unsigned char*)realloc(finalalpha_single,evtnrows * sizeof(unsigned char));
		finalenergy_single = (float*)realloc(finalenergy_single,sizeof(float)*evtnrows);
		finalpi_single = (int*)realloc(finalpi_single,sizeof(int)*evtnrows);

		finalevttime_double  = (double*)realloc(finalevttime_double,evtnrows * sizeof(double));
		finalcztseccnt_double  = (double*)realloc(finalcztseccnt_double,evtnrows * sizeof(double));
		finalpha_double  = (unsigned short*)realloc(finalpha_double,evtnrows * sizeof(unsigned short));
		finalcztntick_double  = (unsigned short*)realloc(finalcztntick_double,evtnrows * sizeof(unsigned short));
		finalveto_double  = (unsigned short*)realloc(finalveto_double,evtnrows * sizeof(unsigned short));
		finaldetid_double  = (unsigned char*)realloc(finaldetid_double,evtnrows * sizeof(unsigned char));
		finalpixid_double  =(unsigned char*) realloc(finalpixid_double,evtnrows * sizeof(unsigned char));
		finaldetx_double  = (unsigned char*)realloc(finaldetx_double,evtnrows * sizeof(unsigned char));
		finaldety_double  = (unsigned char*)realloc(finaldety_double,evtnrows * sizeof(unsigned char));
		finalalpha_double = (unsigned char*)realloc(finalalpha_double,evtnrows * sizeof(unsigned char));
		finalenergy_double = (float*)realloc(finalenergy_double,sizeof(float)*evtnrows);
		finalpi_double = (int*)realloc(finalpi_double,sizeof(int)*evtnrows);

		
		frow      = 1;
		felem     = 1;
		doublenull = 0.;
		intnull = 0;
		bytenull = 0;
		floatnull=0.0;
		
		fits_read_col(fevt, TDOUBLE, 1, frow, felem, evtnrows, &doublenull, evttime,&anynull, &status);
		fits_read_col(fevt, TDOUBLE, 2, frow, felem, evtnrows, &doublenull, evtcztseccnt, &anynull, &status);         
		fits_read_col(fevt, TUSHORT, 3, frow, felem, evtnrows, &intnull, evtcztntick, &anynull, &status); 
		fits_read_col(fevt, TUSHORT, 4, frow, felem, evtnrows, &intnull, evtpha,&anynull, &status);
		fits_read_col(fevt, TBYTE, 5, frow, felem,evtnrows, &bytenull, evtdetid, &anynull, &status);  
		fits_read_col(fevt, TBYTE, 6, frow, felem,evtnrows, &bytenull, evtpixid, &anynull, &status);        
		fits_read_col(fevt, TBYTE, 7, frow, felem, evtnrows, &bytenull, evtdetx, &anynull, &status);   
		fits_read_col(fevt, TBYTE, 8, frow, felem,evtnrows, &bytenull, evtdety,&anynull, &status);
		fits_read_col(fevt, TUSHORT, 9, frow, felem,evtnrows, &intnull, evtveto,&anynull, &status);
		fits_read_col(fevt, TBYTE, 10, frow, felem, evtnrows, &bytenull, evtalpha,&anynull, &status); 
		fits_read_col(fevt, TINT, 11, frow, felem, evtnrows, &intnull, evtpi,&anynull, &status);
		fits_read_col(fevt, TFLOAT, 12, frow, felem, evtnrows, &floatnull, evtenergy,&anynull, &status);
		
		l1=0;
		l2=0;
		//printf("%ld\n",evtnrows);
		for(i=0;i<evtnrows;i++)
		{
			//if(evttime[i+1]<=evttime[i])printf("%lf\t%lf\n",evttime[i],evttime[i+1]);
			//printf("%d\t%d\n",qid,i);
			double_flag = 0; 
			for(j=i+1; j<evtnrows;j++)
			{
				if( evttime[j]>=evttime[i] && evttime[j]<= evttime[i]+double_time )
				{
					double_flag = 1;
					break;
				}
				if(evttime[j] > evttime[i]+1.0){break;}
			}
			
			if(double_flag==1)
			{
				if(redund_flag[i]==0 && redund_flag[j] == 0)
				{
					finalevttime_double[l1]=evttime[i];
					finalcztseccnt_double[l1]  = evtcztseccnt[i];
					finalpha_double[l1]  = evtpha[i];
					finalcztntick_double[l1]  = evtcztntick[i];
					finalveto_double[l1]  = evtveto[i];
					finaldetid_double[l1]  = evtdetid[i];
					finalpixid_double[l1]  = evtpixid[i];
					finaldetx_double[l1]  = evtdetx[i];
					finaldety_double[l1]  = evtdety[i];
					finalalpha_double[l1] = evtalpha[i];
					finalenergy_double[l1] = evtenergy[i];
					finalpi_double[l1] = evtpi[i];
					l1++;
					
					finalevttime_double[l1]=evttime[j];
					finalcztseccnt_double[l1]  = evtcztseccnt[j];
					finalpha_double[l1]  = evtpha[j];
					finalcztntick_double[l1]  = evtcztntick[j];
					finalveto_double[l1]  = evtveto[j];
					finaldetid_double[l1]  = evtdetid[j];
					finalpixid_double[l1]  = evtpixid[j];
					finaldetx_double[l1]  = evtdetx[j];
					finaldety_double[l1]  = evtdety[j];
					finalalpha_double[l1] = evtalpha[j];
					finalenergy_double[l1] = evtenergy[j];
					finalpi_double[l1] = evtpi[j];
					l1++;		
					redund_flag[i]=1;
					redund_flag[j]=1;	
				}
			}
			
			else
			{
				if(redund_flag[i]==0)
				{
					finalevttime_single[l2]=evttime[i];
					finalcztseccnt_single[l2]  = evtcztseccnt[i];
					finalpha_single[l2]  = evtpha[i];
					finalcztntick_single[l2]  = evtcztntick[i];
					finalveto_single[l2]  = evtveto[i];
					finaldetid_single[l2]  = evtdetid[i];
					finalpixid_single[l2]  = evtpixid[i];
					finaldetx_single[l2]  = evtdetx[i];
					finaldety_single[l2]  = evtdety[i];
					finalalpha_single[l2] = evtalpha[i];
					finalenergy_single[l2] = evtenergy[i];
					finalpi_single[l2] = evtpi[i];
					l2++;
				}
			}
		}
		//printf("%d\t%d\t%d\t%d\n",l1,l2,l1+l2,evtnrows);
		//printf("%d\t%d\t%lf\n",sizeof(finalevttime_single),l2, finalevttime_single[l2-1]);
		double tot_exposure;
		fits_read_key(fevt,TDOUBLE,"EXPOSURE",&tot_exposure,NULL, &status);
		
		
		writeEvent(finalevttime_single,finalcztseccnt_single,finalcztntick_single,finalpha_single,finaldetid_single,finalpixid_single,finaldetx_single,finaldety_single,finalveto_single,finalalpha_single,finalpi_single,finalenergy_single,outfile_single,l2,qid+2,tot_exposure);
		writeEvent(finalevttime_double,finalcztseccnt_double,finalcztntick_double,finalpha_double,finaldetid_double,finalpixid_double,finaldetx_double,finaldety_double,finalveto_double,finalalpha_double,finalpi_double,finalenergy_double,outfile_double,l1,qid+2,tot_exposure);
		free(redund_flag);	
			
	}//end of qid
	modifyEventHeaderParams(outfile_single);
	modifyEventHeaderParams(outfile_double);

	if ( fits_close_file(fevt, &status) )       
	        printerror( status );

	//free memory
	free(evttime);free(evtcztseccnt);free(evtcztntick);free(evtveto);free(evtpha);
	free(evtdetid);free(evtpixid);free(evtdetx);free(evtdety);free(evtalpha);free(evtenergy);free(evtpi);
	free(finalcztntick_single);free(finalveto_single);free(finalevttime_single);free(finalcztseccnt_single);
        free(finalpha_single);
	free(finaldetid_single);free(finalpixid_single);free(finaldetx_single);free(finaldety_single);
	free(finalalpha_single);free(finalenergy_single);free(finalpi_single);
	free(finalcztntick_double);free(finalveto_double);free(finalevttime_double);free(finalcztseccnt_double);
	free(finalpha_double);free(finaldetid_double);free(finalpixid_double);free(finaldetx_double);
	free(finaldety_double);free(finalalpha_double);free(finalenergy_double);free(finalpi_double);
		printf("\nEVENT SEPERATION COMPLETED SUCCESSFULLY.\n");
	return;	

}

//creating event file 
void createEventFile(char *outputfile,char *eventfile)
{
	
	fitsfile *fptrOut=NULL,*fptrevt=NULL;      
	int status=0, hdutype,anynull,tfields=12,i,hdunum=2;
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

	// Copy Exposure
	if(fits_movnam_hdu(fptrevt, BINARY_TBL, "EXPOSURE", 0, &status)) 
		printerror( status );         

	if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		printerror( status );

	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );

	if ( fits_close_file(fptrevt, &status) )       
	        printerror( status );
	return;
}

void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i,tstarti,tstopi;
	double tstart,tstop,start[4],stop[4],largest,smallest,exposure,tstartf,tstopf;
	status=0;
	char *comment; 

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

	if ( fits_movabs_hdu(fptrOut, 1, NULL, &status) ) 
	      	printerror( status );
	tstarti=(int)smallest;
	tstopi=(int)largest;

	tstartf=smallest-tstarti;
	tstopf=largest-tstopi;
	exposure=largest-smallest;
	printf("Time : %f\t%f\t%d\t%d\t%10f\t%10f\n",smallest,largest,tstarti,tstopi,tstartf,tstopf);
	
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

void printerror( int status)
{
	if (status)
	{
		fits_report_error(stderr, status); 	
		exit( status );    	
    	}
  	return;
}


void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum,double exposure)
{
	fitsfile *fptrOut;       
	int status, hdutype,intnull;
	long frow, felem,nrows;
	unsigned char bytenull;
	int tstarti,tstopi;
	double doublenull,tstart,tstop,tstartf,tstopf;
	status=0;

	
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
	//printf("Writing %d...\n",hdunum-2);	

	frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;
	
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

int read_input_parameters(int r, char *infile,char *outfile_single,char *outfile_double,int *clobber,int *history)
{
	int status=0,hdutype=0;
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
	
	r=PILGetFname("outfile_single",outfile_single);
	r=PILGetFname("outfile_double",outfile_double);
	//---------------------------------------------------------Overwrite files and remember history------------------------------------------------------------
	r=PILGetBool("clobber", clobber);
	r=PILGetBool("history", history);
	//-------------------------------------------------------------End of input paramters----------------------------------------------------------------------	
	PILClose(r);	//Closing PIL file
	
	return(status);
}


int display_input_parameters(char *infile, char *outfile_single,char *outfile_double,int clobber,int history)
{
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf("                                                    CZTEVENTSEPERATION PARAMETERS \n");
    	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Input file              		: %s\n",infile);
	printf(" Single event output file  		: %s\n",outfile_single);
	printf(" Double event output file  		: %s\n",outfile_double);
    	printf(" Clobber				: %d\n",clobber);
    	printf(" History				: %d\n",history);
    	printf("----------------------------------------------------------------------------------------------------------------------------\n\n");
    	return 0;
}

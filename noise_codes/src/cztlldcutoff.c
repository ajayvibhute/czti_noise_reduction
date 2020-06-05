/* CZTLLDCUTOFF
 * Author	: Ajay Ratheesh, TIFR Mumbai
		
 * Date 	: 18-08-2017
 * This code is used to put an artificial LLD for pixels and Detectors
 */

#include<stdio.h>
#include<math.h>
#include<fitsio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>

//Following are the parameters
#define default_delta_thresh 3.0
#define BUFFSIZE 1024
void printerror( int status);
void processpixLLDcutoff();
void createEventFile(char *outputfile,char *eventfile);
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int bufsize,int hdunum);
void modifyEventHeaderParams(char *outputfile);
int read_input_parameters(int r, char *infile,char *caldb_lld_file,char *detthresholdfile,char *pixthresholdfile,char *outfile,int *clobber,int *history);
int display_input_parameters(char *infile,char *caldb_lld_file,char *detthresholdfile,char *pixthresholdfile,char *outfile,int clobber,int history);

char *CZTNOISECLEAN,infile[BUFFSIZE],parfilename[BUFFSIZE], caldb_lld_file[BUFFSIZE],outfile[BUFFSIZE],detthresholdfile[BUFFSIZE],pixthresholdfile[BUFFSIZE];
int clobber=1, history=1;

int main(int argc,char *argv[])
{

	/*if(argc==5)
        {
		processpixLLDcutoff(argv[1],argv[2],argv[3],argv[4]);
	}
	if(argc<5)
	{
		printf("Enter all command line arguments\n1:Event file name\n2.CALDB lld file\n3.Detid threshold file\n4.Pixid threshold file\n");
		exit(-1);
	}*/
	CZTNOISECLEAN = getenv ("CZTNOISECLEAN");
	if(CZTNOISECLEAN==NULL)
	{
		printf("CZTNOISECLEAN Variable is not set\n");
		exit(0);
	}
	strcpy(parfilename,CZTNOISECLEAN);
	strcat(parfilename,"/paramfiles/cztlldcutoff.par");
	PILSetModuleName(parfilename);
	//Read all inputs through Parameter Interface Library
	int r=PILInit(argc,argv);
	read_input_parameters(r, infile, caldb_lld_file, detthresholdfile, pixthresholdfile, outfile, &clobber, &history);
	display_input_parameters(infile,caldb_lld_file,detthresholdfile,pixthresholdfile, outfile,clobber,history);
	processpixLLDcutoff();

	return 0;

}


void processpixLLDcutoff()
{
	int *evt_flag;
	FILE *logfile;
	long evtnrows, lldnrows;
	int quad1[16384],quad2[16384];
	FILE *f1, *f2;
	long frow, felem;
	fitsfile *fevt, *flld;	
	int status=0,hdutype=0,size=1;
	int qid,i,j,k,*pilld;
	double *evttime,*evtcztseccnt,*finalevttime,*finalcztseccnt,doublenull;
	unsigned short *evtpha,*finalpha,*evtcztntick,*evtveto,*finalcztntick,*finalveto;
	int intnull,anynull,*evtpi,*finalpi;
	unsigned char *evtdetid,*evtpixid,*evtdetx,*evtdety,*evtalpha,*finaldetid,*finalpixid,*finaldetx,*finaldety,*finalalpha,bytenull;
	unsigned char *pixidlld,*detidlld, detid_input[64], detx_input[16384],dety_input[16384],detx_in,dety_in,*caldb_detx,*caldb_dety;
	float *evtenergy,floatnull,*finalenergy;
	float LLD[64][64], LLD_thresh[64][64],det_thresh_input[64],pix_thresh_input[16384];
	char outlogfile[1000];
	int det,pix,row,counter1=0,counter2=0;


	f1=fopen(detthresholdfile,"rt"); 
	f2=fopen(pixthresholdfile,"rt");

		evttime=(double*)malloc(sizeof(double)*evtnrows);
		evtdetid=(unsigned char*)malloc(sizeof(unsigned char)*size);
		evtpixid=(unsigned char*)malloc(sizeof(unsigned char)*size);
		evtdetx=(unsigned char*)malloc(sizeof(unsigned char)*size);
		evtdety=(unsigned char*)malloc(sizeof(unsigned char)*size);
		evtcztseccnt  = (double*)malloc(size * sizeof(double));
		evtpha  = (unsigned short*)malloc(size * sizeof(unsigned short));
		evtcztntick  = (unsigned short*)malloc(size * sizeof(unsigned short));
		evtveto  = (unsigned short*)malloc(size * sizeof(unsigned short));
		evtalpha = (unsigned char*)malloc(size * sizeof(unsigned char));
		evtenergy = (float*)malloc(sizeof(float)*size);
		evtpi = (int*)malloc(sizeof(int)*size);
		evt_flag=(int*)malloc(sizeof(int)*(size));

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

		detidlld=(unsigned char*)malloc(sizeof(unsigned char)*size);
		pixidlld=(unsigned char*)malloc(sizeof(unsigned char)*size);
		caldb_detx  = (unsigned char*)malloc(size * sizeof(unsigned char));
		caldb_dety  = (unsigned char*)malloc(size * sizeof(unsigned char));
		pilld = (int*)malloc(sizeof(int)*size);

	while ((fscanf(f1,"%d",&quad1[counter1])!=EOF))
 	{
		fscanf(f1,"%c",&detid_input[counter1]);
		fscanf(f1,"%f",&det_thresh_input[counter1]);
		counter1++;
	}
	while ((fscanf(f2,"%d",&quad2[counter2]))!=EOF)
 	{
		fscanf(f2,"%c",&detx_input[counter2]);
		fscanf(f2,"%c",&dety_input[counter2]);
		fscanf(f2,"%f",&pix_thresh_input[counter2]);
		counter2++;
	}



	fits_open_file(&fevt,infile,READONLY,&status);
		if(status) { printerror( status );}

	fits_open_file(&flld,caldb_lld_file,READONLY,&status);
		if(status) { printerror( status );}

	char* tempevt = calloc(strlen(infile)+1, sizeof(char));
	strcpy(tempevt, infile);

	char *file = strtok(infile, "."); 
	//sprintf(outfile, "%s_lldcut.evt",file);
	if(clobber==1) remove(outfile);
	createEventFile(outfile,tempevt);

	sprintf(outlogfile, "%s_lldcut.log",file);
	remove(outlogfile);
	logfile=fopen(outlogfile,"w");

	fprintf(logfile,"Inputs for LLD cut off\n1.Event file : %s\n2.LLD file : %s\n3.Det Threshold file : %s\n3.Pix Threshold file : %s\n\n",tempevt,caldb_lld_file,detthresholdfile,pixthresholdfile);
	fprintf(logfile,"Important thresholds ....\n1.default_delta_thresh : %f\n\n",default_delta_thresh); 

        frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;
	floatnull=0.0;

	for(qid=0;qid<4;qid++)
	{
		fprintf(logfile,"Quad %d is processing...\n",qid);
		fits_movabs_hdu(fevt, qid+2, &hdutype, &status);
		fits_movabs_hdu(flld, qid+2, &hdutype, &status);

		fits_get_num_rows(fevt, &evtnrows, &status);
		fits_get_num_rows(flld, &lldnrows, &status);


		if(evtnrows==0)
   		 {
        		 printf("No events in the event file\n");
       		 	 continue;
		 }
    
  	        if(lldnrows==0)
                {
	   		printf("No events in the CALDB lld file\n");
       			continue;
	   	}

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
		evt_flag=(int*)realloc(evt_flag,sizeof(int)*(evtnrows));

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


		detidlld=(unsigned char*)realloc(detidlld,sizeof(unsigned char)*lldnrows);
		pixidlld=(unsigned char*)realloc(pixidlld,sizeof(unsigned char)*lldnrows);
		caldb_detx  = (unsigned char*)realloc(caldb_detx,lldnrows * sizeof(unsigned char));
		caldb_dety  = (unsigned char*)realloc(caldb_dety,lldnrows * sizeof(unsigned char));
		pilld = (int*)realloc(pilld,sizeof(int)*lldnrows);

		
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

		fits_read_col(flld, TBYTE, 1, 1, 1, lldnrows, NULL, detidlld, NULL, &status);
		fits_read_col(flld, TBYTE, 2, 1, 1, lldnrows, NULL, pixidlld, NULL, &status);

		for(i=0;i<64;i++)
		{
			
			for(j=0;j<64;j++)
			{
				LLD_thresh[i][j] = default_delta_thresh;
			}

		}
		int colLLD[5];
		for(det=0;det<16;det++)  // For each detector in the quadrant
        	{
			
			for (pix=0;pix<256;pix++)   // For each pixel in the module
            		{
				row=det*256+pix; 
				fits_read_col(flld,TINT,3,row+1,1,5,NULL,colLLD,NULL,&status);
				pilld[row]=colLLD[1];

				caldb_detx[row]=((detidlld[row]%4)*16)+(pixidlld[row]%16);
				caldb_dety[row]=((detidlld[row]/4)*16)+(pixidlld[row]/16);
			
				if(qid==0 || qid==3)
					caldb_dety[row]=63-caldb_dety[row];
				else
					caldb_detx[row]=63-caldb_detx[row];
				LLD[caldb_detx[row]][caldb_dety[row]]=(pilld[row]*0.5)+5.0;
			}

		}



		for(i=0;i<counter1;i++)
		{
				
			if(qid==quad1[i])
			{
				for(j=0;j<256;j++)
				{
					detx_in = ((detid_input[i]%4)*16)+(j%16);
					dety_in = ((detid_input[i]/4)*16)+(j/16);
					if(qid==0 || qid==3)
					dety_in=63-dety_in;
					else
					detx_in=63-detx_in;


					LLD_thresh[detx_in][dety_in] = det_thresh_input[i];//convert detid and pixid to detx dety
				}
			}
		}
		for(i=0;i<counter2;i++)
		{
			if(qid==quad1[i])
			{
      				LLD_thresh[detx_input[i]][dety_input[i]] = pix_thresh_input[i];//convert detid and pixid to detx dety
			}
		}




		/*for(i=0;i<lldnrows;i++)
		{
			LLD[detx conversion][dety conversion] = conversion of PI to energy;//convert detid and pixid to detx dety
		}*/

		
		for(i=0;i<evtnrows;i++)
		{
			
			evt_flag[i] = 0;
			
		}
		for(i=0;i<evtnrows;i++)
		{
			
			if(evtenergy[i]<(LLD[evtdetx[i]][evtdety[i]]+LLD_thresh[evtdetx[i]][evtdety[i]]))
			{
				evt_flag[i] = 1;
			}
		}
		

		int l=0;
//sum=0;
		fprintf(logfile,"Events removed after applying LLD cut off\nIndex\tTime\tDETX\tDETY\n");
		for(i=0;i<evtnrows;i++)
		{
			if(evt_flag[i]==0)
			{
				finalevttime[l]=evttime[i];
				finalcztseccnt[l]  = evtcztseccnt[i];
				finalpha[l]  = evtpha[i];
				finalcztntick[l]  = evtcztntick[i];
				finalveto[l]  = evtveto[i];
				finaldetid[l]  = evtdetid[i];
				finalpixid[l]  = evtpixid[i];
				finaldetx[l]  = evtdetx[i];
				finaldety[l]  = evtdety[i];
				finalalpha[l] = evtalpha[i];
				finalenergy[l] = evtenergy[i];
				finalpi[l] = evtpi[i];
				l++;
			}
			else if(evt_flag[i]==1)
					fprintf(logfile,"%d\t%f\t%u\t%u\n",i,evttime[i],evtdetx[i],evtdety[i]);
		
		}
		//printf("true\n");
		writeEvent(finalevttime,finalcztseccnt,finalcztntick,finalpha,finaldetid,finalpixid,finaldetx,finaldety,finalveto,finalalpha,finalpi,finalenergy,outfile,l,qid+2);
	
	}//qid end
	modifyEventHeaderParams(outfile);

	free(pilld);free(evttime);free(evtcztseccnt);free(finalevttime);free(finalcztseccnt);
	free(evtpha);free(finalpha);free(evtcztntick);free(evtveto);free(finalcztntick);free(finalveto);
	free(evtpi);free(finalpi);free(evtdetid);free(evtpixid);free(evtdetx);free(evtdety);free(evtalpha);
	free(finaldetid);free(finalpixid);free(finaldetx);free(finaldety);free(finalalpha);
	free(pixidlld);free(detidlld);free(caldb_detx);free(caldb_dety);
	free(evtenergy);free(finalenergy);free(evt_flag);

	printf("\nPIXEL LLD CUTOFF COMPLETED SUCCESSFULLY.\n");

}

//creating event file 
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

void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype,intnull;
	long frow, felem,nrows;
	int tstarti,tstopi;
	double tstart,tstop,tstartf,tstopf,exposure,doublenull;
	unsigned char bytenull;
	status=0;

	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	frow      = 1;//firstindex;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;

	tstart=evttime[0];
	tstop=evttime[writesize-1];

	tstarti=(int)tstart;
	tstopi=(int)tstop;

	tstartf=tstart-tstarti;
	tstopf=tstop-tstopi;
	exposure=tstop-tstart;

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

void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i,tstarti,tstopi;
	double tstart,tstop,start[4],stop[4],largest,smallest,exposure,tstartf,tstopf;
	status=0;
	char *comment; 

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
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


int read_input_parameters(int r, char *infile,char *caldb_lld_file,char *detthresholdfile,char *pixthresholdfile,char *outfile,int *clobber,int *history)
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
	
	r=PILGetFname("caldb_lld_file",caldb_lld_file);
	if (fits_open_file(&fptr, caldb_lld_file, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,caldb_lld_file);
		printerror( status );
	}
	fits_close_file(fptr,&status);

	r=PILGetFname("detthresholdfile",detthresholdfile);
	if((fp = fopen(detthresholdfile, "r")) == NULL) {
		 printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,detthresholdfile);
	}

	r=PILGetFname("pixthresholdfile",pixthresholdfile);
	if((fp = fopen(pixthresholdfile, "r")) == NULL) {
		 printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,pixthresholdfile);
	}

	r=PILGetFname("outfile",outfile);
	//---------------------------------------------------------Overwrite files and remember history------------------------------------------------------------
	r=PILGetBool("clobber", clobber);
	r=PILGetBool("history", history);
	//-------------------------------------------------------------End of input paramters----------------------------------------------------------------------	
	PILClose(r);	//Closing PIL file
	
	return(status);
}


int display_input_parameters(char *infile,char *caldb_lld_file,char *detthresholdfile,char *pixthresholdfile,char *outfile,int clobber,int history)
{
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf("                                                    CZTLLDCUTOFF PARAMETERS \n");
    	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Input file              		: %s\n",infile);
	printf(" CALDB lld file  			: %s\n",caldb_lld_file);
	printf(" Detector threshold file  		: %s\n",detthresholdfile);
	printf(" Pixel threshold file  			: %s\n",pixthresholdfile);
	printf(" Output file  				: %s\n",outfile);
    	printf(" Clobber				: %d\n",clobber);
    	printf(" History				: %d\n",history);
    	printf("----------------------------------------------------------------------------------------------------------------------------\n\n");
    	return 0;
}


/*
btigen.c
Ajay Vibhute, 26 Aug 2017.
code to generate bti files
inputs:- input event file
output:- bti file
//find all events in a second, if it crosses 3024, end of the second and +20microsecond to last event in that time bin
*/
#include<stdio.h>
#include<fitsio.h>
#include<string.h>
#include<math.h>
#define NUMQUADRANT 4
#define BUFSIZE 4096
int main(int argc,char*argv[])
{
	char infile[BUFSIZE], outfile[BUFSIZE],bunchfile[BUFSIZE],extname[100];
	fitsfile *fptr,*fpout,*fpbunch;	
	int status=0,hdunum=0,hdutype=0,timecol=0,cztseccntcol=0,i=0,lasttime=0,bincount=0,tfields=2,btirows=0;
	int *bunchlc,bunchnumbins;	
	long nrows=0,bunchrows=0;
	double *time,*cztseccnt,*bunchtime,bunchstart;
	char *ttype[] = { "TSTART", "TSTOP"};
    	char *tform[] = { "1D",     "1D" };
    	char *tunit[] = { "s",      "s"};


	if(argc<4)
	{
		printf("Enter all command line arguments\n 1: Input event file name\n2. Bunch file \n3. Output BTI file name\n");
		exit(-1);
	}
	strcpy(infile,argv[1]);
	strcpy(bunchfile,argv[2]);
	strcpy(outfile,argv[3]);
	printf("Input file name:%s\n",infile);
	printf("Output file name:%s\n",outfile);	
	if ( fits_open_file(&fptr, infile, READONLY, &status) )
	{
		printf("(%s:%d): Unable to open the input file\n",__FILE__,__LINE__);
		exit(0);
	}

	if ( fits_open_file(&fpbunch, bunchfile, READONLY, &status) )
	{
		printf("(%s:%d): Unable to open the input file\n",__FILE__,__LINE__);
		exit(0);
	}
	remove(outfile);
    	if(fits_create_file(&fpout,outfile,&status))
   	{
		printf("(%s:%d): Unable to create the output file\n",__FILE__,__LINE__);
        	fits_report_error(stderr,status);
        	return (EXIT_FAILURE);
    	}

	for(hdunum=3;hdunum<4;hdunum++)
	{
	        sprintf(extname,"Q%d",hdunum);
        	if(fits_create_tbl(fpout,BINARY_TBL,0,tfields,ttype,tform,tunit,extname,&status))
		{
			printf("(%s:%d): Unable to create the hdu\n",__FILE__,__LINE__);
		
        		fits_report_error(stderr,status);  
			return (EXIT_FAILURE); 

		}

		if ( fits_movabs_hdu(fpbunch, hdunum+2, &hdutype, &status) )
		{
			printf("(%s:%d): Unable to move to hdu no %d of bunch file\n",__FILE__,__LINE__,hdunum);
			exit(-1);
		}
		if(fits_get_num_rows(fpbunch,&bunchrows,&status) )
		{
	
			printf("Error (%s:%d): Error while getting number of rows in bunch file\n",__FILE__,__LINE__);
			exit(-1);
		}
		bunchtime=(double*)malloc(sizeof(double)*bunchrows);
		if(bunchtime==NULL)
		{
			printf("%s:%d: Unable to allocate memory\n",__FILE__,__LINE__);
			exit(-1);
		}

		if(fits_read_col(fpbunch, TDOUBLE, 1,1, 1, bunchrows, NULL, bunchtime,NULL, &status))
		{
			printf("Error (%s:%d): Error while reading time column from bunch file\n",__FILE__,__LINE__);
			exit(-1);
		}

		if ( fits_movabs_hdu(fptr, hdunum+2, &hdutype, &status) )
		{
			printf("(%s:%d): Unable to move to hdu no %d\n",__FILE__,__LINE__,hdunum);
			exit(-1);
		}
		if (hdutype != BINARY_TBL)
		{
			printf("Error (%s:%d): %d th extension supposed to be binary table, but it is not. Exiting...\n",__FILE__,__LINE__,hdunum);
			exit(-1);
		}
		if(fits_get_num_rows(fptr,&nrows,&status) )
		{
	
			printf("Error (%s:%d): Error while getting number of rows\n",__FILE__,__LINE__);
			exit(-1);
		}
		time=malloc(sizeof(double)*nrows);
		cztseccnt=malloc(sizeof(double)*nrows);
		if(time==NULL || cztseccnt==NULL)
		{
			printf("Error (%s:%d): Unable to allocate the memory\n",__FILE__,__LINE__);
			exit(-1);
		}
		if(fits_get_colnum(fptr,CASEINSEN,"Time",&timecol,&status))
		{
			printf("Error (%s:%d): Error while getting time column number\n",__FILE__,__LINE__);
			exit(-1);
		}

		if(fits_get_colnum(fptr,CASEINSEN,"CZTSECCNT",&cztseccntcol,&status))
		{
			printf("Error (%s:%d): Error while getting time column number\n",__FILE__,__LINE__);
			exit(-1);
		}
		if(fits_read_col(fptr, TDOUBLE, timecol,1, 1, nrows, NULL, time,NULL, &status))
		{
			printf("Error (%s:%d): Error while reading time column \n",__FILE__,__LINE__);
			exit(-1);
		}

		if(fits_read_col(fptr, TDOUBLE, cztseccntcol,1, 1, nrows, NULL, cztseccnt,NULL, &status))
		{
			printf("Error (%s:%d): Error while reading time column \n",__FILE__,__LINE__);
			exit(-1);
		}
		lasttime=(int)cztseccnt[0];
		bincount=0;
		long binno=0,numbins=cztseccnt[nrows-1]-cztseccnt[0];
		int *lc=malloc(sizeof(int)*nrows);
		int tmp=0,tmod=0,lasttmod=0,btiindex=0,k=0,br=0,bcount=0;
		double *btistart,*btistop,tmptime=0;
		for (i=0;i<nrows;i++)
			lc[i]=0;
		btirows=0;
		btistart=(double*)malloc(sizeof(double)*1025);
		btistop=(double*)malloc(sizeof(double)*1025);
		br=0;
		int lastbrindex=0,tmpcount=0,tmpbcount=0;
		for(i=1;i<nrows;i++)
		{	
			if(tmptime!=time[i])
			{
				for(br=lastbrindex;bunchtime[br]<=time[i];br++)
				{

					if(bunchtime[br]==time[i])
					{
						bincount++;
						lastbrindex=br+1;
						tmptime=time[i];
						break;

					}
				}
			}
			
			if(lasttime==(long)cztseccnt[i])
			{
				
				bincount++;
			}
			else
			{

				bincount=0;
				bcount=0;
				tmpcount=0;
			}

			if(bincount>3072)
			{
				while(lasttime==(long)cztseccnt[i])
				{
					i++;
				}
				i--;
			//	printf("Bad End %f\t%ld\t%d\t%d\n",cztseccnt[i],(long)cztseccnt[i+1],lasttime,bincount);
			}
			lasttime=(long)cztseccnt[i];
			
			tmod=(long)cztseccnt[i]%100;
			if(lasttmod!=tmod && tmod==0)
			{
				btistart[btiindex]=time[i];
				btistop[btiindex]=time[i]+0.4;
				btiindex++;
			}
			lasttmod=tmod;
			if(btiindex>=1024)
			{

				
				fits_write_col(fpout,TDOUBLE,1,btirows+1,1,btiindex,btistart,&status);
				fits_write_col(fpout,TDOUBLE,2,btirows+1,1,btiindex,btistop,&status);
				btirows+=btiindex;
				btiindex=0;

			}

		}
		if(btiindex!=0)
		{	

			fits_write_col(fpout,TDOUBLE,1,btirows+1,1,btiindex,btistart,&status);
			fits_write_col(fpout,TDOUBLE,2,btirows+1,1,btiindex,btistop,&status);
			btirows+=btiindex;
			btiindex=0;
			
		}
	//	printf("%d\t%d\n",hdunum,btirows);	
	}
    	if(fits_close_file(fpout,&status))
	{
		printf("(%s:%d): Unable to close the output file\n",__FILE__,__LINE__);
        	fits_report_error(stderr,status);  
		return (EXIT_FAILURE); 

	}
	
    	if(fits_close_file(fpbunch,&status))
	{
		printf("(%s:%d): Unable to close the bunch file file\n",__FILE__,__LINE__);
        	fits_report_error(stderr,status);  
		return (EXIT_FAILURE); 

	}
}

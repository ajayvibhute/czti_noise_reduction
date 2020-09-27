/* CZTSUPERBUNCHCLEAN
 * Author	: Ajay Ratheesh, TIFR Mumbai
			: Mayuri Shinde, IUCAA Pune
			: Ajay Vibhute, IUCAA Pune
 * Date 	: 14-08-2017
 * 
 * Updated by Ajay Ratheesh on 20th Dec 2017
 */

#include<stdio.h>
#include<math.h>
#include<fitsio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>
#include "pil.h"
//Following are the parameters
/*#define heavy_bun_thresh_count_10ms 4
#define heavy_bun_thresh_count_1ms 4
#define bun_size_thresh 200
#define dph_time 0.1
#define thresh_hot 0.707
#define allowable_hot_thresh 3.0*/
#define BUFFSIZE 1024

void printerror( int status);
void processSuperBunchClean();
void createEventFile(char *eventfile);
void createLivetimeFile();
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,int bufsize,int hdunum,double exposure);
void writeLivetime(double *time,double *livetime,int livetimesize,int hdunum);
void livetimeGeneration(double *live_time_UT,int live_counter,char *outputfile,int hdunum,double *lttime,double *fracexp,int livetimenrows);
void modifyEventHeaderParams(char *outputfile);
void writeBTItime(double *time_1,double *time_2,int BTIsize,int hdunum);
void createBTIFile();
void writegtiextension(double *gtistart, double *gtistop, double gtinrows, char extention_name[20], double tstart, double tstop, double exposure);
//~ void writegtiextension(double *gtistart, double *gtistop, double gtinrows, int hdunum, double tstart, double tstop, double exposure);
//void livetimeGeneration(double *live_time_UT,int live_counter,char *outputfile,int hdunum,double *lttime,double *fracexp,int livetimenrows,double tstart,double tstop);

int display_input_parameters(char *infile, char *inbunchfile, char *thresholdfile,char *inlivetimefile,char *outlivetimefile,char *outbtifile,char *outfile,int clobber,int history);
int read_input_parameters(int r,char *infile, char *inbunchfile, char *thresholdfile,char *inlivetimefile,char *outlivetimefile,char *outbtifile,char *outfile,int *clobber,int *history);

int heavy_bun_thresh_count_10ms,heavy_bun_thresh_count_1ms,bun_size_thresh;
float dph_time,thresh_hot,allowable_hot_thresh;

char *CZTNOISECLEAN,infile[BUFFSIZE],parfilename[BUFFSIZE], inbunchfile[BUFFSIZE],inlivetimefile[BUFFSIZE],outlivetimefile[BUFFSIZE],outbtifile[BUFFSIZE],outfile[BUFFSIZE],thresholdfile[BUFFSIZE];
int clobber=1, history=1;

int main(int argc,char *argv[])
{

	/*if(argc==5)
        {
		processSuperBunchClean(argv[1],argv[2],argv[3],argv[4]);
	}
	if(argc!=5)
	{
		printf("Enter all command line arguments\n1:Event file name\n2.Bunch file name\n3.Livetime file\n4.Orbit number\n");
		exit(-1);
	}*/
	CZTNOISECLEAN = getenv ("CZTNOISECLEAN");
	if(CZTNOISECLEAN==NULL)
	{
		printf("CZTNOISECLEAN Variable is not set\n");
		exit(0);
	}
	strcpy(parfilename,CZTNOISECLEAN);
	strcat(parfilename,"/paramfiles/cztsuperbunchclean.par");
	PILSetModuleName(parfilename);
	//Read all inputs through Parameter Interface Library
	int r=PILInit(argc,argv);
	//Read all input files & all the selection parameters
	read_input_parameters(r,infile,inbunchfile,thresholdfile,inlivetimefile,outlivetimefile,outbtifile,outfile,&clobber,&history);
	display_input_parameters(infile,inbunchfile,thresholdfile,inlivetimefile,outlivetimefile,outbtifile,outfile,clobber,history);
	processSuperBunchClean();
	
	return 0;

}

void processSuperBunchClean()
{
	int tmpindex=0,tmpindex1=0;
	unsigned int i,j,k,m,e,n,p;//,sum;
	unsigned int qid;
	//int evttime_col,evtdetx_col,evtdety_col,evtPI_col,evtdetid_col,evtpixid_col,evttime_index;
	int buntime_col,bunsize_col,size=1,buntime_df1_col,buntime_df2_col;
	double *buntime,*buntime_heavy,*buntime_real;
	int *bunsize_real;
	//double bun;
	int *evt_flag,*bun_flag;
	int *bunsize;
	long evtnrows, bunnrows,livetimenrows;
	int status=0,hdutype=0;
	int magnitude_temp[4096],detx_temp[4096],dety_temp[4096];
	int detx_close_bunch[3072],dety_close_bunch[3072];
	long frow, felem;
	double *evttime,*evtcztseccnt,*finalevttime,*finalcztseccnt;
	unsigned short *evtpha,*finalpha,*evtcztntick,*evtveto,*finalcztntick,*finalveto;
	int *evtpi,*finalpi;
	unsigned char *evtdetid,*evtpixid,*evtdetx,*evtdety,*evtalpha,*finaldetid,*finalpixid,*finaldetx,*finaldety,*finalalpha;
	float *evtenergy,*finalenergy;
	double *live_time_UT;
	int live_counter=0,livetimesize,BTIcounter;
	unsigned char *buntime_dfs,*buntime_dsl;
	double *buntime_del1,*buntime_del2;

	unsigned int  counter2,counter3, counter4;
	//int  buntime_heavy_index;
	unsigned int dph[64][64], dph_hot[64][64];
//	float distance;
	float cij;
	float total_hotness =0.0 ,n_pair =0.0,n_hot_pix = 0.0,allowable_hot,d;
	char outtxtfile[1000];
	//char outfile[1000],outlivetimefile[1000],outbtifile[1000];;
	double *lttime,*fracexp;
	int BTI_counter=0;
	int BTItimesize;
	char extname[20]; 
	double *BTITSTART, *BTITSTOP;

		evttime=(double*)malloc(sizeof(double)*size);
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

	
		buntime=(double*)malloc(sizeof(double)*size);
		buntime_heavy=(double*)malloc(sizeof(double)*size);
		buntime_real =(double*)malloc(sizeof(double)*size);
		buntime_dfs=(unsigned char*)malloc(sizeof(unsigned char)*size);
		buntime_dsl=(unsigned char*)malloc(sizeof(unsigned char)*size);
		bunsize_real =(int*)malloc(sizeof(int)*size);
		bunsize=(int*)malloc(sizeof(int)*size);
		bun_flag=(int*)malloc(sizeof(int)*size);
		buntime_del1=(double*)malloc(sizeof(double)*size);
		buntime_del2=(double*)malloc(sizeof(double)*size); 
		
		lttime=(double*)malloc(sizeof(double)*size);
		fracexp=(double*)malloc(sizeof(double)*size);

	FILE *thrfile,*f1;
	thrfile = fopen(thresholdfile, "r");

        if (thrfile == NULL) {
                printf("Couldn't open the threshold file.");
                return;
        }

    const size_t label_size = 300;
	char* label = malloc(label_size);

	while (fgets(label, label_size, thrfile) != NULL)  
	{
		char *val;
		strtok_r (label, " ", &val);
		if (strcmp(label,"heavy_bun_thresh_count_10ms") == 0)
			heavy_bun_thresh_count_10ms=atoi(val);
		else if (strcmp(label,"heavy_bun_thresh_count_1ms") == 0)
			heavy_bun_thresh_count_1ms=atoi(val);
		else if (strcmp(label,"bun_size_thresh") == 0)
			bun_size_thresh=atoi(val);
		else if (strcmp(label,"dph_time") == 0)
			dph_time=atof(val);
		else if (strcmp(label,"thresh_hot") == 0)
			thresh_hot=atof(val);	
		else if (strcmp(label,"allowable_hot_thresh") == 0)
			allowable_hot_thresh=atof(val);	 	
	}
	free(label);
	char* tempevt = calloc(strlen(infile)+1, sizeof(char));
	strcpy(tempevt, infile);

	printf("Input event file is %s\n",infile);  
	fitsfile *fevt, *fbun,*flivetime;	
	fits_open_file(&fevt,infile,READONLY,&status);
	if(status) { printerror( status );}

	fits_open_file(&fbun,inbunchfile,READONLY,&status);
	if(status) { printerror( status );}

	fits_open_file(&flivetime,inlivetimefile,READONLY,&status);
	if(status) { printerror( status );}

	char *file = strtok(infile, "."); 
	
	//file = strtok(file, "/");
	//sprintf(outfile, "%s_sbc.evt",file);
	sprintf(outtxtfile, "%s_sbc.txt",file);
	printf("The output log file is  %s\n",outtxtfile);
	//sprintf(outlivetimefile, "%s_sbc_livetime.fits",file);
	//sprintf(outbtifile, "%s_sbc_bti.fits",file);
	f1=fopen(outtxtfile,"a");      

	if(clobber==1)
	{
		remove(outfile);
		remove(outlivetimefile);
		remove(outbtifile);
	}
	createBTIFile();
	createEventFile(tempevt);
	createLivetimeFile();
 
	for(qid=0;qid<4;qid++)
	{   
		BTI_counter=0;
		printf("%d started----------------------------------------------------\n",qid);
		live_counter=0;
		fits_movabs_hdu(fevt, qid+2, &hdutype, &status);
		fits_movabs_hdu(fbun, qid+2, &hdutype, &status);
		fits_movabs_hdu(flivetime, qid+2, &hdutype, &status);

		fits_get_num_rows(fevt, &evtnrows, &status);
		fits_get_num_rows(fbun, &bunnrows, &status);
		fits_get_num_rows(flivetime, &livetimenrows, &status);
	
		if(evtnrows==0)
	    {
		printf("No events in the event file\n");
		continue;
		}
	    
		if(bunnrows==0)
		{
		   printf("No events in the Bunch file\n");
		   continue;
		}
	    
	    /*(fits_get_colnum(fevt,CASEINSEN,"TIME",&evttime_col,&status);
	    fits_get_colnum(fevt,CASEINSEN,"DETID",&evtdetid_col,&status);		
	    fits_get_colnum(fevt,CASEINSEN,"PIXID",&evtpixid_col,&status);
	    fits_get_colnum(fevt,CASEINSEN,"DETX",&evtdetx_col,&status);
	    fits_get_colnum(fevt,CASEINSEN,"DETY",&evtdety_col,&status);*/
	    
	    fits_get_colnum(fbun,CASEINSEN,"TIME",&buntime_col,&status);
	    fits_get_colnum(fbun,CASEINSEN,"NUMEVENT",&bunsize_col,&status);
		fits_get_colnum(fbun,CASEINSEN,"Time_dfs",&buntime_df1_col,&status);
		fits_get_colnum(fbun,CASEINSEN,"Time_dsl",&buntime_df2_col,&status);
	 
	    
	    double tstart,tstop;
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

	
		buntime=(double*)realloc(buntime,sizeof(double)*bunnrows);
		buntime_heavy=(double*)realloc(buntime_heavy,sizeof(double)*bunnrows);
		
		buntime_dfs=(unsigned char*)realloc(buntime_dfs,sizeof(unsigned char)*bunnrows);
		buntime_dsl=(unsigned char*)realloc(buntime_dsl,sizeof(unsigned char)*bunnrows);
		
		buntime_del1=(double*)realloc(buntime_del1,sizeof(double)*bunnrows);
		buntime_del2=(double*)realloc(buntime_del2,sizeof(double)*bunnrows);
		
	
		buntime_real = (double*)realloc(buntime_real,sizeof(double)*bunnrows);
	//	bunsize_real = bunsize=(int*)realloc(bunsize_real,sizeof(int)*bunnrows); valgrind is giving segfault
		bunsize_real = (int*)realloc(bunsize_real,sizeof(int)*bunnrows);
		bunsize=(int*)realloc(bunsize,sizeof(int)*bunnrows);

		bun_flag=(int*)realloc(bun_flag,sizeof(int)*bunnrows);

		lttime=(double*)realloc(lttime,sizeof(double)*livetimenrows);
		fracexp=(double*)realloc(fracexp,sizeof(double)*livetimenrows+1);
	
		frow      = 1;
		felem     = 1;
		//doublenull = 0.;
		//intnull = 0;
		//bytenull = 0;
		//floatnull=0.0;

		fits_read_col(fevt, TDOUBLE, 1, frow, felem, evtnrows, NULL, evttime,NULL, &status);
		fits_read_col(fevt, TDOUBLE, 2, frow, felem, evtnrows, NULL, evtcztseccnt, NULL, &status);         
		fits_read_col(fevt, TUSHORT, 3, frow, felem, evtnrows, NULL, evtcztntick, NULL, &status); 
		fits_read_col(fevt, TUSHORT, 4, frow, felem, evtnrows, NULL, evtpha,NULL, &status);
		fits_read_col(fevt, TBYTE, 5, frow, felem,evtnrows, NULL, evtdetid, NULL, &status);  
		fits_read_col(fevt, TBYTE, 6, frow, felem,evtnrows, NULL, evtpixid, NULL, &status);        
		fits_read_col(fevt, TBYTE, 7, frow, felem, evtnrows, NULL, evtdetx, NULL, &status);   
		fits_read_col(fevt, TBYTE, 8, frow, felem,evtnrows, NULL, evtdety,NULL, &status);
		fits_read_col(fevt, TUSHORT, 9, frow, felem,evtnrows, NULL, evtveto,NULL, &status);
		fits_read_col(fevt, TBYTE, 10, frow, felem, evtnrows, NULL, evtalpha,NULL, &status); 
		fits_read_col(fevt, TINT, 11, frow, felem, evtnrows, NULL, evtpi,NULL, &status);
		fits_read_col(fevt, TFLOAT, 12, frow, felem, evtnrows, NULL, evtenergy,NULL, &status);
		fits_read_key(fevt,TDOUBLE,"TSTART",&tstart,NULL, &status);
		fits_read_key(fevt,TDOUBLE,"TSTOP",&tstop,NULL, &status);
	
		fits_read_col(fbun, TDOUBLE, buntime_col, 1, 1, bunnrows, NULL, buntime,NULL, &status);
		fits_read_col(fbun, TINT, bunsize_col, 1, 1, bunnrows, NULL, bunsize,NULL, &status);
		fits_read_col(fbun, TBYTE, buntime_df1_col, 1, 1, bunnrows, NULL, buntime_dfs,NULL, &status);
		fits_read_col(fbun, TBYTE, buntime_df2_col, 1, 1, bunnrows, NULL, buntime_dsl,NULL, &status);

		fits_read_col(flivetime, TDOUBLE, 1, 1, 1, livetimenrows, NULL, lttime,NULL, &status);
		fits_read_col(flivetime, TDOUBLE, 2, 1, 1, livetimenrows, NULL, fracexp,NULL, &status);

		livetimesize=((int)(evttime[evtnrows-1])-(int)(evttime[0])+1)*10;
		live_time_UT=(double*)malloc(sizeof(double)*livetimesize);

		BTItimesize=((int)(evttime[evtnrows-1])-(int)(evttime[0])+1)*10;
		
		BTITSTART=(double*)malloc(sizeof(double)*BTItimesize);
		BTITSTOP=(double*)malloc(sizeof(double)*BTItimesize);

		for(i=0;i<evtnrows;i++){evt_flag[i]=0;}
		int heavy_bunch_length=0;
		long real_bunch_length=0;
		//double kk=0.0;
		
		
		
		
		
		// for(i=0;i<bunnrows;i++)
		// {
		// 	buntime_del1[i] = buntime[i] - 20.0*(double)buntime_dfs[i]/1000000.0;
		// 	buntime_del2[i] = buntime[i] + 20.0*(double)buntime_dsl[i]/1000000.0;
		// }
		//printf("Bunch bunnrows %ld\n",bunnrows);
		//printf("%lf %lf \n",tstart,tstop);
		
		i=0;
		while(i < bunnrows)
		{
			//printf("entered while\n");
			if(buntime[i]<(tstart-1.0)) {i++;continue;}
			else if(buntime[i]>tstop){break;}
			else{(int*)realloc(bunsize,sizeof(int)*bunnrows);}
			
			buntime_real[real_bunch_length] = buntime[i];
			bunsize_real[real_bunch_length] = bunsize[i];
			
			if(bunsize[i]!=63)
			{
				
				if(bunsize_real[real_bunch_length] > bun_size_thresh)
				{ 
					buntime_heavy[heavy_bunch_length++] = buntime_real[real_bunch_length];
				}
					
				
				real_bunch_length++;
				i++;
				
				continue;
			}
			else
			{
				//~ if( (buntime[i]> 209474728.004126+71283.7) && (buntime[i]< 209474728.004126+71283.7+2.0)) 
				//~ {
					//~ printf("%d %lf %lf \n",bunsize[i],(buntime[i] + 20.0*(double)buntime_dsl[i]/1000000.0), buntime[i+1] );
					
				//~ }
				
				for(j=i+1;j<bunnrows;j++)
				{
					//printf("%d %d %d %d\n",i,j,bunnrows,real_bunch_length);
					//printf("%lf %lf\n",20.0*(double)buntime_dfs[j]/1000000.0,20.0*(double)buntime_dsl[i]/1000000.0);
					//printf("%lf \n",((buntime[j] - 20.0*(double)buntime_dfs[j]/1000000.0)-(buntime[i] + 20.0*(double)buntime_dsl[i]/1000000.0)));
					if( ((buntime[j] - 20.0*(double)buntime_dfs[j]/1000000.0)-(buntime[i] + 20.0*(double)buntime_dsl[i]/1000000.0)) <= 30.0/1000000.0)
					{
						//printf("T\n");
						i=j;
						buntime_real[real_bunch_length] = (buntime[j] + 20.0*(double)buntime_dsl[j]/1000000.0);
						bunsize_real[real_bunch_length] = bunsize_real[real_bunch_length]+bunsize[j];
			
					}
					else
					{
						i++;
						if(j==bunnrows)break;
						if(bunsize_real[real_bunch_length] > bun_size_thresh)
						{ 
							buntime_heavy[heavy_bunch_length++] = buntime_real[real_bunch_length];
							//kk++;
						}
						
						
						real_bunch_length++;
						//printf("%d %d %d %d\n",i,j,bunnrows,real_bunch_length);
						break;
						
					}
				}
				//~ if( (buntime[i]> 209474728.004126+71283.7) && (buntime[i]< 209474728.004126+71283.7+2.0)) 
				//~ {
					//~ printf("%d %lf \n",bunsize_real[real_bunch_length-1], buntime_heavy[real_bunch_length-1]);
					
				//~ }
				
				
				
			} 
			if(i==(bunnrows-1))break;
			//printf("%lf %d\n",buntime_real[real_bunch_length-1],bunsize_real[real_bunch_length-1]);
			//printf("%ld %d\n",i,bunnrows);
			//printf("%d\n",kk);
			
		}
		printf("Initial number of bunches = %ld\n",bunnrows);
		printf("Bunch length correction Completed-----total bunches = %ld\n",real_bunch_length);
		printf("Number of super bunches = %ld\n",heavy_bunch_length);
		
	


		tmpindex = 0;
		tmpindex1 =0;
		//printf("%d\n",heavy_bunch_length);
		for(i=0;i<heavy_bunch_length;i++)
		{
			//printf("%d\t%d\n",qid,i);
			//if(bun_flag[i]==1)
			//{
				total_hotness =0.0;
				n_pair =0.0;
				n_hot_pix = 0.0;
	
	
				//initializing DPH array
				for(k=0;k<64;k++)
				{   
					for(m=0;m<64;m++)
					{	//printf("%d\t%d\t%d\t%d\n",qid,i,k,m );
						dph[k][m]=0;
						dph_hot[k][m]=0;
					}
		
				}
				
	
				counter2 = 0;
				for(j=tmpindex;j<evtnrows;j++)
				{
	
					if((int)evttime[j] > (int)buntime_heavy[i]+1){break;}
					//printf("%d\t%d\t%d\n",qid,i,j);
					//printf("%d\t%d\n",i,j);
					if((evttime[j] >= buntime_heavy[i]) && (evttime[j] <= buntime_heavy[i]+dph_time))
					{	
						if(counter2==0){tmpindex = j;}
						//if(counter4==0){evttime_index = j;}
						//time_evt_close_bunch[counter4] = evttime[i];
						detx_close_bunch[counter2] = evtdetx[j];
						dety_close_bunch[counter2] = evtdety[j];
						dph[evtdetx[j]][evtdety[j]]++;
						counter2++;
				    }
				}
	
				//printf("HBL loop %d %d\n",i,heavy_bunch_length);
				if(counter2==0){continue;}
	
				//Comment : Ajay R on 02/05/2020
				//for(j=0;j<counter2;j++)
				//{
					//printf("%d\n",j);
					//printf("%d\t%d\n",detx_close_bunch[j],dety_close_bunch[j]);
					//dph[detx_close_bunch[j]][dety_close_bunch[j]]++;
		
				//}
	
				//Getting magnitude and coordinates of non-zero dph values
				counter3=0;
				for(k=0;k<64;k++)
				{
					for(m=0;m<64;m++)
					{
						//printf("%d\t%d\n",k,m);
						if(dph[k][m]>0)
						{	
				
							magnitude_temp[counter3] = dph[k][m];
							detx_temp[counter3] = k;
							dety_temp[counter3] = m;
							counter3++;
						}
			
					}
		
				}
							//printf("%d\n",counter3);
				

				for(k=0;k<counter3;k++)
				{
					for(m=k+1;m<counter3;m++)
					{
						//if(i==68998)	printf("%d %d %d %d %d %d %d %d \n",k,m,detx_temp[m],detx_temp[k],dety_temp[m],dety_temp[k],dph[detx_temp[k]][dety_temp[k]],dph[detx_temp[m]][dety_temp[m]]);
			
						//if(qid==2 && i==7799){ printf("%d\t%d\t%d\n",m,k,counter3);}
					//	if(i==68998)		//printf("DPH values \n");
					
						d = sqrt(pow((detx_temp[m]-detx_temp[k]),2)+pow((dety_temp[m]-dety_temp[k]),2));
						if((dph[detx_temp[m]][dety_temp[m]]==1) && (dph[detx_temp[k]][dety_temp[k]]==1)  && (d>1.5))//to avoid isolated genuine double events
						{
							counter4=0;
							for(n=detx_temp[m]-1;n<=detx_temp[m]+1;n++)
							{ 
								if(n<0||n>63){continue;}
								for(p=dety_temp[m]-1;p<=dety_temp[m]+1;p++)
								{
									if(p<0||p>63){continue;}
									if(dph[n][p]>0){counter4++;}
						
								}
							}
							for(n=detx_temp[k]-1;n<=detx_temp[k]+1;n++)
							{ 
								if(n<0||n>63){continue;}
								for(p=dety_temp[k]-1;p<=dety_temp[k]+1;p++)
								{
									if(p<0||p>63){continue;}
									if(dph[n][p]>0){counter4++;}
						
								}
							}
						}
									
						//if(qid==2 && i==7799){ printf("%d\t%d\t%d\n",m,k,counter3);}
			
						else{counter4 =10;} //randomly assigned 10
					//	if(i==68998)	printf("2\n");			
				
						if(counter4>2)
						{
							d = sqrt(pow((detx_temp[m]-detx_temp[k]),2)+pow((dety_temp[m]-dety_temp[k]),2));
							cij = magnitude_temp[k]*magnitude_temp[m]/d; 
							
							//if(i==68998)
							//{	
							//	printf("hotness values %d %d %lf %lf %d %d %d %d\n",k,m,d,cij,detx_temp[m],detx_temp[k],dety_temp[m],dety_temp[k]);		
							//	printf("hotness values %lf %lf %d %d %d %d\n",d,cij,detx_temp[m],detx_temp[k],dety_temp[m],dety_temp[k]);			
							//}
							
							if(cij > thresh_hot)
							{
								total_hotness = total_hotness+cij;
								n_pair++;
								dph_hot[detx_temp[k]][dety_temp[k]]=1;
								dph_hot[detx_temp[m]][dety_temp[m]]=1;

							}
						}
					//printf("%d\t%d\n",m,k);

				 	}
				}	
				//printf("HBL loop 2 %d %d\n",i,heavy_bunch_length);	
	
	
				//printf("%d\n",counter3);
	
				for(k=0;k<64;k++)
				{
					for(m=0;m<64;m++)
					{//printf("%d\t%d\n",k,m);
						n_hot_pix=n_hot_pix+dph_hot[k][m];
					}
				}
				
			
				
				if (n_hot_pix!=0.0){
				allowable_hot = total_hotness/n_hot_pix;}
				else allowable_hot =0.0;
				//need to change(orbit no)
				//fprintf(f1,"%d\t%f\t%f\t%f\t%f\t%s\n",bunch_size,allowable_hot,total_hotness,n_hot_pix,n_pair,orbitno);
				
				fprintf(f1,"%f\t%f\t%f\t%f\t%d\n",allowable_hot,total_hotness,n_hot_pix,n_pair,qid);
				
				
				//~ if( (buntime_heavy[i]> 209474728.004126+71283.7) && (buntime_heavy[i]< 209474728.004126+71283.7+2.0)) 
				//~ {
					//~ printf("%lf\t%f\t%f\t%f\t%f\t%d\n",buntime_heavy[i],allowable_hot,total_hotness,n_hot_pix,n_pair,qid);
					
				//~ }
				
				
				if (allowable_hot>allowable_hot_thresh)
				{	
					//printf("ENTERED FLAGGED\n");
					live_time_UT[live_counter++] = buntime_heavy[i];
					BTITSTART[BTI_counter] = buntime_heavy[i];
					BTITSTOP[BTI_counter]  = buntime_heavy[i]+dph_time;
					BTI_counter++;

					//printf("%lf\t%lf\t%d\n",buntime_heavy[i],live_time_UT[live_counter-1],live_counter);
					for(e=tmpindex;e<evtnrows;e++)
					{
						if((int)evttime[e] > (int)buntime_heavy[i]+1){break;}
						
						//if((int)evttime[e]-5 < (int)buntime_heavy[i]){tmpindex1=e;}
						//tmpindex1=e;
						
						//if( (evttime[e]> 209474728.004126+71283.7) && (evttime[e]< 209474728.004126+71283.7+0.12))
						//printf("%d\t%lf\n",e,evttime[e],buntime_heavy[i]);
						
						
						if((evttime[e] >= buntime_heavy[i]) && (evttime[e] <= buntime_heavy[i]+dph_time))
						{	

							evt_flag[e] = 1;
						}
					}
				}
				
				//~ if( (buntime_heavy[i]> 209474728.004126+71283.7) && (buntime_heavy[i]< 209474728.004126+71283.7+2.0)) 
				//~ {
					//~ printf("%lf\t%f\t%f\t%f\t%f\t%d\n",buntime_heavy[i],allowable_hot,total_hotness,n_hot_pix,n_pair,qid);
					
				//~ }
				
				

			//}
		}
		printf("Main Process Over\n");
		
		//livetimeGeneration(live_time_UT,live_counter,outlivetimefile,qid+2,lttime,fracexp,livetimenrows);
		writeBTItime(BTITSTART,BTITSTOP,BTI_counter,qid+2);
		int l=0;
	
		long qgtinrows;//new_qgtinrows=0;
		
		//double tstart,tstop;
		double exp_loss_due_bunches=0.0;
		//fits_read_key(fevt,TDOUBLE,"TSTART",&tstart,NULL, &status);
		//fits_read_key(fevt,TDOUBLE,"TSTOP",&tstop,NULL, &status);
		
		
		//~ fits_movabs_hdu(fevt, qid+10, &hdutype, &status);
		sprintf(extname, "Q%d_GTI",qid);
		fits_movnam_hdu(fevt, BINARY_TBL, extname, 0, &status);
		
		
		
		fits_get_num_rows(fevt, &qgtinrows, &status);
		double *BTITSTART_merge,*BTITSTOP_merge;
		double *GTITSTART_merge,*GTITSTOP_merge;
		double totime=(double)livetimenrows*(live_time_UT[1]-live_time_UT[0]);
		double *gtitstart,*gtitstop,*new_gtitstart,*new_gtitstop,*gtitstart_sel,*gtitstop_sel;
		
		
		
		gtitstart=(double*)malloc(sizeof(double)*qgtinrows);
		gtitstop=(double*)malloc(sizeof(double)*qgtinrows);
		gtitstart_sel=(double*)malloc(sizeof(double)*qgtinrows);
		gtitstop_sel =(double*)malloc(sizeof(double)*qgtinrows);
		
		
		new_gtitstart=(double*)malloc(sizeof(double)*(qgtinrows+(int)totime*10));
		new_gtitstop=(double*)malloc(sizeof(double)*(qgtinrows+(int)totime*10));
		
		frow      = 1;
		felem     = 1;
		//doublenull = 0.;
		//intnull = 0;
		//bytenull = 0;
		//floatnull=0.0;

		fits_read_col(fevt, TDOUBLE, 1, frow, felem, qgtinrows, NULL, gtitstart,NULL, &status);
		fits_read_col(fevt, TDOUBLE, 2, frow, felem, qgtinrows, NULL, gtitstop, NULL, &status);  
		
		
		
		BTITSTART_merge=(double*)malloc(sizeof(double)*BTItimesize);
		BTITSTOP_merge=(double*)malloc(sizeof(double)*BTItimesize);
		
		GTITSTART_merge=(double*)malloc(sizeof(double)*(BTItimesize+1));
		GTITSTOP_merge=(double*)malloc(sizeof(double)*(BTItimesize+1));
		
		
		
		//To select GTI for the particular file
		long gti_counter=0;
		for(i=0;i<qgtinrows;i++)
		{
			if( gtitstart[i] >= gtitstop[i])
			{
				printf("GTI ERROR PRESENT in Quadrant %d\n",qid);
			}
			
			
			if(((gtitstart[i]>=tstart)&&(gtitstart[i]<tstop))&&((gtitstop[i]>tstart)&&(gtitstop[i]<=tstop)))
			{
			
				gtitstart_sel[gti_counter]= gtitstart[i];
				gtitstop_sel[gti_counter]= gtitstop[i];
				gti_counter++;
			}
			else if(((gtitstart[i]>=tstart)&&(gtitstart[i]<tstop)) && !((gtitstop[i]>tstart)&&(gtitstop[i]<=tstop)))
			{
			
				gtitstart_sel[gti_counter]= gtitstart[i];
				gtitstop_sel[gti_counter]= tstop;
				gti_counter++;
			}
			else if(!((gtitstart[i]>=tstart)&&(gtitstart[i]<tstop)) && ((gtitstop[i]>tstart)&&(gtitstop[i]<=tstop)))
			{
				
				gtitstart_sel[gti_counter]= tstart;
				gtitstop_sel[gti_counter]= gtitstop[i];
				gti_counter++;
			}
			else if(((gtitstart[i]<=tstart) && (gtitstop[i]>=tstop)))
			{
				gtitstart_sel[gti_counter]= tstart;
				gtitstop_sel[gti_counter]= tstop;
				gti_counter++;
				
			}
			
		}
		
		
		double exposure=0.0;
		for(i=0;i<gti_counter;i++)
		{
			exposure = exposure+(gtitstop_sel[i]-gtitstart_sel[i]);
		
		}
		printf("Initial Exposure = %lf\n",exposure);
		
		//To merge the overlaps in Bad Time Intervals
		int ii=0;
		double BTI_temp=0;
		//int BTI_next=0;
		
		//printf("BTI_counter---------%d\n",BTI_counter);
		//fits_read_key(fevt,TDOUBLE,"TSTART",&GTITSTART_merge[0],NULL, &status);

		int temp_i;
	
		
		GTITSTART_merge[ii]=tstart;
		for(i=0;i<BTI_counter;i++)
		{
			temp_i=i;
			BTITSTART_merge[ii] = BTITSTART[i];
			BTITSTOP_merge[ii] = BTITSTOP[i];
			BTI_temp = BTITSTOP[i];
			for(j=i+1;j<BTI_counter;j++)
			{
			
				if(BTITSTART[j] <= BTI_temp)
				{
					//printf("TRUE____________\n");
					BTITSTOP_merge[ii] = BTITSTOP[j];
					BTI_temp = BTITSTOP[j];
					temp_i = j+1;
				}
				else
					break;
				
			}
			i=temp_i;
			GTITSTOP_merge[ii]=BTITSTART_merge[ii];
			GTITSTART_merge[ii+1]=BTITSTOP_merge[ii];
			
			ii++;
		
		}
		
		
		//fits_read_key(fevt,TDOUBLE,"TSTOP",&GTITSTOP_merge[ii],NULL, &status);
		GTITSTOP_merge[ii]=tstop;
		
		if(BTI_counter==0){ii++;}
	
		//printf("%d\t%d\n",BTI_counter,ii);
		//To test overlap
		
		//double large_value = GTITSTOP_merge[0],small_value = GTITSTART_merge[0];
		//printf("%lf\t%lf\t%lf\n",small_value,large_value,large_value-small_value);
		//for(i=1;i<ii;i++)
		//{
		//	printf("%lf\n",GTITSTOP_merge[i]-GTITSTART_merge[i]);
		//	if(GTITSTART_merge[i]<large_value) {printf("OVERLAP PRESENT---%lf\t%lf\n",GTITSTART_merge[i],large_value);}
		//	
		//	if(large_value < GTITSTOP_merge[i]){large_value=GTITSTOP_merge[i];}
		//	
		//
		//
		//
		//}
		//printf("%lf\t%lf\t%lf\n",small_value,large_value,large_value-small_value);
	
	
	
	
		//To get the corresponding livetime
		//double fracexp_new[livetimenrows], lttime_new[livetimenrows];
		//double BTI_exclude_TSTART[livetimenrows],BTI_exclude_TSTOP[livetimenrows];
		
		double *fracexp_new,*lttime_new,*BTI_exclude_TSTART,*BTI_exclude_TSTOP;
		fracexp_new=(double*)calloc(livetimenrows,sizeof(double));
		lttime_new=(double*)calloc(livetimenrows,sizeof(double));
		BTI_exclude_TSTART=(double*)calloc(livetimenrows,sizeof(double));
		BTI_exclude_TSTOP=(double*)calloc(livetimenrows,sizeof(double));
		
		
		long livetimenrows_new=0;
		double livetime_resolution = (lttime[1]-lttime[0])/2.0;	
		long BTI_saturation_no=0;
		
		
		
	
		for(i=0;i<livetimenrows;i++)
		{
			if( lttime[i]>=(double)(long)tstart && lttime[i]<(double)(long)tstop )
			{
				lttime_new[livetimenrows_new] = lttime[i];
				fracexp_new[livetimenrows_new] = fracexp[i];
				livetimenrows_new++;
				//printf("%lf\n",(livetime_resolution*2-fracexp[i]));
				exp_loss_due_bunches = exp_loss_due_bunches + (livetime_resolution*2-fracexp[i]);
				
				if(lttime[i]-tstart < livetime_resolution/2.0)
				{
					fracexp_new[livetimenrows_new] = fracexp[i]-(tstart-(lttime[i]-livetime_resolution));
				}
				
				if(tstop-lttime[i] < livetime_resolution/2.0)
				{
					fracexp_new[livetimenrows_new] = fracexp[i]-((lttime[i]+livetime_resolution)-tstop);
				}
				
				/*if(fracexp[i]<0.2)
				{
					BTI_exclude_TSTART[BTI_saturation_no] = lttime[i]-livetime_resolution;
					BTI_exclude_TSTOP[BTI_saturation_no]  = lttime[i]-livetime_resolution;
					BTI_saturation_no++;
					fracexp[i]=1.0;
				
				}*/
			}
			}
			//printf("Exposure loss due to bunches initial = %lf\n",exp_loss_due_bunches);
		
		//writeLivetime(lttime_new,fracexp_new,livetimenrows_new,outlivetimefile,qid+2);
		
		/*double GTI_saturation_TSTART[livetimenrows+1],GTI_saturation_TSTOP[livetimenrows+1];
		
		GTI_saturation_TSTART[0] = tstart;
		for(i=0;i<BTI_saturation_no;i++)
		{
			GTI_saturation_TSTOP[i] = BTI_exclude_TSTART[i];
			GTI_saturation_TSTART[i+1] = BTI_exclude_TSTOP[i];
			
			
		}
		GTI_saturation_TSTOP[BTI_saturation_no] = tstop;*/

		
		
		
		
		//To merge two GTIs
		double *gtitstart_f,*gtitstop_f;
		long gtinrows_f=0,gtinrows_est=gti_counter+ii;
		
		gtitstart_f=(double*)calloc(gtinrows_est,sizeof(double));
		gtitstop_f=(double*)calloc(gtinrows_est,sizeof(double));
		
		
		//printf("GTI____%ld %d\n",gti_counter,ii);
		
		//printf("GTI duration = %d\n",gti_counter);
		for(i=0;i<gti_counter;i++)
		{
			if(gtitstart_sel[i] > gtitstop_sel[i])
			{
				printf("GTI ERROR PRESENT in Quadrant %d\n",qid);
			}
			
			if (ii==0){
					gtitstart_f[gtinrows_f] =  gtitstart_sel[i];
					gtitstop_f[gtinrows_f] =  GTITSTOP_merge[j];
					gtinrows_f++;}
			
			
			for(j=0;j<ii;j++)
			{
				
				if(GTITSTART_merge[j] > GTITSTOP_merge[j])
				{
					printf("GTI ERROR PRESENT in Quadrant %d\n",qid);
				}
				
				
				
				
				
				
				if(((GTITSTART_merge[j]>=gtitstart_sel[i])&&(GTITSTART_merge[j]<gtitstop_sel[i]))&&((GTITSTOP_merge[j]>=gtitstart_sel[i])&&(GTITSTOP_merge[j]<gtitstop_sel[i])))
				{
					//printf("TRUE_m\n");
					gtitstart_f[gtinrows_f] =  GTITSTART_merge[j];
					gtitstop_f[gtinrows_f] =  GTITSTOP_merge[j];
					gtinrows_f++;
				}
				
				else if((!(GTITSTART_merge[j]>=gtitstart_sel[i])&&(GTITSTART_merge[j]<gtitstop_sel[i]))&&((GTITSTOP_merge[j]>=gtitstart_sel[i])&&(GTITSTOP_merge[j]<gtitstop_sel[i])))
				{
					//printf("TRUE_f\n");
					gtitstart_f[gtinrows_f] =  gtitstart_sel[i];
					gtitstop_f[gtinrows_f] =  GTITSTOP_merge[j];
					gtinrows_f++;
				}
				
				else if(((GTITSTART_merge[j]>=gtitstart_sel[i])&&(GTITSTART_merge[j]<gtitstop_sel[i]))&& !((GTITSTOP_merge[j]>=gtitstart_sel[i])&&(GTITSTOP_merge[j]<gtitstop_sel[i])))
				{
					//printf("TRUE_l\n");
					gtitstart_f[gtinrows_f] =  GTITSTART_merge[j];
					gtitstop_f[gtinrows_f] =  gtitstop_sel[i];
					gtinrows_f++;
				}
				else if((GTITSTART_merge[j]<=gtitstart_sel[i])&& (GTITSTOP_merge[j]>=gtitstop_sel[i]))
				{
					//printf("TRUE_l\n");
					gtitstart_f[gtinrows_f] =  gtitstart_sel[i];
					gtitstop_f[gtinrows_f] = gtitstop_sel[i];
					gtinrows_f++;
				}
				
				
			}
			
		}
		printf("GTI duration = %ld\n",gtinrows_f);
		
		//To test overlap
		
		//double large_value = gtitstop_f[0],small_value = gtitstart_f[0];
		////printf("%lf\t%lf\t%lf\n",small_value,large_value,large_value-small_value);
		//for(i=1;i<gtinrows_f;i++)
		//{
		//	printf("%lf\n",gtitstop_f[i]-gtitstart_f[i]);
		//	if(gtitstart_f[i]<large_value) {printf("OVERLAP PRESENT");}
		//	
		//	if(large_value < gtitstop_f[i]){large_value=gtitstop_f[i];}
		//	
		//
		//
		//
		//}
		//printf("%lf\t%lf\t%lf\n",small_value,large_value,large_value-small_value);
		
		
		
		
		//To generate livetime
		double livfrac = 0.0;
		
		//printf("-----------------------------------------%ld\t%ld\n",livetimenrows_new,gtinrows_f);
		
		for(i=0;i<livetimenrows_new;i++)
		{
		
			livfrac = 0.0;
			for(j=0;j<gtinrows_f;j++)
			{
				if(((gtitstart_f[j] >= lttime_new[i]-livetime_resolution) && (gtitstart_f[j] < lttime_new[i]+livetime_resolution) )  && ((gtitstop_f[j]>= lttime_new[i]-livetime_resolution) && (gtitstop_f[j] < lttime_new[i]+livetime_resolution) ))
				{
					livfrac = livfrac+ gtitstop_f[j]-gtitstart_f[j];
			
				
						
				}
				else if(((gtitstart_f[j]>= lttime_new[i]-livetime_resolution) && (gtitstart_f[j] < lttime_new[i]+livetime_resolution) )  && !((gtitstop_f[j]>= lttime_new[i]-livetime_resolution) && (gtitstop_f[j] < lttime_new[i]+livetime_resolution) ))
				{
					livfrac = livfrac+ lttime_new[i]+livetime_resolution-gtitstart_f[j];
					//printf("%lf\t%lf\n",(lttime_new[i]+livetime_resolution),gtitstart_f[j]);
					
				}
				else if(!((gtitstart_f[j]>= lttime_new[i]-livetime_resolution) && (gtitstart_f[j] < lttime_new[i]+livetime_resolution) )  && ((gtitstop_f[j]>= lttime_new[i]-livetime_resolution) && (gtitstop_f[j] < lttime_new[i]+livetime_resolution) ))
				{
					livfrac = livfrac+  gtitstop_f[j]-(lttime_new[i]-livetime_resolution);
					//printf("%lf\t%lf\n",(lttime_new[i]-livetime_resolution),gtitstop_f[j]);
					
				}
				else if((gtitstart_f[j] <= lttime_new[i]-livetime_resolution) && (gtitstop_f[j] >= lttime_new[i]+livetime_resolution) )  
				{
					livfrac = livfrac+ 2.0*livetime_resolution;
					
					
				}
			}
			
			//printf("%lf\n",livfrac);
			
			//double test=livfrac/(2.0*livetime_resolution)-(1.0-fracexp_new[i])*(2.0*livetime_resolution);
			//printf("%lf\t%lf\t%lf\t%lf\n",livfrac,fracexp_new[i],1-fracexp_new[i]/(2.0*livetime_resolution),test);
			
			if (livfrac==0.0){fracexp_new[i]=livetime_resolution*2;}
			//double test_value = fracexp_new[i];
			else
			{//fracexp_new[i] = livfrac/(2.0*livetime_resolution)-(1-fracexp_new[i])*(2.0*livetime_resolution);
			double x = fracexp_new[i];
			//exp_loss_due_bunches = exp_loss_due_bunches + (2.0*livetime_resolution-fracexp_new[i]);
			fracexp_new[i] = livfrac-(2.0*livetime_resolution-fracexp_new[i])*(livfrac/(2.0*livetime_resolution));
			
			//printf("%lf\t%lf\n",fracexp_new[i],livfrac);
			}
			
			
			
			if(fracexp_new[i]<0.1*livetime_resolution*2)
			{
					//printf("%lf\t%lf\t%lf\n",fracexp_new[i],livfrac,test_value);
					BTI_exclude_TSTART[BTI_saturation_no] = lttime_new[i]-livetime_resolution;
					BTI_exclude_TSTOP[BTI_saturation_no]  = lttime_new[i]+livetime_resolution;
					BTI_saturation_no++;
					fracexp_new[i]=livetime_resolution*2.0;
				
			}
			
		}
		
		//printf("%ld\t%ld\n",livetimenrows,livetimenrows_new);
		/*for(i=0;i<livetimenrows_new;i++)
		{
			printf("%lf\t%lf\n",lttime[i],fracexp_new[i]);
		}*/
		//printf("BTI_saturation_no: %ld\n",BTI_saturation_no);
		
		
		
		//double GTI_saturation_TSTART[livetimenrows+1],GTI_saturation_TSTOP[livetimenrows+1];
		
		double *GTI_saturation_TSTART,*GTI_saturation_TSTOP;
		
		GTI_saturation_TSTART=(double*)calloc((BTI_saturation_no+1),sizeof(double));
		GTI_saturation_TSTOP=(double*)calloc((BTI_saturation_no+1),sizeof(double));
		
		
		//printf("RUN UNTILL HERE\n");
		GTI_saturation_TSTART[0] = tstart;
		int temp_ind=0,temp_ii=0;
		
		if(BTI_exclude_TSTART[0]<tstart && BTI_exclude_TSTOP[0]>tstart)
		{
				GTI_saturation_TSTART[0]=BTI_exclude_TSTOP[0];
				temp_ii++;
		}
		
		
		//printf(" ERROR BEGINING  %lf %lf %lf", tstart,BTI_exclude_TSTART[0],BTI_exclude_TSTOP[0]);
		for(i=temp_ii;i<BTI_saturation_no;i++)
		{
			
			if(BTI_exclude_TSTART[i] == BTI_exclude_TSTOP[i])
			{
				printf("GTI ERROR PRESENT in Quadrant %d\n",qid);
				//printf("%lf \t %lf \n",gtitstart_f[j],gtitstop_f[j]);
			}
			
		
			for(j=temp_ind;j<evtnrows;j++)
			{
					
				if(evttime[j] >= BTI_exclude_TSTART[i] && evttime[j] < BTI_exclude_TSTOP[i])
				{
					evt_flag[j] = 1;
					temp_ind = j;
				}
				
			}
			GTI_saturation_TSTOP[i-temp_ii] = BTI_exclude_TSTART[i];
			GTI_saturation_TSTART[i+1-temp_ii] = BTI_exclude_TSTOP[i];
			
			//printf("%lf %lf %lf %lf %d \n",BTI_exclude_TSTART[i],BTI_exclude_TSTOP[i],GTI_saturation_TSTART[i],GTI_saturation_TSTOP[i],i);
			
		}
		GTI_saturation_TSTOP[BTI_saturation_no-temp_ii] = tstop;
		
		if(BTI_exclude_TSTART[BTI_saturation_no-1]<tstop && BTI_exclude_TSTOP[BTI_saturation_no-1]>tstop)
		{
				temp_ii--;
		}
		
		
		
		long GTI_saturation_no= BTI_saturation_no-temp_ii+1;
		
		
		
		//To merge final GTIs
		double *gtitstart_final,*gtitstop_final;
		long gtinrows_final=0,gtinrows_est_final=gtinrows_f+BTI_saturation_no+1;
		
		gtitstart_final=(double*)calloc(gtinrows_est_final,sizeof(double));
		gtitstop_final=(double*)calloc(gtinrows_est_final,sizeof(double));

		//printf("GTI_saturation = %ld\n",GTI_saturation_no);
		
		for(i=0;i<GTI_saturation_no;i++)
		{
			
			if(GTI_saturation_TSTART[i] == GTI_saturation_TSTOP[i])
			{
				
				//printf("ERROR---1----GTI  Tstart >= Tstop---------\n");
				//printf("%lf \t %lf %d\n",GTI_saturation_TSTART[i],GTI_saturation_TSTOP[i],i);
				continue;
			}
			if(GTI_saturation_TSTART[i] > GTI_saturation_TSTOP[i])
			{
				
				printf("GTI ERROR PRESENT in Quadrant %d\n",qid);
				//printf("%lf \t %lf %d\n",GTI_saturation_TSTART[i],GTI_saturation_TSTOP[i],i);
				//continue;
			}

			
			
			
			for(j=0;j<gtinrows_f;j++)
			{
				if(gtitstart_f[j] > gtitstop_f[j])
				{
					printf("GTI ERROR PRESENT in Quadrant %d\n",qid);
					//printf("%lf \t %lf \n",gtitstart_f[j],gtitstop_f[j]);
				}
				
				if(((gtitstart_f[j]>=GTI_saturation_TSTART[i])&&(gtitstart_f[j]<GTI_saturation_TSTOP[i]))&&((gtitstop_f[j]>=GTI_saturation_TSTART[i])&&(gtitstop_f[j]<GTI_saturation_TSTOP[i])))
				{
					gtitstart_final[gtinrows_final] = gtitstart_f[j];
					gtitstop_final[gtinrows_final] =  gtitstop_f[j];
			
					
					gtinrows_final++;
					
				}
				
				else if(!((gtitstart_f[j]>=GTI_saturation_TSTART[i])&&(gtitstart_f[j]<GTI_saturation_TSTOP[i]))&&((gtitstop_f[j]>=GTI_saturation_TSTART[i])&&(gtitstop_f[j]<GTI_saturation_TSTOP[i])))
				{
			
					gtitstart_final[gtinrows_final] = GTI_saturation_TSTART[i];
					gtitstop_final[gtinrows_final] =  gtitstop_f[j];
					
				
					
					
					gtinrows_final++;
				}
				
				else if(((gtitstart_f[j]>=GTI_saturation_TSTART[i])&&(gtitstart_f[j]<GTI_saturation_TSTOP[i]))&&!((gtitstop_f[j]>=GTI_saturation_TSTART[i])&&(gtitstop_f[j]<GTI_saturation_TSTOP[i])))
				{
			
					gtitstart_final[gtinrows_final] = gtitstart_f[j];
					gtitstop_final[gtinrows_final] =  GTI_saturation_TSTOP[i];
					
					
				
					gtinrows_final++;
				}
				else if((gtitstart_f[j]<=GTI_saturation_TSTART[i]) && (gtitstop_f[j]>=GTI_saturation_TSTOP[i]))
				{
				
				
					gtitstart_final[gtinrows_final] =  GTI_saturation_TSTART[i];
					gtitstop_final[gtinrows_final] = GTI_saturation_TSTOP[i];
				

					gtinrows_final++;
				}
				
				
			}
			
		}
		printf("GTI correction Over\n");
		//for(i=0;i<gtinrows_final;i++)
		//{
		//	printf("%lf\t%lf\t%lf\n",gtitstart_final[i],gtitstop_final[i],gtinrows_final);
		//}
		exposure=0.0;
		
	
		
		//printf("%ld %ld \n",gtinrows_final,gtinrows_est_final);
		for(i=0;i<gtinrows_final;i++)
		{
			if(gtitstart_final[i] >= gtitstop_final[i])
			{
				printf("ERROR---LAST----GTI  Tstart >= Tstop---------\n");
				printf("%lf \t %lf \n",gtitstart_final[i],gtitstop_final[i]);
			}
			exposure = exposure+(gtitstop_final[i]-gtitstart_final[i]);
			
		}
		
		
		
		//exposure = exposure - exp_loss_due_bunches;
		
		printf("Exposure loss due to bunches = %lf\n",exp_loss_due_bunches);
		printf("Final Exposure = %lf\n",exposure);
		writeLivetime(lttime_new,fracexp_new,livetimenrows_new,qid+2);
		//livetimeGeneration(live_time_UT,live_counter,outlivetimefile,qid+2,lttime_new,fracexp_new,livetimenrows_new);
		
		
		
		//~ writegtiextension(gtitstart_final, gtitstop_final,gtinrows_final,qid+10,tstart,tstop,exposure);
		writegtiextension(gtitstart_final, gtitstop_final,gtinrows_final,extname,tstart,tstop,exposure);
		
		
		
		
		
		/*
		
		for(i=0;i<BTI_counter;i++)
		{
			for(j=0;j<qgtinrows;j++)
			{
				if((BTITSTART[i] >= gtitstart[j] && BTITSTART[i] <= gtitstop[j]) && (BTITSTOP[i] >= gtitstart[j] && BTITSTOP[i] <= gtitstop[j]))
				{
					new_gtitstart[new_qgtinrows] = gtitstart[j];
					new_gtitstop[new_qgtinrows] = BTITSTART[i];
					new_qgtinrows++;
					new_gtitstart[new_qgtinrows] = BTITSTOP[i];
					new_gtitstop[new_qgtinrows] = gtitstop[j];
					new_qgtinrows++;
					break;
				}
				
				
				
				
			}
			
		}
			*/
		
		for(i=0;i<evtnrows;i++)
		{
				//sum =sum+evt_flag[i];
				//fprintf(f1,"%d\t%lf\t%d\n",qid,evttime[i],evt_flag[i]);
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
		
		}
		
		
	        writeEvent(finalevttime,finalcztseccnt,finalcztntick,finalpha,finaldetid,finalpixid,finaldetx,finaldety,finalveto,finalalpha,finalpi,finalenergy,l,qid+2,exposure);
		printf("%d end\n",qid);
	free(live_time_UT);free(BTITSTART);free(BTITSTOP);
	 	free(gtitstart_final);free(gtitstop_final);
		free(GTI_saturation_TSTART);free(GTI_saturation_TSTOP);
		free(gtitstart_f);free(gtitstop_f);
		free(fracexp_new);free(lttime_new);free(BTI_exclude_TSTART);free(BTI_exclude_TSTOP);
		free(BTITSTART_merge);free(BTITSTOP_merge);
		free(GTITSTART_merge);free(GTITSTOP_merge);
		free(gtitstart);free(gtitstop);free(new_gtitstart);
		free(new_gtitstop);free(gtitstart_sel);free(gtitstop_sel);
		


	}//qid end
	modifyEventHeaderParams(outfile);
	//modifyEventHeaderParams(outlivetimefile);
	//free memory
	free(evttime);free(evtcztseccnt);free(finalevttime);free(finalcztseccnt);
	free(evtcztntick);free(evtveto);free(finalcztntick);free(finalveto);free(evtpha);free(finalpha);
	free(evtdetid);free(evtpixid);free(evtdetx);free(evtdety);free(evtalpha);free(finaldetid);free(finalpixid);
	free(finaldetx);free(finaldety);free(finalalpha);
	free(evtenergy);free(finalenergy);
	free(evtpi);free(finalpi);

	free(buntime);free(buntime_heavy);free(buntime_real);
	free(buntime_dfs);free(buntime_dsl);
	free(bunsize);free(bun_flag);free(buntime_del1);free(buntime_del2); 
	free(lttime);free(fracexp);free(evt_flag);free(bunsize_real);
	
	fclose(thrfile);fclose(f1);
	free(tempevt);	

	fits_close_file(fevt,&status);
	if(status) { printerror( status );}

	fits_close_file(fbun,&status);
	if(status) { printerror( status );}

	fits_close_file(flivetime,&status);
	if(status) { printerror( status );}

	printf("\nSUPER BUNCH CLEAN COMPLETED SUCCESSFULLY.\n");

	return;
}

/*
void livetimeGeneration(double *live_time_UT,int live_counter,char *outputfile,int hdunum,double *lttime,double *fracexp,long livetimenrows)
{
	int i,j,temp=0;
	double init_fracexp[livetimenrows];
	
	double fracexp_new[livetimenrows], lttime_new[livetimenrows];
	long livetimenrows_new=0;
	
	double livetime_resolution = (lttime[1]-lttime[0]);
	
	
	//To get the corresponding livetime
	for(i=0;i<livetimenrows;i++)
	{
		if( lttime[i]>=(double)(long)tstart && lttime[i]<(double)(long)tstop )
		{
			lttime_new[livetimenrows_new] = lttime[i];
			fracexp_new[livetimenrows_new] = fracexp[i];
			livetimenrows_new++;
		}
		
		
	}
	
	for(i=0;i<livetimenrows_new;i++)
	{
		
		
	}
	

	for(i=0;i<livetimenrows;i++)
	{
		if(fracexp[i]<=0.1){continue;}
		init_fracexp[i]=fracexp[i];
		for(j=temp;j<live_counter;j++)
		{
			if(live_time_UT[j]>lttime[i]-0.5 && live_time_UT[j]<=lttime[i]+0.5)
			{
				if((int)live_time_UT[j]==(int)(live_time_UT[j]+dph_time))
					fracexp[i]-=dph_time;
				else
				{
					fracexp[i]-=lttime[i]+0.5-live_time_UT[j];
					if(i<livetimenrows-1)
						fracexp[i+1]-=(live_time_UT[j]+dph_time)-(lttime[i+1]-0.5);
				}
				temp=j;
			}
			if(live_time_UT[j]-lttime[i]>100.0)
				break;
		}
	}
	for(i=0;i<livetimenrows;i++)
	{
		if(fracexp[i] != init_fracexp[i])
		{		
			printf("%f\t%f\t%f\n",lttime[i],fracexp[i],init_fracexp[i]);
		}	
	}
	writeLivetime(lttime,fracexp,livetimenrows,outputfile,hdunum);
		
}
*/



void livetimeGeneration(double *live_time_UT,int live_counter,char *outputfile,int hdunum,double *lttime,double *fracexp,int livetimenrows)
{
	int i,j,temp=0;
	double init_fracexp[livetimenrows];

	for(i=0;i<livetimenrows;i++)
	{
		if(fracexp[i]<=0.1){continue;}
		init_fracexp[i]=fracexp[i];
		for(j=temp;j<live_counter;j++)
		{
			if(live_time_UT[j]>lttime[i]-0.5 && live_time_UT[j]<=lttime[i]+0.5)
			{
				if((int)live_time_UT[j]==(int)(live_time_UT[j]+dph_time))
					fracexp[i]-=dph_time;
				else
				{
					fracexp[i]-=lttime[i]+0.5-live_time_UT[j];
					if(i<livetimenrows-1)
						fracexp[i+1]-=(live_time_UT[j]+dph_time)-(lttime[i+1]-0.5);
				}
				temp=j;
			}
			if(live_time_UT[j]-lttime[i]>100.0)
				break;
		}
	}
	for(i=0;i<livetimenrows;i++)
	{
		if(fracexp[i] != init_fracexp[i])
		{		
			//printf("%f\t%f\t%f\n",lttime[i],fracexp[i],init_fracexp[i]);
		}	
	}
	writeLivetime(lttime,fracexp,livetimenrows,hdunum);
		
}





//creating event file 
void createEventFile(char *eventfile)
{
	
	fitsfile *fptrOut,*fptrevt;      
	int status, hdutype,tfields=12,i,hdunum=2;
	long frow, felem;
	int mjdrefi=55197,mjdreff=0,equinox=2000;
	float ra_pnt,dec_pnt;
	double timedel,telapse;
	char object[20],obs_id[20],obs_mode[20],date_obs[20],time_obs[20],date_end[20],time_end[20],date[20],creator[20],filename[70],checksum[20],datasum[20],chksumcomm[50],datasumcomm[50];
      
	char extname[20]; 
	          
	
	char *ttype[] = { "TIME", "CZTSECCNT","CZTNTICK","PHA","DetID","pixID","DETX","DETY","veto","alpha","PI","ENERGY"};
	char *tform[] = { "D","D","I","I","B","B","B","B","I","B","I","E"};
	char *tunit[] = {"s","s","micro-sec","counts","","","","","counts","counts","",""};
       
	status=0;
        if (fits_create_file(&fptrOut, outfile, &status))
	       	 printerror( status );       

	//if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
	  //      printerror( status );
	
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
	//~ if(fits_movnam_hdu(fptrevt, BINARY_TBL, "GTI", 0, &status)) 
		//~ printerror( status );         

	//~ if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
		//~ printerror( status );
		
		
		
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
	for(i=0;i<4;i++)
	{
		
		status=0;
		sprintf(extname, "Q%d_GTI",i);
		
		
		//~ if(fits_movabs_hdu(fptrevt,i+8+hdunum, 0, &status)) 
			//~ printerror( status );
			
		if(fits_movnam_hdu(fptrevt, BINARY_TBL, extname, 0, &status)) 
			printerror( status );   
		
		//if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        	//printerror( status );
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields_gti, ttype_gti, tform_gti,tunit_gti, extname, &status) )
			printerror( status );
			
		//~ if ( fits_movabs_hdu(fptrOut, i+10, &hdutype, &status) ) 
	         	//~ printerror( status );
	         	
	    if(fits_movnam_hdu(fptrOut, BINARY_TBL, extname, 0, &status)) 
			printerror( status );  
	         	
	         	
		//if(fits_copy_hdu(fptrevt, fptrOut, 0, &status))
			//printerror( status );

	}
	

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
void createLivetimeFile()
{
	fitsfile *fptrOut, *fptrlive;    
	int status, hdutype;
	double livetimeres;
	long frow, felem;
      	int tfields,i;
	char extname[20];
	
	tfields=2;
	char *ttype[] = {"TIME","FRACEXP"};
	char *tform[] = { "D","E"};
	char *tunit[] = {"",""};
	        
	int hdunum=2;
	
	status=0;
        if (fits_create_file(&fptrOut, outlivetimefile, &status))
	       	 printerror( status );       

	//if ( fits_open_file(&fptrOut, outlivetimefile, READWRITE, &status) ) 
	 //        printerror( status );
	         
	if ( fits_open_file(&fptrlive,inlivetimefile, READONLY, &status) ) 
	         printerror( status );    
	

	
	
	//fits_write_key(fptrOut,TDOUBLE,"HIERARCH LV_BINSIZE",&livetimeres,"Livetime binsize", &status);
	         
	       

    for(i=0;i<4;i++)
    {       
		sprintf(extname, "Q%d",i);
		if ( fits_movabs_hdu(fptrlive, i+hdunum, &hdutype, &status) ) 
	         	//printerror( status );
		status=0;
		fits_read_key(fptrlive,TDOUBLE,"HIERARCH LV_BINSIZE",&livetimeres,NULL, &status);
		
		
		
		if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields, ttype, tform,tunit, extname, &status) )
			printerror( status );
		if ( fits_movabs_hdu(fptrOut, i+hdunum, &hdutype, &status) ) 
					printerror( status );
		fits_write_key(fptrOut,TDOUBLE,"HIERARCH LV_BINSIZE",&livetimeres,"Livetime binsize", &status);
	}
	if ( fits_movabs_hdu(fptrOut, 1, &hdutype, &status) ) 
					printerror( status );
	fits_write_key(fptrOut,TDOUBLE,"HIERARCH LV_BINSIZE",&livetimeres,"Livetime binsize", &status);
	
	
	if ( fits_close_file(fptrlive, &status) )       
	       printerror( status );
	       
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;

}
/*
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
	
    if(fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

	if(fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         printerror( status );

 
	if ( fits_create_tbl( fptrOut, BINARY_TBL, 0, tfields, ttype, tform,tunit, extname, &status) )
		printerror( status );
	if ( fits_movabs_hdu(fptrOut, i+hdunum, &hdutype, &status) ) 
	       	printerror( status );
	
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;

}*/


void writegtiextension(double *gtistart, double *gtistop, double gtinrows, char extention_name[20], double tstart, double tstop, double exposure)
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
	if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
	        printerror( status );
	        
	if (   fits_movnam_hdu(fptrOut, BINARY_TBL, extention_name, 0, &status) )
			printerror( status );
	//~ if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	
	 

	status=0;
	printf("Final Exposure = %lf\n",exposure);
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




void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,int writesize,int hdunum,double exposure)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	status=0;
	int tstarti,tstopi;
	double tstart,tstop,tstartf,tstopf;
	//,exposure;

	
	if ( fits_open_file(&fptrOut, outfile, READWRITE, &status) ) 
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

void writeLivetime(double *time,double *livetime,int livetimesize,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	double tstart,tstop;
	status=0;
	
	if ( fits_open_file(&fptrOut, outlivetimefile, READWRITE, &status) ) 
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
        fits_write_col(fptrOut, TDOUBLE, 2,frow, felem, livetimesize, livetime,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}


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






/*
void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i,tstarti,tstopi;
	double tstart,tstop,start[4],stop[4],largest,smallest,exposure,tstartf,tstopf;
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

	if ( fits_movabs_hdu(fptrOut, 1, NULL, &status) ) 
	      	printerror( status );
	tstarti=(int)smallest;
	tstopi=(int)largest;

	tstartf=smallest-tstarti;
	tstopf=largest-tstopi;
	exposure=largest-smallest;
	
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


void createBTIFile()
{
	fitsfile *fptrOut;      
	int status, hdutype;
	long frow, felem;
      	int tfields,i;
	char extname[20];
	
	tfields=2;
	char *ttype[] = {"TSTART","TSTOP"};
	char *tform[] = { "D","D"};
	char *tunit[] = {"",""};
	        
	int hdunum=2;
	
	status=0;
	
        if (fits_create_file(&fptrOut, outbtifile, &status))
	       	 printerror( status );       

//	if ( fits_open_file(&fptrOut, outbtifile, READWRITE, &status) ) 
//	         printerror( status );

    	for(i=0;i<=3;i++)
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

void writeBTItime(double *time_1,double *time_2,int BTIsize,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype;
	long frow, felem;
	status=0;
	
	if ( fits_open_file(&fptrOut, outbtifile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	frow      = 1;
	felem     = 1;

    fits_write_col(fptrOut, TDOUBLE, 1, frow, felem, BTIsize,time_1,&status);
    fits_write_col(fptrOut, TDOUBLE, 2,frow, felem, BTIsize, time_2,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}

int read_input_parameters(int r,char *infile, char *inbunchfile, char *thresholdfile,char *inlivetimefile,char *outlivetimefile,char *outbtifile,char *outfile,int *clobber,int *history)
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
	
	r=PILGetFname("inbunchfile",inbunchfile);
	if (fits_open_file(&fptr, inbunchfile, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,inbunchfile);
		printerror( status );
	}
	fits_close_file(fptr,&status);

	r=PILGetFname("thresholdfile",thresholdfile);
	if((fp = fopen(thresholdfile, "r")) == NULL) {
		 printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,thresholdfile);
	}
	fclose(fp);
	
	r=PILGetFname("inlivetimefile",inlivetimefile);
	if (fits_open_file(&fptr, inlivetimefile, READONLY, &status))
	{
		printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,inlivetimefile);
		printerror( status );
	}
	fits_close_file(fptr,&status);
	
	r=PILGetFname("outlivetimefile",outlivetimefile);
	r=PILGetFname("outbtifile",outbtifile);
	r=PILGetFname("outfile",outfile);
	
	//---------------------------------------------------------Overwrite files and remember history------------------------------------------------------------
	r=PILGetBool("clobber", clobber);
	r=PILGetBool("history", history);
	//-------------------------------------------------------------End of input paramters----------------------------------------------------------------------	
	PILClose(r);	//Closing PIL file
	
	return(status);
}

int display_input_parameters(char *infile, char *inbunchfile, char *thresholdfile,char *inlivetimefile,char *outlivetimefile,char *outbtifile,char *outfile,int clobber,int history)
{
	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf("                                                    CZTSUPERBUNCHCLEAN PARAMETERS \n");
    	printf("----------------------------------------------------------------------------------------------------------------------------\n");
    	printf(" Input Event file              		: %s\n",infile);
	printf(" Input Bunch file             		: %s\n",inbunchfile);
    	printf(" Input Livetime file	              	: %s\n",inlivetimefile);
	printf(" Threshold file	              		: %s\n",thresholdfile);
	printf(" Output Livetime file	              	: %s\n",outlivetimefile);
	printf(" Output BTI file	              	: %s\n",outbtifile);
	printf(" Output Event file 			: %s\n",outfile);
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









        
        
        
        
        
        
        
      

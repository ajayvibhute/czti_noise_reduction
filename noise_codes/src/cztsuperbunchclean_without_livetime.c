/* CZTSUPERBUNCHCLEAN
 * Author	: Ajay Ratheesh, TIFR Mumbai
 * Date 	: 14-08-2017
 */

#include<stdio.h>
#include<math.h>
#include<fitsio.h>
#include<stdlib.h>
#include<string.h>
#include<float.h>
#define BUFF_SIZE 2048

//Following are the parameters
#define heavy_bun_thresh_count_10ms 2
#define heavy_bun_thresh_count_1ms 1
#define bun_size_thresh 62
#define dph_time 0.1
#define thresh_hot 0.707
#define allowable_hot_thresh 3.0


FILE *f1;
void printerror( int status);
void processSuperBunchClean(char *input_evtfile,char *input_bunfile,char *orbitno);
void createEventFile(char *outputfile,char *eventfile);
//void createLivetimeFile(char *outputfile);
void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int bufsize,int hdunum);
void writeLivetime(double *time,double *livetime,int livetimesize,char *outputfile,int hdunum);
void livetimeGeneration(double *live_time_UT,int live_counter,char *outputfile,int hdunum,double *lttime,double *fracexp,int livetimenrows);
void modifyEventHeaderParams(char *outputfile);
void writeBTItime(double *time_1,double *time_2,int BTIsize,char *outputfile,int hdunum);
void createBTIFile(char *outputfile);


int main(int argc,char *argv[])
{

	if(argc==4)
        {
		processSuperBunchClean(argv[1],argv[2],argv[3]);
		printf("\nSUPER BUNCH CLEAN COMPLETED SUCCESSFULLY.\n");
	}
	else
		printf("Please give event file,bunch file names and orbit number.");
	
	
	return 0;

}

//void processSuperBunchClean(char *input_evtfile,char *input_bunfile,char *livetimefile,char *orbitno)
void processSuperBunchClean(char *input_evtfile,char *input_bunfile,char *orbitno)
{
	int tmpindex=0,tmpindex1=0;
	int BTI_counter=0;
	unsigned int i,j,k,m,e,n,p,sum;
	unsigned int qid;
	int evttime_col,evtdetx_col,evtdety_col,evtPI_col,evtdetid_col,evtpixid_col,buntime_col,bunsize_col,evttime_index;
	double *buntime,*buntime_heavy;
	double bun;
	int *evt_flag,*bun_flag;
	int *bunsize;
	long evtnrows, bunnrows,livetimenrows;
	int status=0,hdutype=0;
	int magnitude_temp[4096],detx_temp[4096],dety_temp[4096];
	int detx_close_bunch[3072],dety_close_bunch[3072];
	long frow, felem;
	double *evttime,*evtcztseccnt,*finalevttime,*finalcztseccnt,doublenull;
	unsigned short *evtpha,*finalpha,*evtcztntick,*evtveto,*finalcztntick,*finalveto;
	int intnull,anynull,*evtpi,*finalpi;
	unsigned char *evtdetid,*evtpixid,*evtdetx,*evtdety,*evtalpha,*finaldetid,*finalpixid,*finaldetx,*finaldety,*finalalpha,bytenull;
	float *evtenergy,floatnull,*finalenergy;
	int BTItimesize;
	double *BTITSTART, *BTITSTOP;
	
	
	int live_counter=0,livetimesize, BTIcounter;

	unsigned int counter1, counter2,counter3, counter4,counter5, counter6,counter7, counter8,counter9,flag;
	int  buntime_heavy_index;
	unsigned int dph[64][64], dph_hot[64][64];
	float distance;
	float cij;
	float total_hotness =0.0 ,n_pair =0.0,n_hot_pix = 0.0,allowable_hot,d;

	char outtxtfile[100];
	char outbtifile[100];
	char outevtfile[100];
	double *lttime,*fracexp;

	char* tempevt = calloc(strlen(input_evtfile)+1, sizeof(char));
	strcpy(tempevt, input_evtfile);

	printf("Input event file is %s\n",input_evtfile);  
	fitsfile *fevt, *fbun,*flivetime;	
	fits_open_file(&fevt,input_evtfile,READONLY,&status);
	if(status) { printerror( status );}

	fits_open_file(&fbun,input_bunfile,READONLY,&status);
	if(status) { printerror( status );}

	//fits_open_file(&flivetime,livetimefile,READONLY,&status);
	//if(status) { printerror( status );}

	char *file = strtok(input_evtfile, "."); 
	sprintf(outevtfile, "%s_sbc.evt",file);
	sprintf(outtxtfile, "%s_sbc.txt",file);
	//sprintf(outlivetimefile, "%s_sbc_livetime.fits",file);
	sprintf(outbtifile, "%s_sbc_bti.fits",file);
	f1=fopen(outtxtfile,"a");      

	remove(outevtfile);
	remove(outbtifile);
	createEventFile(outevtfile,tempevt);
	//createLivetimeFile(outlivetimefile);
	createBTIFile(outbtifile);

	frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;
	floatnull=0.0;
 
	for(qid=0;qid<4;qid++)
	{   
		//live_counter=0;
		BTI_counter=0;
		fits_movabs_hdu(fevt, qid+2, &hdutype, &status);
		fits_movabs_hdu(fbun, qid+2, &hdutype, &status);
		//fits_movabs_hdu(flivetime, qid+2, &hdutype, &status);

		fits_get_num_rows(fevt, &evtnrows, &status);
		fits_get_num_rows(fbun, &bunnrows, &status);
		//fits_get_num_rows(flivetime, &livetimenrows, &status);
	
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
	    
	    fits_get_colnum(fevt,CASEINSEN,"TIME",&evttime_col,&status);
	    fits_get_colnum(fevt,CASEINSEN,"DETID",&evtdetid_col,&status);		
	    fits_get_colnum(fevt,CASEINSEN,"PIXID",&evtpixid_col,&status);
	    fits_get_colnum(fevt,CASEINSEN,"DETX",&evtdetx_col,&status);
	    fits_get_colnum(fevt,CASEINSEN,"DETY",&evtdety_col,&status);
	    
	    fits_get_colnum(fbun,CASEINSEN,"TIME",&buntime_col,&status);
	    fits_get_colnum(fbun,CASEINSEN,"NUMEVENT",&bunsize_col,&status);
	    
	    evttime=(double*)malloc(sizeof(double)*evtnrows);
		evtdetid=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
		evtpixid=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
		evtdetx=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
		evtdety=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
		evtcztseccnt  = (double*)malloc(evtnrows * sizeof(double));
		evtpha  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
		evtcztntick  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
		evtveto  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
		evtalpha = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
		evtenergy = (float*)malloc(sizeof(float)*evtnrows);
		evtpi = (int*)malloc(sizeof(int)*evtnrows);
		evt_flag=(int*)malloc(sizeof(int)*(evtnrows));

		finalevttime  = (double*)malloc(evtnrows * sizeof(double));
		finalcztseccnt  = (double*)malloc(evtnrows * sizeof(double));
		finalpha  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
		finalcztntick  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
		finalveto  = (unsigned short*)malloc(evtnrows * sizeof(unsigned short));
		finaldetid  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
		finalpixid  =(unsigned char*) malloc(evtnrows * sizeof(unsigned char));
		finaldetx  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
		finaldety  = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
		finalalpha = (unsigned char*)malloc(evtnrows * sizeof(unsigned char));
		finalenergy = (float*)malloc(sizeof(float)*evtnrows);
		finalpi = (int*)malloc(sizeof(int)*evtnrows);

	
		buntime=(double*)malloc(sizeof(double)*bunnrows);
		buntime_heavy=(double*)malloc(sizeof(double)*bunnrows);
		bunsize=(int*)malloc(sizeof(int)*bunnrows);
		bun_flag=(int*)malloc(sizeof(int)*bunnrows);

		//lttime=(double*)malloc(sizeof(double)*livetimenrows);
		//fracexp=(double*)malloc(sizeof(double)*livetimenrows);
	
		

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
		fits_read_col(fevt, TINT, 11, frow, felem, evtnrows, &floatnull, evtpi,&anynull, &status);
		fits_read_col(fevt, TFLOAT, 12, frow, felem, evtnrows, &floatnull, evtenergy,&anynull, &status);
	
	
		fits_read_col(fbun, TDOUBLE, buntime_col, 1, 1, bunnrows, NULL, buntime,NULL, &status);
		fits_read_col(fbun, TINT, bunsize_col, 1, 1, bunnrows, NULL, bunsize,NULL, &status);

		//fits_read_col(flivetime, TDOUBLE, 1, 1, 1, livetimenrows, NULL, lttime,NULL, &status);
		//fits_read_col(flivetime, TDOUBLE, 2, 1, 1, livetimenrows, NULL, fracexp,NULL, &status);

		BTItimesize=((int)(evttime[evtnrows-1])-(int)(evttime[0])+1)*10;
		
		//live_time_UT=(double*)malloc(sizeof(double)*BTItimesize);
		BTITSTART=(double*)malloc(sizeof(double)*BTItimesize);
		BTITSTOP=(double*)malloc(sizeof(double)*BTItimesize);

		for(i=0;i<evtnrows;i++){evt_flag[i]=0;}
		int heavy_bunch_length=0;
		for(i=0;i < bunnrows;i++)
		{
			if(bunsize[i] > bun_size_thresh)
			{ 
				buntime_heavy[heavy_bunch_length++] = buntime[i];
			}
		}

		bun_flag=(int*)malloc(sizeof(int)*heavy_bunch_length);
		for(i=0;i<heavy_bunch_length;i++){bun_flag[i]=0;}
		//printf("%lf\t%lf\n",buntime[bunnrows-1],buntime[0]);

		tmpindex=0;
		for(bun=buntime[0];bun<buntime[bunnrows-1];bun=bun+0.001)
		{
			//printf("%lf\n",bun);

			counter1=0;
			//printf("%lf\n",bun);
			for(j=tmpindex;j<heavy_bunch_length;j++)
			{
				if((int)buntime_heavy[j] > (int)(bun+1)){break;}
		
				if((buntime_heavy[j]>bun) && (buntime_heavy[j]<bun+0.001))
				{	
					tmpindex = j;
					if(counter1==0)
					{
						buntime_heavy_index = j;
					}
					if(buntime_heavy[buntime_heavy_index]>=buntime_heavy[j]);
					{
					buntime_heavy_index = j;
					}
					counter1++;
			
				}
				//if(counter1==0)
					//	tmpindex=0;
					//else
					//	tmpindex=buntime_heavy_index-1;
			
			
			}
	
			if(counter1 >= heavy_bun_thresh_count_1ms)
			{	//printf("True-%d\t%d\n",buntime_heavy_index,qid);
				bun_flag[buntime_heavy_index]=1;
		    }
	
		}


		tmpindex=0;
		for(bun=buntime[0];bun<buntime[bunnrows-1];bun=bun+0.01)
		{
			//printf("%lf\n",bun);

			counter1=0;
			//printf("%lf\n",bun);
			for(j=tmpindex;j<heavy_bunch_length;j++)
			{
				if((int)buntime_heavy[j] > (int)(bun+1)){break;}
		
				if((buntime_heavy[j]>bun) && (buntime_heavy[j]<bun+0.01))
				{	
					tmpindex = j;
					if(counter1==0)
					{
						buntime_heavy_index = j;
				
					}
					if(buntime_heavy[buntime_heavy_index]>=buntime_heavy[j]);
					{
					buntime_heavy_index = j;
					}
					counter1++;
			
				}
				//if(counter1==0)
					//	tmpindex=0;
					//else
					//	tmpindex=buntime_heavy_index-1;
			
			
			}
	
			if(counter1 >= heavy_bun_thresh_count_1ms)
			{	//printf("True-%d\t%d\n",buntime_heavy_index,qid);
				bun_flag[buntime_heavy_index]=1;
		    	}
	
		}

		/*
		for(i=0;i<heavy_bunch_length;i++)
			{
				fprintf(f1,"%d\t%lf\t%d\n",qid,buntime_heavy[i],bun_flag[i]);	
			}
		*/
		tmpindex = 0;
		tmpindex1 =0;
		printf("%d\n",heavy_bunch_length);
		for(i=0;i<heavy_bunch_length;i++)
		{
			//printf("%d\t%d\n",qid,i);
			if(bun_flag[i]==1)
			{
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
					if((evttime[j] > buntime_heavy[i]) && (evttime[j] <= buntime_heavy[i]+dph_time))
					{	
						tmpindex = j;
						//if(counter4==0){evttime_index = j;}
						//time_evt_close_bunch[counter4] = evttime[i];
						detx_close_bunch[counter2] = evtdetx[j];
						dety_close_bunch[counter2] = evtdety[j];
						counter2++;
				    }
				}
	

				if(counter2==0){continue;}
	
	
	
				for(j=0;j<counter2;j++)
				{
					//printf("%d\n",j);
					//printf("%d\t%d\n",detx_close_bunch[j],dety_close_bunch[j]);
					dph[detx_close_bunch[j]][dety_close_bunch[j]]++;
		
				}
	
	
	
	
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
						//printf("%d\t%d\n",k,m);
			
						//if(qid==2 && i==7799){ printf("%d\t%d\t%d\n",m,k,counter3);}
			
						if((dph[detx_temp[m]][dety_temp[m]]==1) && (dph[detx_temp[k]][dety_temp[k]]==1))//to avoid isolated genuine double events
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
			
				
						if(counter4>2)
						{
							d = sqrt(pow((detx_temp[m]-detx_temp[k]),2)+pow((dety_temp[m]-dety_temp[k]),2));
							cij = magnitude_temp[k]*magnitude_temp[m]/d;
		
							if(cij > thresh_hot){total_hotness = total_hotness+cij;n_pair++;dph_hot[k][m]=1;}
						}
					//printf("%d\t%d\n",m,k);

				 	}
				}	
		
	
				//printf("%d\n",counter3);
	
				for(k=0;k<64;k++)
				{
					for(m=0;m<64;m++)
					{//printf("%d\t%d\n",k,m);
						n_hot_pix=n_hot_pix+dph_hot[k][m];
					}
				}
	
				allowable_hot = total_hotness/n_hot_pix;
				//need to change(orbit no)
				//fprintf(f1,"%d\t%f\t%f\t%f\t%f\t%s\n",bunch_size,allowable_hot,total_hotness,n_hot_pix,n_pair,orbitno);
				fprintf(f1,"%f\t%f\t%f\t%f\t%s\t%d\n",allowable_hot,total_hotness,n_hot_pix,n_pair,orbitno,qid);
		
			
		
				if (allowable_hot>allowable_hot_thresh)
				{	
					BTITSTART[BTI_counter] = buntime_heavy[i];
					BTITSTOP[BTI_counter]  = buntime_heavy[i]+0.1;
					BTI_counter++;
					
					//live_time_UT[live_counter++] = buntime_heavy[i];
					//printf("%lf\t%lf\t%d\n",buntime_heavy[i],live_time_UT[live_counter-1],live_counter);
					for(e=tmpindex1;e<evtnrows;e++)
					{
						if((int)evttime[e] > (int)buntime_heavy[i]+1){break;}
						if((evttime[e] > buntime_heavy[i]) && (evttime[e] <= buntime_heavy[i]+dph_time))
						{
							tmpindex1=e;
							evt_flag[e] = 1;
						}
					}
				}
			}
		}
		
		//livetimeGeneration(live_time_UT,live_counter,outlivetimefile,qid+2,lttime,fracexp,livetimenrows);
		//BTIGeneration(live_time_UT,live_counter,outlivetimefile,qid+2);
		
		writeBTItime(BTITSTART,BTITSTOP,BTI_counter,outbtifile,qid+2);
		
		int l=0;
		//sum=0;

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
	        writeEvent(finalevttime,finalcztseccnt,finalcztntick,finalpha,finaldetid,finalpixid,finaldetx,finaldety,finalveto,finalalpha,finalpi,finalenergy,outevtfile,l,qid+2);
	}//qid end
	modifyEventHeaderParams(outevtfile);

	return;
}

/*
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
					fracexp[i+1]=(live_time_UT[j]+dph_time)-(lttime[i+1]-0.5);
				}
				temp=j;
			}
			if(live_time_UT[j]-lttime[i]>100.0)
				break;
		}
	}
		//writeLivetime(,,,outputfile,hdunum);
		
		//writeLivetime(lttime,fracexp,livetimenrows,outputfile,hdunum);
		
}*/



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
	
        if (fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         printerror( status );

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

}*/




void createBTIFile(char *outputfile)
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
	
        if (fits_create_file(&fptrOut, outputfile, &status))
	       	 printerror( status );       

	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	         printerror( status );

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


void writeEvent(double *evttime,double *cztseccnt,unsigned short *cztntick,unsigned short *pha,unsigned char *detid,unsigned char *pixid,unsigned char*detx,unsigned char *dety,unsigned short *veto,unsigned char *alpha,int *pi,float *energy,char *outputfile,int writesize,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype,intnull;
	long frow, felem,nrows;
	double doublenull;
	unsigned char bytenull;
	status=0;
	int tstarti,tstopi;
	double tstart,tstop,tstartf,tstopf,exposure;

	
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


/*
void writeLivetime(double *time,double *livetime,int livetimesize,char *outputfile,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype,intnull;
	long frow, felem,nrows;
	double doublenull;
	unsigned char bytenull;
	status=0;
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;

        fits_write_col(fptrOut, TDOUBLE, 1, frow, felem,livetimesize ,time,&status);
        fits_write_col(fptrOut, TDOUBLE, 2,frow, felem, livetimesize, livetime,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}*/

void writeBTItime(double *time_1,double *time_2,int BTIsize,char *outputfile,int hdunum)
{
	fitsfile *fptrOut;       
	int status, hdutype,intnull;
	long frow, felem,nrows;
	double doublenull;
	unsigned char bytenull;
	status=0;
	
	if ( fits_open_file(&fptrOut, outputfile, READWRITE, &status) ) 
	        printerror( status );
	if ( fits_movabs_hdu(fptrOut, hdunum, &hdutype, &status) ) 
	      	printerror( status );

	frow      = 1;
	felem     = 1;
	doublenull = 0.;
	intnull = 0;
	bytenull = 0;

    fits_write_col(fptrOut, TDOUBLE, 1, frow, felem, BTIsize,time_1,&status);
    fits_write_col(fptrOut, TDOUBLE, 2,frow, felem, BTIsize, time_2,&status);
	
	if ( fits_close_file(fptrOut, &status) )       
	        printerror( status );
	return;
}


void modifyEventHeaderParams(char *outputfile)
{
	
	fitsfile *fptrOut;  
	int status,i;
	double tstart,tstop,start[4],stop[4],largest,smallest;
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

	fits_update_key(fptrOut, TDOUBLE,"TSTART",&smallest,"Start time of observation",&status);
	fits_update_key(fptrOut, TDOUBLE,"TSTOP",&largest,"Stop time of observation",&status);
	
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









        
        
        
        
        
        
        
      

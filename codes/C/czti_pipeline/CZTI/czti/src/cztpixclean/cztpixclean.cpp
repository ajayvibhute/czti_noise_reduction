/*cztpixclean.c

This routine identifies the noisy pixels and removes the events from 
those pixels. It also looks for pixels or modules having higher counts 
during T second period and ignores them for that duration.

Mithun N P S 
25/11/15

Edits:

**29-30/11/15
* Converted to cpp and PIL added 
* Badpixel file reading and writing added (needs further updation)
* Exposure extension modification added

08/12/15
*Corrected exposure fraction for cztclock, and changed 
to fractional instead of total exposure time

*/

#include "cztpixclean.h"

//Constructor
cztpixclean::cztpixclean(){
	strcpy(modulename, "cztpixclean_v");
	strcat(modulename, VERSION);
}

//Destructor
cztpixclean::~cztpixclean(){}


int cztpixclean::read(int argc,char **argv)
{
    int status=0;

    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error Initializing PIL***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("par_infile",infile))){
        LOG(ERROR)<<"***Error reading input event file"<<infile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("par_inlivetimefile",inlivetimefile))){
        LOG(ERROR)<<"***Error reading input livetime file"<<inlivetimefile<<"***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("par_writedblevt", &writedblevt))) {
        LOG(ERROR)<<"***Error reading input writedblevt(create dblevt or not)"<<writedblevt<<"***";
		return (EXIT_FAILURE);
    }

    if(PIL_OK!=(status=PILGetFname("par_outfile1",outfile1))){
        LOG(ERROR)<<"***Error reading output single event file"<<outfile1<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("par_outfile2",outfile2))){
        LOG(ERROR)<<"***Error reading output double event file"<<outfile2<<"***";
        return status;
    }
	
    if(PIL_OK!=(status=PILGetFname("par_outlivetimefile",outlivetimefile))){
        LOG(ERROR)<<"***Error reading output livetime file"<<outlivetimefile<<"***";
        return status;
    }


    /*if(PIL_OK!=(status=PILGetFname("par_caldb_badpix",caldb_badpix))){
        LOG(ERROR)<<"***Error reading caldb badpixel file"<<caldb_badpix<<"***";
        return status;
    }*/

    if(PIL_OK!=(status=PILGetFname("par_badpixfile",badpixfile))){
        LOG(ERROR)<<"***Error reading output badpixel file"<<badpixfile<<"***";
        return status;
    }

    if (PIL_OK != (status=PILGetInt("par_nsigma", &nsigma))){
        LOG(ERROR) << "***Error reading nsigma for bad pixel selection.***";
        return  EXIT_FAILURE;
    }   


    if (PIL_OK != (status=PILGetReal4("par_det_tbinsize", &det_tbinsize))){
        LOG(ERROR) << "***Error reading detector lc time binsize.***";
        return  EXIT_FAILURE;
    } 

    if (PIL_OK != (status=PILGetReal4("par_pix_tbinsize", &pix_tbinsize))){
        LOG(ERROR) << "***Error reading pixel lc time binsize.***";
        return  EXIT_FAILURE;
    }

    if (PIL_OK != (status=PILGetInt("par_det_count_thresh", &det_count_thresh))){
        LOG(ERROR) << "***Error reading detector count threshold.***";
        return  EXIT_FAILURE;
    }
    

	if (PIL_OK != (status=PILGetInt("par_pix_count_thresh", &pix_count_thresh))){
        LOG(ERROR) << "***Error reading pixel count threshold.***";
        return  EXIT_FAILURE;
    }

    PILClose(status);
    return (EXIT_SUCCESS);

}
void cztpixclean::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                           CZTPIXCLEAN PARAMETERS                     ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"Modulename            : "<<modulename;
    LOG(INFO)<<"Input Event file      : "<<infile;         //input event data file 
    LOG(INFO)<<"Input livetime file   : "<<inlivetimefile;
    LOG(INFO)<<"Output Event file     : "<<outfile1;        //output event file
    LOG(INFO)<<"Output livetime file  : "<<outlivetimefile;        //output event file
    LOG(INFO)<<"Output badpix file    : "<<badpixfile;
    if(clobber==YES)
    LOG(INFO)<<"Clobber               : YES";
    else
    LOG(INFO)<<"Clobber               : NO";
    if(history==YES)
    LOG(INFO)<<"History               : YES";
    else
    LOG(INFO)<<"History               : NO";    
	LOG(INFO)<<"---------------------------------------------------------------------------";
}

int cztpixclean::read(char *infile,char *inlivetimefile,int writedblevt,char *outfile1, char *outfile2,char *outlivetimefile, char *badpixfile, int nsigma,float det_tbinsize,float pix_tbinsize,int det_count_thresh,int pix_count_thresh)
{
    strcpy(this->infile,infile);
    strcpy(this->outfile1,outfile1);
    strcpy(this->outfile2,outfile2);
	//strcpy(this->caldb_badpix,caldb_badpix);
	strcpy(this->badpixfile,badpixfile);	
    strcpy(this->inlivetimefile,inlivetimefile);
    strcpy(this->outlivetimefile,outlivetimefile);
	this->nsigma=nsigma;
	this->det_tbinsize=det_tbinsize;
    this->pix_tbinsize=pix_tbinsize;
    this->det_count_thresh=det_count_thresh;
    this->pix_count_thresh=pix_count_thresh;
	this->writedblevt=writedblevt;

}



int cztpixclean::cztpixcleanProcess()
{
    int status=0,hdutype=0;
    int qid;
    int quadstart=0,quadend=3;
    int ntick_col,cztsec_col,detid_col,pixid_col,pi_col,detx_col,dety_col,UT_col;
    long nrows,i,j;
    int naxis1,k;
    double exposure_time[4];
    char gtitype[20];

    double **pixel_exposure;
    double **exposure_loss;
    int **allquad_pixflag;
    double *inlivetime, *outlivetime;
    double *timearray,livetime_binsize;
    long ntbins;
	double btstart,btstop;

    pixel_exposure=(double**)malloc(sizeof(double*)*4);
    exposure_loss=(double**)malloc(sizeof(double*)*4);
    allquad_pixflag=(int**)malloc(sizeof(int*)*4);

    for(i=0;i<4;i++)
    {
        pixel_exposure[i]=(double*)malloc(sizeof(double)*4096);
        exposure_loss[i]=(double*)malloc(sizeof(double)*4096);
        allquad_pixflag[i]=(int*)malloc(sizeof(int)*4096);
        for(j=0;j<4096;j++) 
        {   
            allquad_pixflag[i][j]=0;
            exposure_loss[i][j]=0.;
            pixel_exposure[i][j]=0.;
        }
    }

    fitsfile *fptr, *fout, *finlv;

    fits_open_file(&fptr, infile, READONLY, &status);
    if (status)
    {
    printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,infile);
    fits_report_error(stderr, status);
    return (EXIT_FAILURE);
    }

    //Get the gtitype
    
    fits_read_key(fptr,TSTRING,"GTITYPE",gtitype,NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	//Read the Pixel exposure

    fits_movnam_hdu(fptr, BINARY_TBL, "EXPOSURE", 0, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_open_file(&finlv, inlivetimefile, READONLY, &status);
    if (status)
    {
    printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,infile);
    fits_report_error(stderr, status);
    return (EXIT_FAILURE);
    }


	for(i=0;i<4;i++)
	{
	    fits_read_col(fptr, TDOUBLE, i+1, 1, 1, nrows, NULL, pixel_exposure[i],NULL, &status);
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
	}
	
    fits_read_key(finlv,TDOUBLE,"LV_BINSIZE",&livetime_binsize,NULL,&status);//added by Mayuri,16th Dec 2017
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }//added by Mayuri,16th Dec 2017
	
    //create livetime output file
    remove(outlivetimefile);
    status=create_livetime_file(outlivetimefile,livetime_binsize);//added by Mayuri,16th Dec 2017
    //status=create_livetime_file(outlivetimefile);//commented by Mayuri,16th Dec 2017

    if(status)
    {
        LOG(ERROR)<<"Unable to create output live time file";
        return(EXIT_FAILURE);
    }

    char gtiextnam[20];

    for(qid=quadstart;qid<=quadend;qid++)
    {
        LOG(INFO)<<"Quad "<<qid<<" processing started"; 

        if(strcasecmp(gtitype,"common")==0)
            sprintf(gtiextnam,"%s","GTI");

        else if (strcasecmp(gtitype,"quad")==0)
            sprintf(gtiextnam,"Q%d_GTI",qid);
        else
        {
            LOG(ERROR)<<"Unknown GTI type. exiting..";
            return EXIT_FAILURE;
        }

        //Read the live time for this quadrant
        

        try{
        fits_movabs_hdu(finlv, qid+2, &hdutype, &status);

        fits_get_num_rows(finlv, &ntbins, &status);
        }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }

        inlivetime=(double*)malloc(sizeof(double)*ntbins);
        outlivetime=(double*)malloc(sizeof(double)*ntbins);		
        timearray=(double*)malloc(sizeof(double)*ntbins);                           

        try{
        fits_read_col(finlv, TDOUBLE, 1, 1, 1, ntbins, NULL, timearray,NULL, &status);
        fits_read_col(finlv, TDOUBLE, 2, 1, 1, ntbins, NULL, inlivetime,NULL, &status);
        }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }
  

		for(i=0;i<ntbins;i++) outlivetime[i]=inlivetime[i];

	    double livetime_binsize;

		if(ntbins>=2) livetime_binsize=timearray[1]-timearray[0];
		else 
		{
			LOG(WARNING)<<"Unable to get live time binsize from livetime file. Setting to default value 1";
			livetime_binsize=1.0;
		}

		LOG(INFO)<<"TIME BIN SIZE FOR LIVETIME IS "<<livetime_binsize;


        fits_movnam_hdu(fptr, BINARY_TBL, gtiextnam, 0, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        status=updatelivetime_gti(fptr,timearray, outlivetime,ntbins,livetime_binsize);

        if(status)
        {
            LOG(ERROR)<<"Unable to update live time with GTI";
            return(EXIT_FAILURE);
        }


        // Move to quadrant data

	    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        //Get the exposure time from header
        fits_read_key(fptr,TDOUBLE,"EXPOSURE",&exposure_time[qid],NULL,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        //Multiply exposure time with fraction to get pixel exposure times    
        for(i=0;i<4096;i++)
            pixel_exposure[qid][i]*=exposure_time[qid];

        fits_get_num_rows(fptr, &nrows, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        if(nrows==0)
        {
            LOG(INFO)<<"************NO EVENTS IN QUADRANT "<<qid<<" **************";
            continue;
        }

		LOG(INFO)<<"Events in input file is: "<<nrows;
        
		fits_get_colnum(fptr,CASEINSEN,"TIME",&UT_col,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_get_colnum(fptr,CASEINSEN,"DETID",&detid_col,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_get_colnum(fptr,CASEINSEN,"PIXID",&pixid_col,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_get_colnum(fptr,CASEINSEN,"DETX",&detx_col,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_get_colnum(fptr,CASEINSEN,"DETY",&dety_col,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_get_colnum(fptr,CASEINSEN,"PI",&pi_col,&status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        double cztstart,cztstop;
        double *UTtime;
        int *detid,*pixid,*PI;
        char *evt_flag;

        UTtime=(double*)malloc(sizeof(double)*nrows);
        detid=(int*)malloc(sizeof(int)*nrows);
        pixid=(int*)malloc(sizeof(int)*nrows);
        PI=(int*)malloc(sizeof(int)*nrows);
        evt_flag=(char*)malloc(sizeof(char)*(nrows));

        fits_read_col(fptr, TDOUBLE, UT_col, 1, 1, nrows, NULL, UTtime,NULL, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_read_col(fptr, TINT, detid_col, 1, 1, nrows, NULL, detid,NULL, &status);
        fits_read_col(fptr, TINT, pixid_col, 1, 1, nrows, NULL, pixid,NULL, &status);
        fits_read_col(fptr, TINT, pi_col, 1, 1, nrows, NULL, PI,NULL, &status);

        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        long tstart=(long)(UTtime[0]);
        long tstop=(long)(UTtime[nrows-1])+1;

        float *detlc,*pixlc;

		// Identify single and double events and flag them
		for(i=0;i<nrows-1;i++) 
		{
			if(UTtime[i+1]-UTtime[i] > 2.5e-5)
				evt_flag[i]=1;
			else
			{	
				evt_flag[i]=2;   //is double event
				evt_flag[i+1]=2;
				i++;
			}
		}
		if(evt_flag[nrows-1]!=2) evt_flag[nrows-1]=1;


        detlc=(float*)malloc(sizeof(float)*16);
        pixlc=(float*)malloc(sizeof(float)*4096);

		// Get DPH and find noisy and dead pixels and flag them
		long dph[4096];
		int pixflag[4096];
		int lld[4096];

		//Read the CALDB badpix file
		int badpixExtnum,lldExtnum;	
    	double head_tstart,head_tstop;  
	    char quadnam[20];	
		string caldb_badpix;
		string caldb_lld;

    	fits_read_key(fptr, TDOUBLE, "TSTART", &head_tstart, NULL, &status);
    	report_error(status, (string) "Error in reading keyword TSTART");

    	fits_read_key(fptr, TDOUBLE, "TSTOP", &head_tstop, NULL, &status);
    	report_error(status, (string) "Error in reading keyword TSTOP");
		
    	sprintf(quadnam,"QUADRANT%d",qid);
		
		if(QueryCaldb("ASTROSAT","CZTI",(string)quadnam,"BADPIX",head_tstart,head_tstop,caldb_badpix,badpixExtnum))
    	{
            LOG(ERROR) << "Not able to get CALDB BADPIX file";
        	return (EXIT_FAILURE);
    	}
		else
			LOG(INFO)<<"CALDB badpix file "<<caldb_badpix<<" is used";

        if(QueryCaldb("ASTROSAT","CZTI",(string)quadnam,"LLD",head_tstart,head_tstop,caldb_lld,lldExtnum))
        {
            LOG(ERROR) << "Not able to get CALDB LLD file";
            return (EXIT_FAILURE);
        }
        else
            LOG(INFO)<<"CALDB LLD file "<<caldb_lld<<" is used";

		read_lldfile(caldb_lld,qid,lld);

		read_badpixfile(caldb_badpix,qid, pixflag);

		for(i=0;i<4096;i++) dph[i]=0;
		for(i=0;i<nrows;i++) {dph[detid[i]*256+pixid[i]]++;}
		
		// Find bad pixels by iterative nsigma clipping
		findbadpix(dph,pixflag,nsigma);

        int ngoodpix;

		//Set the pixel exposure to zero for bad pixels
		for(i=0;i<4096;i++)
		{
			if(pixflag[i]>=3)
				pixel_exposure[qid][i]=0.;
            else
                ngoodpix++;            
		}

		long nevt=0,nsec=0;
		long prev_det_tbin=-1,prev_pix_tbin=-1;
		long start_evt,stop_evt;
        long det_tbin,pix_tbin;

		for(i=0;i<nrows;i++) 
		{
			//  Flag noisy pixel events
			if(pixflag[detid[i]*256+pixid[i]]>=3) evt_flag[i]=0;

			//Flag events below threshold PI channel of each module
			//if(PI[i]<lld[detid[i]*256+pixid[i]]) evt_flag[i]=0;		

			//Check for noisy modules in time scale of det_tbinsize		
			if(evt_flag[i]==1||evt_flag[i]==2)
			{
				det_tbin=(int)((UTtime[i]-tstart)/(det_tbinsize));
				
				if(det_tbin==prev_det_tbin)
				{
					detlc[detid[i]]++;
				}
				else //The prev_det_tbin is completed
				{
					//Mark i-1 as stop event of prev_tbin
					stop_evt=i-1;

					//Flag events between start and stop of previous bin if needed
					if(prev_det_tbin>=0)
					{
						for(j=start_evt;j<=stop_evt;j++)
						{
						if(detlc[detid[j]] > det_count_thresh)
							evt_flag[j]=0;				
						}

						//Subtract the effective exposure time of the detector
						for(j=0;j<16;j++)
						{
							if(detlc[j]>det_count_thresh)
							{
								for(k=0;k<256;k++)
								{
									if(pixflag[j*256+k]<3)	
										exposure_loss[qid][j*256+k]+=det_tbinsize;
								}
							 
							 //LOG(INFO)<<"PIX DT "<<(double)tstart+(det_tbin+1)*det_tbinsize-(tstart+det_tbin*det_tbinsize);
							 btstart=(double)tstart+(double)det_tbin*det_tbinsize;
							 btstop=(double)tstart+(double)(det_tbin+1)*det_tbinsize;	 
							 re_updatelivetime(btstart,btstop,inlivetime,outlivetime,(long)(timearray[0]-livetime_binsize/2.0),1.0/16.0,0,ntbins,livetime_binsize);
							}
							
						}
                        

					}

					//reinit detlc and pixlc
					for(k=0;k<16;k++) detlc[k]=0;
							
					//Mark this as start event of this bin
					start_evt=i;

					//Add the current event in detlc and pixlc
					detlc[detid[i]]++;	
				}

				//Update the previous tbin
				prev_det_tbin=det_tbin;
			
			}
		} //End of Loop over rows of evt file for detector threshold


		// Now at pixel level
        for(i=0;i<nrows;i++)
        {

            //Check for noisy pixels in time scale of pix_tbinsize 
            if(evt_flag[i]==1||evt_flag[i]==2)
            {
                pix_tbin=(int)((UTtime[i]-tstart)/(pix_tbinsize));

                if(pix_tbin==prev_pix_tbin)
                {
					pixlc[detid[i]*256+pixid[i]]++;
                }
                else //The prev_det_tbin is completed
                {
                    //Mark i-1 as stop event of prev_tbin
                    stop_evt=i-1;

                    //Flag events between start and stop of previous bin if needed
                    if(prev_pix_tbin>=0)
                    {
                        for(j=start_evt;j<=stop_evt;j++)
                        {
                        if(pixlc[detid[j]*256+pixid[j]] > pix_count_thresh)
                            evt_flag[j]=0;
                        }

                        //Subtract the effective exposure time of the pixel
                        for(j=0;j<4096;j++)
                        {
                            if(pixlc[j]>pix_count_thresh)
			    {
                                    exposure_loss[qid][j]+=pix_tbinsize;
				    
									btstart=(double)tstart+(double)pix_tbin*pix_tbinsize;
									btstop=(double)tstart+(double)(pix_tbin+1)*pix_tbinsize;	
									re_updatelivetime(btstart,btstop,inlivetime,outlivetime,(long)(timearray[0]-livetime_binsize/2.0),1.0/(float)(ngoodpix),0,ntbins,livetime_binsize);
			    }
                        }

						
                    }

                    //reinit detlc and pixlc
                    for(k=0;k<4096;k++) pixlc[k]=0;

                    //Mark this as start event of this bin
                    start_evt=i;

                    //Add the current event in pixlc
                	pixlc[detid[i]*256+pixid[i]]++;
				}
                //Update the previous tbin
                prev_pix_tbin=pix_tbin;
            }

        } //End of Loop over rows of evt file

	
		// Make sure that for double events companion of bad event is also flagged bad
		for(i=1;i<nrows-1;i++)
		{
			if(evt_flag[i]==2)
			{
				if((evt_flag[i+1]!=2&&UTtime[i+1]-UTtime[i]<3.0e-5)||(evt_flag[i-1]!=2&&UTtime[i]-UTtime[i-1]<3.0e-5))
					evt_flag[i]=0;
			}
		}
		
		if((evt_flag[0]==2)&&(evt_flag[1]!=2)) evt_flag[0]=0;
		if(evt_flag[nrows-1]==2&&evt_flag[nrows-2]!=2) evt_flag[nrows-1]=0;


		//Count good events
		
		long good=0;
		for(i=0;i<nrows;i++) {if(evt_flag[i]==1) good++;}
        
        LOG(INFO)<<"Good events :"<<good;

        good=0;
        for(i=0;i<nrows;i++) {if(evt_flag[i]==2) good++;}
        LOG(INFO)<<"Double events :"<<good/2;

		//Write the badpix flag to allquad array
		for(j=0;j<4096;j++) allquad_pixflag[qid][j]=pixflag[j];

        //Compute the new fractional exposures of pixels
        for(j=0;j<4096;j++)
            pixel_exposure[qid][j]=(pixel_exposure[qid][j]-exposure_loss[qid][j])/(exposure_time[qid]);

		//UPDATE LIVETIME WITH GTI
/*
        fits_movnam_hdu(fptr, BINARY_TBL, gtiextnam, 0, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        status=updatelivetime_gti(fptr,timearray, outlivetime,ntbins,livetime_binsize);

        if(status)
        {
            LOG(ERROR)<<"Unable to update live time with GTI";
            return(EXIT_FAILURE);
        }
*/

        //Write the new livetime file
        writelivetime(timearray,outlivetime,ntbins,qid,outlivetimefile);


		//Create new evt file with good events only
		remove_events(infile,writedblevt,outfile1,outfile2,evt_flag,qid,pixel_exposure);

		LOG(INFO)<<"Done writing event file for this quadrant";

		free(UTtime);
		free(detid);
		free(pixid);
		free(PI);
		free(evt_flag);	
		free(detlc);
		free(pixlc);
		free(inlivetime);	
       	free(outlivetime);
		free(timearray);
        
		LOG(INFO)<<"Memory freed for this quadrant";
		}

	//Write all quadrant bad pixel file
	status= write_badpix(badpixfile,allquad_pixflag);
	if (status){return (EXIT_FAILURE);}

    free(exposure_loss);
    free(pixel_exposure);

	//UPDATING KEYWORDS
    updateKeywords(outfile1, modulename);
	return(EXIT_SUCCESS);
}

int remove_events(char* infile,int writedblevt,char* outfile1,char *outfile2,char *evt_flag,int qid, double **pixel_exposure)
{

    int status=0,hdutype=0;
    int j;
    int quadstart=0,quadend=0;
    int ntick_col,cztsec_col,detid_col,pixid_col,pi_col;
    long nrows,i;
    int naxis1;

    fitsfile *fptr, *fout1,*fout2;

    fits_open_file(&fptr, infile, READONLY, &status);
    if (status)
    {
    printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,infile);
    fits_report_error(stderr, status);
    return (EXIT_FAILURE);
    }


	// Only if qid is zero create file and copy primary header
	if(qid==0)
	{
    remove(outfile1);
    fits_create_file(&fout1,outfile1,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    remove(outfile2);
	if(writedblevt==1){
    fits_create_file(&fout2,outfile2,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
	}
    fits_copy_header(fptr, fout1,  &status);
	if(writedblevt==1){
	fits_copy_header(fptr, fout2,  &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
	}
	}
	else
	{
		fits_open_file(&fout1, outfile1, READWRITE, &status);
		 if(writedblevt==1){	
                	fits_open_file(&fout2, outfile2, READWRITE, &status);
	         		if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
		}
	}	

    LOG(INFO)<<"Quad "<<qid<<" started writing output files";

    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_copy_header(fptr, fout1,  &status);
	if(writedblevt==1){
    fits_copy_header(fptr, fout2,  &status);
	}
    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_read_key(fptr,TINT,"NAXIS1",&naxis1,NULL, &status);

    unsigned char *tblrow;
    tblrow=(unsigned char*)malloc(sizeof(unsigned char)*naxis1);

	long nevnt1=0,nevnt2=0;	
    for(i=0;i<nrows;i++)
    {
		if(evt_flag[i]==1)
		{
            fits_read_tblbytes(fptr, i+1, 1,naxis1, tblrow, &status);
			if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
            nevnt1+=1;
            fits_write_tblbytes (fout1, nevnt1, 1,naxis1, tblrow, &status);			
			if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }	
		}	
	    else if(evt_flag[i]==2)
        {
            fits_read_tblbytes(fptr, i+1, 1,naxis1, tblrow, &status);
            if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
            nevnt2+=1;
			if(writedblevt==1){	
            fits_write_tblbytes (fout2, nevnt2, 1,naxis1, tblrow, &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }   
	   		}		
        }  

	}

    fits_update_key(fout1, TLONG, "NAXIS2",&nevnt1, NULL, &status);
	if(writedblevt==1){	
    fits_update_key(fout2, TLONG, "NAXIS2",&nevnt2, NULL, &status);
	}

//    printf("Total events old %ld new single %ld double %ld\n",nrows,nevnt1,nevnt2);


    if(qid==3)
    {

        //Copy the remaining extensions and write the exposure

        int num_hdus;

        fits_get_num_hdus(fptr, &num_hdus,&status);

        for(i=6;i<=num_hdus;i++)
        {
        fits_movabs_hdu(fptr, i, &hdutype, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_copy_hdu(fptr, fout1, 0,&status);
        if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
		if(writedblevt==1){
        fits_copy_hdu(fptr, fout2, 0,&status);
        if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
        }
		}

	    // Modify exposure for single events
		fits_movnam_hdu(fout1, BINARY_TBL, "EXPOSURE", 0, &status);
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    	fits_get_num_rows(fout1, &nrows, &status);
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
		
        for(i=0;i<4;i++)
        {
            fits_write_col(fout1, TDOUBLE, i+1, 1,1, 4096, pixel_exposure[i], &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
        }

        // Modify exposure for double events
		if(writedblevt==1){	 
        fits_movnam_hdu(fout2, BINARY_TBL, "EXPOSURE", 0, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        fits_get_num_rows(fout2, &nrows, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        for(i=0;i<4;i++)
        {
            fits_write_col(fout2, TDOUBLE, i+1, 1,1, 4096, pixel_exposure[i], &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
        }


        LOG(INFO)<<"Exposure extension modified";
		}
	}
fits_close_file(fptr,&status);
if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

//close output file

fits_close_file(fout1,&status);
if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
if(writedblevt==1){
fits_close_file(fout2,&status);
if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
}
}

int findbadpix(long *dph,int *pixflag,int nsig)
{
	int i;
	double avg=100,prev_avg=0,stdev;

	for(i=0;i<4096;i++) 
	{
		if(dph[i]==0)
			pixflag[i]=4;
	}

	while (avg!=prev_avg)
	{
		avg=compute_mean(dph,pixflag);
		prev_avg=avg; 
		stdev=compute_stdev(dph,pixflag,avg);
		
		for(i=0;i<4096;i++)
		{
			if(pixflag[i]<3&&dph[i]>avg+nsig*stdev)
				pixflag[i]=3;
		}
		avg=compute_mean(dph,pixflag);
	}

}

double compute_mean(long *dph,int *pixflag)
{
	int i,npix=0;
	double avg=0.;

    for(i=0;i<4096;i++)
    {
        if(pixflag[i]<3) 
        {
            avg+=dph[i];
            npix++;
        }
    }
	
	avg/=(double)npix;
	return(avg);
}

double compute_stdev(long *dph,int *pixflag,double avg)
{
    int i,npix=0;
    double stdev=0.;

    for(i=0;i<4096;i++)
    {   
        if(pixflag[i]<3) 
        {   
            stdev+=(dph[i]-avg)*(dph[i]-avg);
            npix++;
        }   
    }   
    
	stdev/=npix;
	stdev=sqrt(stdev);

    return(stdev);
}

int read_lldfile(string caldb_lld,int qid,int *lld)
{
    fitsfile *fptr;
    int status=0,hdutype=0;
    int lld_col,i;
    long nrows;

    fits_open_file(&fptr, (char *)caldb_lld.c_str(), READONLY, &status);
    if (status)
    {
    printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,(char*)caldb_lld.c_str());
    fits_report_error(stderr, status);
    return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_colnum(fptr,CASEINSEN,"LLD",&lld_col,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	for(i=0;i<4096;i++)
	{

    	fits_read_col(fptr, TINT, lld_col, i+1, 1, 1, NULL, &lld[i],NULL, &status);
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	}

    fits_close_file(fptr,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    return(EXIT_SUCCESS);

}


int read_badpixfile(string caldb_badpix,int qid, int *pixflag)
{
	fitsfile *fptr;
	 int status=0,hdutype=0;
	int pixflag_col;
	long nrows;

    fits_open_file(&fptr, (char *)caldb_badpix.c_str(), READONLY, &status);
    if (status)
    {
    printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,(char*)caldb_badpix.c_str());
    fits_report_error(stderr, status);
    return (EXIT_FAILURE);
    }

   	fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
	
   	fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_colnum(fptr,CASEINSEN,"PIX_FLAG",&pixflag_col,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	fits_read_col(fptr, TINT, pixflag_col, 1, 1, nrows, NULL, pixflag,NULL, &status);	
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	fits_close_file(fptr,&status);
	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	return(EXIT_SUCCESS);
}

int write_badpix(char *fname,int **flg)
{
    int status=0,i;
    fitsfile *fout;

    int det,pix,row;

	remove(fname);

    // Create file
    fits_create_file(&fout,fname,&status);
    if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    int bitpix=16;
    int naxis=0;

    status = create_primaryfits(fout);
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }

    int quad;

    for (quad=0;quad<4;quad++)
    {
    //Quadrant number as string
    char string_quad[15];
    sprintf(string_quad,"%d",quad);
    char detname[10]="QUADRANT";
    strcat(detname,string_quad);
    char extname[10]="Q";
    strcat(extname,string_quad);


    int tfields=3;

    char *ttype[] = {"DETID", "PIXID", "PIX_FLAG"};
    char *tform[] = {"B","B","B"};

    // Create binary table extension
    fits_create_tbl(fout, BINARY_TBL, 0,  tfields,ttype,tform,NULL, extname, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    // Write keywords to header

    char as[1024]="ASTROSAT";
	char czt[1024]="CZTI";
	char ogip[1024]="OGIP";
	char ogip1[1024]="OGIP 1.0";

      // Write keywords to header
    fits_write_key(fout,TSTRING,"TELESCOP",as,"Name of mission/satellite",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
    fits_write_key(fout,TSTRING,"INSTRUME",czt,"Name of Instrument/detector",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    if(detname!=NULL)
    {
    fits_write_key(fout,TSTRING,"DETNAM",detname,"Quadrant ID",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
    }

    fits_write_key(fout,TSTRING,"HDUCLASS",ogip,"OGIP Standard",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
    fits_write_key(fout,TSTRING,"LONGSTRN",ogip1,"The OGIP Long String Convention may be used",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

/*    fits_write_key(fout,TSTRING,"ORIGIN","CZTI POC","Origin of FITS file",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"VERSION","1","Version number",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"CONTENT","CZTI Badpixel list","File content",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
    fits_write_key(fout,TSTRING,"COMMENT",NULL,NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"COMMENT","Pixel Quality flags",NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"COMMENT","Good-0; Spectrscopically bad-1",NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"COMMENT","Flickering-2; Noisy-3; Dead/Inactive-4",NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"COMMENT",NULL,NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

*/

    // Write the data to FITS File  
    for (det=0;det<16;det++)
    {
        for (pix=0;pix<256;pix++)
        {
            row=det*256+pix;

            // Write detid and pixid
            fits_write_col(fout, TBYTE, 1, row+1,1, 1, &det, &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

            fits_write_col(fout, TBYTE, 2, row+1,1, 1, &pix, &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

            fits_write_col(fout, TBYTE, 3, row+1,1, 1, &flg[quad][row], &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        }
    }


    fits_write_chksum(fout, &status);
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }

    }

    //Close the fits file
    fits_close_file(fout,&status);
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }


    return 0;
}

int create_primaryfits(fitsfile *fout)
{
    int status=0;

	char as[1024]="ASTROSAT";
	char czt[1024]="CZTI";
    int bitpix=16;
    int naxis=0;

    // Create Primary extension  
    fits_create_img(fout,bitpix,naxis,NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    // Add necessary keyword in the header

    fits_write_key(fout,TSTRING,"TELESCOP",as,"Name of mission/satellite",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TSTRING,"INSTRUME",czt,"Name of Instrument/detector",&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_chksum(fout, &status);
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }

    return 0;
}


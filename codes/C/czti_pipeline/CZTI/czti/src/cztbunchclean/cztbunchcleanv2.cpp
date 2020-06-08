/* cztbunchclean.c

This routine identifies the bunches in CZTI data and generates 
event file without bunch events. It also removes events just after the 
bunch based on certain parameters

Mithun N P S 
24/11/15

Edits:

** 27/11/15
* Fractional exposure calculation for pixels

*Write exposure as extension of fits file
*Changed to cpp and included PIL

06/12/15 
*Corrected exposure fraction for cztclock, and changed 
to fractional instead of total exposure time

12/12/15
* Changed all calculations to UT instead of 
 czttime. Only for removing events from first 
 400ms of 100th sec czttime is used. 

* Now computes live time of each second bin during the observation 
  period and writes of livetime file

18/02/16

*Changes made to process the event data after CZTI on-board software
 patch done on 17/02/16

01/04/16

* Memory leaks fixed

08/06/20

* Now writes an output bunchfile and an optional output bunch event file: 
  For data with on-board bunchclean, the output bunch file has a copy of the 
  list bunches from the input file barring the bunches that fall within first 
  400ms of 100s tick of CZT clock. For the data without on-board bunch clean, 
  information of the bunches identified are written out in the output bunchfile. 
  The output buncheventfile will have the events that are part of the bunch 
  (except those within 400ms of 100s cztclock tick), removed from the output 
  event file. 


*/

#include "cztbunchcleanv2.h"

//Constructor
cztbunchclean::cztbunchclean(){
    strcpy(modulename, "cztbunchclean_v");
	strcat(modulename,VERSION);
}

//Destructor
cztbunchclean::~cztbunchclean(){}

int cztbunchclean::read(int argc,char **argv)
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

    if(PIL_OK!=(status=PILGetFname("par_bunchfile",bunchfile))){
        LOG(ERROR)<<"***Error reading input bunch file"<<bunchfile<<"***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("par_outfile",outfile))){
        LOG(ERROR)<<"***Error reading output event file"<<outfile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("par_livetimefile",livetimefile))){
        LOG(ERROR)<<"***Error reading output livetime file"<<livetimefile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("par_outbunchfile",outbunchfile))){
        LOG(ERROR)<<"***Error reading output bunch file"<<outbunchfile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetString("par_outbunchevtfile",outbunchevtfile))){
        LOG(ERROR)<<"***Error reading output bunch event file"<<outbunchevtfile<<"***";        
        return status;
    }

    if(strcasecmp(outbunchevtfile, "-") == 0)
        outbunchevtstatus=0;
    else
        outbunchevtstatus=1;
    
    if (PIL_OK != (status=PILGetInt("par_bunchdeftime", &bunchdeftime))){
        LOG(ERROR) << "***Error reading Bunch definition time.***";
        return  EXIT_FAILURE;
    }

    if (PIL_OK != (status=PILGetInt("par_bunch_length_thresh", &bunch_length_thresh))){
        LOG(ERROR) << "***Error reading Bunch length threshold.***";
        return  EXIT_FAILURE;
    }

    if (PIL_OK != (status = PILGetReal4("par_skipT1", &skipT1))) {
        LOG(ERROR) << "***Error reading SkipT1***";
        return status;
    }

    if (PIL_OK != (status = PILGetReal4("par_skipT2", &skipT2))) {
        LOG(ERROR) << "***Error reading SkipT2***";
        return status;
    }

    if (PIL_OK != (status = PILGetReal4("par_skipT3", &skipT3))) {
        LOG(ERROR) << "***Error reading SkipT3***";
        return status;
    }

    if (PIL_OK != (status = PILGetReal("par_livetime_binsize", &livetime_binsize))) {
        LOG(ERROR) << "***Error reading livetime_binsize***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("clobber", &clobber))) {
        LOG(ERROR) << "***Error Reading clobber:" << clobber << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("history", &history))) {
        LOG(ERROR) << "***Error Reading history parameter" << history << "***";
        return status;
    }


	PILClose(status);
    return (EXIT_SUCCESS);
}


int cztbunchclean::read(char *infile,char *bunchfile,char *outfile,char * livetimefile,int bunchdeftime,int bunch_length_thresh,float skipT1 ,float skipT2,float skipT3,double livetime_binsize, int clobber, int history){
    strcpy(this->infile,infile);
    strcpy(this->bunchfile,bunchfile);
	strcpy(this->outfile,outfile);
    strcpy(this->livetimefile,livetimefile);
	this->bunchdeftime=bunchdeftime;
	this->bunch_length_thresh=bunch_length_thresh;
    this->skipT1=skipT1;
    this->skipT2=skipT2;
    this->skipT3=skipT3;
	this->livetime_binsize=livetime_binsize;
    this->clobber = clobber;
    this->history = history;
    return (EXIT_SUCCESS);			
}

void cztbunchclean::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                          CZTBUNCHCLEAN PARAMETERS                            ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
	LOG(INFO)<<  "Modulename                    : "<<modulename;
    LOG(INFO) << "Input Event file              : " << infile;
    LOG(INFO) << "Bunch file                    : " << bunchfile;
    LOG(INFO) << "Output event file             : " << outfile;
    LOG(INFO) << "Output livetime file          : " << livetimefile;
    LOG(INFO) << "Bunch definition time (us)    : " << bunchdeftime;
    LOG(INFO) << "Bunch length threshold        : " << bunch_length_thresh;
    LOG(INFO) << "Time parameter SkipT1 (sec)   : " <<skipT1;
    LOG(INFO) << "Time parameter SkipT2 (sec)   : " <<skipT2;
    LOG(INFO) << "Time parameter SkipT3 (sec)   : " <<skipT3;
    LOG(INFO) << "Livetime binsize              : " << livetime_binsize;
    if(clobber==YES)
    LOG(INFO) << "Clobber                       : YES";
    else
    LOG(INFO) << "Clobber                       : NO";
    if(history==YES)
    LOG(INFO) << "History                       : YES";
    else
    LOG(INFO) << "History                       : NO";

    LOG(INFO)<<"----------------------------------------------------------------------------";

}	

int cztbunchclean::cztbunchcleanProcess()
{
    int status=0,hdutype=0;
    int qid;
    int quadstart=0,quadend=3;
    int cztsec_col,detid_col,time_col,rownum_col,numevent_col,det1_col,det2_col,det3_col,det4_col;
    int pixid_col,DetId_fevt_col,PixId_fevt_col,DetId_sevt_col,PixId_sevt_col,DetId_tevt_col,PixId_tevt_col,Time_dsl_col,Time_dfs_col,bunch_time_col;
    long nrows,i,j,k,l,nrows_gti,nrows_bunch;
    int ncols_bunch;
    double *livetime, *timearray;
    long tstarti,tstopi,num_timebins;

    double clock_conversion=0.999998125;//(UT/CZT) 

    fitsfile *fptr,*fgti,*fbunch,*foutbunch;

    //FILE EXISTENCE CHECK AND UNLINKING IF IT DOES
    if (clobber == YES) {
        if (FileExists(outfile)) {
            LOG(INFO) << outfile << "  :FileExists.. Replacing the old file";
            if (unlink(outfile) != 0) {
                LOG(ERROR) << "***Error in deleting " << outfile << "***";
                return (EXIT_FAILURE);
            } 
        }
        if (FileExists(livetimefile)) {
            LOG(INFO) << livetimefile << "  :FileExists.. Replacing the old file";
            if (unlink(livetimefile) != 0) {
                LOG(ERROR) << "***Error in deleting " << livetimefile << "***";
                return (EXIT_FAILURE);
            }
        }
        if (FileExists(outbunchfile)) {
            LOG(INFO) <<  outbunchfile << "  :FileExists.. Replacing the old file";
            if (unlink( outbunchfile) != 0) {
                LOG(ERROR) << "***Error in deleting " <<  outbunchfile << "***";
                return (EXIT_FAILURE);
            }
        }

        if (outbunchevtstatus && FileExists(outbunchevtfile)) {
            LOG(INFO) <<  outbunchevtfile << "  :FileExists.. Replacing the old file";
            if (unlink( outbunchevtfile) != 0) {
                LOG(ERROR) << "***Error in deleting " <<  outbunchevtfile << "***";
                return (EXIT_FAILURE);
            }
        }

    }
    else {
        if (FileExists(outfile)||FileExists(livetimefile)||FileExists(outbunchfile)||(outbunchevtstatus && FileExists(outbunchevtfile))) {
            LOG(ERROR) << "***Output file already exists***";
            LOG(ERROR) << "Use clobber=yes for overwriting the file";
            return (EXIT_FAILURE);
        }
    }

    fits_open_file(&fptr, infile, READONLY, &status);
    if (status)
    {
    	printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,infile);
    	fits_report_error(stderr, status);
    	return (EXIT_FAILURE);
    }

    fits_open_file(&fbunch, bunchfile, READONLY, &status);
    if (status)
    {
    	printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,bunchfile);
    	fits_report_error(stderr, status);
    	return (EXIT_FAILURE);
    }


    fits_create_file(&foutbunch,outbunchfile,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_copy_header(fbunch, foutbunch,  &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }


//    status=create_livetime_file(livetimefile);//commented by Mayuri,16th Dec 2017
    status=create_livetime_file(livetimefile,livetime_binsize);//added by Mayuri,16th Dec 2017

    if(status){
    LOG(ERROR) <<"Error in creating livetime file";
    return (EXIT_FAILURE);
    }

	//Read the Quad GTI extensions and compute the total exposure time for each quadrant

    double **pixel_exposure;
    pixel_exposure=(double**)malloc(sizeof(double*)*4);
    for(i=0;i<4;i++) pixel_exposure[i]=(double*)malloc(sizeof(double)*4096);
    for(qid=0;qid<4;qid++){for(j=0;j<4096;j++) pixel_exposure[qid][j]=0;}
    double allquad_exp_time[4];

    bunchfileAddcolStatus=0; 

    for(qid=quadstart;qid<=quadend;qid++)
    {
        LOG(INFO)<<"Processing Quadrant"<<qid;

        char gtihdunam[100];

        sprintf(gtihdunam,"Q%d_GTI",qid);

        fits_movnam_hdu(fptr, BINARY_TBL, gtihdunam, 0, &status);
        if (status) {
            LOG(ERROR) <<"Error in moving to hdu "<<gtihdunam;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
   
        fits_get_num_rows(fptr, &nrows_gti, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

        double gti_start[nrows_gti],gti_stop[nrows_gti];
        allquad_exp_time[qid]=0;

        try{
        fits_read_col(fptr, TDOUBLE, 1, 1, 1, nrows_gti, NULL, gti_start,NULL, &status);

        fits_read_col(fptr, TDOUBLE, 2, 1, 1, nrows_gti, NULL, gti_stop,NULL, &status);
        } catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }

        
        // Go to quadrant data hdu

        try{

        fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
        fits_get_num_rows(fptr, &nrows, &status);

        fits_movabs_hdu(fbunch, qid+2, &hdutype, &status);
        fits_get_num_rows(fbunch, &nrows_bunch, &status);
        fits_get_num_cols(fbunch, &ncols_bunch, &status);        

        // copy hdu of bunchfile
        fits_copy_header(fbunch, foutbunch,  &status);

        } catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }


        if(nrows==0)
        {
            LOG(WARNING)<<"No events in quadrant "<<qid;
            continue;
        }
        
        if(nrows_bunch==0)
            LOG(INFO)<<"The event file doesn't have on-board bunch clean";
        else
            LOG(INFO)<<"*****This file is with on-board bunch clean*****";
        

        try{
        fits_get_colnum(fptr,CASEINSEN,"TIME",&time_col,&status);
        fits_get_colnum(fptr,CASEINSEN,"CZTSECCNT",&cztsec_col,&status);
        fits_get_colnum(fptr,CASEINSEN,"DETID",&detid_col,&status);
        fits_get_colnum(fptr,CASEINSEN,"PIXID",&pixid_col,&status);


        fits_get_colnum(fbunch,CASEINSEN,"Time",&bunch_time_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"Time_dfs",&Time_dfs_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"Time_dsl",&Time_dsl_col,&status);

        fits_get_colnum(fbunch,CASEINSEN,"ROWNUM",&rownum_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"NUMEVENT",&numevent_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"detid1",&det1_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"detid2",&det2_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"detid3",&det3_col,&status);
        fits_get_colnum(fbunch,CASEINSEN,"detid4",&det4_col,&status);

        } catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }


        if (ncols_bunch==15){
        
            fits_get_colnum(fbunch,CASEINSEN,"DetId_fevt",&DetId_fevt_col,&status);
            fits_get_colnum(fbunch,CASEINSEN,"DetId_sevt",&DetId_sevt_col,&status);
            fits_get_colnum(fbunch,CASEINSEN,"DetId_levt",&DetId_tevt_col,&status);
            fits_get_colnum(fbunch,CASEINSEN,"PixId_fevt",&PixId_fevt_col,&status);
            fits_get_colnum(fbunch,CASEINSEN,"PixId_sevt",&PixId_sevt_col,&status);
            fits_get_colnum(fbunch,CASEINSEN,"PixId_levt",&PixId_tevt_col,&status);
        
            bunchfileAddcolStatus=1;
        }
         
        num_timebins=((long)(((long)(gti_stop[nrows_gti-1])-(long)(gti_start[0]))/livetime_binsize)+1);

        tstarti=(long)(gti_start[0]);

        livetime=(double*)malloc(sizeof(double)*num_timebins);
        timearray=(double*)malloc(sizeof(double)*num_timebins);

        for(j=0;j<num_timebins;j++) 
        {
            timearray[j]=tstarti+j*livetime_binsize+0.5*livetime_binsize;
            livetime[j]=livetime_binsize;
        }

        double *UTtime,*czttime;
        int *detid, *pixid;
        int *evt_index;
        int *numevent_bunch;
        char *evt_flag;
        long *bunch_row;
        char **bunch_det;
        int *bunch_flag;

        UTtime=(double*)malloc(sizeof(double)*nrows);
        czttime=(double*)malloc(sizeof(double)*nrows);
        detid=(int*)malloc(sizeof(int)*nrows);
        pixid=(int*)malloc(sizeof(int)*nrows);        
        evt_index=(int*)malloc(sizeof(int)*nrows);
        evt_flag=(char*)malloc(sizeof(char)*nrows);

        for(i=0;i<nrows;i++)
        {
            evt_flag[i]=0;
            evt_index[i]=0;
        }

        if(nrows_bunch!=0) 
        {
            bunch_row=(long*)malloc(sizeof(long)*nrows_bunch);
            numevent_bunch=(int*)malloc(sizeof(int)*nrows_bunch);
            bunch_det=(char**)malloc(sizeof(char*)*4);
            for(j=0;j<4;j++) bunch_det[j]=(char*)malloc(sizeof(char)*nrows_bunch);
            bunch_flag=(int*)malloc(sizeof(int)*nrows_bunch);
            for (i=0;i<nrows_bunch;i++) bunch_flag[i]=0;
        }


	    //vectors to store bunch information if nrows_bunch==0
    	vector <double> bunch_time;
    	vector <long>evt_row_num;
    	vector <unsigned char> time_dfs;//Time difference between first event and second event in a bunch
    	vector <unsigned char>time_dsl;//Time difference between second event and last event in a bunch
    	vector <unsigned char>num_bunchevents;//Number of events in a bunch minus one.
    	vector <unsigned char>detid1;//Detector ID of event 3
    	vector <unsigned char>detid2;//Detector ID of event 4
    	vector <unsigned char>detid3;//Detector ID of event 5
    	vector <unsigned char>detid4;//Detector ID of event 6
    	vector <unsigned char>DetId_fevt;//Detector ID of first event
    	vector <unsigned char>PixId_fevt;//Pixel ID for first event
    	vector <unsigned char>DetId_sevt;//Detector ID of second event
    	vector <unsigned char>PixId_sevt;//Pixel ID for second event
    	vector <unsigned char>DetId_tevt;//Detector ID of last event
    	vector <unsigned char>PixId_tevt;//Pixel ID for last event


        for(j=0;j<nrows_gti;j++)
            allquad_exp_time[qid]+=gti_stop[j]-gti_start[j];

        try{

            fits_read_col(fptr, TDOUBLE, time_col, 1, 1, nrows, NULL, UTtime,NULL, &status);
            fits_read_col(fptr, TDOUBLE, cztsec_col, 1, 1, nrows, NULL, czttime,NULL, &status);
            fits_read_col(fptr, TINT, detid_col, 1, 1, nrows, NULL, detid,NULL, &status);
            fits_read_col(fptr, TINT, pixid_col, 1, 1, nrows, NULL, pixid,NULL, &status);
        
            //Read the bunch file
            if(nrows_bunch!=0) 
            {  
                fits_read_col(fbunch, TLONG, rownum_col, 1, 1, nrows_bunch, NULL, bunch_row,NULL, &status);
                fits_read_col(fbunch, TINT, numevent_col, 1, 1, nrows_bunch, NULL,numevent_bunch,NULL, &status);            
                fits_read_col(fbunch, TBYTE, det1_col, 1, 1, nrows_bunch, NULL,bunch_det[0],NULL, &status);
                fits_read_col(fbunch, TBYTE, det2_col, 1, 1, nrows_bunch, NULL,bunch_det[1],NULL, &status);
                fits_read_col(fbunch, TBYTE, det3_col, 1, 1, nrows_bunch, NULL,bunch_det[2],NULL, &status);
                fits_read_col(fbunch, TBYTE, det4_col, 1, 1, nrows_bunch, NULL,bunch_det[3],NULL, &status);
            }
        } catch(ErrorHandler errHandler){
            logError(errHandler);
            return EXIT_FAILURE;
        }

        /*
        FILE *f1;
        char lvfile[1024];
        sprintf(lvfile,"livetime_Q%d.txt",qid);
        f1=fopen(lvfile,"w");

        for(i=0;i<num_timebins;i++)
            fprintf(f1,"%lf\n",livetime[i]);

        fclose(f1);
        */  

        int delta;
        int nbunch_evt=1;
        double bunch_tstart,bunch_tstop;
        unsigned int bunch_time_dfs,bunch_time_dsl;
        int tmod100;

        LOG(INFO)<<"Number of events in quadrant "<<qid<<" : "<<nrows;

		// Identify the bunches and assign evt_index 
        double tot_loss=0.; 
        int row_index=0;

        for(i=0;i<nrows-1;i++)
        {
            tmod100=(long)(czttime[i])%100;
			//printf("%ld\t %d\t %lf\n",i,tmod100,UTtime[i]);
            
            if(!((tmod100==0)&&((czttime[i]-(long)czttime[i])<=0.400))) // Ignore first 400ms of every 100s
            {

            // If delta is too large (in case of data gap) set it to 1e7 
            if((UTtime[i+1]-UTtime[i])<10)
                delta=(int)((UTtime[i+1]-UTtime[i])*1e6);
            else
                delta=1e7;

            if(nrows_bunch==0)
            {
                if(delta <= bunchdeftime && i!=nrows-2)// To stop at event before the last 
                {
                    nbunch_evt++;
                    
                    if(nbunch_evt==2){
                        bunch_tstart=UTtime[i];
                        bunch_time_dfs=(unsigned int)(round(delta/20.0));
                    }

                    // when we realize that it is actually a bunch, push the info to the vectors for bunchfile
                    
                    if(nbunch_evt==3) {
                        bunch_time.push_back(UTtime[i]);    // Time of second event of bunch     
                        evt_row_num.push_back(i-1);         // event_row at start of bunch
                        time_dfs.push_back((unsigned char)bunch_time_dfs);
                        DetId_fevt.push_back(detid[i-1]);
                        PixId_fevt.push_back(pixid[i-1]);                        
                        DetId_sevt.push_back(detid[i]);
                        PixId_sevt.push_back(pixid[i]);                        
                    }
                    if(nbunch_evt==4) detid1.push_back(detid[i]);
                    if(nbunch_evt==5) detid2.push_back(detid[i]);
                    if(nbunch_evt==6) detid3.push_back(detid[i]);
                    if(nbunch_evt==7) detid4.push_back(detid[i]);

                }
                else
                {
                    bunch_tstop=UTtime[i];

                    // if it is actually a bunch, push the info to the vectors for bunchfile
                    if(nbunch_evt>=3) {

                        bunch_time_dsl=(unsigned int)(round((bunch_tstop-bunch_tstart)*1e6/20.0))-bunch_time_dfs;
                        if(bunch_time_dsl < 256) time_dsl.push_back((unsigned char)bunch_time_dsl);
                        else time_dsl.push_back(255);
                        
                        if(nbunch_evt<256) num_bunchevents.push_back(nbunch_evt-1); // No of events-1
                        else num_bunchevents.push_back(255);     
                        
                        DetId_tevt.push_back(detid[i]);
                        PixId_tevt.push_back(pixid[i]);
                        if(nbunch_evt<7) detid4.push_back(99);
                        if(nbunch_evt<6) detid3.push_back(99);
                        if(nbunch_evt<5) detid2.push_back(99);
                        if(nbunch_evt<4) detid1.push_back(99);
                    }
        
				    for(j=i;j>i-nbunch_evt;j--)
                        evt_index[j]=nbunch_evt;
                
				    nbunch_evt=1;

                }
            }
            else // If bunch file has entries
            {
                if(row_index<nrows_bunch && i==bunch_row[row_index])
                {
                    nbunch_evt=numevent_bunch[row_index]+1;
                    //printf("%d\t%ld\t%d\t%d\n",qid,i,nbunch_evt,row_index);
                    row_index++;
                    for(j=i;j<=i+2;j++)
                        evt_index[j]=nbunch_evt;
                    
                    i+=2; //Skip NEXT TWO ROWS, AS they are part of this bunch
                    nbunch_evt=1;
                }
                else
                {
                    if(delta <= bunchdeftime && i!=nrows-2)
                    {
                        nbunch_evt++;
                        //if(nbunch_evt>2) LOG(INFO)<<"SEE "<<nbunch_evt<<"\n";
                    }
                    else
                    {
                        for(j=i;j>i-nbunch_evt;j--)
                            evt_index[j]=nbunch_evt;

                        nbunch_evt=1;
                    }                        
                }
            }

            }
			else
            {
                if(row_index<nrows_bunch && i==bunch_row[row_index]) 
                {
                    bunch_flag[row_index]=1;  // Bad bunch (within 400 ms of 100s tick)
                    row_index++; 
                    for(j=i;j<=i+2;j++) evt_index[j]=0;
                    i+=2;
                }
                else
				    evt_index[i]=0;
            }
	    }

        if(nrows_bunch!=0) {
            long bad_bunch=0;
            LOG(INFO)<<"**row_index "<<row_index<<" "<<nrows_bunch;
            for(i=0;i<nrows_bunch;i++){
                if(bunch_flag[i]==1) bad_bunch++;
            }
            LOG(INFO)<<"Bunches within 400ms of 100s ticks "<<bad_bunch;
        }

        if(nrows_bunch==0) LOG(INFO)<<"**Bunches identified "<<bunch_time.size()<<"  "<<detid4.size()<<"\n";

    	long tot_evt=0;
	    for(i=0;i<nrows-1;i++)
	    {
            //printf("%ld\t%d\n",i,evt_index[i]);
		    if(evt_index[i]==1||evt_index[i]==2) tot_evt++;
            //else if (evt_index[i]!=0 && qid==3) printf("%ld\t%d\n",i,evt_index[i]);
	    }

	    LOG(INFO)<<"Number of events after cleaning stage 1: "<<tot_evt;

		// Write out the bunch index and id to optional output file
        //FILE *f1;
        //char index_file[1024];
        //sprintf(index_file,"%s%d%s","event_index_Q",qid,".txt");
        //f1=fopen(index_file,"w");
        //for(i=0;i<nrows-1;i++) fprintf(f1,"%ld\t%ld\t%d\n",i,evt_bunch_id[i],evt_index[i]);
        //fclose(f1);
		

		// Removing bunches and post bunch events and generate event file

        double prev_bunch_tstop=0.,prev_bunch_tstart=0.;
        int prev_bunch_module=-1;
        int prev_bunch_length=0;

        double prev2_btstart=0,prev2_btstop=0;
        double prev2_btstartmod=0,prev2_btstopmod=0;
        int prev2_bunch_module=-1;
//      double prev2_bunch_tstop=0.,prev2_bunch_tstart=0.;
//      int prev2_bunch_length=0;

        double delta_tstop;
        long currentGtiRow=0;

        double currentgtiStart=0;
        double CurrentgtiStop=0;
  
        long bunch_count=0;	
        for(i=0;i<nrows-1;)
        {
            
			if(evt_index[i]==1) // If the event is single (not in bunch)
            {
				
				delta_tstop=UTtime[i]-prev_bunch_tstop; //diff between current event and last bunch end time
				
				//Condition for event selection post bunch	
				if(((delta_tstop>skipT1)&&(prev_bunch_length<bunch_length_thresh)&&((detid[i]!=prev_bunch_module)||(delta_tstop>skipT2)))||(delta_tstop>skipT3))
                {
					// This single event is selected
                    evt_flag[i]=1;
                }
				
				i++; // Increment i by one
            }
			else if(evt_index[i]==2) // If it is double event
			{
				//Condition for event selection post bunch
				
				delta_tstop=UTtime[i]-prev_bunch_tstop;
					
				if(((delta_tstop>skipT1)&&(prev_bunch_length<bunch_length_thresh)&&((detid[i]!=prev_bunch_module)||(delta_tstop>skipT2)))||(delta_tstop>skipT3))
                {
					//This double event is selected
					evt_flag[i]=2;
					evt_flag[i+1]=2;
				}
                    
				i+=2; //Increment i by two	
			}
            else if(evt_index[i] > 2) // If the event is part of bunch (>2 events within bunchdeftime)
            {
				
				evt_flag[i]=10; // Set the first event of bunch

                //Set the prev2 parameters (the previous bunch parameters)
//                prev2_bunch_tstart=prev_bunch_tstart;
//                prev2_bunch_tstop= prev_bunch_tstop;
//                prev2_bunch_length=prev_bunch_length;
//               prev2_bunch_module=prev_bunch_module;

				// Get start and stop time of bunch and length of bunch
                prev_bunch_tstart=UTtime[i];
                if(nrows_bunch==0) prev_bunch_tstop=UTtime[i+evt_index[i]-1];
                else prev_bunch_tstop=UTtime[i+2];
                prev_bunch_length=evt_index[i];

				//Module histogram of events in bunch
                int modhist[16];
                for(k=0;k<16;k++) modhist[k]=0;

                if(nrows_bunch==0)
                {
                    for(j=i;j<i+evt_index[i];j++) //loop over events in the same bunch
                        modhist[detid[j]]++;
                }
                else
                {
                    for(j=i;j<=i+2;j++)
                        modhist[detid[j]]++;
                    
                    if(prev_bunch_length>3)
                    {
                        
                        //Add histogram with the other detector details from bunch file
                        for(k=0;k<4;k++)
                        {
                            if(prev_bunch_length-3 >= k) modhist[bunch_det[k][bunch_count]]++;
                        }
                    }
                    
                    bunch_count++;         
                }
				//Find module with maximum counts
                prev_bunch_module=-1;
                int maxhist=0;
                for(k=0;k<16;k++){if(modhist[k]>maxhist) prev_bunch_module=k;}

				if(nrows_bunch==0) i+=evt_index[i]; // Increment i by no. of events in bunch 
                else i+=3;

				//Compute the loss of exposure of modules and update livetime 

                
                double btstart=0,btstop=0;
                double btstartmod=0,btstopmod=0;
                double btstartomod=0,btstopomod=0;

                if(prev_bunch_length<bunch_length_thresh)
                {
                    btstart=prev_bunch_tstart;
                    btstop=prev_bunch_tstop+skipT1;
                    if(skipT2>skipT1)
                    {
                        btstartmod=prev_bunch_tstop+skipT1;
                        btstopmod=prev_bunch_tstop+skipT2;
                    }
                    else
                    {
                        btstartmod=0;
                        btstopmod=0;
                    }
                }
                else
                {
                    btstart=prev_bunch_tstart;
                   
                    if(skipT1>skipT3)
                        btstop=prev_bunch_tstop+skipT1;
                    else
                        btstop=prev_bunch_tstop+skipT3;

                    btstopmod=0;
                    btstartmod=0;
                }                

                //printf("%lf\t%lf\t%lf\n",btstop-btstart,btstopmod-btstartmod,btstopomod-btstartmod);

               
                if(btstopmod==0) //In this case removal of time is at quadrant level only
                { 
                    if(btstart<prev2_btstop)
                        btstart=prev2_btstop;

                    
//                    if(btstop<currentgtiStop)
//                    {
//                        btstop=currentgtiStop;
//                        //increment to next gti?
//                    }
                   
                    if(prev2_btstopmod!=0)
                    { 
                        if(btstart<prev2_btstop)
                        {
                            btstart=prev2_btstopmod;
                            btstartomod=prev2_btstop;
                            
                            if(btstop<prev2_btstopmod)
                            {
                                btstopomod=btstop;
                                btstop=prev2_btstopmod;
                            }
                            else
                                btstopomod=prev2_btstopmod;

                        }
                        else if(btstart<prev2_btstopmod)
                        {
                            btstart=prev2_btstopmod;
                            btstartomod=prev2_btstartmod;

                            if(btstop<prev2_btstopmod)
                            {
                                btstop=prev2_btstopmod;
                                btstopomod=btstop;
                            }
                            else
                                btstopomod=prev2_btstopmod;
                        }
                              
                    }

                    updatelivetime(btstart,btstop,livetime,tstarti,1.0,0,num_timebins,livetime_binsize);
                    updatelivetime(btstartomod,btstopomod,livetime,tstarti,15.0/16.0,0,num_timebins,livetime_binsize);
       
		    if(btstop>btstart){ 
                    for(k=0;k<4096;k++)
                       pixel_exposure[qid][k]+=(btstop-btstart);
		    }

                    if(btstopomod!=0&&btstopomod>btstartomod)
                    {
                        for(k=0;k<16;k++)
                        {
                            if(k!=prev2_bunch_module)
                            {
                                for(l=0;l<256;l++)
                                    pixel_exposure[qid][k*256+l]+=btstopomod-btstartomod;
                            }
                        }   
                    }

                }
                else 
                {
                   
                    if(prev2_btstopmod==0)    
                    {
                        if(btstart<prev2_btstop)
                        {
                            btstart=prev2_btstop;
                            if(btstop<prev2_btstop)
                                btstop=prev2_btstop;
                        }
                        if(btstartmod<prev2_btstop)
                        {
                            btstartmod=prev2_btstop;
                            if(btstopmod<prev2_btstop)
                                btstopmod=prev2_btstop; 
                        }
                    }
                    else
                    {
                        
                        if(btstart<prev2_btstop)
                        {
                            btstart=btstop;
                            btstop=btstop;
                            btstartomod=prev2_btstop;
                            btstopomod=btstop;

                            if(prev_bunch_module==prev2_bunch_module)
                                btstartmod=prev2_btstopmod;
                            else
                            {//Do not change btstrtmod and stopmod 
                            }

                        }
                        else if(btstart<prev2_btstopmod)
                        {
                            if(btstartmod<prev2_btstopmod)
                            {
                                btstartomod=btstart;
                                btstopomod=btstop;

                                if(prev_bunch_module==prev2_bunch_module)
                                    btstartmod=prev2_btstopmod;
                                
                                btstart=btstop;
                            }
                            else
                            {
                                btstartomod=btstart;
                                btstopomod=prev2_btstopmod;
                                btstart=prev2_btstopmod;
                            }
                        
                        }
                        
                    }
                 
                    updatelivetime(btstart,btstop,livetime,tstarti,1.0,0,num_timebins,livetime_binsize);
                    updatelivetime(btstartmod,btstopmod,livetime,tstarti,1.0/16.0,0,num_timebins,livetime_binsize);
                    updatelivetime(btstartomod,btstopomod,livetime,tstarti,15.0/16.0,0,num_timebins,livetime_binsize);

		    if(btstop>btstart){
                    for(k=0;k<4096;k++)
                        pixel_exposure[qid][k]+=(btstop-btstart);
		    }
                    //printf("%lf\t%lf\t%lf\n",btstop-btstart,btstopmod-btstartmod,btstopomod-btstartmod);
                    if(btstopmod>btstartmod){
		    for(l=0;l<256;l++)
                        pixel_exposure[qid][prev_bunch_module*256+l]+=btstopmod-btstartmod;
		    }
                    if(btstopomod!=0&&btstopomod>btstartomod)
                    {
                        for(k=0;k<16;k++)
                        {
                            if(k!=prev2_bunch_module)
                            {
                                for(l=0;l<256;l++)
                                    pixel_exposure[qid][k*256+l]+=btstopomod-btstartomod;
                            }
                        }
                    }


                }

                //Update prev2 params
                prev2_btstart=btstart;
                prev2_btstop=btstop;
                prev2_btstartmod=btstartmod;
                prev2_btstopmod=btstopmod;
                prev2_bunch_module=prev_bunch_module;

            }// End of if statement for bunch events
			else if(evt_index[i]==0)// In case of 100s ignored events
			{
				evt_flag[i]=0;
				i++;
			}


        }// Loop on nrows for cleaning events


	tot_evt=0;

        for(i=0;i<nrows-1;i++)
        {
                if(evt_flag[i]==1||evt_flag[i]==2) tot_evt++;
        }

	LOG(INFO)<<"Number of events after cleaning stage 2 : "<<tot_evt;

        //Fractional exposure computation
       
        for(j=0;j<4096;j++)
            pixel_exposure[qid][j]=(allquad_exp_time[qid]-pixel_exposure[qid][j])/(allquad_exp_time[qid]);            
    
    //Write out the livetime to file

    writelivetime(timearray,livetime,num_timebins,qid,livetimefile);


    //Write new event file with single and double events and one entry for each bunch
    remove_events(infile,outfile,evt_flag,qid,pixel_exposure);

	//If the optional bunch event file name is set, write out the bunch events
	if(outbunchevtstatus) write_bunch_events(infile,outbunchevtfile,evt_index,qid);


	// Write the bunch information to bunch file if it is empty (for cases with no on-board bunchclean)
    // otherwise write out the bunches with bunc_flag!=1
	if(nrows_bunch==0)
	{
	    long out_nrows_bunch=bunch_time.size();


		try{
            fits_write_col(foutbunch, TDOUBLE,  bunch_time_col, 1,1,out_nrows_bunch, bunch_time.data(), &status);
            fits_write_col(foutbunch, TBYTE, Time_dfs_col, 1,1,out_nrows_bunch, time_dfs.data(), &status);
            fits_write_col(foutbunch, TBYTE, Time_dsl_col, 1,1,out_nrows_bunch, time_dsl.data(), &status);
            fits_write_col(foutbunch, TBYTE, numevent_col, 1,1,out_nrows_bunch, num_bunchevents.data(), &status);
            fits_write_col(foutbunch, TBYTE, det1_col, 1,1,out_nrows_bunch, detid1.data(), &status);
            fits_write_col(foutbunch, TBYTE, det2_col, 1,1,out_nrows_bunch, detid2.data(), &status);
            fits_write_col(foutbunch, TBYTE, det3_col, 1,1,out_nrows_bunch, detid3.data(), &status);
            fits_write_col(foutbunch, TBYTE, det4_col, 1,1,out_nrows_bunch, detid4.data(), &status);
            fits_write_col(foutbunch, TLONG, rownum_col, 1,1,out_nrows_bunch, evt_row_num.data(), &status);
            if(bunchfileAddcolStatus)
            {
                fits_write_col(foutbunch, TBYTE, DetId_fevt_col, 1,1,out_nrows_bunch, DetId_fevt.data(), &status);
                fits_write_col(foutbunch, TBYTE, PixId_fevt_col, 1,1,out_nrows_bunch, PixId_fevt.data(), &status);
                fits_write_col(foutbunch, TBYTE, DetId_sevt_col, 1,1,out_nrows_bunch, DetId_sevt.data(), &status);
                fits_write_col(foutbunch, TBYTE, PixId_sevt_col, 1,1,out_nrows_bunch, PixId_sevt.data(), &status);
                fits_write_col(foutbunch, TBYTE, DetId_tevt_col, 1,1,out_nrows_bunch, DetId_tevt.data(), &status);
                fits_write_col(foutbunch, TBYTE, PixId_tevt_col, 1,1,out_nrows_bunch, PixId_tevt.data(), &status);
            }
        } 
        catch(ErrorHandler errHandler){
            logError(errHandler);
            return EXIT_FAILURE;
        }

	}
    else
    {

		int naxis1;
        fits_read_key(fbunch,TINT,"NAXIS1",&naxis1,NULL, &status);
 
        unsigned char *tblrow;
        tblrow=(unsigned char*)malloc(sizeof(unsigned char)*naxis1);
 
        long noutbunch=0;
 
        for(i=0;i<nrows_bunch;i++)
        {
            if(bunch_flag[i]==0)
            {
                fits_read_tblbytes(fbunch, i+1, 1,naxis1, tblrow, &status);
                if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
                noutbunch+=1;
                fits_write_tblbytes (foutbunch, noutbunch, 1,naxis1, tblrow, &status);
                if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
            }
         
        }
 
        free(tblrow);
        fits_update_key(foutbunch, TLONG, "NAXIS2",&noutbunch, NULL, &status);

	    LOG(INFO)<<"Bunches in output bunchfile: "<<noutbunch;
        
    }


    //Free allocated memory

	free(UTtime);
	free(czttime);
	free(detid);
	free(evt_index);
	free(evt_flag);
    	free(timearray);
    	free(livetime);

        if(nrows_bunch!=0)
	{
		free(bunch_row);
		free(numevent_bunch);
		for(j=0;j<4;j++) free(bunch_det[j]);
		free(bunch_det);
	}

    }//Quadrant loop ends here

    //Free allocated memory

    for(i=0;i<4;i++) free(pixel_exposure[i]);
    free(pixel_exposure);

    fits_close_file(fptr,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_close_file(fbunch,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_close_file(foutbunch,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    //writing history to all HDUs of output file
    if (history == YES) {
        vector<string> vhistory;
        getHistory(vhistory);
        if(writeHistory(outfile, vhistory)){ //writes history to each HDU of output filE
            LOG(ERROR) << "Error in writing History to output event file.";
        }
    }

   if(updateKeywords(outfile, modulename)){
        LOG(WARNING) << "Error in updating keywords.";
    }


	return (EXIT_SUCCESS);
}

int cztbunchclean::getHistory(vector<string> &vhistory) {
    //char *user = getlogin();
    strcpy(modulename, "cztbunchclean_v");
    strcat(modulename,VERSION);

    char *user = getenv("USER");
	vhistory.push_back("Module run by " + (string) user);
    vhistory.push_back("Parameter List START for " + (string) modulename);
    vhistory.push_back("P1 infile=" + (string) infile);
    vhistory.push_back("P2 outfile=" + (string) outfile);
    char bunchdeftime_str[25],bunch_length_thresh_str[25],skipT1_str[25];
    char skipT2_str[25],skipT3_str[25];
    sprintf(bunchdeftime_str,"%d",bunchdeftime);
    sprintf(bunch_length_thresh_str,"%d",bunch_length_thresh);
    sprintf(skipT1_str,"%f",skipT1);
    sprintf(skipT2_str,"%f",skipT2);
    sprintf(skipT3_str,"%f",skipT3);
    vhistory.push_back("P3 bunchdeftime=" + (string) bunchdeftime_str);
    vhistory.push_back("P4 bunch_length_thresh=" + (string) bunch_length_thresh_str );
    vhistory.push_back("P5 skipT1=" + (string) skipT1_str );
    vhistory.push_back("P6 skipT2=" + (string) skipT2_str );
    vhistory.push_back("P7 skipT3=" + (string) skipT3_str );    
    if (clobber == YES)
        vhistory.push_back("P8 clobber=yes");
    else
        vhistory.push_back("P8 clobber=no");
    if (history == YES)
        vhistory.push_back("P9 history=yes");
    else
        vhistory.push_back("P9 history=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

int remove_events(char* infile,char* outfile,char *evt_flag,int qid,double **pixel_exposure)
{

    int status=0,hdutype=0;
    int pi_col;
    long nrows,i;
    int naxis1;
	int maxPI=511;

    fitsfile *fptr, *fout;

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
    fits_create_file(&fout,outfile,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_copy_header(fptr, fout,  &status);
    }
    else
    {
        fits_open_file(&fout, outfile, READWRITE, &status);
    }

    LOG(INFO)<<"Quad "<<qid<<": began writing to new event file";

    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }


    fits_copy_header(fptr, fout,  &status);

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_read_key(fptr,TINT,"NAXIS1",&naxis1,NULL, &status);

    unsigned char *tblrow;
    tblrow=(unsigned char*)malloc(sizeof(unsigned char)*naxis1);

    long nevnt=0;
    
	for(i=0;i<nrows-1;i++)
    {
        if(evt_flag[i]==1||evt_flag[i]==2)
        {
            fits_read_tblbytes(fptr, i+1, 1,naxis1, tblrow, &status);
            if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
            nevnt+=1;
            fits_write_tblbytes (fout, nevnt, 1,naxis1, tblrow, &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
        }
        
    }

	free(tblrow);
    fits_update_key(fout, TLONG, "NAXIS2",&nevnt, NULL, &status);
    
    LOG(INFO)<<"Events after cleaning: "<<nevnt;


	if(qid==3)
	{

		//Copy the remaining extensions and write the exposure
	
		int num_hdus;

		fits_get_num_hdus(fptr, &num_hdus,&status);

		for(i=6;i<=num_hdus;i++)
		{
	    fits_movabs_hdu(fptr, i, &hdutype, &status);
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
		
		fits_copy_hdu(fptr, fout, 0,&status);
		if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }	
		}

		char extname[1024]="EXPOSURE";

	    int tfields=4;

	    char *ttype[] = {"EXPOSURE_Q0","EXPOSURE_Q1","EXPOSURE_Q2","EXPOSURE_Q3"};
	    char *tform[] = {"D","D","D","D"};

	    fits_create_tbl(fout, BINARY_TBL, 0,  tfields,ttype,tform,NULL, extname, &status);
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

		for(i=0;i<4;i++)
		{
            fits_write_col(fout, TDOUBLE, i+1, 1,1, 4096, pixel_exposure[i], &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
		}
		
		LOG(INFO)<<"Exposure extension appended to the output event file";
				
	}

fits_close_file(fptr,&status);
if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

//close output file

fits_close_file(fout,&status);
if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

return(0);
}

int write_bunch_events(char* infile,char* outfile,int *evt_index,int qid)
{

    int status=0,hdutype=0;
    int pi_col;
    long nrows,i;
    int naxis1;

    fitsfile *fptr, *fout;

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
    fits_create_file(&fout,outfile,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_copy_header(fptr, fout,  &status);
    }
    else
    {
        fits_open_file(&fout, outfile, READWRITE, &status);
    }

    LOG(INFO)<<"Quad "<<qid<<": began writing to new event file";

    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }


    fits_copy_header(fptr, fout,  &status);

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_read_key(fptr,TINT,"NAXIS1",&naxis1,NULL, &status);

    unsigned char *tblrow;
    tblrow=(unsigned char*)malloc(sizeof(unsigned char)*naxis1);

    long nevnt=0;

    for(i=0;i<nrows;i++)
    {
        if(evt_index[i]>2)
        {
            fits_read_tblbytes(fptr, i+1, 1,naxis1, tblrow, &status);
            if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }
            nevnt+=1;
            fits_write_tblbytes (fout, nevnt, 1,naxis1, tblrow, &status);
            if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
        }

    }

    free(tblrow);
    fits_update_key(fout, TLONG, "NAXIS2",&nevnt, NULL, &status);

    LOG(INFO)<<"Bunch events for this quadrant: "<<nevnt;	

	fits_close_file(fptr,&status);
	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	//close output file

	fits_close_file(fout,&status);
	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	return(0);
}


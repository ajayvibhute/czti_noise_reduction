/*
 * cztdatasel.cpp
 *
 * Rewritten based on the existing code. 
 * Edited to take five extension GTI and select
 * events based on common or quad gtitype.
 * 
 *
 * Mithun NPS (07/12/15)

 Edits:

 22/04/16
 
 * Changed tstart,tstop keywords of primary header to min{tstart(quad)} and max{tstop{quad}}

 * */

#include"cztdatasel.h"
#include"stdio.h"

using namespace std;

cztdatasel::cztdatasel(){
    strcpy(modulename,"cztdatasel_v");
    strcat(modulename,VERSION);
}

int cztdatasel::read(int argc, char** argv){
    int status=0;
    
    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error initializing PIL***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("infile",infile))){
        LOG(ERROR)<<"***Error reading input filename***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("gtifile",gtifile))){
        LOG(ERROR)<<"***Error reading GTI filename***";
        return status;
    }

/*    
    if(PIL_OK!=(status=PILGetString("extname",gtiextname))){
        LOG(ERROR)<<"***Error reading GTI extension name***";      
        return status;
    }//edited by TANUL just to include the GTI extension in the specified filename
    commented by Mithun (07/12/2015)
*/

    //Added to ask which GTITYPE is needed
    if(PIL_OK!=(status=PILGetString("gtitype",gtitype))){
        LOG(ERROR)<<"***Error reading GTI Type***";      
        return status;
    }    

    if(PIL_OK!=(status=PILGetFname("outfile",outfile))){
        LOG(ERROR)<<"***Error reading output filename***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetBool("clobber",&clobber))){
        LOG(ERROR)<<"***Error Reading clobber***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetBool("history",&history))){
        LOG(ERROR)<<"***Error Reading history parameter***";
        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

int cztdatasel::read(char* infile, char* gtitype, char* outfile, int clobber, int history){
    strcpy(this->infile,infile);
    strcpy(this->gtitype,gtitype);
    strcpy(this->outfile,outfile);
    this->clobber=clobber;
    this->history=history;
    return (EXIT_SUCCESS);
}

void cztdatasel::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                           CZTDATASEL PARAMETERS                     ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"Modulename          : "<<modulename;
    LOG(INFO)<<"Input Event file    : "<<infile;         //input event data file 
    LOG(INFO)<<"GTI file            : "<<gtifile;
    LOG(INFO)<<"GTITYPE             : "<<gtitype;
    LOG(INFO)<<"Output Event file   : "<<outfile;        //output event file
    if(clobber==YES)
    LOG(INFO)<<"Clobber             : YES";
    else
    LOG(INFO)<<"Clobber             : NO";
    if(history==YES)
    LOG(INFO)<<"History             : YES";
    else
    LOG(INFO)<<"History             : NO";    
	LOG(INFO)<<"---------------------------------------------------------------------------";
}

int cztdatasel::cztdataselProcess(){

     int status=0,i;
     fitsfile *fin,*fout;
     double tstart,tstop;


     if(strcmp(outfile,infile)==0){  
        LOG(ERROR)<<"outfile is same as infile.";
        return (EXIT_FAILURE);
        }
//     else if(strcmp(outfile," ")==0 || strcmp(outfile,"-")==0){
//         strcpy(outfile,infile);
//     }
     else {                           //infile and outfile are different
        if(FileExists(outfile)){
            if(clobber==YES){
                    if(unlink(outfile)!=0){
                        //cerr<<"Error in deleting "<<.outfile; 
                        //return (EXIT_FAILURE);
                    }
            }
        else{
            LOG(ERROR)<<""<<outfile<<" already exists";
            LOG(ERROR)<<"Use clobber=yes for overwriting the file";
            return (EXIT_FAILURE);
           }
        }
     
         //copying input to output file
         fits_open_file(&fin,infile,READWRITE,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<infile<<"***";
             return (EXIT_FAILURE);
         }
         fits_create_file(&fout,outfile,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in creating file "<<outfile<<"***";
             return (EXIT_FAILURE);
         }

         //Changed to copy only till TEMP extension 

         fits_movnam_hdu(fin, BINARY_TBL, "TEMP", 0, &status);
         if (status) {
            LOG(ERROR) <<"Error in moving to TEMP hdu  of input event file";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
         }
         
         fits_copy_file(fin,fout,1,1,0,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in copying file***";
             return (EXIT_FAILURE);
         }
         fits_close_file(fout,&status); if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
         fits_close_file(fin,&status);  if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
     }

     // Infile copied to outfile

     fitsfile *fgti,*fevt;

     char gtifilter_string[2000];
     
     //Opening the gti file
     fits_open_file(&fgti,gtifile,READONLY,&status);
     if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<gtifile<<"***";
             return (EXIT_FAILURE);
     }
          
     //changed on 9 sept 2017 by Mayuri     
     //opening event file 
     fits_open_file(&fevt,outfile,READWRITE,&status);   
     if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<outfile<<"***";
             return (EXIT_FAILURE);
     }

      

/*Assuming that event file has four extensions with four quad data.
 Also that GTI file has five extensions one common GTI and one for 
each quadrant GTI. Names of extensions are assumed to be
GTI,Q0_GTI etc.
gtitype decides the filtering condition
*/

    char extname[20],ext[5];
    int hdunum;
    double exp_time,avg_exp=0;
    cztHeaderParam headKey;
	double tstartmin=0,tstopmax=0;

    fits_update_key(fevt, TSTRING, "GTITYPE",&gtitype, NULL, &status);

    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in Updating header keywords";
        return (EXIT_FAILURE);
    }

    if(strcasecmp(gtitype,"COMMON")==0)
    {
        strcpy(extname,"GTI");
        fits_movnam_hdu(fgti, BINARY_TBL, extname, 0, &status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in moving to EXTNAME -"<<extname<< " in GTI file***";
            return (EXIT_FAILURE);    
        }
        
        fits_copy_hdu(fgti, fevt, 0,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in copying "<<extname<< " in GTI file to evt file";
            return (EXIT_FAILURE);
        }        

        fits_get_hdu_num(fgti, &hdunum);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in getting hdunumber";
            return (EXIT_FAILURE);
        }

        computeExptime(fgti,&exp_time,&tstart,&tstop);
//	cout<<"tstart is changed";
//	 cout<<tstart<<"\t"<<tstop;
        fits_movabs_hdu(fevt,1, NULL, &status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in moving to Primary of evt file *******";
            return (EXIT_FAILURE);
        }

		fits_update_key(fevt, TDOUBLE, "EXPOSURE",&exp_time, NULL, &status);
		headKey.writeTimekey(tstart,tstop,fevt);

        sprintf(gtifilter_string,"gtifilter(\"%s[%d]\")",gtifile,hdunum-1);

        for(i=0;i<4;i++)
        {

        sprintf(extname,"Q%d",i);

        fits_movnam_hdu(fevt, BINARY_TBL, extname, 0, &status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in moving to EXTNAME -"<<extname<< " in event file***";
            return (EXIT_FAILURE);
        }

        fits_select_rows(fevt,fevt,(char *)gtifilter_string,&status);
        if(status) { fits_report_error(stderr,status); return status; }

        fits_update_key(fevt, TDOUBLE, "EXPOSURE",&exp_time, NULL, &status);   
       
        fits_update_key(fevt, TSTRING, "GTITYPE",&gtitype, NULL, &status);
 
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in Updating header keywords";
            return (EXIT_FAILURE);
        }
	headKey.writeTimekey(tstart,tstop,fevt);


        }        
    }
    else if(strcasecmp(gtitype,"QUAD")==0)
    {

        for(i=0;i<4;i++)
        {
            sprintf(extname,"Q%d_GTI",i);
	    //sprintf(ext,"Q%d",i); //changed on 9 sept 2017 by Mayuri
	    status=0;
            fits_movnam_hdu(fgti, BINARY_TBL, extname, 0, &status);
            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in moving to EXTNAME -"<<extname<< " in GTI file***";
                return (EXIT_FAILURE);
            }

            fits_copy_hdu(fgti, fevt, 0,&status);
            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in copying "<<extname<< " in GTI file to evt file";
                return (EXIT_FAILURE);
            }

            fits_get_hdu_num(fgti, &hdunum);
            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in getting hdunumber";
                return (EXIT_FAILURE);
            }
	    status=0;
		
            computeExptime(fgti,&exp_time,&tstart,&tstop);
	    
	   // cout<<"QUAD tstart is changed"<<"\n";
	   // cout<<tstart<<"\t"<<tstop<<"\n";
			if(i==0) 
			{
				tstartmin=tstart;
				tstopmax=tstop;
			}
			else
			{
				if(tstart<tstartmin) tstartmin=tstart;
				if(tstop>tstopmax) tstopmax=tstop;	
			}

			avg_exp+=exp_time;

            sprintf(gtifilter_string,"gtifilter(\"%s[%d]\")",gtifile,hdunum-1);
        
            sprintf(extname,"Q%d",i);

            fits_movnam_hdu(fevt, BINARY_TBL, extname, 0, &status);
            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in moving to EXTNAME -"<<extname<< " in event file***";
                return (EXIT_FAILURE);
            }

            fits_select_rows(fevt,fevt,(char *)gtifilter_string,&status);
            if(status) { fits_report_error(stderr,status); return status; }
	    LOG(INFO)<<extname<<" Exposure : "<<exp_time;
            fits_update_key(fevt, TDOUBLE, "EXPOSURE",&exp_time, NULL, &status);
             
            fits_update_key(fevt, TSTRING, "GTITYPE",&gtitype, NULL, &status);

            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in Updating header keywords";
                return (EXIT_FAILURE);
            }

       	    headKey.writeTimekey(tstart,tstop,fevt);

        }

            fits_movabs_hdu(fevt,1, NULL, &status);
            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in moving to EXTNAME -"<<extname<< " in GTI file***";
                return (EXIT_FAILURE);
            }

			//Write average exposure and min{tstart} and max{tstop} to Primary header	
			avg_exp/=4.0;

			fits_update_key(fevt, TDOUBLE, "EXPOSURE",&avg_exp, NULL, &status);
			headKey.writeTimekey(tstartmin,tstopmax,fevt);	

    }
    else
    {
        LOG(ERROR)<<"***Invalid gtitype "<<gtitype<<"***";
        return (EXIT_FAILURE);
    }


    fits_close_file(fgti,&status);  if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    //Copy the EXPOSURE extension from input file to output file
    // New addition
    fits_open_file(&fin,infile,READWRITE,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in opening file "<<infile<<"***";
        return (EXIT_FAILURE);
    }

    fits_movnam_hdu(fin, BINARY_TBL, "EXPOSURE", 0, &status);
    if (status) {
       LOG(ERROR) <<"Error in moving to EXPOSURE hdu  of input event file";
       fits_report_error(stderr, status);
       return (EXIT_FAILURE);
    }

    fits_copy_file(fin,fevt,0,1,0,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in copying file***";
        return (EXIT_FAILURE);
    }

    fits_close_file(fin,&status);  if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    //Copying exposure extension done

    int numhdu;

    fits_get_num_hdus(fevt,&numhdu,&status);
    if(status) { fits_report_error(stderr,status); return status; }

    //writing history to all HDUs of output file
    if(history==YES){
        vector<string> historystr;
        getHistory(historystr);
        for(int i=1;i<=numhdu;i++){
             fits_movabs_hdu(fevt,i,NULL,&status);
             if(status) { fits_report_error(stderr,status); return status; }
             for(int j=0;j<historystr.size();j++){
                 fits_write_history(fevt,historystr[j].c_str(),&status);
                 if(status) { fits_report_error(stderr,status); return status; }
             }
        }
    }

    fits_close_file(fevt,&status);
    if(status) { fits_report_error(stderr,status); return status; }


// Below is the old source code of cztdatasel (commented by Mithun)

/*Assumption (OLD,NOT APPLICABLE NOW)
 * Quadrant 0 data in event file is at HDU 2
 * Quadrant 1 data in event file is at HDU 3
 * Quadrant 2 data in event file is at HDU 4
 * Quadrant 3 data in event file is at HDU 5
 * In case of four extensions in GTI the GTI extensions must start from HDU 2 and 
 * each extension header must have quadrant id
 * 
 * GTI file must have either one GTI extension valid for all four quadrants or
 * must have one GTI extension each for one quadrant 
 * Extension name for GTI must be GTI
 */

/*    
     fitsfile *fgti,*fevt;
     char errstring[2000];
     
     //opening event file 
     fits_open_file(&fevt,outfile,READWRITE,&status);   //infile is copied to outfile
     if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<outfile<<"***";
             return (EXIT_FAILURE);
         }
 
    
     //removing already existing GTI extensions in event file
     removeGTIext(fevt);
     
     char gtifilter_string[2000];
     
     fits_open_file(&fgti,gtifile,READONLY,&status);
     if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<gtifile<<"***";
             return (EXIT_FAILURE);
     }
     
     int numhdu=0;
     fits_get_num_hdus(fgti,&numhdu,&status);
     if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in getting number of HDUs in GTI file***";
             return (EXIT_FAILURE);
     }

     char extname[FLEN_VALUE];
     bool flag=true;                   //flag to show that GTI file contains all four extensions for four quadrants
     int quadid;
*/     
      //GTI file contains one extension
     /*fits_read_key(fgti,TSTRING,"EXTNAME",ext_name,NULL,&status);
     fits_movnam_hdu(fgti, 1, &extname, 0, &status);
     if(status){
              fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in moving to EXTNAME -"<<extname<< " in GTI file***";
             return (EXIT_FAILURE);    
     }*/

/*     //edited by TANUL
     fits_movnam_hdu(fgti,BINARY_TBL,gtiextname,0,&status);
     if(status){
         fits_report_error(stderr, status);
         LOG(ERROR)<<"***Error in moving to the specified GTI extension in GTI file***";
         return(EXIT_FAILURE);
     }
     
     sprintf(gtifilter_string,"gtifilter(\"%s\")",gtifile);
     LOG(INFO)<<"Filter :"<<gtifilter_string;
         
     fits_copy_hdu(fgti,fevt,0,&status);
     if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error copying HDU***";
        return (EXIT_FAILURE);
    }

    fits_close_file(fgti,&status);
    if(status){
         fits_report_error(stderr,status);
         LOG(ERROR)<<"***Error closing GTI file***";
         return (EXIT_FAILURE);
     }
*/

     //edited by TANUL
     
     //earlier GTI reading module that used to check the HDU NUMBER. if its 2 and extension name is GTI, it appends the GTI extension to the event file. 
     //the same module contains code to incorporate multiple GTI extensions..
     /*if(numhdu==2){
         fits_movabs_hdu(fgti,2,NULL,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in moving to 2nd HDU in gti file***";
             return (EXIT_FAILURE);
         }
         
         fits_read_key(fgti,TSTRING,"EXTNAME",extname,NULL,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error reading EXTNAME from GTI file***";
             return (EXIT_FAILURE);
         }
         
         if(strcasecmp(extname,"GTI")!=0){
             LOG(ERROR)<<"***Error: GTI extension not found in GTI file***";
             return (EXIT_FAILURE);
         }         
         sprintf(gtifilter_string,"gtifilter(\"%s\")",.gtifile);
         LOG(INFO)<<"Filter :"<<gtifilter_string;
         
        fits_copy_hdu(fgti,fevt,0,&status);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error copying HDU***";
             return (EXIT_FAILURE);
         }
     }
     else{
         //case where GTI file has more than one GTI extension
         for(int i=2;i<=numhdu;i++){
             fits_movabs_hdu(fgti,i,NULL,&status);
             if(status){
                 sprintf(errstring,"***Error moving to HDU %d of %s - cztdatasel_process()***",i,.gtifile);
                 fits_report_error(stderr,status);
                 LOG(ERROR)<<errstring;
                 return (EXIT_FAILURE);
             }
             fits_read_key(fgti,TSTRING,"EXTNAME",extname,NULL,&status);
             if(status){
                 sprintf(errstring,"***Error reading EXTNAME from GTI file");
                 fits_report_error(stderr,status);
                 LOG(ERROR)<<errstring;
                 return (EXIT_FAILURE);
             }
             if(strcasecmp(extname,"GTI")==0){
                fits_read_key(fgti,TINT,"QUADID",&quadid,NULL,&status);
                if(status){
                     fits_report_error(stderr,status);
                     LOG(ERROR)<<"***Error reading QUADID from GTI file***";
                     return (EXIT_FAILURE);
                 }
                if(quadid==i-2)  { 
                        fits_copy_hdu(fgti,fevt,0,&status);
                        flag=flag & true; 
                        continue;
                }
                else flag=flag & false;
                if(quadid>3)     break;
             }
         }
     } 

     
     if(flag==false){
         LOG(ERROR)<<"***GTI file is not in proper format***";
         return (EXIT_FAILURE);
     }
*/     
   
/* 
    //loop for filtering
    long nrows;             //for number of rows in event data extension of event file
    
    //loop for HDU 2 to 5 for event file that contains event data for each quadrant
    for(int i=2;i<=5;i++){
        fits_movabs_hdu(fevt,i,NULL,&status);
        if(status) { fits_report_error(stderr,status); return status; }
        fits_get_num_rows(fevt,&nrows,&status);
        if(status) { fits_report_error(stderr,status); return status; }
        if(nrows<=0) continue;
        if(numhdu==2){                                //If GTI file contains one GTI extension
            LOG(INFO)<<"GTI with one extension";
            fits_select_rows(fevt,fevt,(char *)gtifilter_string,&status);
            if(status) { fits_report_error(stderr,status); return status; }
            continue;
        }
        else if(flag){                              //If GTI file file contains extension for each quadrant separately
            LOG(INFO)<<"GTI with four extensions";
            sprintf(gtifilter_string,"gtifilter(\"%s[%d]\")",gtifile,i);
            //cout<<"Applying filter :"<<gtifilter_string;
            fits_select_rows(fevt,fevt,(char *)gtifilter_string,&status);
            if(status) { fits_report_error(stderr,status); return status; }
            continue;
        }
        else{
            return (EXIT_FAILURE);
        }
    }

*/

/*
    //Updating Header Time
    //UPDATING KEYWORDS
    if (updateHdrTime(outfile, "Q0", "TIME")) {
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q0";
        return EXIT_FAILURE;
    }
    if (updateHdrTime(outfile, "Q1", "TIME")) {
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q1";
        return EXIT_FAILURE;
    }
    if (updateHdrTime(outfile, "Q2", "TIME")) {
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q2";
        return EXIT_FAILURE;
    }
    if (updateHdrTime(outfile, "Q3", "TIME")) {
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q3";
        return EXIT_FAILURE;
    }
*/  
    //UPDATING KEYWORDS
    updateKeywords(outfile, modulename);  
    return (EXIT_SUCCESS);
}

//Compute exposure time 
int computeExptime(fitsfile *fgti,double *exp_time, double *tstart, double *tstop)
{
    int status=0;
    long nrows,j;

    fits_get_num_rows(fgti, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    double *gti_start,*gti_stop,exposure_time=0.;
    gti_start=(double*)malloc(sizeof(double)*nrows);
    gti_stop=(double*)malloc(sizeof(double)*nrows);


    fits_read_col(fgti, TDOUBLE, 1, 1, 1, nrows, NULL,gti_start,NULL, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
    fits_read_col(fgti, TDOUBLE, 2, 1, 1, nrows, NULL, gti_stop,NULL, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    for(j=0;j<nrows;j++) exposure_time+=gti_stop[j]-gti_start[j];

    *exp_time=exposure_time;
    *tstart=gti_start[0];
    *tstop=gti_stop[nrows-1];
   // cout<<"TSTART "<<gti_start[0]<<" TSTOP "<<gti_start[nrows-1];  

    return(EXIT_SUCCESS);
}


//function to remove GTI extensions from the file
int removeGTIext(fitsfile *fptr){
    int status=0;
    fits_movabs_hdu(fptr,1,NULL,&status);
    if(status) { fits_report_error(stderr,status); return status; }
    while(fits_movnam_hdu(fptr,BINARY_TBL,"GTI",0,&status)==0){
        fits_delete_hdu(fptr,NULL,&status);
        if(status) { fits_report_error(stderr,status); return status; }
    }
    return (EXIT_SUCCESS);
}

//Function definition for creating history for the module
int cztdatasel::getHistory(vector<string> &vhistory){
   // char *user=getlogin();//commented by Mayuri,5 Oct 2017
    strcpy(modulename,"cztdatasel_v");
    strcat(modulename,VERSION);
    char *user=getenv("USER");
	vhistory.push_back("Module run by "+(string)user);
    vhistory.push_back("Parameter List START for "+(string)modulename);
    vhistory.push_back("P1 infile="+(string)infile);
    vhistory.push_back("P2 gtifile="+(string)gtifile);
    vhistory.push_back("P3 outfile="+(string)outfile);
    vhistory.push_back("P4 gtitype="+(string)gtitype);    
    if(clobber==YES) 
        vhistory.push_back("P5 clobber=yes");
    else
        vhistory.push_back("P5 clobber=no");
    if(history==YES)
        vhistory.push_back("P6 history=yes");
    else
        vhistory.push_back("P6 history=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

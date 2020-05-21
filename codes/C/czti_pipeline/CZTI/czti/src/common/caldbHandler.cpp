
#include "caldbHandler.h"
#include "ExpMap.h"


using namespace std;

/* CALDB file query routine
*   Rakesh Khanna (12/12/15)
*/
int QueryCaldb(string telescope,string instrument,string detname,string codename,double tstart,double tstop,string &caldb_filename,int &extnum)
{

    string caldbpath;

    if(getenv("CALDB") == NULL)
    {
        LOG(ERROR)<<"CALDB environment variable not set. Exiting..";
        return(EXIT_FAILURE);
    }

    caldbpath=getenv("CALDB");

    string cifpath=caldbpath+"/data/as1/czti/caldb.indx";
    
    fitsfile *fptr;
    int status=0, hdutype=0;
    long nrows,i;
    int tel_col=1,instr_col=2,detnam_col=3,filter_col=4,cal_dev_col=5,cal_dir_col=6;
    int cal_file_col=7,cal_clas_col=8,cal_dtyp_col=9,cal_cnam_col=10,cal_cbd_col=11;
    int cal_xno_col=12,cal_vsd_col=13,cal_vst_col=14,ref_time_col=15,cal_qual_col=16;
    int cal_date_col=17,cal_desc_col=18;
    char *cif_telescope,*cif_instrument;
    char *cif_detname,cif_filtname[10],*cif_codename; 
    char *cif_caldir,*cif_calfile;
    int cif_calextnum,cif_calquality;
    double cif_reftime,start_reftime,end_reftime;
    double mjdreftime=0,maxRefTime=-1;    
    long maxRefTime_index=-1;

    cif_telescope=(char*)malloc(sizeof(char *)*10);
    if(cif_telescope==NULL){printf("Error allocating Memory\n");exit(0);}

    cif_instrument=(char*)malloc(sizeof(char *)*10);
    if(cif_instrument==NULL){printf("Error allocating Memory\n");exit(0);}

    cif_detname=(char*)malloc(sizeof(char *)*20);
    if(cif_detname==NULL){printf("Error allocating Memory\n");exit(0);}

    cif_codename=(char*)malloc(sizeof(char *)*20);
    if(cif_codename==NULL){printf("Error allocating Memory\n");exit(0);}

    cif_caldir=(char*)malloc(sizeof(char *)*70);
    if(cif_caldir==NULL){printf("Error allocating Memory\n");exit(0);}

    cif_calfile=(char*)malloc(sizeof(char *)*40);
    if(cif_calfile==NULL){printf("Error allocating Memory\n");exit(0);}

    mjdreftime=(tstart/86400.0)+55197.0;

    //Open CIF file and move to 1st extension
    fits_open_file(&fptr,cifpath.c_str(),READONLY,&status);
    if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}    

    fits_movabs_hdu(fptr,2,&hdutype,&status);
    if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    if(nrows==0)
    { 
        cout<<"CIF file seems to be empty...\n";
        return(EXIT_FAILURE);
    }

    int row_flag[nrows];
    for(i=0;i<nrows;i++) row_flag[i]=0;

    fits_read_col_str(fptr,tel_col,1,1,1,NULL,&cif_telescope,NULL,&status);
    if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}    

    fits_read_col_str(fptr,instr_col,1,1,1,NULL,&cif_instrument,NULL,&status);
    if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

    if(strcmp(cif_telescope,telescope.c_str())|| strcmp(cif_instrument,instrument.c_str())) 
    {
        cout<<"Mission/Instrument name is not matching\n";
        return(EXIT_FAILURE);
    }
   
     for(i=0;i<nrows;i++)
        {
        fits_read_col_str(fptr,detnam_col,i+1,1,1,NULL,&cif_detname,NULL,&status);
        if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}    
        
        fits_read_col_str(fptr,cal_cnam_col,i+1,1,1,NULL,&cif_codename,NULL,&status);
        if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}
    
        fits_read_col(fptr,TINT,cal_xno_col,i+1,1,1,NULL,&cif_calextnum,NULL,&status);
        if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}
        
        fits_read_col(fptr,TINT,cal_qual_col,i+1,1,1,NULL,&cif_calquality,NULL,&status);
        if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

        fits_read_col(fptr,TDOUBLE,ref_time_col,i+1,1,1,NULL,&cif_reftime,NULL,&status);
        if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}


	if(((!strcasecmp(cif_codename,codename.c_str()))||(!strcasecmp(codename.c_str(),"-")))
                && (!strcmp(cif_detname,detname.c_str())||!strcmp(detname.c_str(),"-"))
                && cif_calquality==0&&mjdreftime>cif_reftime) 
        {         
                row_flag[i]=1;
        }

        }
    
        //Select the latest of the data sets
        for(i=0;i<nrows;i++)
        {
             if(row_flag[i]==1)
             {
        
                fits_read_col(fptr,TDOUBLE,ref_time_col,i+1,1,1,NULL,&cif_reftime,NULL,&status);
                if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

                if(cif_reftime>maxRefTime)
                { 
                    maxRefTime=cif_reftime;
                    maxRefTime_index=i;   
                }
             }
        }

 
        // Obtaining required caldb file and directory
        if(maxRefTime_index!=-1)
        {
            fits_read_col_str(fptr,cal_dir_col,maxRefTime_index+1,1,1,NULL,&cif_caldir,NULL,&status);
            if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

            fits_read_col_str(fptr,cal_file_col,maxRefTime_index+1,1,1,NULL,&cif_calfile,NULL,&status);
            if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

            //Appending file path
            caldb_filename=caldbpath+'/'+cif_caldir+'/'+cif_calfile;
        }
        else
        {
            LOG(ERROR)<<"No caldb file found with given input conditions";
            return(EXIT_FAILURE);
        }
            //Close the CIF file

    fits_close_file(fptr,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); } 

    //Free allocated memory
    free(cif_telescope);
    free(cif_instrument);
    free(cif_detname);
    free(cif_codename);
    free(cif_caldir);
    free(cif_calfile);

   return(EXIT_SUCCESS);
}

/* Query to get the cif file name

*/
int getcif(string telescope,string instrument,string cif)
{

    if(getenv("CALDB") == NULL)
    {
        LOG(ERROR)<<"CALDB environment variable not set. Exiting..";
        return(EXIT_FAILURE);
    }

    if(getenv("CALDBCONFIG") == NULL)
    {
        cout<<"CALDBCONFIG environment variable not set... Quitting..\n";
        return(EXIT_FAILURE);
    }
    
    
    cout<<getenv("CALDB")<<"\n";
    cout<<getenv("CALDBCONFIG")<<"\n";

    cif="/home/mithun/bin/heasoft/caldb/data/as1/czti/caldb.indx";

    return(EXIT_SUCCESS);
}


/** CALDB HANDLER **/

void CaldbHandler::generate_caldb_filepaths(string basepath) {
    string CALDBbasepath = "";
    if (basepath == "$CALDBnew" || basepath == "-") {
        char* caldbenvpath = getenv("CALDBnew");
        CALDBbasepath = (string) caldbenvpath;
    } else {
        CALDBbasepath = basepath;
    }
    this->basepath = CALDBbasepath;
    CaldbFiles["TELDEFQ0"] = (this->basepath + "/bcf/"  + "AS1czt020150211v01.teldef");
    CaldbFiles["TELDEFQ1"] = this->basepath + "/bcf/"  + "AS1czt120150211v01.teldef";
    CaldbFiles["TELDEFQ2"] = this->basepath + "/bcf/"  + "AS1czt220150211v01.teldef";
    CaldbFiles["TELDEFQ3"] = this->basepath + "/bcf/"  + "AS1czt320150211v01.teldef";
    CaldbFiles["BADPIX"] = this->basepath + "/bcf/"  + "AS1cztbadpix20150526v01.fits";
    CaldbFiles["CAMERAGEO"] = this->basepath + "/bcf/"  + "AS1cztcamerageo20150629v01.fits";
    CaldbFiles["DETMAP"] = this->basepath + "/bcf/"  + "AS1cztdetmap20150526v01.fits";
    CaldbFiles["EBOUNDS"] = this->basepath + "/bcf/"  + "AS1cztebounds20150526v01.fits";
    CaldbFiles["EFFAREA"] = this->basepath + "/bcf/"  + "AS1czteff_area20150617v01.fits";
    CaldbFiles["GAINOFFSET"] = this->basepath + "/bcf/"  + "AS1cztgain20150526v01.fits";
    CaldbFiles["COMPMASK"] = this->basepath + "/bcf/"  + "AS1cztmask_pattern20150526v01.fits";
    CaldbFiles["MASK"] = this->basepath + "/bcf/"  + "CZTIMask64x64.fits";
    CaldbFiles["RESPONSE"] = this->basepath + "/bcf/"  + "AS1cztresp_par20150617v01.fits";
    CaldbFiles["CATALOG"] = this->basepath + "/bcf/"  + "BAT_58m_catalog_100924_new.fits.gz";
}

//Getters
string CaldbHandler::get_caldb_file_path(string filetype){
    ErrorHandler errHandler;
    string filepath="";
    if(CaldbFiles.find(filetype) == CaldbFiles.end()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = INVALID_CALDB_FILETYPE;
        errHandler.errorMsg = "Invalid CALDB filetype: " + filetype;
        throw errHandler;
    } else {
        filepath = CaldbFiles[filetype];
    }
    return filepath;
}
/** CALDB HANDLER END **/
/**Teldef class functions to read and store teldef information**/
Teldef::Teldef() {
}


int Teldef::read_teldef_file(string teldefFilename) {
    int status=0;
    
    fitsfile *fptr;  
    // reading fits file
    fits_open_file(&fptr, teldefFilename.c_str(), READONLY, &status );
    if (status){
        LOG(ERROR) << "*** Error in opening teldef file: " << teldefFilename << " ***";
        fits_report_error(stderr, status);
        return(EXIT_FAILURE);
    }
    
    //temporary variables to store header information from teldef file which is of type character array.
    char tempftype[MAX_KEYWORD_SIZE];
    char tempcrd0[MAX_KEYWORD_SIZE];
    char tempcrd1[MAX_KEYWORD_SIZE];
    char tempcrd2[MAX_KEYWORD_SIZE];
    char temprawXcol[MAX_KEYWORD_SIZE];
    char temprawYcol[MAX_KEYWORD_SIZE];
    char temprawUnit[MAX_KEYWORD_SIZE];
    char tempdetXcol[MAX_KEYWORD_SIZE];
    char tempdetYcol[MAX_KEYWORD_SIZE];
    char tempdetUnit[MAX_KEYWORD_SIZE];
    char tempskyXcol[MAX_KEYWORD_SIZE];
    char tempskyYcol[MAX_KEYWORD_SIZE];
    char tempskyUnit[MAX_KEYWORD_SIZE];
    char tempskyFrom[MAX_KEYWORD_SIZE];
    
    // reading teldef keywords from header file
    fits_read_key(fptr, TSTRING,  "CCNM0001", tempftype, NULL, &status );
    report_error(status, (string) "Error in reading keyword CCNM001");
    fileType = (string) tempftype;
    
    fits_read_key(fptr, TSTRING,  "COORD0", tempcrd0, NULL, &status );
    report_error(status, (string) "Error in reading keyword COORD0");
    coord0 = (string) tempcrd0;
    
    fits_read_key(fptr, TSTRING,  "COORD1", tempcrd1, NULL, &status );
    report_error(status, (string) "Error in reading keyword COORD1");
    coord1 = (string) tempcrd1;
    
    fits_read_key(fptr, TSTRING,  "COORD2", tempcrd2, NULL, &status );
    report_error(status, (string) "Error in reading keyword COORD0");
    coord2 = (string) tempcrd2;
   
    fits_read_key(fptr, TSTRING,  "RAW_XCOL", temprawXcol, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_XCOL");
    rawXcol = (string) temprawXcol;
    
    fits_read_key(fptr, TFLOAT,  "RAW_XSIZ", &rawXsize, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_XSIZ");

    fits_read_key(fptr, TFLOAT,  "RAW_XSCL", &rawXscl, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_XSCL");

    fits_read_key(fptr, TFLOAT,  "RAWXPIX1", &rawXpix1, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAWXPIX1");    
     
    fits_read_key(fptr, TSTRING,  "RAW_YCOL", temprawYcol, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_XCOL");
    rawYcol = (string) temprawXcol;
    
    fits_read_key(fptr, TFLOAT,  "RAW_YSIZ", &rawYsize, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_YSIZ");

    fits_read_key(fptr, TFLOAT,  "RAW_YSCL", &rawYscl, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_YSCL");

    fits_read_key(fptr, TFLOAT,  "RAWYPIX1", &rawYpix1, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAWYPIX1");
    
    fits_read_key(fptr, TSTRING,  "RAW_UNIT", temprawUnit, NULL, &status );
    report_error(status, (string) "Error in reading keyword RAW_UNIT");
    rawUnit = (string) temprawUnit;
    

    fits_read_key(fptr, TSTRING,  "DET_XCOL", tempdetXcol, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_XCOL");
    detXcol = (string) tempdetXcol;
    
    fits_read_key(fptr, TFLOAT,  "DET_XSIZ", &detXsize, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_XSIZ");

    fits_read_key(fptr, TFLOAT,  "DET_XSCL", &detXscl, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_XSCL");

    fits_read_key(fptr, TFLOAT,  "DETXPIX1", &detXpix1, NULL, &status );
    report_error(status, (string) "Error in reading keyword DETXPIX1");    
     
    fits_read_key(fptr, TSTRING,  "DET_YCOL", tempdetYcol, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_XCOL");
    detYcol = (string) tempdetXcol;
    
    fits_read_key(fptr, TFLOAT,  "DET_YSIZ", &detYsize, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_YSIZ");

    fits_read_key(fptr, TFLOAT,  "DET_YSCL", &detYscl, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_YSCL");

    fits_read_key(fptr, TFLOAT,  "DETYPIX1", &detYpix1, NULL, &status );
    report_error(status, (string) "Error in reading keyword DETYPIX1");
    
    fits_read_key(fptr, TSTRING,  "DET_UNIT", tempdetUnit, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_UNIT");
    detUnit = (string) tempdetUnit;
    
    fits_read_key(fptr, TFLOAT,  "COE_X_A", &coeXA, NULL, &status );
    report_error(status, (string) "Error in reading keyword COE_X_A");
    
    fits_read_key(fptr, TFLOAT,  "COE_X_B", &coeXB, NULL, &status );
    report_error(status, (string) "Error in reading keyword COE_X_B");
    
    fits_read_key(fptr, TFLOAT,  "COE_X_C", &coeXC, NULL, &status );
    report_error(status, (string) "Error in reading keyword COE_X_C");
    
    fits_read_key(fptr, TFLOAT,  "COE_Y_A", &coeYA, NULL, &status );
    report_error(status, (string) "Error in reading keyword COE_Y_A");
    
    fits_read_key(fptr, TFLOAT,  "COE_Y_B", &coeYB, NULL, &status );
    report_error(status, (string) "Error in reading keyword COE_Y_B");
    
    fits_read_key(fptr, TFLOAT,  "COE_Y_C", &coeYC, NULL, &status );
    report_error(status, (string) "Error in reading keyword COE_Y_C");
    
    fits_read_key(fptr, TFLOAT,  "DET_XOFF", &detXoff, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_XOFF");    

    fits_read_key(fptr, TFLOAT,  "DET_YOFF", &detYoff, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_YOFF");
    
    fits_read_key(fptr, TFLOAT,  "DETXFLIP", &detXflip, NULL, &status );
    report_error(status, (string) "Error in reading keyword DETXFLIP"); 

    fits_read_key(fptr, TFLOAT,  "DETYFLIP", &detYflip, NULL, &status );
    report_error(status, (string) "Error in reading keyword DETYFLIP");
    
    fits_read_key(fptr, TFLOAT,  "DET_SCAL", &detScal, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_SCAL");
    
    fits_read_key(fptr, TFLOAT,  "DET_ROTD", &detRotd, NULL, &status );
    report_error(status, (string) "Error in reading keyword DET_ROTD");
    

    fits_read_key(fptr, TSTRING,  "SKY_XCOL", tempskyXcol, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKY_XCOL");
    skyXcol = (string) tempskyXcol;
    
    fits_read_key(fptr, TFLOAT,  "SKY_XSIZ", &skyXsize, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKY_XSIZ");

    fits_read_key(fptr, TFLOAT,  "SKYXPIX1", &skyXpix1, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKYXPIX1");    
     
    fits_read_key(fptr, TSTRING,  "SKY_YCOL", tempskyYcol, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKY_XCOL");
    skyYcol = (string) tempskyXcol;
    
    fits_read_key(fptr, TFLOAT,  "SKY_YSIZ", &skyYsize, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKY_YSIZ");

    fits_read_key(fptr, TFLOAT,  "SKYYPIX1", &skyYpix1, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKYYPIX1");
    
    fits_read_key(fptr, TSTRING,  "SKY_UNIT", tempskyUnit, NULL, &status );
    report_error(status, (string) "Error in reading keyword SKY_UNIT");
    skyUnit = (string) tempskyUnit;

    fits_read_key(fptr, TSTRING, "SKY_FROM", tempskyFrom, NULL, &status);
    report_error(status, (string) "Error in reading keyword SKY_FROM");
    skyFrom = (string) tempskyFrom;

    fits_read_key(fptr, TFLOAT, "ALIGNM11", &alignM11, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM11");

    fits_read_key(fptr, TFLOAT, "ALIGNM12", &alignM12, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM12");

    fits_read_key(fptr, TFLOAT, "ALIGNM13", &alignM13, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM13");

    fits_read_key(fptr, TFLOAT, "ALIGNM21", &alignM21, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM21");

    fits_read_key(fptr, TFLOAT, "ALIGNM22", &alignM22, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM22");

    fits_read_key(fptr, TFLOAT, "ALIGNM23", &alignM23, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM23");

    fits_read_key(fptr, TFLOAT, "ALIGNM31", &alignM31, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM31");

    fits_read_key(fptr, TFLOAT, "ALIGNM32", &alignM32, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM32");

    fits_read_key(fptr, TFLOAT, "ALIGNM33", &alignM33, NULL, &status);
    report_error(status, (string) "Error in reading keyword ALIGNM33");
    
    fits_read_key(fptr, TFLOAT, "FOCALLEN", &focalLength, NULL, &status);
    report_error(status, (string) "Error in reading keyword FOCALLEN");
    
    fits_read_key(fptr, TFLOAT, "OPTAXISY", &optAxisY, NULL, &status);
    report_error(status, (string) "Error in reading keyword OPTAXISY");
    
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "ERROR IN CLOSING TELDEF FILE.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    return status;
            
}

int Teldef::display(){

    LOG(INFO) << "TELDEF FILE KEYWORDS";
    LOG(INFO) << "------------------------";
    LOG(INFO) << "File type: " << fileType;
    LOG(INFO) << "Coord0   : " << coord0;
    LOG(INFO) << "Coord1   : " << coord1;
    LOG(INFO) << "Coord2   : " << coord2;
    
    LOG(INFO) << "rawXcol  : " << rawXcol;
    LOG(INFO) << "rawXsize : " << rawXsize;
    LOG(INFO) << "rawXpix  : " << rawXpix1;
    LOG(INFO) << "rawXscl  : " << rawXscl;
    LOG(INFO) << "rawYcol  : " << rawYcol;
    LOG(INFO) << "rawYsize : " << rawYsize;
    LOG(INFO) << "rawYpix  : " << rawYpix1;
    LOG(INFO) << "rawYscl  : " << rawYscl;
    LOG(INFO) << "rawUnit  : " << rawUnit;

    LOG(INFO) << "detXcol  : " << detXcol;
    LOG(INFO) << "detXsize : " << detXsize;
    LOG(INFO) << "detXpix  : " << detXpix1;
    LOG(INFO) << "detXscl  : " << detXscl;
    LOG(INFO) << "detYcol  : " << detYcol;
    LOG(INFO) << "detYsize : " << detYsize;
    LOG(INFO) << "detYpix  : " << detYpix1;
    LOG(INFO) << "detYscl  : " << detYscl;
    LOG(INFO) << "detUnit  : " << detUnit;
    
    LOG(INFO) << "coeXA    : " << coeXA;
    LOG(INFO) << "coeXB    : " << coeXB;
    LOG(INFO) << "coeXC    : " << coeXC;
    LOG(INFO) << "coeYA    : " << coeYA;
    LOG(INFO) << "coeYB    : " << coeYB;
    LOG(INFO) << "coeYC    : " << coeYC;

    LOG(INFO) << "detXoff  : " << detXoff;
    LOG(INFO) << "detYoff  : " << detYoff;
    LOG(INFO) << "detXflip : " << detXflip;
    LOG(INFO) << "detYflip : " << detYflip;
    LOG(INFO) << "detScal  : " << detScal;
    LOG(INFO) << "detRotd  : " << detRotd;
    
    LOG(INFO) << "skyXcol  : " << skyXcol;
    LOG(INFO) << "skyXsize : " << skyXsize;
    LOG(INFO) << "skyXpix  : " << skyXpix1;
    LOG(INFO) << "skyYcol  : " << skyYcol;
    LOG(INFO) << "skyYsize : " << skyYsize;
    LOG(INFO) << "skyYpix  : " << skyYpix1;
    LOG(INFO) << "skyUnit  : " << skyUnit;
    LOG(INFO) << "skyFrom  : " << skyFrom;
    
    LOG(INFO) << "alignM11 : " << alignM11;
    LOG(INFO) << "alignM12 : " << alignM12;
    LOG(INFO) << "alignM13 : " << alignM13;
    LOG(INFO) << "alignM21 : " << alignM21;
    LOG(INFO) << "alignM22 : " << alignM22;
    LOG(INFO) << "alignM23 : " << alignM23;
    LOG(INFO) << "alignM31 : " << alignM31;
    LOG(INFO) << "alignM32 : " << alignM32;
    LOG(INFO) << "alignM33 : " << alignM33;
    
    LOG(INFO) << "focalLen : " << focalLength;
    LOG(INFO) << "optAxisY : " << optAxisY;


    
    
    
}


/** Ebounds class functions to read and store energy bounds **/
Ebounds::Ebounds(){
}

int Ebounds::read_ebounds_file(string eboundsFilename){
    int status=0;   // status variable
    int i,j=0;  // counter variables
    int colnum=0;
    long nrows=0;
    string errorMsg="";
    fitsfile *fptr; // Pointer to EBOUNDS CALDB FILE.
    unsigned short *channelArray;
    float *eminArray;
    float *emaxArray;
 
    // reading fits file
    fits_open_file(&fptr, eboundsFilename.c_str(), READONLY, &status );
    errorMsg = "*** Error in opening ebounds file: " + eboundsFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_movnam_hdu(fptr, BINARY_TBL, "EBOUNDS", 0, &status);
    errorMsg = "Error in reading EBOUNDS extension in EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in ebounds extension of EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "CHANNEL", &colnum, &status);
    errorMsg = "Error in getting column number of CHANNEL column in EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    channelArray = new unsigned short[nrows];
    for(i=0; i<nrows; i++){
        channelArray[i] =0;
    }
    fits_read_col(fptr, TSHORT, colnum, 1, 1, nrows, NULL, channelArray, NULL, &status);
    errorMsg = "Error in reading CHANNEL column of ebounds extension of EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "E_MIN", &colnum, &status);
    errorMsg = "Error in getting column number of E_MIN column in EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;} 
    eminArray = new float[nrows];
    for(i=0; i<nrows; i++){
        eminArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, eminArray, NULL, &status);
    errorMsg = "Error in reading E_MIN column of ebounds extension of EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    
    fits_get_colnum(fptr, CASEINSEN, "E_MAX", &colnum, &status);
    errorMsg = "Error in getting column number of E_MAX column in EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;} 
    emaxArray = new float[nrows];
    for(i=0; i<nrows; i++){
        emaxArray[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows, NULL, emaxArray, NULL, &status);
    errorMsg = "Error in reading E_MAX column of ebounds extension of EBOUNDS file: " + eboundsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}   
    
    fits_close_file(fptr, &status);
    errorMsg = "*** Error in closing EBOUNDS file: " + eboundsFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    // fits file has been read and closed
    
    
    // Assigning data read from ebounds file to class variables
    for(i=0; i<nrows; i++){
        channel.push_back(channelArray[i]);
        eMin.push_back(eminArray[i]);
        eMax.push_back(emaxArray[i]);
    }
    // Data assigned to class variables
   
    //deleting dynamically declared variables;
    delete[] channelArray, eminArray, emaxArray;
    
    return status;
}

unsigned short Ebounds::find_channel_from_energy(float energy, int& status){
    int i=0; //counter variable
    status=0;
//    if (energy > eMax[eMax.size()-1]) {
//        status = EXIT_FAILURE;
//        return status;
//    }
    
	if(energy < eMin[0])
		return channel[0];
	else if (energy >= eMax[get_nrows()-1])
		return channel[get_nrows()-1];
	else
	{
		for (i=0; i<channel.size(); i++){
        	if(energy>= eMin[i] && energy<eMax[i])
            	return channel[i];
    	}    
    }
}

float Ebounds::find_energy_from_channel(unsigned short channelNo){
    string errorMsg="";
    int ichannel=0; //index on channel
    bool channelFlag=false;
    float energy=0.0;
    
    for(ichannel=0; ichannel<channel.size(); ichannel++){
        if(channel[ichannel] == channelNo){
            channelFlag=true;
            energy = (eMin[ichannel] + eMax[ichannel])/2.0; //calculating the mean of emin and emax
        }
    }
    if (channelFlag == false) {
        errorMsg = "Channel number " + itoa(channelNo) + " is not present in ebounds file.";
        throw errorMsg;

    }
    return energy;
}

float Ebounds::get_energy(int rowno){
    float energy=0.0;
    string errorMsg="";
    if(rowno<0 || rowno>channel.size()){
        errorMsg = "Row number should lie in the range 0-" + itoa(channel.size()-1);
        throw errorMsg;
    }
    else{
        energy = (eMin[rowno] + eMax[rowno])/2.0; //calculating the mean of emin and emax 
    }
    return energy;
}
unsigned int Ebounds::get_channel(int rowno){
    unsigned int PI;
    string errorMsg="";
    if(rowno<0 || rowno>channel.size()){
        errorMsg = "Row number should lie in the range 0-" + itoa(channel.size()-1);
        throw errorMsg;
    }
    else{
        PI = channel[rowno]; 
    }
    return PI;
}

int Ebounds::display(){
    int status=0;
    int minChannel=0, maxChannel=0;
    float minEmin=0.0, maxEmin=0.0;
    float minEmax=0.0, maxEmax=0.0;
    
    compute_min_and_max(channel.begin(), channel.end(), minChannel, maxChannel);
    compute_min_and_max(eMin.begin(), eMin.end(), minEmin, maxEmin);
    compute_min_and_max(eMax.begin(), eMax.end(), minEmax, maxEmax);
    LOG(INFO) << ">>> Statistics: EBOUNDS FILE <<<";

    LOG(INFO) << "Number of channels        :" << channel.size();
    LOG(INFO) << "Minimum channel number [Energy Range]:" << minChannel << " [" << minEmin << "-" << minEmax <<"]";
    LOG(INFO) << "Maximum channel number [Energy Range]:" << maxChannel << " [" << maxEmin << "-" << maxEmax <<"]";

    return status;

}

/****************GAIN CLASS***********************/
GainOffset::GainOffset(){
    
}
int GainOffset::read_gainoffset_file(string gainsFilename, string extname){
    int status=0; //status variable
    int i,j=0; //counter variables
    int colnum=0;
    long nrows=0;
    string errorMsg="";
    fitsfile *fptr; //Pointer to GAINS CALDB FILE;
    char *detidArray;
    char *pixidArray;
    float *go00Array;
    float *go05Array;
    float *go10Array;
    float *go15Array;
    float *go20Array;
    
    //reading fits file
    fits_open_file(&fptr, gainsFilename.c_str(), READONLY, &status);
    errorMsg = "*** Error in opening GAINS file: " + gainsFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_movnam_hdu(fptr, BINARY_TBL, (char *) extname.c_str(), 0, &status);
    errorMsg = "Error in reading " + extname+ " extension in GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in " + extname + " extension of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "DETID", &colnum, &status);
    errorMsg = "Error in getting column number of DETID column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    detidArray = new char[nrows];
    for(i=0; i<nrows; i++){
        detidArray[i] =0;
    }
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detidArray, NULL, &status);
    errorMsg = "Error in reading DETID column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_colnum(fptr, CASEINSEN, "PIXID", &colnum, &status);
    errorMsg = "Error in getting column number of PIXID column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    pixidArray = new char[nrows];
    for(i=0; i<nrows; i++){
        pixidArray[i] =0;
    }
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, pixidArray, NULL, &status);
    errorMsg = "Error in reading PIXID column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    //reading Gain Offset Array for temperature 0
    fits_get_colnum(fptr, CASEINSEN, "GAIN_00", &colnum, &status);
    errorMsg = "Error in getting column number of GAIN_00 column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    go00Array = new float[nrows*2];
    for (i = 0; i < nrows*2; i++) {
        go00Array[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows*2, NULL, go00Array, NULL, &status);
    errorMsg = "Error in reading GAIN_00 column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    //reading Gain Offset Array for temperature 05
    fits_get_colnum(fptr, CASEINSEN, "GAIN_05", &colnum, &status);
    errorMsg = "Error in getting column number of GAIN_05 column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    go05Array = new float[nrows*2];
    for (i = 0; i < nrows*2; i++) {
        go05Array[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows*2, NULL, go05Array, NULL, &status);
    errorMsg = "Error in reading GAIN_05 column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    //reading Gain Offset Array for temperature 10
    fits_get_colnum(fptr, CASEINSEN, "GAIN_10", &colnum, &status);
    errorMsg = "Error in getting column number of GAIN_10 column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    go10Array = new float[nrows*2];
    for (i = 0; i < nrows*2; i++) {
        go10Array[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows*2, NULL, go10Array, NULL, &status);
    errorMsg = "Error in reading GAIN_10 column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    //reading Gain Offset Array for temperature 15
    fits_get_colnum(fptr, CASEINSEN, "GAIN_15", &colnum, &status);
    errorMsg = "Error in getting column number of GAIN_15 column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    go15Array = new float[nrows*2];
    for (i = 0; i < nrows*2; i++) {
        go15Array[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows*2, NULL, go15Array, NULL, &status);
    errorMsg = "Error in reading GAIN_15 column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    //reading Gain Offset Array for temperature 20
    fits_get_colnum(fptr, CASEINSEN, "GAIN_20", &colnum, &status);
    errorMsg = "Error in getting column number of GAIN_20 column in extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    go20Array = new float[nrows*2];
    for (i = 0; i < nrows*2; i++) {
        go20Array[i] = 0.0;
    }
    fits_read_col(fptr, TFLOAT, colnum, 1, 1, nrows*2, NULL, go20Array, NULL, &status);
    errorMsg = "Error in reading GAIN_20 column of extension " + extname + " of GAINS file: " + gainsFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_close_file(fptr, &status);
    errorMsg = "*** Error in closing GAINS file: " + gainsFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    // fits file has been read and closed
    
    // Assigning data read to class variables
    for(i=0; i<nrows; i++){
        detID.push_back(detidArray[i]);
        pixID.push_back(pixidArray[i]);
        gain00.push_back(go00Array[i*2]);
        gain05.push_back(go05Array[i*2]);
        gain10.push_back(go10Array[i*2]);
        gain15.push_back(go15Array[i*2]);
        gain20.push_back(go20Array[i*2]);
        offset00.push_back(go00Array[2*i + 1]);
        offset05.push_back(go05Array[2*i + 1]);
        offset10.push_back(go10Array[2*i + 1]);
        offset15.push_back(go15Array[2*i + 1]);
        offset20.push_back(go20Array[2*i + 1]);

    }   // Data assigned to class variables
    
    //deleting dynamically declared variables
    delete[] detidArray, pixidArray;
    delete [] go00Array, go05Array, go10Array, go15Array, go20Array;
    
    return EXIT_SUCCESS;
}

int GainOffset::display(){
    int status=0;
    int i=0; //counter variable
    long nrows= detID.size();
    
    //printing all gains and offsets on terminal
    LOG(INFO)<< "Displaying Gain Offset file.";
    cout << setw(10) << "DETID" << setw(10) << "PIXID" << setw(10) << "G00" << setw(10) << "OFF00" 
            << setw(10) << "G05" << setw(10) << "OFF05" << setw(10) << "G10" << setw(10) << "OFF10" 
            << setw(10) << "G15" << setw(10) << "OFF15" << setw(10) << "G20" << setw(10) << "OFF20" << endl;
    cout << setw(10) << "-----" << setw(10) << "-----" << setw(10) << "---" << setw(10) << "-----" 
            << setw(10) << "---" << setw(10) << "-----" << setw(10) << "---" << setw(10) << "-----" 
            << setw(10) << "---" << setw(10) << "-----" << setw(10) << "---" << setw(10) << "-----" << endl;
    for(i=0; i<nrows; i++){
        cout<<setw(10)<<((int)detID[i])<<setw(10)<<((int)pixID[i])<<setw(10)<<gain00[i]<<setw(10)<<offset00[i]
                <<setw(10)<<gain05[i]<<setw(10)<<offset05[i]<<setw(10)<<gain10[i]<<setw(10)<<offset10[i]
                <<setw(10)<<gain15[i]<<setw(10)<<offset15[i]<<setw(10)<<gain20[i]<<setw(10)<<offset20[i]<<endl;
    }
    

    return status;
    
}

vector<float> GainOffset::get_gain(int temperature, int& status) {
    vector<float> gain;
    switch (temperature) {
        case 1:
            gain = gain00;
            break;
        case 2:
            gain = gain05;
            break;
        case 3:
            gain = gain10;
            break;
        case 4:
            gain = gain15;
            break;
        case 5:
            gain = gain20;
            break;
        default:
            status = EXIT_FAILURE;
            LOG(ERROR) << "No gain values for Temperature: T" << temperature;
    }
    return gain;
}
vector<float> GainOffset::get_offset(int temperature, int& status) {
    vector<float> offset;
    switch (temperature) {
        case 1:
            offset = offset00;
            break;
        case 2:
            offset = offset05;
            break;
        case 3:
            offset = offset10;
            break;
        case 4:
            offset = offset15;
            break;
        case 5:
            offset = offset20;
            break;
        default:
            status = EXIT_FAILURE;
            LOG(ERROR) << "No offset values for Temperature: T" << temperature;
    }
    return offset;
}

/***************RESPONSE PAR CLASS******************/
ResponsePar::ResponsePar() {

}

int ResponsePar::read_response_par_file(string responseParsFilename, string extname){
   int status=0; //status variable;
   int i,j=0;   //counter variable;
   int colnum=0;
   long nrows=0;
   string errorMsg="";
   fitsfile *fptr; //Pointer to RESPONSE_PAR CALDB File
   char *detidArray;
   char *pixidArray;
   char *detxArray;
   char *detyArray;
   float *mResArray;
   float *cResArray;
   
   // reading fits file
   fits_open_file(&fptr, responseParsFilename.c_str(), READONLY, &status);
   errorMsg = "*** Error in RESPONSE file: " + responseParsFilename + " ***";
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}


   fits_movnam_hdu(fptr, BINARY_TBL, (char *) extname.c_str(), 0, &status);
   errorMsg = "Error in reading " + extname + " extension in RESPONSE file: " + responseParsFilename;
   if (report_error(status, errorMsg)) {return EXIT_FAILURE;}
   
   fits_get_num_rows(fptr, &nrows, &status);
   errorMsg = "Error in getting number of rows in " + extname + " extension of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

   fits_get_colnum(fptr, CASEINSEN, "DETID", &colnum, &status);
   errorMsg = "Error in getting column number of DETID column in extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
   detidArray = new char[nrows];
   for(i=0; i<nrows; i++){
       detidArray[i] =0;
   }
   fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detidArray, NULL, &status);
   errorMsg = "Error in reading DETID column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   
   fits_get_colnum(fptr, CASEINSEN, "PIXID", &colnum, &status);
   errorMsg = "Error in getting column number of PIXID column in extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
   pixidArray = new char[nrows];
   for(i=0; i<nrows; i++){
       pixidArray[i] =0;
   }
   fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, pixidArray, NULL, &status);
   errorMsg = "Error in reading PIXID column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   fits_get_colnum(fptr, CASEINSEN, "DETX", &colnum, &status);
   errorMsg = "Error in getting column number of DETX column in extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   detxArray = new char[nrows];
   for (i = 0; i < nrows; i++) {
       detxArray[i] = 0;
   }
   fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detxArray, NULL, &status);
   errorMsg = "Error in reading DETX column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   
   fits_get_colnum(fptr, CASEINSEN, "DETY", &colnum, &status);
   errorMsg = "Error in getting column number of DETY column in extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   detyArray = new char[nrows];
   for (i = 0; i < nrows; i++) {
       detyArray[i] = 0;
   }
   fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detyArray, NULL, &status);
   errorMsg = "Error in reading DETY column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
   fits_get_colnum(fptr, CASEINSEN, "M_RES", &colnum, &status);
   errorMsg = "Error in reading M_RES column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}  
   mResArray = new float[nrows];
   for (i=0; i<nrows; i++){
       mResArray[i]=0;
   }
   fits_read_col(fptr, TFLOAT, colnum, 1,1, nrows, NULL, mResArray, NULL, &status);
   errorMsg = "Error in reading M_RES column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   
   fits_get_colnum(fptr, CASEINSEN, "C_RES", &colnum, &status);
   errorMsg = "Error in reading C_RES column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}  
   cResArray = new float[nrows];
   for (i=0; i<nrows; i++){
       cResArray[i]=0;
   }
   fits_read_col(fptr, TFLOAT, colnum, 1,1, nrows, NULL, cResArray, NULL, &status);
   errorMsg = "Error in reading C_RES column of extension " + extname + " of RESPONSE file: " + responseParsFilename;
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   
   fits_close_file(fptr, &status);
   errorMsg = "*** Error in closing RESPONSE file: " + responseParsFilename + " ***";
   if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   // fits file has been read

   

   // Assigning data read to class variables
   for(i=0; i<nrows; i++){  
        detID.push_back(detidArray[i]);
        pixID.push_back(pixidArray[i]);
        detX.push_back(detxArray[i]);
        detY.push_back(detyArray[i]);
        mRes.push_back(mResArray[i]);
        cRes.push_back(cResArray[i]);
   } // Data assigned to class variables
   
  
    // deleting dynamically declared variables
   delete[] detidArray, pixidArray, detxArray, detyArray;
   delete[] mResArray, cResArray;
    
   return status;
}

int ResponsePar::display() {
    int status=0;
    int i=0; //counter variable
    vector <unsigned char> uniqueDetID;
    uniqueDetID = detID;
    if (uniqueDetID.size()) {
        vector <unsigned char>::iterator pd = unique(uniqueDetID.begin(), uniqueDetID.end());
        uniqueDetID.erase(pd, uniqueDetID.end());
        LOG(INFO) << ">>> UNIQUE DETECTOR IDS FOUND in RESPONSE PAR FILE <<<";
        for (i = 0; i < uniqueDetID.size(); i++) {
            LOG(INFO) << (int) uniqueDetID[i];
        }
    }
    return status;    
}
/************* EFFECTIVE AREA CLASS ******************/

EffArea::EffArea() {

}

int EffArea::read_effarea_file(string effareaFilename, string extname){
    int status=0;
    int i,j=0;
    int colnum=0;
    long nrows=0;
    string errorMsg="";
    fitsfile *fptr; //Pointer to EFFAREA CALDB file
    char *detidArray;
    char *pixidArray;
    char *detxArray;
    char *detyArray;
    float *areaArray;

    //Clearing class vectors
    detID.clear();
    pixID.clear();
    detX.clear();
    detY.clear();
    area.clear();
    
    //reading fits file
    fits_open_file(&fptr, effareaFilename.c_str(), READONLY, &status);
    errorMsg = "*** Error in opening EFFAREA file: " + effareaFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_movnam_hdu(fptr, BINARY_TBL, (char *) extname.c_str(), 0, &status);
    errorMsg = "Error in reading " + extname+ " extension in EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in " + extname + " extension of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    
    fits_get_colnum(fptr, CASEINSEN, "DETID", &colnum, &status);
    errorMsg = "Error in getting column number of DETID column in extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    detidArray = new char[nrows];
    for(i=0; i<nrows; i++){
        detidArray[i] =0;
    }
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detidArray, NULL, &status);
    errorMsg = "Error in reading DETID column of extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_colnum(fptr, CASEINSEN, "PIXID", &colnum, &status);
    errorMsg = "Error in getting column number of PIXID column in extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    pixidArray = new char[nrows];
    for(i=0; i<nrows; i++){
        pixidArray[i] =0;
    }
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, pixidArray, NULL, &status);
    errorMsg = "Error in reading PIXID column of extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "DETX", &colnum, &status);
    errorMsg = "Error in getting column number of DETX column in extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    detxArray = new char[nrows];
    for (i = 0; i < nrows; i++) {
        detxArray[i] = 0;
    }
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detxArray, NULL, &status);
    errorMsg = "Error in reading DETX column of extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}


    fits_get_colnum(fptr, CASEINSEN, "DETY", &colnum, &status);
    errorMsg = "Error in getting column number of DETY column in extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    detyArray = new char[nrows];
    for (i = 0; i < nrows; i++) {
        detyArray[i] = 0;
    }

    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detyArray, NULL, &status);
    errorMsg = "Error in reading DETY column of extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}


    fits_get_colnum(fptr, CASEINSEN, "AREA", &colnum, &status);
    errorMsg = "Error in getting column number of AREA column in extension " + extname + " of EFFAREA file: " + effareaFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    int areaSize=300;
    vector <float> tempArea;
    areaArray = new float[areaSize];
    for (i = 0; i < areaSize; i++) {
        areaArray[i] = 0.0;
    }        
    for (i=0; i<nrows; i++){
        tempArea.clear();
        fits_read_col(fptr, TFLOAT, colnum, i+1,1, areaSize, NULL, areaArray, NULL, &status);
        errorMsg = "Error in reading AREA column of extension " + extname + " of EFFAREA file: " + effareaFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
        for(j=0; j<300; j++){
            tempArea.push_back(areaArray[j]);
        }
        area.push_back(tempArea);
    }
    
    fits_close_file(fptr, &status);
    errorMsg = "*** Error in closing EFFAREA file: " + effareaFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    // fits file has been read and closed
    
    
    // Assigning data read to class variables
    for(i=0; i< nrows; i++){
        detID.push_back(detidArray[i]);
        pixID.push_back(pixidArray[i]);
        detX.push_back(detxArray[i]);
        detY.push_back(detyArray[i]);
        //NOTE: area vector has already been assigned values while reading from
        //      fits file.
    }
    
    //deleting dynamically declared arrays
    delete[] detidArray, pixidArray, detxArray, detyArray;
    delete [] areaArray;
    return status;
}


float EffArea::add_effarea(vector <float> area1d, float eMin, float eMax) {
    int i=0; //counter variable
    int iMin = floor(eMin)-1;// Minimum Index
    int iMax = floor(eMax)-1;// Maximum Index
    float sumEffArea=0.0;
    for(i=iMin; i<=iMax; i++){
        sumEffArea= sumEffArea + area1d[i];
    }
    
    return sumEffArea;
}

int EffArea::calculate_normalized_effarea(float eMin, float eMax){
    int status=0;
    long tempIndex=0;
    long i,ix,iy=0; //counter variable
    long nrows=0;
    vector <float> tempNormEffArea;
    float sumEffArea=0.0;
    maxEffArea=0.0;
    
    nrows=detID.size();
    if(nrows!=PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD){
        LOG(ERROR)<<"Number of rows in effective area file : "<< nrows;
        LOG(ERROR)<<"Expected number of rows : " << PIXELS_PER_ROW_QUAD*PIXELS_PER_COL_QUAD;
        return EXIT_FAILURE;
    }
//    if(floor(eMax)-floor(eMin)==0){
//        LOG(ERROR)<<"Minimum Energy Bin Size should be greater than or equal to 1";
//        return EXIT_FAILURE;
//    }
    
    //Initializing normalized effective area 2D vector
    tempNormEffArea.clear();
    for(ix=0; ix<PIXELS_PER_ROW_QUAD; ix++){
        tempNormEffArea.push_back(0.0);
    }
    normalizedEffArea.clear();
    for(iy=0; iy<PIXELS_PER_COL_QUAD; iy++){
        normalizedEffArea.push_back(tempNormEffArea);
    }
    
    LOG(INFO)<< "Size of normalize effective area 2D vector : [YxX]" 
            << normalizedEffArea.size() << "x" << tempNormEffArea.size(); 
    //Normalized effective area initialized.
    
    for(i=0; i<nrows; i++){
        sumEffArea = add_effarea(area[i], eMin, eMax);
        if(sumEffArea>maxEffArea){
            maxEffArea=sumEffArea;
            tempIndex=i;
        }
        normalizedEffArea[detY[i]][detX[i]]=sumEffArea;
    }
    LOG(INFO)<<"Maximum Effective Area is:"<< setprecision(20) <<maxEffArea;
    //LOG(INFO)<<"Index is :"<< tempIndex;
    for(ix=0; ix<PIXELS_PER_ROW_QUAD; ix++){
        for(iy=0; iy<PIXELS_PER_COL_QUAD; iy++){
            normalizedEffArea[ix][iy]=normalizedEffArea[ix][iy]/maxEffArea;
        }
    }
    return EXIT_SUCCESS;
}

int EffArea::calculate_full_normalized_effarea(string effareaFilename, float eMin, float eMax) {
    int status=0;
    int i=0;
    char extname[FLEN_VALUE];
    vector <float> tempNormEffArea;
    tempNormEffArea.resize(NPIXELS_X, 0.0);
    fullNormalizedEffArea.resize(NPIXELS_Y, tempNormEffArea);
    for(i=0; i<NUMQUAD; i++){
        quadToHDU(i, extname);
        if(read_effarea_file(effareaFilename, (string)extname)){
            LOG(ERROR) << "Error in reading HDU " << extname << " of effective area file " << effareaFilename;
            return EXIT_FAILURE;
        }
        if(calculate_normalized_effarea(eMin, eMax)){
            LOG(ERROR)<< "Error in calculating normalized effective area for HDU " << extname;
            return EXIT_FAILURE;
        }
        if(rearrange_quads(&fullNormalizedEffArea, normalizedEffArea, i)){
            LOG(ERROR) << "Error in arranging individual quadrant to generate "
                    << " full normalized effective area.";
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
    
}

vector<float> EffArea::get_effective_area(unsigned char pixx, unsigned char pixy, vector<float> ekeV){
    vector<float> energies;
    int ienergy = 0;
    int iekeV=0; //index on user specified energy
    int ipix=0;
    int pixrecord=-1;
    float effarea1=0.0, effarea2=0.0; //y1,y2
    float energy1=0.0, energy2=0.0; //x1,x2
    vector <float> effarea;
    string errorMsg="";
    float minekeV=0.0; //minimum of user specified energy
    float maxekeV=0.0; //maximum of user specified energy
    int nenergy = area[0].size();
    energies.resize(nenergy, 0.0);
    effarea.resize(ekeV.size(), 0.0);
    
    for(ienergy=0; ienergy<nenergy; ienergy++){
        energies[ienergy] = 1.5 + float(ienergy);
    }
    
    for(ipix=0; ipix<(int)get_nrows(); ipix++){
        if(pixx==detX[ipix] && pixy==detY[ipix]){
            pixrecord=ipix; //record number of pixel for which ekeV has to be evaluated.
            break;
        }
    }
    
    //Getting minimum and maximum energy values in user supplied vector
    compute_min_and_max(ekeV.begin(), ekeV.end(), minekeV, maxekeV);
    if (minekeV < energies[0] || maxekeV > energies[nenergy - 1]) {
        errorMsg = "Effective area can only be evaluated in energy range " + itoa(energies[0]) + " - " + itoa(energies[nenergy - 1]);
        throw errorMsg;
    }
    
    if(pixrecord==-1){
        errorMsg = "Effective area value not available for" + itoa((int)pixx) + "," + itoa((int)pixy);
        throw errorMsg;
    }
    else {
        for (iekeV = 0; iekeV < ekeV.size(); iekeV++) {
            for (ienergy = 1; ienergy < nenergy; ienergy++) {
                if (ekeV[iekeV] <= energies[ienergy] && ekeV[iekeV] >= energies[ienergy - 1]) {
                    effarea1 = area[pixrecord][ienergy - 1];
                    effarea2 = area[pixrecord][ienergy];
                    energy1 = energies[ienergy - 1];
                    energy2 = energies[ienergy];
                    effarea[iekeV] = effarea1 + ((effarea2 - effarea1) / (energy2 - energy1)) * (ekeV[iekeV]-energy1);
                    //DLOG(INFO) << energy1 << " " << energy2 << ", " << effarea1 << ", " << effarea2 <<  " Energy(keV): " << ekeV[iekeV] << " Effective area: " << effarea[iekeV];
                }
                
            }
        }
    }
    
    return effarea;
}

void EffArea::reset() {
    detID.clear();
    pixID.clear();
    detX.clear();
    detY.clear();
    normalizedEffArea.clear();
    fullNormalizedEffArea.clear();
}


/************** CZTI CAMERA GEOMETRY CLASS**************/

CztiCameraGeometry::CztiCameraGeometry() {

}

int CztiCameraGeometry::define_camera_geometry(){
    int isurf=0; //counter variable
    int status=0; //status variable
    int nsmax=100; //maximum number of surfaces
    int nsurf=62; //number of surfaces in the model
    
    // All dimensions are in cm
    // Surfaces parallel to yz plane; 29
    // Collimator slats: 12
    // Collimator housing walls, upper 36cm: 8
    // Collimator housing walls, lower 4 cm: 8
    // Middle division in calibration housing: 1
    
    //initializing surface variables
    for (isurf = 0; isurf <= nsurf; isurf++) {
        zmin.push_back(0.0);
        zmax.push_back(0.0);
        xminz.push_back(0.0);
        xminc.push_back(0.0);
        xmaxz.push_back(0.0);
        xmaxc.push_back(0.0);
        yminz.push_back(0.0);
        yminc.push_back(0.0);
        ymaxz.push_back(0.0);
        ymaxc.push_back(0.0);
        xnp.push_back(0.0);
        ynp.push_back(0.0);
        znp.push_back(0.0);
        dnp.push_back(0.0);
        thickAl.push_back(0.0);
        thickTa.push_back(0.0);
        openfrac.push_back(0.0);
        indexSurf.push_back(isurf);
    }

    //Disabling all other surfaces which are not part of the model
    for (isurf = 1; isurf <= nsurf; isurf++) {
        disable.push_back(0);
    }
    for (isurf = nsurf + 1; isurf <= 100; isurf++) {
        disable.push_back(1);
    }
    
    for(isurf=1; isurf<=29; isurf++) {
        xnp[isurf] = 1.0;
        ynp[isurf] = 0.0;
        znp[isurf] = 0.0;
        openfrac[isurf] = 0.0;
        xminz[isurf] = 0.0;
        xmaxz[isurf] = 0.0;
        yminz[isurf] = 0.0;
        ymaxz[isurf] = 0.0;
    }

    for (isurf = 1; isurf <= 12; isurf++) {
        thickAl[isurf] = 0.04;
        thickTa[isurf] = 0.007;
        zmin[isurf] = 7.7;
        zmax[isurf] = 47.7;
    }

    for (isurf = 1; isurf <= 3; isurf++) {
        yminc[isurf] = -18.25;
        ymaxc[isurf] = -1.8;
        dnp[isurf] = -18.325 + isurf * 4.15;
        xminc[isurf] = dnp[isurf] - 0.0235;
        xmaxc[isurf] = dnp[isurf] + 0.0235;
    }

    for (isurf = 4; isurf <= 6; isurf++) {
        yminc[isurf] = -18.25;
        ymaxc[isurf] = -1.8;
        dnp[isurf] = 1.725 + (isurf - 3)*4.15;
        xminc[isurf] = dnp[isurf] - 0.0235;
        xmaxc[isurf] = dnp[isurf] + 0.0235;
    }

    for (isurf = 7; isurf <= 12; isurf++) {
        yminc[isurf] = 1.8;
        ymaxc[isurf] = 18.25;
        dnp[isurf] = dnp[isurf - 6];
        xminc[isurf] = dnp[isurf] - 0.0235;
        xmaxc[isurf] = dnp[isurf] + 0.0235;
    }

    //Walls of Collimator Housing

    //Top Parts
    for (isurf = 13; isurf <= 20; isurf++) {
        zmin[isurf] = 11.7;
        zmax[isurf] = 47.7;
        thickAl[isurf] = 0.0;
        thickTa[isurf] = 0.007;
    }

    for (isurf = 13; isurf <= 16; isurf++) {
        yminc[isurf] = -18.25;
        ymaxc[isurf] = -1.8;
    }

    for (isurf = 17; isurf <= 20; isurf++) {
        yminc[isurf] = 1.8;
        ymaxc[isurf] = 18.25;
    }
    
    for (isurf = 13; isurf <= 14; isurf++) {
 	   dnp[isurf]=-18.25+16.45*(isurf-13);
	   dnp[29-isurf]=-dnp[isurf];
	   dnp[isurf+4]=dnp[isurf];
	   dnp[33-isurf]=-dnp[isurf];       
    }
    
    for (isurf = 13; isurf <= 20; isurf++) {
	   xminc[isurf]=dnp[isurf]-0.035;
	   xmaxc[isurf]=dnp[isurf]+0.035;        
    }

    //Bottom Parts
    for (isurf = 21; isurf <= 28; isurf++) {
        zmin[isurf] = 7.7;
        zmax[isurf] = 11.7;
        thickAl[isurf] = 0.0;
        thickTa[isurf] = 0.01;
        dnp[isurf] = dnp[isurf - 8];
        xminc[isurf] = xminc[isurf - 8];
        xmaxc[isurf] = xmaxc[isurf - 8];
        yminc[isurf] = yminc[isurf - 8];
        ymaxc[isurf] = ymaxc[isurf - 8];
    }

    // Middle Partition of Calibration Housing
    isurf = 29;
    chsm = 0.3239552;
    chsc = 20.744455;
    zmin[isurf] = 0.0;
    zmax[isurf] = 7.7;
    yminz[isurf] = chsm;
    yminc[isurf] = -chsc;
    ymaxz[isurf] = -chsm;
    ymaxc[isurf] = chsc;
    dnp[isurf] = 0.0;
    xminc[isurf] = -0.035;
    xmaxc[isurf] = 0.035;
    thickAl[isurf] = 0.0;
    thickTa[isurf] = 0.007;

    // Surfaces parallel to the zx plane (28)
    // Collimator slats: 12
    // Collimator Housing Walls, upper 36cm: 8
    // Collimator housing walls, lower 4 cm: 8
    for (isurf = 30; isurf <= 57; isurf++) {
        xnp[isurf] = 0.0;
        ynp[isurf] = 1.0;
        znp[isurf] = 0.0;
        openfrac[isurf] = 0.0;
        xminz[isurf] = 0.0;
        xmaxz[isurf] = 0.0;
        yminz[isurf] = 0.0;
        ymaxz[isurf] = 0.0;
    }

    // Collimator slats
    for (isurf = 30; isurf <= 41; isurf++) {
        thickAl[isurf] = 0.04;
        thickTa[isurf] = 0.007;
        zmin[isurf] = 7.7;
        zmax[isurf] = 47.7;
    }

    for (isurf = 30; isurf <= 32; isurf++) {
        xminc[isurf] = -18.25;
        xmaxc[isurf] = -1.8;
        dnp[isurf] = -18.325 + (isurf - 29)*4.15;
        yminc[isurf] = dnp[isurf] - 0.0235;
        ymaxc[isurf] = dnp[isurf] + 0.0235;
    }

    for (isurf = 33; isurf <= 35; isurf++) {
        xminc[isurf] = -18.25;
        xmaxc[isurf] = -1.8;
        dnp[isurf] = 1.725 + (isurf - 32)*4.15;
        yminc[isurf] = dnp[isurf] - 0.0235;
        ymaxc[isurf] = dnp[isurf] + 0.0235;
    }

    for (isurf = 36; isurf <= 41; isurf++) {
        xminc[isurf] = 1.8;
        xmaxc[isurf] = 18.25;
        dnp[isurf] = dnp[isurf - 6];
        yminc[isurf] = dnp[isurf] - 0.0235;
        ymaxc[isurf] = dnp[isurf] + 0.0235;
    }

    //Walls of Collimator Housing
    for (isurf = 42; isurf <= 49; isurf++) {
        zmin[isurf] = 11.7;
        zmax[isurf] = 47.7;
        thickAl[isurf] = 0.0;
        thickTa[isurf] = 0.007;
    }

    for (isurf = 42; isurf <= 45; isurf++) {
        xminc[isurf] = -18.25;
        xmaxc[isurf] = -1.8;
    }

    for (isurf = 46; isurf <= 49; isurf++) {
        xminc[isurf] = 1.8;
        xmaxc[isurf] = 18.25;
    }

    for (isurf = 42; isurf <= 43; isurf++) {
        dnp[isurf] = -18.25 + 16.45 * (isurf - 42);
        dnp[87 - isurf] = -dnp[isurf];
        dnp[isurf + 4]= dnp[isurf];
        dnp[91 - isurf] = -dnp[isurf];
    }

    for (isurf = 42; isurf <= 49; isurf++) {
        yminc[isurf] = dnp[isurf] - 0.035;
        ymaxc[isurf] = dnp[isurf] + 0.035;
    }

    // Bottom parts
    for (isurf = 50; isurf <= 57; isurf++) {
        zmin[isurf] = 7.7;
        zmax[isurf] = 11.7;
        thickAl[isurf] = 0.0;
        thickTa[isurf] = 0.01;
        dnp[isurf] = dnp[isurf - 8];
        xminc[isurf] = xminc[isurf - 8];
        xmaxc[isurf] = xmaxc[isurf - 8];
        yminc[isurf] = yminc[isurf - 8];
        ymaxc[isurf] = ymaxc[isurf - 8];
    }

    //Slanted surfaces [calibration housing]
    chcos = 0.3081869;
    chsin = 0.9513258;
    chdcom = 19.73473;
    chsm = 0.3239552;
    chsc = 20.744455;
    for (isurf = 58; isurf <= 62; isurf++) {
        zmin[isurf] = 0.0;
        zmax[isurf] = 7.7;
        znp[isurf] = chcos;
        dnp[isurf] = chdcom;
        thickAl[isurf] = 0.15;
    }

    isurf = 58;
    xnp[isurf] = -chsin;
    ynp[isurf] = 0.0;
    thickTa[isurf] = 0.05;
    openfrac[isurf] = 0.5;
    xminz[isurf] = chsm;
    xminc[isurf] = -chsc - 0.1;
    xmaxz[isurf] = chsm;
    xmaxc[isurf] = -chsc + 0.1;
    yminz[isurf] = chsm;
    yminc[isurf] = -chsc;
    ymaxz[isurf] = -chsm;
    ymaxc[isurf] = chsc;

    isurf = 59;
    xnp[isurf] = chsin;
    ynp[isurf] = 0.0;
    thickTa[isurf] = 0.007;
    openfrac[isurf] = 0.0;
    xminz[isurf] = -chsm;
    xminc[isurf] = chsc - 0.1;
    xmaxz[isurf] = -chsm;
    xmaxc[isurf] = chsc + 0.1;
    yminz[isurf] = chsm;
    yminc[isurf] = -chsc;
    ymaxz[isurf] = -chsm;
    ymaxc[isurf] = chsc;

    isurf = 60;
    xnp[isurf] = 0.0;
    ynp[isurf] = -chsin;
    thickTa[isurf] = 0.02;
    openfrac[isurf] = 0.0;
    xminz[isurf] = chsm;
    xminc[isurf] = -chsc;
    xmaxz[isurf] = -chsm;
    xmaxc[isurf] = chsc;
    yminz[isurf] = chsm;
    yminc[isurf] = -chsc - 0.1;
    ymaxz[isurf] = chsm;
    ymaxc[isurf] = -chsc + 0.1;

    isurf = 61;
    xnp[isurf] = 0.0;
    ynp[isurf] = chsin;
    thickTa[isurf] = 0.02;
    openfrac[isurf] = 0.0;
    xminz[isurf] = chsm;
    xminc[isurf] = -chsc;
    xmaxz[isurf] = -chsm;
    xmaxc[isurf] = chsc;
    yminz[isurf] = -chsm;
    yminc[isurf] = chsc - 0.1;
    ymaxz[isurf] = -chsm;
    ymaxc[isurf] = chsc + 0.1;

    // Additional 0.07mm Tantalum shield on surface 58

    isurf = 62;
    xnp[isurf] = -chsin;
    ynp[isurf] = 0.0;
    thickAl[isurf] = 0.0;
    thickTa[isurf] = 0.007;
    openfrac[isurf] = 0.0;
    xminz[isurf] = chsm;
    xminc[isurf] = -chsc - 0.1;
    xmaxz[isurf] = chsm;
    xmaxc[isurf] = -chsc + 0.1;
    yminz[isurf] = chsm;
    yminc[isurf] = -chsc;
    ymaxz[isurf] = -chsm;
    ymaxc[isurf] = chsc;
    
    return status;
}

int CztiCameraGeometry::create_camera_geometry_file(string cameraGeomFileName, string geomTemplate){
    int status=0; //status variable
    fitsfile *fcztgeom;
    string templateFileName = ""; //CZTI geometry fits file template.
    
    templateFileName = template_full_path_generator(geomTemplate);
    if(templateFileName==""){
        LOG(ERROR)<< "Not able to generate Event file template path.";
        LOG(ERROR)<< "Probably Environment Variables are not declared properly.";
        return(EXIT_FAILURE);
    }
    LOG(INFO) << "Template file used to create CZTI geometry CALDB file: " << templateFileName;
    
    fits_create_template(&fcztgeom, (char*) cameraGeomFileName.c_str(), (char*) templateFileName.c_str(), &status);
    if(status){
        LOG(ERROR) << "Error in creating CZTI camera geometry CALDB file from corresponding template.";
        fits_report_error(stderr, status);
        return(EXIT_FAILURE);
    }
    
    fits_close_file(fcztgeom, &status);
    if (status) {
        LOG(ERROR) << "***Error in closing CALDB CZTI Geometry file : " << cameraGeomFileName << "***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    LOG(INFO) << "CZTI camera geometry CALDB file created.";
    return (EXIT_SUCCESS);
}

int CztiCameraGeometry::write_camera_geometry_file(string cameraGeomFileName, string geomTemplate){
    int status=0; //status variable
    int i,j=0; //counter variables
    int colnum=0; //
    fitsfile *fcztgeom;
    
    if(create_camera_geometry_file(cameraGeomFileName, geomTemplate)){
        LOG(ERROR) << "*** Error in creating camera geometry empty file. ***";
        return (EXIT_FAILURE);
    };
    
    fits_open_file(&fcztgeom, (char*) cameraGeomFileName.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening Camera Geometry File: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(fcztgeom, BINARY_TBL, (char*) "SURFACE_DEFINITION", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to SURFACE_DEFINITION HDU in file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing Index
    fits_get_colnum(fcztgeom, CASEINSEN, "Index", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column Index in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_write_col(fcztgeom, TUSHORT, colnum, 1, 1, indexSurf.size(), indexSurf.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing zmin column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Writing zmin
    fits_get_colnum(fcztgeom, CASEINSEN, "zmin", &colnum, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting column number for column zmin in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, zmin.size(), zmin.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing zmin column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Writing zmax
    fits_get_colnum(fcztgeom, CASEINSEN, "zmax", &colnum, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting column number for column zmax in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, zmax.size(), zmax.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing zmax column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing xminz
    fits_get_colnum(fcztgeom, CASEINSEN, "xminz", &colnum, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting column number for column xminz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, xminz.size(), xminz.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing xminz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Writing xmaxz
    fits_get_colnum(fcztgeom, CASEINSEN, "xmaxz", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xmaxz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, xmaxz.size(), xmaxz.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing xmaxz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing xminc
    fits_get_colnum(fcztgeom, CASEINSEN, "xminc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xminc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, xminc.size(), xminc.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing xminc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing xmaxc
    fits_get_colnum(fcztgeom, CASEINSEN, "xmaxc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xmaxc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, xmaxc.size(), xmaxc.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing xmaxc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing yminz
    fits_get_colnum(fcztgeom, CASEINSEN, "yminz", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column yminz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, yminz.size(), yminz.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing yminz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing yminc
    fits_get_colnum(fcztgeom, CASEINSEN, "yminc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column yminc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, yminc.size(), yminc.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing yminc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing ymaxz
    fits_get_colnum(fcztgeom, CASEINSEN, "ymaxz", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column ymaxz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, ymaxz.size(), ymaxz.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing ymaxz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing ymaxc
    fits_get_colnum(fcztgeom, CASEINSEN, "ymaxc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column ymaxc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, ymaxc.size(), ymaxc.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing ymaxc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing xnp
    fits_get_colnum(fcztgeom, CASEINSEN, "xnp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xnp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, xnp.size(), xnp.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing xnp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing ynp
    fits_get_colnum(fcztgeom, CASEINSEN, "ynp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column ynp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, ynp.size(), ynp.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing ynp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing znp
    fits_get_colnum(fcztgeom, CASEINSEN, "znp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column znp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, znp.size(), znp.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing znp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing dnp
    fits_get_colnum(fcztgeom, CASEINSEN, "dnp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column dnp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, dnp.size(), dnp.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing dnp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing thickAl
    fits_get_colnum(fcztgeom, CASEINSEN, "thickAl", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column thickAl in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, thickAl.size(), thickAl.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing thickAl column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing thickTa
    fits_get_colnum(fcztgeom, CASEINSEN, "thickTa", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column thickTa in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, thickTa.size(), thickTa.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing thickTa column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing openfrac
    fits_get_colnum(fcztgeom, CASEINSEN, "openfrac", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column openfrac in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fcztgeom, TFLOAT, colnum, 1, 1, openfrac.size(), openfrac.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing openfrac column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Deleting the first row which contains 0 index
    fits_delete_rows(fcztgeom, 1, 1, &status);
    if (status) {
        LOG(ERROR) << "Error in deleting first row in SURFACE_DEFINITION extension of "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fcztgeom, &status);
    if (status) {
        LOG(ERROR) << "***Error in closing file: " << cameraGeomFileName << "***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int CztiCameraGeometry::read_camera_geometry_file(string cameraGeomFileName){
    int status=0;
    int i,j=0; //counter variables
    int colnum=0;
    long nrows=0;
    unsigned int* indexArray;
    float *zminArray, *zmaxArray;
    float *xminzArray, *xmincArray, *xmaxzArray, *xmaxcArray;
    float *yminzArray, *ymincArray, *ymaxzArray, *ymaxcArray;
    float *xnpArray, *ynpArray, *znpArray, *dnpArray;
    float *thickAlArray, *thickTaArray, *openfracArray;
    fitsfile *fcztgeom;
    
    fits_open_file(&fcztgeom, (char*) cameraGeomFileName.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(fcztgeom, BINARY_TBL, (char*) "SURFACE_DEFINITION", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to SURFACE_DEFINITION HDU in file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_get_num_rows(fcztgeom, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows in SURFACE_DEFINITION extension of file : "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //initializing all arrays
    indexArray = new unsigned int[nrows];
    zminArray = new float[nrows];
    zmaxArray = new float[nrows];
    xminzArray = new float[nrows];
    xmincArray = new float[nrows];
    xmaxzArray = new float[nrows];
    xmaxcArray = new float[nrows];
    yminzArray = new float[nrows];
    ymincArray = new float[nrows];
    ymaxzArray = new float[nrows];
    ymaxcArray = new float[nrows];
    xnpArray = new float[nrows];
    ynpArray = new float[nrows];
    znpArray = new float[nrows];
    dnpArray = new float[nrows];
    thickAlArray = new float[nrows];
    thickTaArray = new float[nrows];
    openfracArray = new float[nrows];
    
    for(i=0; i<nrows; i++){
        indexArray[i]=0;
        zminArray[i]=0.0;
        zmaxArray[i]=0.0;
        xminzArray[i]=0.0;
        xmincArray[i]=0.0;
        xmaxzArray[i]=0.0;
        xmaxcArray[i]=0.0;
        yminzArray[i]=0.0;
        ymincArray[i]=0.0;
        ymaxzArray[i]=0.0;
        ymaxcArray[i]=0.0;
        xnpArray[i]=0.0;
        ynpArray[i]=0.0;
        znpArray[i]=0.0;
        dnpArray[i]=0.0;
        thickAlArray[i]=0.0;
        thickTaArray[i]=0.0;
        openfracArray[i]=0.0;
    }
    
    //Reading Index
    fits_get_colnum(fcztgeom, CASEINSEN, "Index", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column Index in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_read_col(fcztgeom, TUSHORT, colnum, 1, 1, nrows, NULL, indexArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading Index column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Reading zmin
    fits_get_colnum(fcztgeom, CASEINSEN, "zmin", &colnum, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting column number for column zmin in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, zminArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading zmin column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Reading zmax
    fits_get_colnum(fcztgeom, CASEINSEN, "zmax", &colnum, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting column number for column zmax in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, zmaxArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading zmax column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading xminz
    fits_get_colnum(fcztgeom, CASEINSEN, "xminz", &colnum, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting column number for column xminz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, xminzArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading xminz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Reading xmaxz
    fits_get_colnum(fcztgeom, CASEINSEN, "xmaxz", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xmaxz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, xmaxzArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading xmaxz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading xminc
    fits_get_colnum(fcztgeom, CASEINSEN, "xminc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xminc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, xmincArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading xminc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading xmaxc
    fits_get_colnum(fcztgeom, CASEINSEN, "xmaxc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xmaxc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, xmaxcArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading xmaxc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading yminz
    fits_get_colnum(fcztgeom, CASEINSEN, "yminz", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column yminz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, yminzArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading yminz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading yminc
    fits_get_colnum(fcztgeom, CASEINSEN, "yminc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column yminc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, ymincArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading yminc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading ymaxz
    fits_get_colnum(fcztgeom, CASEINSEN, "ymaxz", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column ymaxz in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, ymaxzArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading ymaxz column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading ymaxc
    fits_get_colnum(fcztgeom, CASEINSEN, "ymaxc", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column ymaxc in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, ymaxcArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading ymaxc column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading xnp
    fits_get_colnum(fcztgeom, CASEINSEN, "xnp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column xnp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, xnpArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading xnp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading ynp
    fits_get_colnum(fcztgeom, CASEINSEN, "ynp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column ynp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, ynpArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading ynp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading znp
    fits_get_colnum(fcztgeom, CASEINSEN, "znp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column znp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, znpArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading znp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading dnp
    fits_get_colnum(fcztgeom, CASEINSEN, "dnp", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column dnp in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, dnpArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading dnp column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading thickAl
    fits_get_colnum(fcztgeom, CASEINSEN, "thickAl", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column thickAl in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, thickAlArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading thickAl column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading thickTa
    fits_get_colnum(fcztgeom, CASEINSEN, "thickTa", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column thickTa in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, thickTaArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading thickTa column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading openfrac
    fits_get_colnum(fcztgeom, CASEINSEN, "openfrac", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for column openfrac in SURFACE_DEFINITION extension of"
                << " file: " << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_col(fcztgeom, TFLOAT, colnum, 1, 1, nrows, NULL, openfracArray, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading openfrac column in SURFACE_DEFINITION extension of file: "
                << cameraGeomFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_close_file(fcztgeom, &status);
    if (status) {
        LOG(ERROR) << "***Error in closing file: " << cameraGeomFileName << "***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //clearing class variables
    indexSurf.clear();
    zmin.clear();
    zmax.clear();
    xminz.clear();
    xmaxz.clear();
    xminc.clear();
    xmaxc.clear();
    yminz.clear();
    ymaxz.clear();
    yminc.clear();
    ymaxc.clear();
    xnp.clear();
    ynp.clear();
    znp.clear();
    dnp.clear();
    thickAl.clear();
    thickTa.clear();
    openfrac.clear();
    
    //Assigning data read to class variables
    for(i=0; i<nrows; i++){
        indexSurf.push_back(indexArray[i]);
        zmin.push_back(zminArray[i]);
        zmax.push_back(zmaxArray[i]);
        xminz.push_back(xminzArray[i]);
        xminc.push_back(xmincArray[i]);
        xmaxc.push_back(xmaxcArray[i]);
        xmaxz.push_back(xmaxzArray[i]);
        yminz.push_back(yminzArray[i]);
        yminc.push_back(ymincArray[i]);
        ymaxc.push_back(ymaxcArray[i]);
        ymaxz.push_back(ymaxzArray[i]);
        xnp.push_back(xnpArray[i]);
        ynp.push_back(ynpArray[i]);
        znp.push_back(znpArray[i]);
        dnp.push_back(dnpArray[i]);
        thickAl.push_back(thickAlArray[i]);
        thickTa.push_back(thickTaArray[i]);
        openfrac.push_back(openfracArray[i]);
    }
    
    delete[] zminArray, zmaxArray;
    delete[] xmincArray, xmaxcArray, xminzArray, xmaxzArray;
    delete[] ymincArray, ymaxcArray, yminzArray, ymaxzArray;
    delete[] xnpArray, ynpArray, znpArray, dnpArray;
    delete[] thickAlArray, thickTaArray;
    delete[] openfracArray;
    return (EXIT_SUCCESS);
    
}

/*****************STAR CATALOG CLASS **************************/

//CATALOG STRUCT
int Catalog::set_catalog_values(string srcname, string cprtname, float ra, float dec, float thetaxr,
        float thetayr, float snr, float ctptra, float ctptdec, float flux, float fluxlo, float fluxhi) {
    sourceNameBAT.push_back(srcname);
    counterpartName.push_back(cprtname);
    RA.push_back(ra);
    DEC.push_back(dec);
    this->thetaxr.push_back(thetaxr);
    this->thetayr.push_back(thetayr);
    SNR.push_back(snr);
    ctptRA.push_back(ctptra);
    ctptDEC.push_back(ctptdec);
    Flux.push_back(flux);
    FluxLO.push_back(fluxlo);
    FluxHI.push_back(fluxhi);
    
    return EXIT_SUCCESS;
}

//CATALOG STRUCT END
int StarCatalog::find_sources_in_FOV(Catalog &catalog, double RA, double DEC, double TWIST, 
                double thetaxr, double thetayr){
    int status=0;
    long nsources=0; //number of sources in current catalog;
    long isource=0; //index variable to iterate over each source;
    double sourceThetax=0.0;
    double sourceThetay=0.0;
    
    nsources = this->catalog.sourceNameBAT.size();
    LOG(INFO) << "Source lying in FOV (thetaxr, thetayr) " << "(" << thetaxr << "," << thetayr << ")";
    for(isource=0; isource<nsources; isource++){
        to_thetaX_thetaY(RA, DEC, TWIST, (double) (this->catalog.RA[isource] * TORAD), 
                (double) (this->catalog.DEC[isource] * TORAD),
                sourceThetax, sourceThetay);
        if(sourceThetax>=(-thetaxr) && sourceThetax<=thetaxr && sourceThetay>=(-thetayr) && sourceThetay<=thetayr){
            LOG(INFO) << "SOURCE: " << this->catalog.sourceNameBAT[isource] <<  
                    " THETAX(rad): " << sourceThetax << ", THETAY(rad): " << sourceThetay;
            catalog.set_catalog_values(this->catalog.sourceNameBAT[isource],
                    this->catalog.counterpartName[isource],
                    this->catalog.RA[isource],
                    this->catalog.DEC[isource], 
                    sourceThetax, sourceThetay,
                    this->catalog.SNR[isource], 
                    this->catalog.ctptRA[isource],
                    this->catalog.ctptDEC[isource],
                    this->catalog.Flux[isource], 
                    this->catalog.FluxLO[isource], 
                    this->catalog.FluxHI[isource]);
        }
    }
    //&& sourceThetax<=thetaxr && sourceThetay>=(-thetayr) && sourceThetay<=-thetayr
    
    
    return status;
}
int StarCatalog::read_catalog_file(string catalogFilename, int catalogExtnum){
    int status=0;
    int hduType=0;
    fitsfile *fcatalog;
    long nrows=0; //number of rows in catalog extension
    
    fits_open_file(&fcatalog, (char*) catalogFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening catalog file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_movabs_hdu(fcatalog, catalogExtnum, &hduType, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU number " << catalogExtnum << 
                " of catalog file " << catalogFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fcatalog, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows in catalog file " << catalogFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading Fits columns
    if(read_fits_string_column(fcatalog, "BAT_NAME", TSTRING, 1, 1, nrows, catalog.sourceNameBAT)){
        LOG(ERROR) << "Error in reading BAT_NAME column of catalog file: " << catalogFilename;
    }
    if(read_fits_string_column(fcatalog, "COUNTERPART_NAME", TSTRING, 1, 1, nrows, catalog.counterpartName)){
        LOG(ERROR) << "Error in reading COUNTERPART_NAME column of catalog file: " << catalogFilename;
    }
    if(read_fits_column(fcatalog, "RA", TFLOAT, 1, 1, nrows, catalog.RA)){
        LOG(ERROR) << "Error in reading RA column of catalog file: " << catalogFilename;
    }
    if(read_fits_column(fcatalog, "DEC", TFLOAT, 1, 1, nrows, catalog.DEC)){
        LOG(ERROR) << "Error in reading DEC column of catalog file: " << catalogFilename;
    }
    if(read_fits_column(fcatalog, "SNR", TFLOAT, 1, 1, nrows, catalog.SNR)){
        LOG(ERROR) << "Error in reading SNR column of catalog file: " << catalogFilename;
    }
    if(read_fits_column(fcatalog, "CTPT_RA", TFLOAT, 1, 1, nrows, catalog.ctptRA)){
        LOG(ERROR) << "Error in reading CTPT_RA column of catalog file: " << catalogFilename;
    }
    if(read_fits_column(fcatalog, "CTPT_DEC", TFLOAT, 1, 1, nrows, catalog.ctptDEC)){
        LOG(ERROR) << "Error in reading CTPT_DEC column of catalog file: " << catalogFilename;
    }
    if (read_fits_column(fcatalog, "FLUX", TFLOAT, 1, 1, nrows, catalog.Flux)) {
        LOG(ERROR) << "Error in reading FLUX column of catalog file: " << catalogFilename;
    }
    if (read_fits_column(fcatalog, "FLUX_LO", TFLOAT, 1, 1, nrows, catalog.FluxLO)) {
        LOG(ERROR) << "Error in reading FLUX_LO column of catalog file: " << catalogFilename;
    }
    if (read_fits_column(fcatalog, "FLUX_HI", TFLOAT, 1, 1, nrows, catalog.FluxHI)) {
        LOG(ERROR) << "Error in reading FLUX_HI column of catalog file: " << catalogFilename;
    }
    
    //closing fits file
    fits_close_file(fcatalog, &status);
    if (status) {
        LOG(ERROR) << "Error in closing fits file " << catalogExtnum;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return status;
}


/*****************STAR CATALOG CLASS END**********************/
/*******************BADPIXEL CLASS**********************/
// This class has been defined in badpixCALDB.h file.
/*******************************************************/



/***********INDEPENDENT FUNCTIONS******************/

string caldb_full_path_generator(string filename){
    string baseCALDBdirPath="";
    string CALDBfileFullPath="";
    try{
        string baseCALDBdirPath = getenv("CALDBnew");

        string CALDBfileFullPath = baseCALDBdirPath + "/bcf/" + filename;
        return CALDBfileFullPath;
    }
    catch (std::exception& e) {
        LOG(ERROR) << "Error: " << e.what();
        return CALDBfileFullPath="";
    }

}

string caldb_filename_generator(string content, string version, string creationTime, string extname, string quadNo){
    string outFilename;
    outFilename = "AS1czt"+quadNo+content+creationTime+version+"."+extname;
    return outFilename;
}

int update_CALDB_keywords(char* filename, char* creator){
    int status=0, numhdu;
    fitsfile *fptr;
    char utcDate[100];
    char utcTime[100];
    
    //Getting current time and formatting it.
    time_t now = time(0);
    tm* gmtm = gmtime(&now);
    sprintf(utcDate, "%d-%02d-%02d", gmtm->tm_year + 1900, gmtm->tm_mon + 1, gmtm->tm_mday);
    sprintf(utcTime, "%02d:%02d:%02d", gmtm->tm_hour, gmtm->tm_min + 1, gmtm->tm_sec);
    
    fits_open_file(&fptr, filename, READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in update_CALDB_keywords()";
        fits_report_error(stderr, status);
        return status;
    }
    fits_get_num_hdus(fptr, &numhdu, &status);
    if (status) {
        LOG(ERROR) << "Error in update_CALDB_keywords()";
        fits_report_error(stderr, status);
        return status;
    }
    for (int i = 1; i <= numhdu; i++) {
        fits_movabs_hdu(fptr, i, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in update_CALDB_keywords()";
            fits_report_error(stderr, status);
            return status;
        }
        fits_write_date(fptr, &status);
        fits_update_key(fptr, TSTRING, "ORIGIN", (char *) ORIGIN, NULL, &status);
        fits_update_key(fptr, TSTRING, "CREATOR", creator, NULL, &status);
        fits_update_key(fptr, TSTRING, "FILENAME", filename, NULL, &status);
        fits_update_key(fptr, TSTRING, "CVSD0001", utcDate, NULL, &status);
        fits_update_key(fptr, TSTRING, "CVST0001", utcTime, NULL, &status);
        fits_write_chksum(fptr, &status);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in update_CALDB_keywords()";
        fits_report_error(stderr, status);
        return status;
    }
    return (EXIT_SUCCESS);    
    
}



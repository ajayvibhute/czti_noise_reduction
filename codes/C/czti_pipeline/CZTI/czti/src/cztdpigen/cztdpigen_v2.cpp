#include "cztdpigen_v2.h"
#include "level2validation.h"
#include "coordinateTransformation.h"

Cztdpigen::Cztdpigen(){
    strcpy(modulename,"cztdpigen_v");
    strcat(modulename,VERSION);
    strcpy(quadsToProcess,"");
    nebins=0;    //will be set in cztdpigen_process
    emin=NULL;  //will be set in cztdpigen_process
    emax=NULL;  //will be set in cztdpigen_process
    ntbins=0;
    ndpi=0;
    tstart=NULL;
    tstop=NULL;
}

int Cztdpigen::read(int argc,char **argv){
    int status=0;
   
    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error Initializing PIL***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("par_infile",infile))){
        LOG(ERROR)<<"***Error Reading input file:"<<infile<<"***";
        return status;
    }
    
    /*if(PIL_OK!=(status=PILGetFname("par_effAreafile",effAreafile))){
        LOG(ERROR)<<"***Error Reading Effective area CALDB file:"<<effAreafile<<"***";
        return status;
    }*/

    if(PIL_OK!=(status=PILGetFname("par_badpixFile",badpixFile))){
        LOG(ERROR)<<"***Error Reading Badpixel file:"<<badpixFile<<"***";
        return status;
    }

    if (PIL_OK != (status = PILGetInt("par_badpixThreshold", &badpixThreshold))) {
        LOG(ERROR) << "***Error reading badpix threshold parameter.***";
        return status;
    }
    
    if (PIL_OK != (status = PILGetFname("par_outDPHfile", outDPHfile))) {
        LOG(ERROR) << "***Error Reading outDPHfile:" << outDPHfile << "***";
        return status;
    }
    
    if (PIL_OK != (status = PILGetFname("par_outDPIfile", outDPIfile))) {
        LOG(ERROR) << "***Error Reading outDPIfile:" << outDPIfile << "***";
        return status;
    }

    do {
        if (PIL_OK != (status = PILGetString("par_quadsToProcess", quadsToProcess))) {
            LOG(ERROR) << "***Error reading quadrant:" << quadsToProcess << "***";
            return status;
        }
    } while (isQuadrantNumbersValid(quadsToProcess));

    /*if (PIL_OK != (status = PILGetFname("par_caldbdir", caldbdir))) {
        LOG(ERROR) << "***Error Reading CALDB directory:" << caldbdir << "***";
        return status;
    }*/
    
    do {

        if (PIL_OK != (status = PILGetString("par_ebins", ebins))) {
            LOG(ERROR) << "***Error Reading ebins:" << ebins << "***";
            return status;
        }
    } while (isBinsValid(ebins));

    do {
        if (PIL_OK != (status = PILGetString("par_timerange", timerange))) {
            LOG(ERROR) << "***Error Reading timerange***";
            return status;
        }
    } while (isBinsValid(timerange));
    
    
    if(PIL_OK!=(status=PILGetBool("par_history",&history))){
        LOG(ERROR)<<"***Error Reading history parameter***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetBool("par_clobber",&clobber))){
        LOG(ERROR)<<"***Error Reading clobber***";
        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

int Cztdpigen::read(char* infile, char* badpixelFile, char* outDPHfile, 
        char* outDPIfile, char *quadsToProcess, int badpixThreshold,
        char* ebins, char *timerange,
        int clobber, int history){
    strcpy(this->infile,infile);
    //strcpy(this->effAreafile,effAreafile);
    strcpy(this->badpixFile, badpixelFile);
    strcpy(this->outDPHfile,outDPHfile);
    strcpy(this->outDPIfile,outDPIfile);
    strcpy(this->quadsToProcess, quadsToProcess);
    strcpy(this->ebins,ebins);
    strcpy(this->timerange,timerange);
    //strcpy(this->caldbdir, caldbdir);
    this->badpixThreshold = badpixThreshold;
    this->clobber=clobber;
    this->history=history;
    return (EXIT_SUCCESS);
}

void Cztdpigen::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                          CZTDPIGEN PARAMETERS                            ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
	LOG(INFO)<<"Modulename                  : "<<modulename;
    LOG(INFO)<<"Input event file            : "<<infile;
    LOG(INFO)<<"Input badpixel file         : "<<badpixFile;
    LOG(INFO)<<"Badpixel Threshold          : "<<badpixThreshold;
    LOG(INFO)<<"Output DPH file             : "<<outDPHfile;
    LOG(INFO)<<"Output DPI file             : "<<outDPIfile;
    //LOG(INFO)<<"Effective Area File       :"<<effAreafile;
    LOG(INFO)<<"Quadrants to be processed   : "<<quadsToProcess;
    LOG(INFO)<<"Time Range                  : "<<timerange;
    LOG(INFO)<<"Energy Bins                 : "<<ebins;
    //LOG(INFO)<<"CALDB directroy path:     :"<<caldbdir;
    LOG(INFO)<<"Clobber                     : "<<clobber;
    LOG(INFO)<<"History                     : "<<history;
    LOG(INFO)<<"----------------------------------------------------------------------------";
}

int Cztdpigen::cztdpigenProcess(){
	/*Added by ajay vibhute, Dec 9 2015 5:27*/
    int timefilter=0,energyfilter=1;
    int status=0; //status variable
    DPI dpih;
    int i,j,k=0; //counter variable
    fitsfile *fptrIn; //pointer to input event file
    fitsfile *fptrDpiOut; //pointer to DPI output file
    long nrows=0; 
    int quadID=0;
    long** dpi;
    long* fulldpi;
    int outputFileIndex=0;
    vector<string> vec_history;
    char extname[FLEN_VALUE];
    string fullCALDBfilePath; //to store full caldb file path for any caldb file
    vector<string> outDPIfilenames; // string vector containing name of output DPI files.
    vector<string> outDPHfilenames;
    EffArea effectiveArea;
    string DPItemplate="dpiTemplate";
    string DPHtemplate="dphTemplate";
    
    int* quadArray; //Array containing quadrant numbers to be processed.
    int tempNoOcuurance=0; // temporary variable to store number of commas in csv input
    int no_quads=0; // stores number of quadrants that has to be processed.
    int quadid=0; //ID of quadrant for which DPI has to be made.
    
    //VERIFYING INPUT EVENT FILE;
    FitsFileInfo fInInfo; //to store info of any opened file.
    int numhdus=0; //to validate any fits file.
    vector<string> vec_hdunames;
    vector<int> vec_hdutypes;
    vector<int> vec_ncols;
    vector<long> vec_nrows;
    vector<int> quadsToProcessVec;

    double head_tstart,head_tstop;    
    //Opening input event file and checking its validity.
    fits_open_file(&fptrIn,infile,READONLY,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"Error in opening energy added input event file:"<<infile;
        return (EXIT_FAILURE);
    }    

    fits_read_key(fptrIn, TDOUBLE, "TSTART", &head_tstart, NULL, &status);
    report_error(status, (string) "Error in reading keyword TSTART");

    fits_read_key(fptrIn, TDOUBLE, "TSTOP", &head_tstop, NULL, &status);
    report_error(status, (string) "Error in reading keyword TSTOP");

    
//    fInInfo.get_fitsinfo(fptrIn); //storing info of input event file in variable fInInfo
//    //Validating input event file
//    numhdus=8;
//    string hdunames_array[]={"Primary","Q0","Q1","Q2","Q3","VETOSPECTRUM","SSM Data", "TEMP"};
//    int hdutypes_array[]={IMAGE_HDU, BINARY_TBL, BINARY_TBL, BINARY_TBL, BINARY_TBL, BINARY_TBL, BINARY_TBL, BINARY_TBL};
//    int ncols_array[]={0,12,12,12,12,3,7,4};    
//    
//    for(i=0; i<numhdus; i++){
//        vec_hdunames.push_back(hdunames_array[i]);
//        vec_hdutypes.push_back(hdutypes_array[i]);
//        vec_ncols.push_back(ncols_array[i]);  
//    }
//    if (fInInfo.validate_fitsfile(numhdus, vec_hdunames, vec_hdutypes, &status, vec_ncols, vec_nrows)) {
//        LOG(ERROR) << "Input Event is not in proper format.";
//        //return (EXIT_FAILURE);
//    }
//    //EVENT FILE VALIDATED

    //Finding number of DPIs to be created
    nebins = get_nbins(ebins);
    ntbins = 1;//get_nbins(timerange);
    ndpi = nebins*ntbins; //number of DPI files to be created.
    if(ndpi<=0) { 
        LOG(ERROR)<<"Number of DPI to be generated is zero/less than zero";
        return (EXIT_FAILURE);
    }
    // Number of DPIs to be created found and is equal to ndpi
    
    // Getting time and energy ranges
    tstart = new double[5];
    tstop = new double[5];
    emin = new double[nebins];
    emax = new double[nebins];
    
    try{
        quadsToProcessVec = get_quadsToProcess(quadsToProcess);
    } catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }
    
    if(getTime(fptrIn,timerange,ntbins,tstart,tstop, quadsToProcessVec)){
        LOG(ERROR)<<"***Error in preparing time bins***";
        return (EXIT_FAILURE);
    }


	/*Added by ajay Vibhute, Dec 9 2015*/
    if(strcmp(timerange,"-")==0)
    {
	    timefilter=0;
	    LOG(INFO)<<"Not applying any time filter";
    }
    else
	    timefilter=1;

/*    
    string caldbBasepath = (string) caldbdir;
    LOG(INFO)<<caldbdir;
*/
    string eboundsFile;	
    int eboundsExtnum;	
		
    if(QueryCaldb("ASTROSAT","CZTI","-","EBOUNDS",head_tstart,head_tstop,eboundsFile,eboundsExtnum))
    {
        LOG(ERROR) << "Not able to get CALDB EBOUNDS file";
        return (EXIT_FAILURE);
    }




    if(getEnergy(ebins,nebins,emin,emax, eboundsFile)){
        LOG(ERROR)<<"***Error in preparing energy bins***";
        return (EXIT_FAILURE);
    }


    //printing number of DPIs to be created & time and energy ranges
    LOG(INFO) << "Number of DPI files to be created: " << ndpi;
    LOG(INFO)<<"Time Bins specified by user:";
    for(i=0;i<ntbins;i++){
        LOG(INFO)<<setprecision(15)<<tstart[i]<<"-"<<setprecision(15)<<tstop[i];
    }
    LOG(INFO)<<"Energy Bins specified by user:";
    for(i=0;i<nebins;i++){
        LOG(INFO)<<emin[i]<<"-"<<emax[i]; 
    }
    
    // Generating Output DPI File Names
    status=getNoutputfiles(ndpi, outDPIfile, outDPHfile, outDPIfilenames, outDPHfilenames);
    if(status){
        LOG(ERROR) << "Error in generating names of output DPI files.";
        return (EXIT_FAILURE);
    }
    LOG(INFO)<<"Output DPI & DPH files to be generated----";
    for(i=0;i<outDPIfilenames.size();i++){
        //LOG(INFO)<<i+1<<" : "<<outDPIfilenames[i] << ", " << outDPHfilenames[i];
        LOG(INFO)<<"File "<<i+1<<" : "<<outDPIfilenames[i];
        LOG(INFO)<<"File "<<i+1<<" : "<<outDPHfilenames[i];
    }
    DLOG(INFO)<<"--------------------------------------------------";    
    // Output DPI file names generated and stored in string vector outDPIfilenames.
    
    // Checking whether output DPI files exist or not
    // If yes then deletes it (for clobber=yes)
    // Otherwise raises an error
    for(i=0; i<outDPIfilenames.size(); i++) {
        if (FileExists((char*) outDPIfilenames[i].c_str())) {
            if (clobber == YES) {
                if (unlink((char*) outDPIfilenames[i].c_str()) != 0) {
                    LOG(ERROR) << "Error in deleting Output DPI File: " << outDPIfilenames[i];
                }
            } else {
                LOG(INFO) << outDPIfilenames[i] << " already exists.";
                LOG(INFO) << "Use clobber=yes for overwriting the file.";
                return (EXIT_FAILURE);
            }
        }
        if (FileExists((char*) outDPHfilenames[i].c_str())) {
            if (clobber == YES) {
                if (unlink((char*) outDPHfilenames[i].c_str()) != 0) {
                    LOG(ERROR) << "Error in deleting Output DPI File: " << outDPHfilenames[i];
                }
            } else {
                LOG(INFO) << outDPHfilenames[i] << " already exists.";
                LOG(INFO) << "Use clobber=yes for overwriting the file.";
                return (EXIT_FAILURE);
            }
        }
    }
    // Output DPI files existence check finished.    

    //Generating Quadrant array containing information regarding quadrants to be processed
    stringFinder(quadsToProcess, ",", 0, &tempNoOcuurance);
    no_quads = tempNoOcuurance + 1;
    quadArray = new int[no_quads];

    if (csvStrToInt(quadsToProcess, ",", quadArray, &no_quads)) {
        LOG(ERROR) << "***Error in converting quadrant string array into integer array***";
    }    
    //Quadrant Array generated.
    
    
    //Generating DPI from event file
    for(i = 0; i<ntbins; i++) {
        for (j = 0; j < nebins; j++) {
            dpih.set_badpixfile((string) badpixFile);
            dpih.set_badpixThreshold(badpixThreshold);
            //opening output fits files
            if (dpih.create_dpi_file((string) outDPHfilenames[outputFileIndex], DPHtemplate, infile)) {
                LOG(ERROR) << "Error in creating DPI file " << outDPHfilenames[outputFileIndex];
                return EXIT_FAILURE;
            }
            if(dpih.create_dpi_file((string) outDPIfilenames[outputFileIndex], DPItemplate, infile)){
                LOG(ERROR)<<"Error in creating DPI file " << outDPIfilenames[outputFileIndex];
                return EXIT_FAILURE;
            }

            LOG(INFO) << "Writing DPI/DPH to file " << outDPIfilenames[outputFileIndex];
            for (k = 0; k < no_quads; k++) {
                quadid = quadArray[k];
                quadToHDU(quadid, extname);
     		/*Modified by ajay Vibhute, Dec 9 2015 5:28*/
			/*Again by Mithun on 11 Jan 2016*/           

        fits_movabs_hdu(fptrIn, k+2, NULL, &status);       
		//dpih.headerParam.quadid=k;	
    	dpih.headerParam.readFromHeader(fptrIn);
		LOG(INFO)<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUAD ID"<<k;
		if(dpih.create_dph2d(infile, quadid, emin[j], emax[j], tstart[k], tstop[k], timefilter,energyfilter)) {
                    LOG(ERROR) << "Error in DPH generation for " << emin[j] << "-" <<
                            emax[j] << " and " << tstart[i] << "-" << tstop[i];
                    return EXIT_FAILURE;
                }
		


		//Get caldb eff_area file
		char quadnam[20];
		int effareaExtnum;
		sprintf(quadnam,"QUADRANT%d",k);

		if(QueryCaldb("ASTROSAT","CZTI",(string)quadnam,"EFF_AREA",head_tstart,head_tstop,effAreafile,effareaExtnum))
    		{
        		LOG(ERROR) << "Not able to get CALDB EFF_AREA file";
        		return (EXIT_FAILURE);
    		}

                if(dpih.create_dpi2d(effAreafile, quadid)){
                    LOG(ERROR) << "Error in DPI generation for " << emin[j] << "-" <<
                            emax[j] << " and " << tstart[i] << "-" << tstop[i]; 
                    return EXIT_FAILURE;
                }
                if (dpih.write_dpih_quad(outDPHfilenames[outputFileIndex], "DPH", extname)) {
                    LOG(ERROR) << "Error in writing DPH file";
                    return (EXIT_FAILURE);
                }
                if (dpih.write_dpih_quad(outDPIfilenames[outputFileIndex], "DPI", extname)) {
                    LOG(ERROR) << "Error in writing DPH file";
                    return (EXIT_FAILURE);
                }
                if(dpih.add_quad_to_fulldpih(quadid, "DPH")){
                    LOG(ERROR) << "Error in placing Quadrant DPH into FULL DPH";
                    return EXIT_FAILURE;
                }
                if(dpih.add_quad_to_fulldpih(quadid, "DPI")){
                    LOG(ERROR) << "Error in placing Quadrant DPI into FULL DPI";
                    return EXIT_FAILURE;
                }
            }
            if (dpih.write_dpih_full((string) outDPHfilenames[outputFileIndex], "DPH", "FULL_DPH")) {
                LOG(ERROR) << "Error in writing FULL DPH";
                return EXIT_FAILURE;
            }
            if (dpih.write_dpih_full((string) outDPIfilenames[outputFileIndex], "DPI", "FULL_DPI")) {
                LOG(ERROR) << "Error in writing FULL DPH";
                return EXIT_FAILURE;
            }

            outputFileIndex++;
        }
    }
    delete[] tstart, tstop, emin, emax, quadArray;
    fits_close_file(fptrIn, &status);
    if (status) {
        LOG(ERROR) << "***Error closing input file***";
        return (EXIT_FAILURE);
    }

    //updating keywords and writing history
    char DPIfilename[FLEN_FILENAME];
    char DPHfilename[FLEN_FILENAME];
    for (int i = 0; i < outDPIfilenames.size(); i++) {
        strcpy(DPIfilename, outDPIfilenames[i].c_str());
        strcpy(DPHfilename, outDPHfilenames[i].c_str());
        if (history == YES) {
            vector<string> vhistory;
            getHistory(vhistory);
            writeHistory(DPIfilename, vhistory);
            writeHistory(DPHfilename, vhistory);
        }
        updateKeywords(DPIfilename, modulename);
	updateKeywords(DPHfilename, modulename);

    }   
    
    return status;
    
}


int getNoutputfiles(int ndpi, char* outDPIfile, char* outDPHfile, vector<string>& outputDPIFileNames, 
        vector<string>& outputDPHFileNames){
    outputDPIFileNames.clear();
    outputDPHFileNames.clear();
    int i=0; //status variable
    int status=0;
    char DPIfilename[FLEN_FILENAME];
    char DPHfilename[FLEN_FILENAME];
    if(ndpi<=0) { return (EXIT_FAILURE); }
    else if(ndpi==1){
        outputDPIFileNames.push_back((string(outDPIfile)));
        outputDPHFileNames.push_back((string(outDPHfile)));
    }
    else {
        for(i=0;i<ndpi;i++){
            sprintf(DPIfilename,"%s_%d",outDPIfile,i);
            sprintf(DPHfilename,"%s_%d",outDPHfile,i);
            outputDPIFileNames.push_back((string)DPIfilename);
            outputDPHFileNames.push_back((string)DPHfilename);
        }
    }
    
    return (EXIT_SUCCESS);
}

int getTime(fitsfile *fptr,char *timerange,int ntbins,double *tstart,double *tstop,
        vector <int> quadsToProcess){
    int status=0;
    EventFileHandler evt;
    if(timerange==NULL){
        LOG(ERROR)<<"Time Range value not correct - getEventTime";
        return (EXIT_FAILURE);
    }
    else if(strcmp(timerange,"-")==0){   //get time from file
        try{
            evt.get_common_time(fptr, quadsToProcess, tstart, tstop);
        } catch (ErrorHandler errHandler){
            logError(errHandler);
            return EXIT_FAILURE;
        }
    }
    else{
         getbins(timerange,tstart,tstop,ntbins);
    }
     return (EXIT_SUCCESS);
}

int getEnergy(char *ebins,int nebins,double *emin,double *emax, string eboundsFile){


/*	
    CaldbHandler caldbhandler;
    string eboundFilePath="";
    
    try{
        caldbhandler.generate_caldb_filepaths(CALDBbasepath);
        eboundFilePath = caldbhandler.get_caldb_file_path("EBOUNDS");
        
    } catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }
*/
	
    Ebounds ebounds;
    int nchannels;   

    if(strcmp(ebins,"-")==0){
        if(ebounds.read_ebounds_file(eboundsFile)){
            LOG(ERROR) << "Error in reading Ebounds file.";
            return EXIT_FAILURE;
        }
        nchannels = ebounds.get_nrows();
        try{
        *emin = ebounds.find_energy_from_channel(0);
        *emax = ebounds.find_energy_from_channel(nchannels-1);
        } catch(string errorMsg){
            LOG(ERROR) << errorMsg;
            return EXIT_FAILURE;
        }
    }
    else{                    //if some ranges are provided
        getbins(ebins,emin,emax,nebins);
    }
    return (EXIT_SUCCESS);
}

int Cztdpigen::getHistory(vector<string> &vhistory){
    //char *user=getlogin();
    strcpy(modulename,"cztdpigen_v");
    strcat(modulename,VERSION);

    char *user = getenv("USER");
	string str="Module run by "+(string)user;
    vhistory.push_back(str);
    vhistory.push_back("Parameter List START for "+(string)modulename);
    vhistory.push_back("P1 infile="+(string)infile);
    vhistory.push_back("P2 effAreafile="+(string)effAreafile);
    vhistory.push_back("P3 badpixelFile="+(string)badpixFile);
    vhistory.push_back("P4 outDPHfile=" + (string) outDPHfile);
    vhistory.push_back("P5 outDPIfile="+(string)outDPIfile);
    vhistory.push_back("P6 quadsToProcess="+(string)quadsToProcess);
    vhistory.push_back("P7 ebins="+(string)ebins);
    vhistory.push_back("P8 timerange="+(string)timerange);
    if(clobber==YES) 
        vhistory.push_back("P9 clobber=yes");
    else
        vhistory.push_back("P9 clobber=no");
    if(history==YES)
        vhistory.push_back("P10 history=yes");
    else
        vhistory.push_back("P10 history=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

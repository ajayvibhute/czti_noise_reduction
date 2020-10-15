#include "cztgaas_v2.h"

Cztgaas::Cztgaas(){
    strcpy(modulename,"cztgaas_v");
    strcat(modulename,VERSION);    
}

int Cztgaas::read(int argc, char** argv) {
    int status = 0;
    
    if (PIL_OK != (status = PILInit(argc, argv))) {
        LOG(ERROR) << "***Error Initializing PIL***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("par_evtfile", eventfile))) {
        LOG(ERROR) << "***Error Reading input event data file:" << eventfile << "***";
        return status;
    }

    /*if (PIL_OK != (status = PILGetFname("par_teldef0", teldeffile0))) {
        LOG(ERROR) << "***Error Reading teldeffile file for quadrant 0:" << teldeffile0 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_teldef1", teldeffile1))) {
        LOG(ERROR) << "***Error Reading teldeffile file for quadrant 1:" << teldeffile1 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_teldef2", teldeffile2))) {
        LOG(ERROR) << "***Error Reading teldeffile file for quadrant 2:" << teldeffile2 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_teldef3", teldeffile3))) {
        LOG(ERROR) << "***Error Reading teldeffile file for quadrant 3:" << teldeffile3 << "***";
        return status;
    }*/

    if (PIL_OK != (status = PILGetFname("par_mkffile", mkffile))) {
        LOG(ERROR) << "***Error Reading attitude file name:" << mkffile << "***";
        return status;
    }
    
    if (PIL_OK != (status = PILGetFname("par_outAspectFile", outAspectFile))) {
        LOG(ERROR) << "***Error Reading output aspect array file name:" << outAspectFile << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("par_history", &history))) {
        LOG(ERROR) << "***Error Reading history parameter" << history << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("par_clobber", &clobber))) {
        LOG(ERROR) << "***Error Reading clobber:" << clobber << "***";
        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

int Cztgaas::read(char *eventfile, char *outAspectFile, char* mkffile,
        int clobber, int history) {
    strcpy(this->eventfile, eventfile);
    //strcpy(this->teldeffile0, teldef0);
    //strcpy(this->teldeffile1, teldef1);
    //strcpy(this->teldeffile2, teldef2);
    //strcpy(this->teldeffile3, teldef3);
    strcpy(this->mkffile, mkffile);
    strcpy(this->outAspectFile, outAspectFile);
    this->clobber = clobber;
    this->history = history;

    return (EXIT_SUCCESS);
}

void Cztgaas::display(){
    LOG(INFO)<<"-----------------------------------------------------------------";
    LOG(INFO)<<"                     CZTGAAS PARAMETERS                            ";
    LOG(INFO)<<"-----------------------------------------------------------------";
    LOG(INFO)<<"Taskname           : "<<modulename;
    LOG(INFO)<<"Event data file    : "<<eventfile;         //input science data file 
    //LOG(INFO)<<"TELDEF file Quad0  : "<<teldeffile0;
    //LOG(INFO)<<"TELDEF file Quad1  : "<<teldeffile1;
    //LOG(INFO)<<"TELDEF file Quad2  : "<<teldeffile2;
    //LOG(INFO)<<"TELDEF file Quad3  : "<<teldeffile3;
    LOG(INFO)<<"MKF file           : "<<mkffile;
    LOG(INFO)<<"Output aspect file : "<<outAspectFile;       //output aspect file
    LOG(INFO)<<"Clobber            : "<<clobber;
    LOG(INFO)<<"History            : "<<history;
    LOG(INFO)<<"-----------------------------------------------------------------";

}

int Cztgaas::cztgaas_process(){
    fitsfile *fptrIn, *fptrOutAs; //Input event file & Output aspect file pointers
    int status=0; //counter variable
    int i, j, k = 0; // counter variables
    int hdunum = 0; //stores current hdu number
    double t=0.0;
    int quadNo =0; //stores quadrant number currently being processed.
    int colnum = 0; //stores column number of fits extension
    int hdutype = 0; //stores hdutype
    int ncols = 0; //store number of columns in current hdu
    long nrows = 0; //store number of rows in current hdu
    string errorMsg = "";
    string extName="";
    //string teldefFilePath="";
    string teldefFileName;
    char tempChar[MAX_KEYWORD_SIZE];
    char tempChar2[MAX_KEYWORD_SIZE];
    char tempOutFilename[MAX_KEYWORD_SIZE];
    long tstarti, tstopi;
    double tstartf, tstopf;
    double tstart, tstop;
    vector<double> vec_tstart, vec_tstop, vec_time;
    vector<float> vec_nX, vec_nY, vec_nZ, vec_nXt, vec_nYt, vec_nZt;
    float nX, nY, nZ, nXt, nYt, nZt=0.0;
    float detX, detY, detZ =0.0;
    float avg_nX, avg_nY, avg_nZ, avg_nXt, avg_nYt, avg_nZt=0.0;
    float RA, DEC, TWIST =0.0;
    float RAdeg, DECdeg, TWISTdeg =0.0;
    Mkf mkf;
    Teldef teldef;


    // Checking whether output average aspect file exists or not
    // If yes then deletes it (for clobber=yes)
    // Otherwise raises an error
    for (i = 0; i < 4; i++) {
        sprintf(tempOutFilename,"");
        sprintf(tempChar2,"");
        strcat(tempOutFilename, outAspectFile);
        sprintf(tempChar2, "_Q%d", i);
        strcat(tempOutFilename, tempChar2);

        if (FileExists(tempOutFilename)) {
            if (clobber == YES) {
                if (unlink(tempOutFilename) != 0) {
                    LOG(ERROR) << "Error in deleting Output Event File: " << tempOutFilename;
                }
            } else {
                LOG(INFO) << tempOutFilename << " already exists.";
                LOG(INFO) << "Use clobber=yes for overwriting the file.";
                return (EXIT_FAILURE);
            }
        }
    }
    // Output average aspect file existence check finished.

    // Reading TSTART & TSTOP keywords for each quadrant from input event file
    //fits_open_file(&fptrIn, eventfile, READWRITE, &status); //commented by mayuri
    fits_open_file(&fptrIn, eventfile, READONLY, &status);
    if(status){
        LOG(ERROR)<<"***Error opening file: "<< eventfile;
        fits_report_error(stderr, status); 
        return (EXIT_FAILURE);
    }

  /*//storing teldefFileNames in a vector
    teldefFileNames.clear();
    teldefFileNames.push_back(teldeffile0);
    teldefFileNames.push_back(teldeffile1);
    teldefFileNames.push_back(teldeffile2);
    teldefFileNames.push_back(teldeffile3);
  *///teldef File names stored in a vector
    
    LOG(INFO) << "Reading mkf file: " << mkffile;
    //reading mkf file
    if(mkf.read_mkf_file(mkffile)){
        LOG(ERROR) << "Error reading mkf file: " << mkffile;
        status = EXIT_FAILURE;
        return status;
    }
    LOG(INFO) << "MKF file hase been read.";

    // To store tstart and tstop for all extensions of EVENT file.
    vec_tstart.clear();
    vec_tstop.clear();
    for(i=0; i<4; i++){
        LOG(INFO) << "Processing for Quadrant " << i;
        quadNo=i;
        extName = quadExtName(quadNo);
        fits_movnam_hdu(fptrIn, BINARY_TBL, (char*) extName.c_str(), 0, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to extension name " << extName << " in event file " << eventfile << ".";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        /*
        fits_read_key(fptrIn, TLONG, "TSTARTI", &tstarti, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading TSTARTI keyword in " << extName << " of event file " << eventfile << ".";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }*/
        fits_read_key(fptrIn, TDOUBLE, "TSTART", &tstart, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading TSTART keyword in " << extName << " of event file " << eventfile << ".";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
/*
        fits_read_key(fptrIn, TLONG, "TSTOPI", &tstopi, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading TSTOPI keyword in " << extName << " of event file " << eventfile << ".";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }*/
        fits_read_key(fptrIn, TDOUBLE, "TSTOP", &tstop, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading TSTOP keyword in " << extName << " of event file " << eventfile << ".";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        //tstart = (double)tstarti + tstartf;
        //tstop = (double) tstopi+ tstopf;
        vec_tstart.push_back(tstart);
        vec_tstop.push_back(tstop);
        
	int teldefExtnum;
	char quadnam[20];
	sprintf(quadnam,"QUADRANT%d",i);

        //teldefFilePath = caldb_full_path_generator(teldefFileNames[i]);
        if(QueryCaldb("ASTROSAT","CZTI",(string)quadnam,"TELDEF",tstart,tstop,teldefFileName,teldefExtnum))
    	{
        	LOG(ERROR) << "Not able to get CALDB EBOUNDS file";
        	return (EXIT_FAILURE);
    	}

	if(teldef.read_teldef_file(teldefFileName)){
            LOG(ERROR) << "Error reading Teldef file: " << teldefFileName;
            status = EXIT_FAILURE;
            return EXIT_FAILURE;
        }

        //clearing vectors
        vec_time.clear();
        vec_nX.clear();
        vec_nXt.clear();
        vec_nY.clear();
        vec_nYt.clear();
        vec_nZ.clear();
        vec_nZt.clear();
        
        int num = (int)(tstop-tstart);
        LOG(INFO) << "TSTOP-TSTART=" << num;
        t=tstart;
        mkf.set_last_index(0);
        
        for(j=0; j<num; j++){
            //DLOG(INFO) << setprecision(12) << t;
            //Central pointing direction
            detX=0; detY=0; detZ=1;
            if(get_inertial_vector(mkf, teldef, t, detX, detY, detZ, nX, nY, nZ)){
                LOG(ERROR) << "Error in getting inertial attitude vector.";
                status = EXIT_FAILURE;
                return EXIT_FAILURE;
            }
            //DLOG(INFO)<< "nX, nY, nZ" << nX << "   " << nY << "   " << nZ;

            //For computation of twist angle
            detX=0; detY=1; detZ=0;
            if(get_inertial_vector(mkf, teldef, t, detX, detY, detZ, nXt, nYt, nZt)){
                LOG(ERROR) << "Error in getting inertial attitude vector.";
                status = EXIT_FAILURE;
                return EXIT_FAILURE;
            }
            //DLOG(INFO) << setprecision(12) << "TIME:" << t << " nX:"<< nX<< " nY:"<< nY<< " nZ:"<< nZ << " nXt:"<< nXt<< " nYt:"<< nYt<< " nZt:"<< nZt;            
            vec_time.push_back(t);
            vec_nX.push_back(nX);
            vec_nXt.push_back(nXt);
            vec_nY.push_back(nY);
            vec_nYt.push_back(nYt);
            vec_nZ.push_back(nZ);
            vec_nZt.push_back(nZt);
            
            t=t+1.0;
        }
        sprintf(tempOutFilename, "");
        sprintf(tempChar2, "");
        strcat(tempOutFilename, outAspectFile);
        sprintf(tempChar2, "_Q%d", i);
        strcat(tempOutFilename, tempChar2);
        LOG(INFO) << "Creating Average Aspect File: " << tempOutFilename;
        //Create aspect file based on template
        if(create_aspect_file((string)tempOutFilename)){
            LOG(ERROR) << "Error in creating aspect file: " << tempOutFilename;
        }
        
        // opening and writing data to aspect file
        fits_open_file(&fptrOutAs, tempOutFilename, READWRITE, &status);
        if (status) {
            LOG(ERROR) << "Error in opening output aspect file " << tempOutFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_movnam_hdu(fptrOutAs, BINARY_TBL, "ASPECT", NULL, &status);
        if (status) {
            LOG(ERROR) <<"Error in moving to ASPECT HDU in " << tempOutFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        // Writing TIME column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "TIME", &colnum, &status);
        errorMsg = "Error in getting column number of TIME column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TDOUBLE, colnum, 1, 1, vec_time.size(), vec_time.data(), &status);
        errorMsg = "Error in writing TIME column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        // Writing Nx column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "Nx", &colnum, &status);
        errorMsg = "Error in getting column number of Nx column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TFLOAT, colnum, 1, 1, vec_nX.size(), vec_nX.data(), &status);
        errorMsg = "Error in writing Nx column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        // Writing Ny column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "Ny", &colnum, &status);
        errorMsg = "Error in getting column number of Ny column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TFLOAT, colnum, 1, 1, vec_nY.size(), vec_nY.data(), &status);
        errorMsg = "Error in writing Ny column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        // Writing Nz column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "Nz", &colnum, &status);
        errorMsg = "Error in getting column number of Nz column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TFLOAT, colnum, 1, 1, vec_nZ.size(), vec_nZ.data(), &status);
        errorMsg = "Error in writing Nz column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        // Writing Nxt column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "Nxt", &colnum, &status);
        errorMsg = "Error in getting column number of Nxt column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TFLOAT, colnum, 1, 1, vec_nXt.size(), vec_nXt.data(), &status);
        errorMsg = "Error in writing Nxt column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        // Writing Nyt column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "Nyt", &colnum, &status);
        errorMsg = "Error in getting column number of Nyt column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TFLOAT, colnum, 1, 1, vec_nYt.size(), vec_nYt.data(), &status);
        errorMsg = "Error in writing Nyt column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        // Writing Nzt column    
        fits_get_colnum(fptrOutAs, CASEINSEN, "Nzt", &colnum, &status);
        errorMsg = "Error in getting column number of Nzt column in Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
        fits_write_col(fptrOutAs, TFLOAT, colnum, 1, 1, vec_nZt.size(), vec_nZt.data(), &status);
        errorMsg = "Error in writing Nzt column in ASPECT extension of Aspect file: " + (string)tempOutFilename;
        if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

        calculate_RA_DEC_Twist(vec_nX, vec_nY, vec_nZ, vec_nXt, vec_nYt, vec_nZt, RA, DEC, TWIST, status);
        LOG(INFO) << "RA: " << RA << " DEC: " << DEC << " TWIST: " << TWIST;

        RAdeg = RA * 180/M_PI; //to degrees
        DECdeg = DEC  * 180/M_PI;
        TWISTdeg = TWIST  * 180/M_PI;
        

        fits_write_key(fptrOutAs, TFLOAT, "RA", &RAdeg, "", &status);
        if (status) {
            LOG(ERROR) << "Error in writing key RA(deg) in average aspect file " << tempOutFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_key(fptrOutAs, TFLOAT, "DEC", &DECdeg, "", &status);
        if (status) {
            LOG(ERROR) << "Error in writing key DEC(deg) in average aspect file " << tempOutFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_key(fptrOutAs, TFLOAT, "TWIST", &TWISTdeg, "", &status);
        if (status) {
            LOG(ERROR) << "Error in writing key TWIST(deg) in average aspect file " << tempOutFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        
        fits_close_file(fptrOutAs, &status);
 

        //writing history to output file
        if (history == YES) {
            vector<string> vhistory;
            get_history(vhistory);
            writeHistory(tempOutFilename, vhistory);
        }

        //updating keywords in output file
        updateKeywords(tempOutFilename, modulename);
        
    }
    
    fits_close_file(fptrIn, &status);

    return status;
}

int create_aspect_file(string aspectFileName){
    fitsfile* fptrAspect;
    int status=0; // status variable
    int i,j=0; //counter variables

    string templateFileName = "";

    templateFileName = template_full_path_generator("aspectTemplate");
    
    if(templateFileName==""){
        LOG(ERROR)<< "Not able to generate Event file template path.";
        LOG(ERROR)<< "Probably Environment Variables are not declared properly.";
        return(EXIT_FAILURE);
    }    
    LOG(INFO) << "Template used: " << templateFileName;
    fits_create_template(&fptrAspect, (char*) aspectFileName.c_str(), (char*) templateFileName.c_str(), &status);
    if (status) {
        LOG(ERROR) << "***Error creating Aspect file: " << aspectFileName << "***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    

    fits_close_file(fptrAspect, &status);
    if (status) {
        LOG(ERROR) << "***Error in closing Aspect File : " << aspectFileName << " ***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    LOG(INFO) << "Output aspect file " << aspectFileName << " created successfully.";
    return (EXIT_SUCCESS);


    
    return status;
}

int Cztgaas::get_history(vector<string> &vhistory){
    //char *user=getlogin();
    strcpy(modulename,"cztgaas_v");
    strcat(modulename,VERSION); 

    char *user = getenv("USER");
	vhistory.push_back("Module run by "+(string)user);
    vhistory.push_back("Parameter List START for "+(string)modulename);
    vhistory.push_back("P1 eventfile="+(string)eventfile);
    /*vhistory.push_back("P2 teldeffile_Q0="+(string)teldeffile0);
    vhistory.push_back("P3 teldeffile_Q1="+(string)teldeffile1);
    vhistory.push_back("P4 teldeffile_Q2="+(string)teldeffile2);
    vhistory.push_back("P5 teldeffile_Q3="+(string)teldeffile3);
    */
    vhistory.push_back("P2 mkfFile="+(string)mkffile);
    vhistory.push_back("P3 outAspectFile="+(string)outAspectFile);
       
    if(clobber==YES) 
        vhistory.push_back("P4 clobber=yes");
    else
        vhistory.push_back("P4 clobber=no");
    
    if(history==YES)
        vhistory.push_back("P5 history=yes");
    else
        vhistory.push_back("P5 history=no");
    
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

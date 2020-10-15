#include "cztflagbadpix.h"

cztflagbadpix::cztflagbadpix() {
    strcpy(modulename, "cztflagbadpix_v");
    strcat(modulename, VERSION);
}

int cztflagbadpix::read(int argc, char** argv) {
    int status=0;
    char badpixfile[PIL_LINESIZE];
    int ifile=0;
    if (PIL_OK != (status = PILInit(argc, argv))) {
        LOG(ERROR) << "***Error Initializing PIL***";
        return status;
    }

    do {
        if (PIL_OK != (status = PILGetInt("nbadpixFiles", &nbadpixFiles))) {
            LOG(ERROR) << "***Error Reading number of badpixfiles***";
            return status;
        }
    }while (nbadpixFiles < 1);    
    
    for(ifile=0; ifile<nbadpixFiles; ifile++) {
        if (PIL_OK != (status = PILGetFname("badpixfile", badpixfile))) {
            LOG(ERROR) << "***Error Reading input badpixfile " << ifile << " : " << badpixfile << "***";
            return status;
        }
        badpixFiles.push_back((string)badpixfile);
    }

    if (PIL_OK != (status = PILGetFname("outfile", outfile))) {
        LOG(ERROR) << "***Error reading output file***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("history", &history))) {
        LOG(ERROR) << "***Error reading history parameter***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("clobber", &clobber))) {
        LOG(ERROR) << "***Error Reading clobber***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("debug", &debug))) {
        LOG(ERROR) << "***Error Reading debug***";
        return status;
    }
    PILClose(status);
    
    
    return (EXIT_SUCCESS);
}

int cztflagbadpix::read(int nbadpixFiles, vector<string> badpixFiles, char* outfile, int debug, int history, int clobber) {
    this->nbadpixFiles = nbadpixFiles;
    this->badpixFiles = badpixFiles;
    strcpy(this->outfile, outfile);
    this->clobber = clobber;
    this->history = history;
    this->debug = debug;
    return (EXIT_SUCCESS);
}

void cztflagbadpix::display() {
    int ifile=0; 
    LOG(INFO) << "----------------------------------------------------------------------------";
    LOG(INFO) << "                          CZTFLAGBADPIX PARAMETERS                          ";
    LOG(INFO) << "----------------------------------------------------------------------------";
    LOG(INFO)<<"Modulename                         : " << modulename;
    LOG(INFO)<<"Number of badpixel files           : " << nbadpixFiles;
    for(ifile=0; ifile<badpixFiles.size(); ifile++){
    LOG(INFO)<<"Bad pixel file " << ifile << "     : " << badpixFiles[ifile];
    }
    LOG(INFO)<<"Output badpixel file               : " << outfile;
    LOG(INFO)<<"Clobber                            : " << clobber;
    LOG(INFO)<<"History                            : " << history;
    LOG(INFO)<<"Debug                              : " << debug;
    LOG(INFO)<<"------------------------------------------------------------------------------";
}

int cztflagbadpix::cztflagbadpix_process() {
    int status=0;
    Badpix badpix; //to read, write and store badpixel file information
    vector <BadpixTable> bptables; //to store bad pixel tables for all 4 quadrants after flagging
                                   //from multiple files.
    if (FileExists(outfile)) {
        if (clobber == YES) {
            unlink(outfile);
        } else {
            LOG(ERROR) << "" << outfile << " already exists";
            LOG(ERROR) << "Use clobber=yes for overwriting the file";
            return (EXIT_FAILURE);
        }
    }
    
    try {
        bptables = flag_badpix_files(badpixFiles);
    }    catch (string s) {
        LOG(ERROR) << s;
        return EXIT_FAILURE;
    }
    try {
        badpix.set_badpix_tables(bptables);
    } catch (string s) {
        LOG(ERROR) << s;
        return EXIT_FAILURE;
    }
    
    if(badpix.create_L2badpix_file((string) outfile, BADPIXTEMPLATE)){
        LOG(ERROR) << "Error in creating empty level-2 badpixel file.";
        return EXIT_FAILURE;
    }
    
    if(badpix.write_L2badpix_file((string) outfile)){
        LOG(ERROR) << "Error in writing level-2 badpixel file.";
        return EXIT_FAILURE;
    }

    if (history == YES) {
        vector<string> vhistory;
        get_history(vhistory);
        writeHistory(outfile, vhistory);
    }
    //UPDATING KEYWORDS
    updateKeywords(outfile, modulename);
    
    return status;
}

int cztflagbadpix::get_history(vector<string>& vhistory) {
    //char *user = getlogin();
    strcpy(modulename, "cztflagbadpix_v");
    strcat(modulename, VERSION);

    char *user = getenv("USER");
	int ifile=0;
    string badpixfilenames="";
    for (ifile = 0; ifile < badpixFiles.size(); ifile++) {
        badpixfilenames += badpixFiles[ifile] + " ";
    }
    string str = "Module run by " + (string) user;
    vhistory.push_back(str);
    vhistory.push_back("Parameter List START for " + (string) modulename);
    vhistory.push_back("P1 nbadpixFiles=" + itoa(nbadpixFiles));
    vhistory.push_back("P2 badpixFiles=" + badpixfilenames);
    vhistory.push_back("P3 outfile=" + (string) outfile);
    if (clobber == YES)
        vhistory.push_back("P4 clobber=yes");
    else
        vhistory.push_back("P4 clobber=no");
    if (history == YES)
        vhistory.push_back("P5 history=yes");
    else
        vhistory.push_back("P5 history=no");
    if (debug == YES)
        vhistory.push_back("P6 debug=yes");
    else
        vhistory.push_back("P6 debug=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}




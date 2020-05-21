#include "cztiopath.h"

FileDirectoryHandler::FileDirectoryHandler() {
    l1pathFlag=false;
    l2pathFlag=false;
    nModes=0;
}

string FileDirectoryHandler::get_l1_file_prefix(string sciencefilepath, string filetype) {
    char *infilepath = strdup((char*) sciencefilepath.c_str());
    char *temp = basename(infilepath); //retrieve onlyfilename
    char *ext = strchr(temp, '.'); //extract name before extension
    string strtemp="";
    if (ext != NULL) strtok(temp, ".");
   
    //changing level1 in fileprefix to level2
    strtemp = (string) temp;
    return strtemp;
}

string FileDirectoryHandler::get_l2_file_prefix(string l1fileprefix) {
    int pos = l1fileprefix.find("level1");
    string l2fileprefix="";
    if (pos > 0 && pos < l1fileprefix.length()) {
        l1fileprefix.replace(pos, string("level1").length(), "level2");
    } else {
        LOG(ERROR) << "Cannot generate level-2 file prefix";
        return EXIT_FAILURE;
    }
    
    return l2fileprefix;
}
void FileDirectoryHandler::generate_l2_filepaths(){
    int imode=0;
    string modename="";
    ErrorHandler errHandler;
    if(l1pathFlag==true){
        l2Paths["INDIR"] = l2basepath + "/" + CZTDIRNAME;
        l2Paths["AUXDIR_IN"] = l2basepath + "/" + "aux";
        l2Paths["AUX1DIR_IN"] = l2Paths["AUXDIR_IN"] + "/" + "aux1";
        l2Paths["AUX2DIR_IN"] = l2Paths["AUXDIR_IN"] + "/" + "aux2";
        l2Paths["AUX3DIR_IN"] = l2Paths["AUXDIR_IN"] + "/" + "aux3";
    
        for(imode=0; imode<TOTALMODES; imode++){
            modename = modePatterns[imode];
            if(modeAvailable[modename]==true){
                l2Paths[modename]=l2Paths["INDIR"]+"/mode"+modename;
                model2FilePaths[modename].scienceFilename=l2Paths[modename];
            }
        }
    }
    
    
}
void FileDirectoryHandler::generate_l1_filepaths(){
    vector <string> tempdirlist;
    vector <string> tempfilelist;
    vector <string> tempSearchedFiles;
    vector <string> filetypes;
    vector <string> extraChecks; //checks to validate file input.
    ErrorHandler errHandler;
    string modename="";
    int ifile=0, idir=0, imode=0, itype=0;
    
    l1Paths["INDIR"]=l1basepath+"/"+CZTDIRNAME;
    l1Paths["AUXDIR_IN"]=l1basepath+"/"+"aux";
    l1Paths["AUX1DIR_IN"]=l1Paths["AUXDIR_IN"]+"/"+"aux1";
    l1Paths["AUX2DIR_IN"]=l1Paths["AUXDIR_IN"]+"/"+"aux2";
    l1Paths["AUX3DIR_IN"]=l1Paths["AUXDIR_IN"]+"/"+"aux3";

    //getting observation ID
    if(find_observation_id() !=NULL) {
        observationID = (string) find_observation_id();
    } else {
        errHandler.severity = errERROR;
        errHandler.errorStatus = INVALID_L1_BASEPATH;
        errHandler.errorMsg = "Error in getting observation ID from level-1 basepath " + l1basepath+ ".";
        throw errHandler;
    }
    
    //setting extraChecks
    extraChecks.push_back(observationID);
    cout << "observation id: " << observationID;
    //Getting file paths of all files in l1 base directory
    try{
        getFiles(l1Paths["INDIR"], &tempfilelist);
    } catch(ErrorHandler errHandler){
        logError(errHandler);
        if(errHandler.errorStatus!=EMPTY_DIRECTORY){
            throw errHandler;
        }
    }
    
    //storing mkf, attitude, orbit, tct & xml[if available] files
    filetypes.clear();
    string filetypesArr[]={MKF_EXT, ORB_EXT, ATT_EXT, TCT_EXT, LBT_EXT};
    filetypes.assign(filetypesArr, filetypesArr+4);
    //MKF
    for (ifile = 0; ifile < filetypes.size(); ifile++) {
        tempSearchedFiles = find_files_in_filelist(tempfilelist, filetypes[ifile], &extraChecks);
        if(tempSearchedFiles.size()==0){
            errHandler.severity = errERROR;
            errHandler.errorStatus = FILE_NOT_FOUND;
            errHandler.errorMsg = filetypes[ifile] + " file not found.";
            throw errHandler;
        } else if(tempSearchedFiles.size()!=1){
            errHandler.severity = errERROR;
            errHandler.errorStatus = MULTIPLE_FILES_FOUND;
            errHandler.errorMsg = "Number of " + filetypes[ifile] + " files found= " + itoa(tempSearchedFiles.size());
            throw errHandler;
        }
        l1Paths[filetypes[ifile]]=tempSearchedFiles[0];
    }
    
    
    //Finding number of mode directories in l1 base directory
    tempdirlist.clear();
    tempfilelist.clear();
    try {
        searchDir(l1Paths["INDIR"], "mode", &tempfilelist, &tempdirlist);
        check_mode_availability(tempdirlist);
        //checking existence of SS mode, otherwise raising error
        if(modeAvailable[modePatterns[SS]]==false){
            errHandler.severity = errERROR;
            errHandler.errorStatus = DIRECTORY_NOT_FOUND;
            errHandler.errorMsg = "modeSS directory not found.";
            throw errHandler;
        }
        //storing corresponding mode file paths.
        string modefiletypeArr[]={EVT_EXT, GTI_EXT, BTI_EXT};
        filetypes.clear();
        filetypes.assign(modefiletypeArr, modefiletypeArr+3);
        for(imode=0; imode<TOTALMODES; imode++){
            modename = modePatterns[imode];
            if(modeAvailable[modename]==false){
                model1FilePaths[modename].scienceFilename="";
                model1FilePaths[modename].btiFilename="";
                model1FilePaths[modename].gtiFilename="";
            } else {
                //getting modeFilePaths for science, bti & gti file.
                //storing mode directory paths
                l1Paths[modename] = l1Paths["INDIR"] + "/mode" + modename;
                extraChecks.clear();
                extraChecks.push_back(observationID);
                extraChecks.push_back(modename);
                tempfilelist.clear();
                tempdirlist.clear();
                searchDir(l1Paths[modename], modename, &tempfilelist, &tempdirlist);
                
                for (itype = 0; itype < filetypes.size(); itype++) {
                    tempSearchedFiles = find_files_in_filelist(tempfilelist, filetypes[itype], &extraChecks);
                    if (tempSearchedFiles.size() != 1) {
                        modeAvailable[modename] = false;
                        model1FilePaths[modename].l1prefix = "";
                        model1FilePaths[modename].scienceFilename = "";
                        model1FilePaths[modename].btiFilename = "";
                        model1FilePaths[modename].gtiFilename = "";
                        break;
                    } else {
                        if (filetypes[itype] == EVT_EXT) {
                            model1FilePaths[modename].scienceFilename = tempSearchedFiles[0];
                        } else if (filetypes[itype] == GTI_EXT) {
                            model1FilePaths[modename].gtiFilename = tempSearchedFiles[0];
                        } else if (filetypes[itype] == BTI_EXT) {
                            model1FilePaths[modename].btiFilename = tempSearchedFiles[0];
                        }
                    }
                }
                if(model1FilePaths[modename].scienceFilename != ""){
                    model1FilePaths[modename].l1prefix = get_l1_file_prefix(model1FilePaths[modename].scienceFilename, "l1");
                }
            }
        }
    } catch(ErrorHandler errHandler){
        throw errHandler;
    }

    l1pathFlag=true; //level-1 file paths generated successfully
}

char* FileDirectoryHandler::find_observation_id(){
    char *level1dir = basename((char*) l1basepath.c_str());
    if (level1dir == NULL) return NULL;
    int len = strlen(level1dir);
    char *temp = new char[len + 1];
    strcpy(temp, level1dir);
    char *t1 = strchr(temp, '_');
    if (t1 == NULL) return NULL;
    char *t2 = strrchr(t1, '_');
    if (t2 == NULL) return NULL;
    int obsidlen = strlen(t1) - 1 - strlen(t2);
    if (obsidlen <= 0) return NULL;
    t1++;
    char *obsid = new char[obsidlen + 1];
    if (obsid == NULL) return NULL;
    strncpy(obsid, t1, obsidlen);
    obsid[obsidlen] = '\0';
    
    this->observationID = (string) obsid;
    delete[] temp;
    return obsid;
}

void FileDirectoryHandler::check_mode_availability(vector<string> modeDirNames) {
    int idir = 0;
    int imodePattern = 0;
    size_t pos;
    modePatterns.resize(16, "");
    modePatterns[M0] = "M0";
    modePatterns[M1] = "M1";
    modePatterns[M2] = "M2";
    modePatterns[M3] = "M3";
    modePatterns[M4] = "M4";
    modePatterns[M5] = "M5";
    modePatterns[M6] = "M6";
    modePatterns[M7] = "M7";
    modePatterns[M8] = "M8";
    modePatterns[M9] = "M9";
    modePatterns[MA] = "MA";
    modePatterns[MB] = "MB";
    modePatterns[MC] = "MC";
    modePatterns[MD] = "MD";
    modePatterns[ME] = "ME";
    modePatterns[SS] = "SS";

    for (idir = 0; idir < modeDirNames.size(); idir++) {
        for (imodePattern = 0; imodePattern < modePatterns.size(); imodePattern++) {
            pos = modeDirNames[idir].find(modePatterns[imodePattern]);

            if (pos != std::string::npos) {
                modeAvailable[modePatterns[imodePattern]] = true;
                nModes++;
            }
        }
    }

}
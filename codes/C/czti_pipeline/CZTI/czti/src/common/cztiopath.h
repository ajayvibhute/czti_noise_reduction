/* 
 * @file  cztiopath.h
 * @author Tanul Gupta
 * @date Created on May 8, 2015, 1:03 PM
 * @brief Create input and output paths of czti level-1 and level-2 directory
 * @details 
 * @version 0.2.0
 */

#ifndef CZTIOPATH_H
#define CZTIOPATH_H

#include "utils.h"
#include "errorHandler.h"
#include <vector>
#include <string>

#define CZTDIRNAME "czti"

typedef string extnames;
const extnames MKF_EXT="mkf";
const extnames ORB_EXT="orb";
const extnames ATT_EXT="att";
const extnames TCT_EXT="tct";
const extnames LBT_EXT="lbt";
const extnames EVT_EXT="fits";
const extnames GTI_EXT="gti";
const extnames BTI_EXT="bti";

const int TOTALMODES=16;
typedef int ModeCodes;
const ModeCodes M0=0;
const ModeCodes M1=1;
const ModeCodes M2=2;
const ModeCodes M3=3;
const ModeCodes M4=4;
const ModeCodes M5=5;
const ModeCodes M6=6;
const ModeCodes M7=7;
const ModeCodes M8=8;
const ModeCodes M9=9;
const ModeCodes MA=10;
const ModeCodes MB=11;
const ModeCodes MC=12;
const ModeCodes MD=13;
const ModeCodes ME=14;
const ModeCodes SS=15;

struct sciencel1Files{
    string l1prefix;
    string scienceFilename;
    string gtiFilename;
    string btiFilename;
};
struct sciencel2Files{
    string l2prefix;
    string scienceFilename;
    string gtiFilename;
    string btiFilename;
};
class FileDirectoryHandler{
public:
    string l1basepath;
    string l2basepath;
    bool l1pathFlag; //set to true if all level-1 file paths generated
    bool l2pathFlag; //set to true if all level-2 file paths generated
    string observationID;
    int nModes; //number of modes for which data is available (set by generate_l1paths)
    vector <string> modePatterns;
    map <string, bool> modeAvailable;
    map <string, sciencel1Files> model1FilePaths;
    map <string, sciencel2Files> model2FilePaths;
    map<string, string> l1Paths; /**1.INDIR
                                  * 2.AUXDIR_IN
                                  * 3.AUX1DIR_IN
                                  * 4.AUX2DIR_IN
                                  * 5.AUX3DIR_IN
                                  * 6.mkf
                                  * 7.orb
                                  * 8.att
                                  * 9.tct
                                  * 10.lbt
                                  */
    map<string, string> l2Paths;

    FileDirectoryHandler();
    void generate_l1_filepaths();
    void generate_l2_filepaths();
    void check_mode_availability(vector<string> modeDirNames);
    char* find_observation_id();
    string get_l1_file_prefix(string sciencefilepath);
    string get_l2_file_prefix(string l1fileprefix);
    
    //setters
    void set_l1basepath(string l1basepath){this->l1basepath=l1basepath;}
    void set_l2basepath(string l2basepath){this->l2basepath=l2basepath;}
    
    //getters
    string get_obsid(){return observationID;}
};
#endif /*CZTIOPATH_H*/
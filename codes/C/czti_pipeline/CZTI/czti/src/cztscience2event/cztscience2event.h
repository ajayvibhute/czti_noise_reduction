/**
 * @file cztscience2event.h
 * @author Tanul Gupta, Preeti Tahlani
 * @date created on July 6, 2012, 5:22 PM; modified on April 7, 2015, 12:39 PM
 *       Major update for optimization on Aug 14, 2015, 9:17 AM.
 * @brief cztscience2event header file
 * @details This class is designed to extract level1 event data from science files.
 *        ____        _          _____                          _
      |  _ \  __ _| |_ __ _  |  ___|__  _ __ _ __ ___   __ _| |_
      | | | |/ _` | __/ _` | | |_ / _ \| '__| '_ ` _ \ / _` | __|
      | |_| | (_| | || (_| | |  _| (_) | |  | | | | | | (_| | |_
      |____/ \__,_|\__\__,_| |_|  \___/|_|  |_| |_| |_|\__,_|\__|

      +-------+--------------+       +-------+--------------+
      | DataID| Description  |       | ModeID| Description  |
      ........|..............|       ........|..............|
      |  0    | Q0 det_data  |       |  0    |Normal        |
      |  1    | Q1 det_data  |       |  1    |VSD           |
      |  2    | Q2 det_data  |       |  2    |2 wevrep      |
      |  3    | Q3 det_data  |       |  3    |VSD + 2wevrep |
      |  4    | Q0 hdr_data  |       |  4    |Fix pkt.      |
      |  5    | Q1 det_data  |       |  5    |Fxpkt+VSD     |
      |  6    | Q2 det_data  |       |  6    |Fxpkt+2wevrep |
      |  7    | Q3 det_data  |       |  7    |Fx+VSD+2wevrp |
      |  8    | Q0 ssm_data  |       |  8    |------        |
      |  9    | Q1 ssm_data  |       |  9    |SAA           |
      |  10   | Q2 ssm_data  |       |  10   |Shadow        |
      |  11   | Q3 ssm_data  |       |  11   |SAA+Shadow    |
      |  12   | 0eepromdata  |       |  12   |mml>1         |
      |  13   | 1eepromdata  |       |  13   |mml>1+SAA     |
      |  14   | 2eepromdata  |       |  14   |mml>1+Shadow  |
      |  15   | 3eepromdata  |       |  15   |mml>1+shdw+SAA|
      |  16   | PE header    |       |       |              |
      +-------+--------------+       +-------+--------------+

 * Event data exists only for Data IDs: 0-3, 8-11 and Mode IDs:0-7
 * This program checks whether data belongs to this range, and then extract data
 * that satisfies this criteria in event file.
 * @version 0.2.0
 */

#ifndef CZTSCIENCE2EVENT_H
#define	CZTSCIENCE2EVENT_H

//Include section
#include <vector>
#include <iomanip>
#include "glog/logging.h"
#include <string>
#include <fstream>
#include <fitsio.h>
#include "cztstring.h"
#include "utils.h"
#include "cztHeaderParam.h"
#include "level2validation.h"
#include "level1handler.h"
#include "l1evtdecode.h"
#include "caldbHandler.h"
#include "macrodef.h"


using namespace std;


/*! 
 * \brief cztscience2event class:
    This class is designed to extract level1 event data from science files.
 * 
    It contains the function required for effective extraction of data based on
    CZTI mode id.
 */
class cztscience2event{

private:
    char modulename[NAMESIZE];
    char infile[PIL_LINESIZE];        /**< Input Science Data file */
    char bunchfile[PIL_LINESIZE];
    char TCTfile[PIL_LINESIZE];
    char outfile[PIL_LINESIZE]; 
    char hdrInfoFile[PIL_LINESIZE];
    char hdrInfoTemplate[PIL_LINESIZE]; 
    char GTIfile[PIL_LINESIZE];        //output gti file
                                       //GTI extension from science data file is extracted and gti file is created                                                                 
    char BTIfile[PIL_LINESIZE];        //output bti file
                                       //BTI extension from science data file is extracted and bti file is created                         
    int applyLUT;              //whether LUT has to be applied to convert 12 bit PHA to 10 bit PHA (for cztscience2event module)
    char LUTfile[PIL_LINESIZE];  //LUTfile to used for conversion of 12 bit PHA to 10 bit PHA (for cztscience2event module)
    int nPackets;              //number of level-1 packets to be read in buffer
    int BigEndian;                    //True if data packet is in bigendian
    int clobber;                      //Overwrite Existing file
    int history;
    int debug; //if yes then debug info is displayed otherwise not.

 public:
    cztscience2event();
    int read(int argc, char **argv);
    int read(char *infile,char *TCTfile,char *outfile,char *hdrInfoFile,char *gtifile,
        char *btifile,int nPackets, int bigendian=YES,
        int clobber=YES,int history=NO, int debug=NO);
    void display();
    
    /**
 * Function for converting science data file from level-1 to event file in level-2
 * @return 
 */
    int cztscience2eventProcess();
    

    
    /**
    * Function to create output GTI and BTI files. It copies the GTI/BTI extensions
    * from input science data file to output GTI/BTI files.  
    * @return GTI and BTI files
    */
    int createGTIBTIfile();
    
   /**
    * Creates history for the module to be written to the output files
    * @param par : Object of cztscience2event_param class
    * @param vhistory : vector ho hold history strings
    * @return 
    */
    int getHistory(vector<string> &vhistory);
    
};

#endif	/* CZTSCIENCE2EVENT_H */


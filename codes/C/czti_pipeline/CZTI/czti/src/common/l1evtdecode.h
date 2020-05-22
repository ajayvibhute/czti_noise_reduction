/* 
 * @file  l1evtdecode.h
 * @author Tanul Gupta
 * @date Created on August 6, 2015, 3:50 PM
 * @brief decodes, extracts and writes l1 data.
 * @details This class decodes, extracts, and writes l1data.
 * @version 0.2.0
 */

#ifndef L1EVTDECODE_H
#define L1EVTDECODE_H

#include <iostream>
#include "glog/logging.h"
#include "utils.h"
#include <string>
#include <vector>
#include <fitsio.h>


using namespace std;

#define PKTSIZE 2048
#define MAXEVENTS_IN_ONE_PACKET 340 
#define VETOSPEC_SIZE 256
#define CZTSPEC_SIZE 512
#define CZTSPECVETO_SIZE 512
#define CZTSPECALPHA_SIZE 512
#define NUM_DET_PER_QUADRANT 16
#define FEBSAA 100           //number of feb headers for 100 second data
#define MAX_PACKETS_FRAME 10
#define MAXFRAMEDATASIZE 1008+1020*(MAX_PACKETS_FRAME-1)

struct FEBheader {
    //----------Word 0 Status word-----------
    //D15-D12 = 0
    bool VetoSpecRange; //D11
    char BaseAdd; //D10-D8
    bool BufferNo; //D7
    bool CmdStatus; //D6
    bool EventReadMode; //D5
    bool CZTstatus; //D4
    char CZTNo; //D3-D0
    //----------Word 1 HK word--------------
    //D15 = 0
    char ChannelNo; //D14-D12
    short ADCoutput; //D11-D0
    //----------Word 2 Command-----------
    unsigned short cmd1sec; //command in 1 second
    //----------Word 3 Alpha count--------
    unsigned short AlphaCount;
    //----------Word 4 Veto count---------
    unsigned short VetoCount;
    //----------Word 5 --------------
    unsigned short CZTcount_lt_ULD; //CZT count < ULD
    //----------Word 6--------------
    unsigned short CZTcount_gt_ULD; //CZT count >=ULD
    //----------Word 7--------------
    unsigned short CZTdataread; //data read from CZT
    //----------Word 8-23 CZT temperatures---------
    unsigned short temperature[NUM_DET_PER_QUADRANT]; //temperature of 16 CZT detectors of 1 quadrant

    //for SAA mode FEB header data
    //----------Word 0 Status word-----------
    //D15-D12 = 0
    bool VetoSpecRange1[FEBSAA]; //D11
    char BaseAdd1[FEBSAA]; //D10-D8
    bool BufferNo1[FEBSAA]; //D7
    bool CmdStatus1[FEBSAA]; //D6
    bool EventReadMode1[FEBSAA]; //D5
    bool CZTstatus1[FEBSAA]; //D4
    char CZTNo1[FEBSAA]; //D3-D0
    //----------Word 1 HK word--------------
    //D15 = 0
    char ChannelNo1[FEBSAA]; //D14-D12
    short ADCoutput1[FEBSAA]; //D11-D0
    //----------Word 2 Command-----------
    unsigned short cmd1sec1[FEBSAA]; //command in 1 second
    //----------Word 3 Alpha count--------
    unsigned short AlphaCount1[FEBSAA];
    //----------Word 4 Veto count---------
    unsigned short VetoCount1[FEBSAA];
    //----------Word 5 --------------
    unsigned short CZTcount_lt_ULD1[FEBSAA]; //CZT count < ULD
    //----------Word 6--------------
    unsigned short CZTcount_gt_ULD1[FEBSAA]; //CZT count >=ULD
    //----------Word 7--------------
    unsigned short CZTdataread1[FEBSAA]; //data read from CZT

    void read(unsigned short *data, int size, int modeid, int dataid);
};

class Packet{
public:
    bool isFirstPacket;
    bool isValid;
    //Header fields
    unsigned short SYNC1;
    unsigned short SYNC2;
    unsigned char mml;         //memory level
    unsigned char pktNo;
    unsigned char dataID;
    unsigned char modeID;
    unsigned short wcnt;
    unsigned char quadID;
    unsigned short *hbtdata; //to store event data
    FEBheader febhdr;
    
    //Frame Keywords for 0th packet
    unsigned short frameno; //frame number   
    unsigned short statusFlag;
    unsigned short writePktNo;
    unsigned short readPktNo;
    int command;
    int secondcount;
    long laxpctime;
    unsigned short dcnt;
    unsigned char errorcount;
    unsigned char errorflag;
    unsigned char bootpageno;    
    
    bool validity();   //Function to check packet validity
 
    Packet();
    ~Packet(); 
        
};

class Frame{
public:
    int isvalid;                    //-1:not frame; 0:invalid frame; 1:valid frame
    unsigned short frameno;          //frame number   
    unsigned short statusFlag;
    unsigned short writePktNo;
    unsigned short readPktNo;
    int command;
    int secondcount;
    double laxpctime;
    unsigned short dcnt;
    unsigned char errorcount;
    unsigned char errorflag;
    unsigned char bootpageno;
    vector <Packet> packets;
    //From packet 0 data
    unsigned char dataID;
    unsigned char modeID;
    
    Frame();
    ~Frame();
    int validity();
    
};

struct GTIvec {
    vector <double> tstart;
    vector <double> tstop;

    int reset();
};

struct BTIvec {
    vector <double> tstart;
    vector <double> tstop;

    int reset();
};

class L1evtHandler{
public:
    string l1Filename;
    fitsfile *fptr;
    unsigned char** pktArray;
    int pktArrayRows;
    GTIvec gti;
    BTIvec bti;
    
    L1evtHandler(string l1Filename);
    int read_l1evt_file_pkts(string extname, long startRowNumber=1, long endRowNumber=-1);
    long get_l1evt_ext_nrows(string extname);
    int read_l1evt_file_gti();
    int read_l1evt_file_bti();
    ~L1evtHandler();
    
    //SETTERS
    /**
     * Clears packet array class object. 
     * @return 
     */
    int clear_pktArray(); 
    
    //GETTERS
    unsigned char** get_pktArray(){return this->pktArray;}
    int get_pktArrayRows(){return pktArrayRows;}
    
    
};



class QuadEvtData{
public:
    vector <double> time;
    vector <double> cztseccnt;
    vector <unsigned short> cztntick;
    vector <unsigned short> pha;
    vector <unsigned char> detid;
    vector <unsigned char> pixid;
    vector <unsigned char> alpha;
    vector <unsigned short> veto;
    vector <unsigned char> detx;
    vector <unsigned char> dety;
    vector <unsigned char> dataID;
    vector <unsigned char> modeID;
    vector <unsigned char> quadID;

    /**
     *Ajay Vibhute, Jan 15 2015. 9:30AM.
     *Added to take of new onboard software which cleans bunches on board.
     */

    //vectors to store bunch information
    vector <double> bunch_time;
    vector <long>evt_row_num;
    vector <unsigned char> time_dfs;//Time difference between first event and second event in a bunch
    vector <unsigned char>time_dsl;//Time difference between second event and last event in a bunch
    vector <unsigned char>num_bunchevents;//Number of events in a bunch minus one.
    vector <unsigned char>detid1;//Detector ID of event 3
    vector <unsigned char>detid2;//Detector ID of event 4
    vector <unsigned char>detid3;//Detector ID of event 5
    vector <unsigned char>detid4;//Detector ID of event 6
    //Ajay Vibhute, 21 May 
    //For new noise reduction code adding additional information
    vector <unsigned char>DetId_fevt;//Detector ID of first event
    vector <unsigned char>PixId_fevt;//Pixel ID for first event
    vector <unsigned char>DetId_sevt;//Detector ID of second event
    vector <unsigned char>PixId_sevt;//Pixel ID for second event
    vector <unsigned char>DetId_tevt;//Detector ID of third event
    vector <unsigned char>PixId_tevt;//Pixel ID for third event
    //SETTERS
    int reset();
    
    //GETTERS
    long get_nrows(){return cztntick.size();}
    
};

class HeaderData{
public:
    vector <double> time;
    vector <double> cztseccnt;
    vector <unsigned char> dataID;
    vector <unsigned char> modeID;
    vector <unsigned short> wordCount;
    vector <unsigned short> frameNo;
    vector <unsigned short> statusFlag;
    vector <unsigned short> writePktNo;
    vector <unsigned short> readPktNo;
    vector <unsigned int> command;
    vector <double> laxpcTime;
    vector <unsigned short> dcnt;
    vector <unsigned char> errorCount;
    vector <unsigned char> errorFlag;
    vector <unsigned char> bootPageNo;
    vector <unsigned char> cztNo;
    vector <unsigned char> cztStatus;
    vector <unsigned char> evtReadMode;
    vector <unsigned char> cmdStatus;
    vector <unsigned char> bufferNo;
    vector <unsigned char> baseAdd;
    vector <unsigned char> vetoSpecRange;
    vector <unsigned char> channelNo;
    vector <unsigned short> ADCoutput;
    vector <unsigned short> cmd1sec;
    vector <unsigned short> alphaCount;
    vector <unsigned short> vetoCount;
    vector <unsigned short> cztCount_lt_uld;
    vector <unsigned short> cztCount_gt_uld;
    vector <unsigned short> cztdataRead;
    
    //SETTERS
    int reset();
};

class VetoSpectrumData{
public:
    vector <bool> vetoSpecRange;
    vector < vector < unsigned short> > vetoSpectrum;
    vector <double> time;
    vector <double> cztseccnt;
    vector <unsigned char> quadID;

    //GETTERS
    long get_nrows(){return time.size();}
    //SETTERS
    int reset();
};

class SSMdata{
public:
    vector <double> time;
    vector <double> cztseccnt;
    vector < vector <unsigned short> > vetoSpec; 
    vector < vector <unsigned short> > cztSpec; 
    vector < vector <unsigned short> > cztSpecVeto; 
    vector < vector <unsigned short> > cztSpecAlpha;
    vector <unsigned char> quadID;

    //GETTERS
    long get_nrows(){return time.size();}    
    //SETTERS
    int reset();
};

class TemperatureData{
public:
    vector <double> time;
    vector <double> cztseccnt;
    vector < vector <float> > temperature;
    vector <unsigned char> quadID;

    int remove_junk_temperature();
    //GETTERS
    long get_nrows(){return time.size();
    }
    //SETTERS
    int reset();
};
//Stores l1 event data for a single qudrant 
class L1evtData{
public:
    bool currentFrameValid; //if current frame is not valid then all subsequent packets corresponding
                            // to that frame is rejected.
    int pktArrayRows;
    unsigned char** pktArray; //packet array containing 2048xpktArrayRows
    int lastFrameIndex;
    int currentFrameIndex;
    unsigned char lastPktNo;
    vector <Frame> l1evtFrames;
    vector <Packet> l1rejectedPackets;
    vector <unsigned short> missingFrames; 
    vector <unsigned short> packetsRejected; //Packets are valid but rejected because no frame
                                             //is present to hold them.
    QuadEvtData quadEventData; //stores event data
    VetoSpectrumData vsd; //veto spectrum data for each frame
    SSMdata ssmd; //SSM data
    TemperatureData tempd; //Temperture Data
    HeaderData hdrd;
    int totalPackets;
    int totalValidPackets;
    int totalInvalidPackets;
    int totalFrames;
    int totalMissingFrames;
    int countP0; //number of packets with 0 packet number in current dataset.
	
	long total_events_extracted; // Added by Mithun (19/05/16) to get correct row number for bunches

    L1evtData();
    
    /**
     * Decodes data and extract frame header, packet header & febhdr information. This is the function
     * which perform basic function of segregation when called in sequence.
     * @param wordArray
     * @param packetIndex
     * @return 
     */
    int decode(unsigned short wordArray[PKTSIZE], int packetIndex, bool displayInfo=false);

    /**
     * segregates data into frames with each packet kept in its respective frame.
     * @return 
     */
    int segregate_data(int &lastpackets,bool displayInfo=false);
    
    /**
     * Extracts event information from segregated frames
     * @return 
     */
    int extract_data(string TCTfilename);
    
    int find_missing_frames();
    
    int get_total_valid_packets();
    
    //Display
    int display_l1data_info();
    int display_l1extracted_data();
    //Getters
    
    //Setters
    int reset();

    
    int set_pktArray(unsigned char**pktArray, int pktArrayRows){
        int i=0, j=0;
        this->pktArrayRows = pktArrayRows;
        this->pktArray = allocateMemory<unsigned char>(this->pktArrayRows,PKTSIZE);
        DLOG(INFO) << "Packet array rows: " << pktArrayRows;
        for(j=0; j<pktArrayRows; j++){
            for(i=0; i<PKTSIZE; i++){
                this->pktArray[j][i]=pktArray[j][i];
            }
        }
        return (EXIT_SUCCESS);
    }
};

#endif /*L1EVTDECODE_H*/





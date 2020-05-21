/**
 * @file decode.h
 * @author Tanul Gupta, Preeti Tahlani
 * @date created on July 6, 2012, 5:22 PM; modified on April 7, 2015, 12:39 PM
 * @brief decode header file
 * @details This class is designed for extracting data from level1 packets
 * @version 0.2.0
 */

#ifndef DECODE_H
#define	DECODE_H

#define PACKETSIZE 2048
#define MAXEVENTS_IN_ONE_PACKET 340 
#define VETOSPEC_SIZE 256
#define NUM_DET_PER_QUADRANT 16
#define FEBSAA 100           //number of feb headers for 100 second data
#define MAX_PACKETS_FRAME 10
#define MAXFRAMEDATASIZE 1008+1020*(MAX_PACKETS_FRAME-1)


/**
 * This function converts character packet array into unsigned short packet array.
 * In this software level-1 data is obtained in packets of character array having
 * a size of 2048 elements. This function converts it into an unsigned short array
 * of 1024 elements.
 * @param data: packet array (Bytes)
 * @param size: size of packet array
 * @param word_data: [OUTPUT] packet array (unsigned short)
 * @param bigendian: set to 1 if word_data to be formed is in bigendian otherwise 0
 */
void toWords(unsigned char *data,long size,unsigned short *word_data,bool bigendian);

struct FEBheader{
    
    //----------Word 0 Status word-----------
    //D15-D12 = 0
    bool VetoSpecRange;          //D11
    char BaseAdd;                //D10-D8
    bool BufferNo;               //D7
    bool CmdStatus;              //D6
    bool EventReadMode;          //D5
    bool CZTstatus;              //D4
    char CZTNo;                  //D3-D0
    //----------Word 1 HK word--------------
    //D15 = 0
    char ChannelNo;              //D14-D12
    short ADCoutput;             //D11-D0
    //----------Word 2 Command-----------
    unsigned short cmd1sec;      //command in 1 second
    //----------Word 3 Alpha count--------
    unsigned short AlphaCount;
    //----------Word 4 Veto count---------
    unsigned short VetoCount;
    //----------Word 5 --------------
    unsigned short CZTcount_lt_ULD;    //CZT count < ULD
    //----------Word 6--------------
    unsigned short CZTcount_gt_ULD;    //CZT count >=ULD
    //----------Word 7--------------
    unsigned short CZTdataread;        //data read from CZT
    //----------Word 8-23 CZT temperatures---------
    unsigned short temperature[NUM_DET_PER_QUADRANT];   //temperature of 16 CZT detectors of 1 quadrant
    
    //for SAA mode FEB header data
      //----------Word 0 Status word-----------
    //D15-D12 = 0
    bool VetoSpecRange1[FEBSAA];          //D11
    char BaseAdd1[FEBSAA];                //D10-D8
    bool BufferNo1[FEBSAA];               //D7
    bool CmdStatus1[FEBSAA];              //D6
    bool EventReadMode1[FEBSAA];          //D5
    bool CZTstatus1[FEBSAA];              //D4
    char CZTNo1[FEBSAA];                  //D3-D0
    //----------Word 1 HK word--------------
    //D15 = 0
    char ChannelNo1[FEBSAA];              //D14-D12
    short ADCoutput1[FEBSAA];             //D11-D0
    //----------Word 2 Command-----------
    unsigned short cmd1sec1[FEBSAA];      //command in 1 second
    //----------Word 3 Alpha count--------
    unsigned short AlphaCount1[FEBSAA];
    //----------Word 4 Veto count---------
    unsigned short VetoCount1[FEBSAA];
    //----------Word 5 --------------
    unsigned short CZTcount_lt_ULD1[FEBSAA];    //CZT count < ULD
    //----------Word 6--------------
    unsigned short CZTcount_gt_ULD1[FEBSAA];    //CZT count >=ULD
    //----------Word 7--------------
    unsigned short CZTdataread1[FEBSAA];        //data read from CZT
    
    void read(unsigned short *data,int size,int modeid,int dataid); 
};

class Packet{
private:
    /**
     * \brief Checks the validity of packet data by ensuring following things:
     * 1. SYNC1=F9A4
     * 2. SYNC2=2BB1
     * 3. DATAid lies between [0,16]
     * 4. mml lies between [0,3]
     * 5. modeid lies between [0,15]
     * 6. wcnt(word count) lies between [0,1020]
     * 7. if packet number is 0, then framecounter lies between [0,65535]
     * @return true:valid packet false:invalid packet
     */
    bool validity();   //Function to check packet validity
   
public:    
    bool isFirstpacket;
    bool isValid;
    //Header Fields
    unsigned short SYNC1,SYNC2;
    unsigned char mml;         //memory level
    unsigned char pktNo;
    unsigned char dataID;
    unsigned char modeID;
    unsigned short wcnt;   
    unsigned short framecounter;
    unsigned short statusFlag;
    unsigned short writePktNo;
    unsigned short readPktNo;
    int command;
    int secondcount;
    long laxpctime;
    long cztTime;
    unsigned short dcnt;
    unsigned char errorcount;
    unsigned char errorflag;
    unsigned char bootpageno;
    unsigned short quadID;
    bool isVSD;                  //whether vetospec disabled
    bool isReduced;    //whether event data in reduced mode
    
    unsigned short *hbtdata; //to store event data
    FEBheader febhdr;
      
    /**
     * \brief Function to decode data and put in respective fields:
     * For each packet following information is extracted:
     * 1. SYNC1
     * 2. SYNC2
     * 3. mml
     * 4. pktNo
     * 5. dataID
     * 6. modeID
     * 7. wcnt (Number of valid data words in the packet)
     * 8. quadID (0,1,2,3 for dataIDs 0-11 else -1)
     * 9. isValid (true if satisfies validity criteria as described in private function validity() of this class| else false)
     * 10.hbtdata: contains event data
     * Now, if pktNo==0 i.e. if the packet is first packet, then these values are extracted as well:
     * 11.framecounter 
     * 12.statusFlag
     * 13.writePktNo
     * 14.readPktNo
     * 15.command
     * 16.errorflag
     * 17.isFirstpacket=true;
     * 18.febhdr
     * @param pkt
     * @param bigendian
     */
    void decode(unsigned char pkt[PACKETSIZE],bool bigendian);
};

/** \brief This class will handle Frames and store important frame data from packets.
 */
class Frame{
public:
    static int count;                //to store total count of frames
    unsigned short framedata[MAXFRAMEDATASIZE];     //buffer to store frame data
    long size;                                     //frame data size, size if framedata
    unsigned short frameno;                       //frame number 
    long startpkt,endpacket;      //wrt to science data file
    unsigned char dataID,modeID,quadID;    //data information
    int secondcount;                       //time count for second in CZTI clock
    long laxpctime;
    long frameindex;                 //index for framedata
    unsigned short temperature[NUM_DET_PER_QUADRANT];   //temperature of 16 CZT detectors of 1 quadrant
    bool isValid;                   //whether frame is valid
    
    
    /**
     * This function set values of following frame values as taken from packet p:
     * 1. dataID
     * 2. frameno :frame counter
     * 3. modeID
     * 4. quadID
     * 5. laxpctime
     * 6. secondcount
     * 7. frameindex: set to 0
     * 8. size: initialize to 0; [It is equal to wcnt as decoded by packet]
     * 9. isValid: 0: Valid; 1:Invalid - frame is valid if the packet is valid
     * @param p: Packet data
     */
    void setValues(Packet &p);     //should be called only for 1st packet
}; 

#endif	/* DECODE_H */


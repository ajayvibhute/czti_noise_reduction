

#include<iostream>
#include <vector>
#include "glog/logging.h"
#include"decode.h"

using namespace std;

void toWords(unsigned char *data,long size,unsigned short *word_data,bool bigendian){
    int i,j,msb,lsb;
    for(i=0,j=0;i<size;j++,i=i+2){
        if(bigendian){
            msb=data[i];
            lsb=data[i+1];
        }
        else{
            lsb=data[i];
            msb=data[i+1];
        }
        word_data[j]=(unsigned short)((msb & 0xff)<<8 | (lsb & 0xff));
    }
}

void FEBheader::read(unsigned short *data,int size,int mode,int dataid){
    DLOG(INFO)<<"READING FEB HEADER WITH DATAID: "<< dataid << " MODE:" << mode;
    int startword,startword2,word_FEB,word_FEB2;
    if(dataid>=4 && dataid<=7){
        startword=16;
        startword2=200;
        word_FEB=24;
        word_FEB2=800;
    }
    else if(dataid>=8 && dataid<=11){
        startword=16;
        startword2=-1;
        word_FEB=24;
        word_FEB2=-1;
    }
    else if(mode>=0 && mode<=7 && dataid>=0 && dataid<=3){
        startword=16;
        startword2=-1;
        word_FEB=8;
        word_FEB2=-1;
    }
    else{
        LOG(INFO) << "Data ID:"<<dataid<<"  Mode ID:"<<mode;
        LOG(INFO)<<"***This combination of mode id and data id is not handled for reading FEB header***";
        return;
    }
    int i=0,j=0;

    if(word_FEB2==800){
        for(i=0,j=startword2;i<FEBSAA;i++,j=j+8){
            VetoSpecRange1[i]=(bool)((data[j]>>11) & (0x1));
            BaseAdd1[i]=(char)((data[j]>>8) & (0x7));
            BufferNo1[i]=(bool)((data[j]>>7) & (0x1));
            CmdStatus1[i]=(bool)((data[j]>>6) & (0x1));
            EventReadMode1[i]=(bool)((data[j]>>5) & (0x1));
            CZTstatus1[i]=(bool)((data[j]>>4) & (0x1));
            CZTNo1[i]=(char)(data[j] & 0xf);

            ChannelNo1[i]=(char)((data[j+1]>>12) & (0x7));
            ADCoutput1[i]=(unsigned short)(data[j+1] & (0xfff));

            cmd1sec1[i]=data[j+2];
            AlphaCount1[i]=data[j+3];
            VetoCount1[i]=data[j+4];
            CZTcount_lt_ULD1[i]=data[j+5];
            CZTcount_gt_ULD1[i]=data[j+6];
            CZTdataread1[i]=data[j+7];
        }

    }
    
    VetoSpecRange=(bool)((data[startword]>>11) & (0x1));
    BaseAdd=(char)((data[startword]>>8) & (0x7));
    BufferNo=(bool)((data[startword]>>7) & (0x1));
    CmdStatus=(bool)((data[startword]>>6) & (0x1));
    EventReadMode=(bool)((data[startword]>>5) & (0x1));
    CZTstatus=(bool)((data[startword]>>4) & (0x1));
    CZTNo=(char)(data[startword] & 0xf);

    ChannelNo=(char)((data[startword+1]>>12) & (0x7));
    ADCoutput=(unsigned short)(data[startword+1] & (0xfff));

    cmd1sec=data[startword+2];
    AlphaCount=data[startword+3];
    VetoCount=data[startword+4];
    CZTcount_lt_ULD=data[startword+5];
    CZTcount_gt_ULD=data[startword+6];
    CZTdataread=data[startword+7];

    if(word_FEB==24){
        for(i=0;i<NUM_DET_PER_QUADRANT;i++){
            temperature[i]=data[startword+8+i];
        }
    }

} 

void Packet::decode(unsigned char pkt[PACKETSIZE],bool bigendian){
    unsigned short data[PACKETSIZE/2];
    toWords(pkt,PACKETSIZE,data,bigendian);
    isFirstpacket=false;
    SYNC1=data[0];
    SYNC2=data[1];
    mml=(unsigned char)((data[2]>>13) & 0x3);
    dataID=(unsigned char)((data[2]>>8) & 0x1f);
    pktNo=(unsigned char)((data[2]>>4) & 0xf);
    modeID=(unsigned char)(data[2] &0xf);                                                                                                                                                                                                                                                                                                                                                              
    wcnt=data[3];
    
    if(pktNo==0){                        // if the packet is first packet
        framecounter=data[4];
        statusFlag=data[5];
        writePktNo=data[6];
        readPktNo=data[7];
        command=(int)(((data[9]<<16) & 0xffff0000) | (data[8] & 0xffff));
        secondcount=(int)(((data[10]<<16) & 0xffff0000) | (data[11] & 0xffff));
        laxpctime = ((unsigned int)((data[14]<<12)&0xff000 |(data[13]>>4)) + ((data[12]>>10)|((data[13]&4)<<6))/1000 + ((data[12] & 0x03ff)/1000000));
        dcnt=data[15] & 0x3fff;
        bootpageno=(data[15]>>14) &0x3;
        errorcount=(data[14]>>12) & 0xf;
        errorflag=(data[14]>>8) & 0xf;
        isFirstpacket=true;
        wcnt=wcnt-12;
     
       }
    quadID=dataID<12 ? dataID%4 : -1;
    isValid=validity(); //Checking the validity of packet
     
    if(isFirstpacket)  febhdr.read(data,PACKETSIZE/2,modeID,dataID);
 
    hbtdata=new unsigned short[wcnt];
    int i,j;
    int datastart;
   
    if(isFirstpacket){
        datastart=16;
    }
    else{
        datastart=4;
    }
    
    if(isValid){
    for(i=0,j=datastart;i<wcnt;i++,j++)   hbtdata[i]=data[j];
    }

}

bool Packet::validity(){
    bool flag=true;
   
    flag=(SYNC1==0xf9a4) ? (flag & true) : (flag & false);
    if(flag==false) { //cerr<<"\nInvalid SYNC1:"<<SYNC1;   
        return flag; }
    
    flag=(SYNC2==0x2bb1) ? (flag & true) : (flag & false); 
    if(flag==false) { //cerr<<"\nInvalid SYNC2:"<<SYNC2;    
        return flag; }
    
    flag=(dataID>=0 && dataID<16) ? (flag & true) : (flag & false);
    if(flag==false) { //cerr<<"\nInvalid Data ID:"<<dataID; 
        return flag; }
    
    flag=(mml>=0 && mml<4) ? (flag &true) : (flag & false);
    if(flag==false) { //cerr<<"\nInvalid memory level:"<<mml; 
        return flag; }
    
    flag=(modeID>=0 && modeID<16) ? (flag & true) : (flag & false);
    if(flag==false) { //cerr<<"\nInvalid Mode ID:"<<modeID; 
        return flag; }
    
    flag=(wcnt>=0 && wcnt<=1020) ? (flag & true) : (flag & false);
    if(flag==false) {   DLOG(INFO)<< "Invalid WCNT:" << wcnt; 
        return flag; }
    
    if(pktNo==0){
       flag=(framecounter>=0 && framecounter<65536) ? (flag & true) : (flag & false);
       if(flag==false) {
           return flag; }

     }
    return flag;
}


//If the first packet is invalid, then full frame is considered invalid
 int Frame::count=0;
 void Frame::setValues(Packet &p){
     this->dataID=p.dataID;
     this->frameno=p.framecounter;
     this->modeID=p.modeID;
     this->quadID=p.quadID;
     this->laxpctime=p.laxpctime;
     this->secondcount=p.secondcount;
     this->frameindex=0;
     this->size=0;
     this->isValid=p.isValid;
     
     if (p.dataID>=4 && p.dataID<=11){
         for(int i=0; i < NUM_DET_PER_QUADRANT; i++){
             this->temperature[i]=p.febhdr.temperature[i];
         }
     }
 }

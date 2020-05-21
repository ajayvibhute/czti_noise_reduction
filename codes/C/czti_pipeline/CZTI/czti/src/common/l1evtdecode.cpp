#include "l1evtdecode.h"
#include "level1handler.h"

using namespace std;

void toWords(unsigned char *data, long size, unsigned short *word_data, bool bigendian) {
    int i, j, msb, lsb;
    for (i = 0, j = 0; i < size; j++, i = i + 2) {
        if (bigendian) {
            msb = data[i];
            lsb = data[i + 1];
        } else {
            lsb = data[i];
            msb = data[i + 1];
        }
        word_data[j] = (unsigned short) ((msb & 0xff) << 8 | (lsb & 0xff));
    }
}

//L1evt CLASS
L1evtHandler::L1evtHandler(string l1Filename){
    int status=0;
    this->l1Filename = l1Filename;
    //opening this file and storing in a fitsfile pointer
    fits_open_file(&fptr, this->l1Filename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level 1 event file: " << l1Filename;
        fits_report_error(stderr, status);
        exit(EXIT_FAILURE);
    }
    
    
}

L1evtHandler::~L1evtHandler(){
    int status=0;
    //closing level-1 event file
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-1 event file: " << l1Filename;
        fits_report_error(stderr, status);
        exit(EXIT_FAILURE);
    }
}

//SETTERS
int L1evtHandler::clear_pktArray() {
    for (long i = 0; i < pktArrayRows; i++) {
        delete [] pktArray[i];
    }
    delete[] pktArray;
    return EXIT_SUCCESS;
}

long L1evtHandler::get_l1evt_ext_nrows(string extname){
    int status=0;
    long nrows=0;
    int colnum=0;
    fits_movnam_hdu(fptr, BINARY_TBL, (char*) extname.c_str(), 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU" << extname << " of Level 1 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_colnum(fptr, CASEINSEN, "DataArray", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for Data Array column.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fptr, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows for Data Array column of HDU " << extname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return nrows;
    
}
int L1evtHandler::read_l1evt_file_pkts(string extname, long startRowNumber, long endRowNumber){
    int status=0;
    int colnum=0;
    int firstelem=1;
    long i=0, j=0;
    long nelements=0;
    long nrows=0;
    unsigned char* dataArray;
    
    //Reading Packet Arrays
    fits_movnam_hdu(fptr, BINARY_TBL, (char*) extname.c_str(), 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU" << extname << " of Level 1 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_colnum(fptr, CASEINSEN, "DataArray", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for Data Array column.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fptr, &nrows, &status);
    if (status) {
        LOG(ERROR) <<"Error in getting number of rows for Data Array column of HDU " << extname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    if(endRowNumber>(nrows)){
        endRowNumber = nrows;
    }
    if (endRowNumber==-1){
        nelements= (nrows - startRowNumber + 1) * PKTSIZE;
        dataArray = new unsigned char [nelements];
        fits_read_col(fptr, TBYTE, colnum, startRowNumber, 1, nelements, NULL, dataArray, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading dataArray column.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    else {
        nelements = (endRowNumber-startRowNumber+1) * PKTSIZE;
        dataArray = new unsigned char [nelements];
        fits_read_col(fptr, TBYTE, colnum, startRowNumber, 1, nelements, NULL, dataArray, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in reading dataArray column.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        } 
    }
    pktArrayRows = nelements/PKTSIZE;
    pktArray= allocateMemory<unsigned char>(pktArrayRows, PKTSIZE);
    for(i=0; i< pktArrayRows; i++ ){
        for(j=0; j<PKTSIZE; j++){
            pktArray[i][j]=dataArray[i*PKTSIZE + j];
        }
    }
    delete[] dataArray;
    //Packet Arrays read
    
    return EXIT_SUCCESS;
}

int L1evtHandler::read_l1evt_file_gti() {
    int status=0;
    string extname="GTI";
    //Reading GTI Arrays
    fits_movnam_hdu(fptr, BINARY_TBL, (char*) extname.c_str(), 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU" << extname << " of Level 1 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    
    
    if(read_fits_column(fptr, "START", TDOUBLE, 1, 1, -1, gti.tstart)){
        LOG(ERROR) << "Error in reading GTI extension of level-1 file " << l1Filename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fptr, "STOP", TDOUBLE, 1, 1, -1, gti.tstop)){
        LOG(ERROR) << "Error in reading GTI extension of level-1 file " << l1Filename;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;

}
int L1evtHandler::read_l1evt_file_bti() {
    int status=0;
    string extname="BTI";
    //Reading GTI Arrays
    fits_movnam_hdu(fptr, BINARY_TBL, (char*) extname.c_str(), 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU" << extname << " of Level 1 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    
    
    if(read_fits_column(fptr, "START", TDOUBLE, 1, 1, -1, bti.tstart)){
        LOG(ERROR) << "Error in reading GTI extension of level-1 file " << l1Filename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fptr, "STOP", TDOUBLE, 1, 1, -1, bti.tstop)){
        LOG(ERROR) << "Error in reading GTI extension of level-1 file " << l1Filename;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;

}

//L1evt END
//Temperature Data Class
int TemperatureData::remove_junk_temperature() {
    int status=0;
    ErrorHandler errHandler;
    int itemp=0;
    int idet=0; //0-15
    float modeTemperature=0.0;
    map<float,long> modemap;
    
    if(temperature.size()<=0){
        errHandler.severity = errWARNING;
        errHandler.errorStatus = WARNING;
        errHandler.errorMsg = "NO TEMPERATURE DATA AVAILABLE.";
        throw errHandler;
    }
    
    for(itemp=0; itemp<temperature.size(); itemp++){
        //If first temperature read from the detector is zero, replace it 
        //with the Mode of other 15 detector temperatures.  
        if (temperature[itemp][0]==0 or temperature[itemp][0]>50){
            //Calculating mode of temperature values to replace temperature of 0th detector
            try{calculate_mode(temperature[itemp], &modeTemperature, &modemap);}
            catch(ErrorHandler errHandler){
                throw errHandler;
            }
            temperature[itemp][0]=modeTemperature;
        }
        
        // If any other temperature other than the first read temperature is 
        // either zero or junk, replace it with the previous temperature.
        for(idet=1; idet<NUM_DET_PER_QUAD; idet++){
           if (temperature[itemp][idet] == 0 or temperature[itemp][0] > 50) {
               temperature[itemp][idet]=temperature[itemp][idet-1];
           }   
        }
    }
    return status;
}

//Temperature Data Class end

//L1EVTDATA CLASS
L1evtData::L1evtData() {
    currentFrameIndex=-1;
    countP0=0;
    currentFrameValid=false;
    pktArrayRows=0;
    lastFrameIndex=0;
    lastPktNo=0;
    l1evtFrames.clear();
    l1rejectedPackets.clear();
    missingFrames.clear();
    packetsRejected.clear();
    totalPackets=0;
    totalValidPackets=0;
    totalInvalidPackets=0;
    totalFrames=0;
    totalMissingFrames=0;
    countP0=0;
}

//** setters
int L1evtData::reset() {

    //clear pktArray
    for (long i = 0; i < pktArrayRows; i++) {
  	delete [] pktArray[i];
    }   
    delete [] pktArray;
    
	//Clear HBT data in packets
    for(long i=0;i<l1evtFrames.size();i++)
    {
	for(long j=0;j<l1evtFrames[i].packets.size();j++)
	    delete [] l1evtFrames[i].packets[j].hbtdata;
   	l1evtFrames[i].packets.clear();
    }

    currentFrameIndex = -1;
    countP0 = 0;
    currentFrameValid = false;
    pktArrayRows = 0;
    lastFrameIndex = 0;
    lastPktNo = 0;
    l1evtFrames.clear();
    l1rejectedPackets.clear();
    missingFrames.clear();
    packetsRejected.clear();
    totalPackets = 0;
    totalValidPackets = 0;
    totalInvalidPackets = 0;
    totalFrames = 0;
    totalMissingFrames = 0;
    countP0 = 0;
    quadEventData.reset();
    vsd.reset();
    ssmd.reset();
    hdrd.reset();
    tempd.reset();
    
	return EXIT_SUCCESS;    
}

//** basic functions
int L1evtData::decode(unsigned short wordArray[PKTSIZE], int packetIndex, bool displayInfo){
    int status=0;
    Packet p;
    Frame f;
    int i=0, j=0;
    int datastart=0;
    p.isFirstPacket=false;
    p.SYNC1=wordArray[0];
    p.SYNC2=wordArray[1];
    p.mml=(unsigned char)((wordArray[2]>>13) & 0x3);
    p.dataID=(unsigned char)((wordArray[2]>>8) & 0x1f);
    p.pktNo=(unsigned char)((wordArray[2]>>4) & 0xf);
    p.modeID=(unsigned char)(wordArray[2] & 0xf);                                                                                                                                                                                                                                                                                                                                                              
    p.wcnt=wordArray[3]; 
    
    if(p.pktNo==0){
        f.frameno=wordArray[4];
        p.frameno=f.frameno;
        f.dataID=p.dataID;
        f.modeID=p.modeID;
        f.statusFlag=wordArray[5];
        f.writePktNo=wordArray[6];
        f.readPktNo=wordArray[7];
        f.command=(int)(((wordArray[9]<<16) & 0xffff0000) | (wordArray[8] & 0xffff));
        f.secondcount=(int)(((wordArray[10]<<16) & 0xffff0000) | (wordArray[11] & 0xffff));

	/*f.secondcount-=2;
	f.secondcount-=((int)(f.secondcount/1048576)*1048576.0);
	*/
	f.laxpctime = ((double)((wordArray[14]<<12)&0xff000 |(wordArray[13]>>4)) + ((double)((wordArray[12]>>10)|((wordArray[13]&15)<<6)))/1000.0 + ((double)((wordArray[12] & 0x03ff))/1000000.0));
        f.dcnt=wordArray[15] & 0x3fff;
        f.bootpageno=(wordArray[15]>>14) &0x3;
        f.errorcount=(wordArray[14]>>12) & 0xf;
        f.errorflag=(wordArray[14]>>8) & 0xf;
        f.isvalid = f.validity();
        p.isFirstPacket=true;
        p.wcnt=p.wcnt-12;
        countP0++;
    }
    p.quadID = p.dataID<12 ? p.dataID%4 : -1;
    p.isValid= p.validity(); //checking the validity of packet
    
    if(p.isValid){
        p.hbtdata = new unsigned short[p.wcnt]; //to store HBT data       
        //Extracting FEBHDR
        if(p.isFirstPacket){
            currentFrameValid = true;
            p.febhdr.read(wordArray, PKTSIZE/2, p.modeID, p.dataID);
            this->currentFrameIndex++;
            //Extracting HBT data for frame number 0
            for (i = 0, j = 16; i < p.wcnt; i++, j++) {
                p.hbtdata[i] = wordArray[j];
            }
            
            //This code snippet rejects packets with non-zero packet number in the beginning of the data
            //i.e. before occurrence of first 0 pkt number.
            //At the same time only valid packets and frames are accounted for.
            f.packets.push_back(p);
            l1evtFrames.push_back(f);
            
        }
        else if(p.isFirstPacket==false){
            //Extracting HBT data for normal packets
            for (i = 0, j = 4; i < p.wcnt; i++, j++) {
                p.hbtdata[i] = wordArray[j];
            }
        }
    }
    else {
	LOG(INFO)<<"************PACKET REJECTION DUE TO SYNC WORDS***************";    
        packetsRejected.push_back(packetIndex); //Invalid packets|frames are pushed back here
        l1rejectedPackets.push_back(p);
    }

    if (l1evtFrames.size() > 0) {
        if (p.isFirstPacket == false &&
                p.pktNo < (l1evtFrames[currentFrameIndex].packets.back()).pktNo) {
            currentFrameValid = false;
        }
    }
    
    if (currentFrameValid && p.isFirstPacket==false) {
        l1evtFrames[currentFrameIndex].packets.push_back(p);
    }

    //Packets that are invalid due to 
    if (currentFrameValid == false) {
	LOG(INFO)<<"Frame rejection";    
        packetsRejected.push_back(packetIndex);
        l1rejectedPackets.push_back(p);
    }
    
    
    //DISPLAY BASIC INFO
    if(displayInfo==true){
    DLOG(INFO) << "Pkt Index: " << packetIndex << " PktNo: " << (int) p.pktNo << 
            " FrameNo.: " << currentFrameIndex << " FValid: " << f.isvalid 
            << " PValid: "<< p.isValid; 
        DLOG(INFO) << "Pkt Index: " << packetIndex << " PacketNo: " << (int)p.pktNo << " DataID: " << (int)p.dataID << "  ModeID:" << (int)p.modeID;
    }
    
    return EXIT_SUCCESS;
}

int L1evtData::segregate_data(int &lastpackets,bool displayInfo) {
    int status=0;
    int i=0;
    int last_frame_pkts=0;
    unsigned short* wordArray;
    bool bigendian=1;
    wordArray = new unsigned short[PKTSIZE];
    l1evtFrames.clear();
    l1rejectedPackets.clear();
    for(i=0; i<pktArrayRows; i++){
        toWords(pktArray[i],PKTSIZE,wordArray, bigendian);
        if(decode(wordArray,i, displayInfo)){
            LOG(ERROR) << "Error in decoding packet " << i  << ".";
            return EXIT_FAILURE;
        }
    }

    LOG(INFO)<<"*************************FRAME INDEX "<<currentFrameIndex<<" SIZE "<<l1evtFrames[currentFrameIndex].packets.size();	

    if (lastpackets==0)
    {
    	last_frame_pkts=l1evtFrames[currentFrameIndex].packets.size();

    	l1evtFrames.pop_back();	
    	currentFrameIndex-=1;
    }
     
    lastpackets=last_frame_pkts;

    totalPackets = pktArrayRows;
    totalInvalidPackets = packetsRejected.size();
    get_total_valid_packets();
    find_missing_frames();
    
    delete[] wordArray;
    return EXIT_SUCCESS;
}

/*
 *Ajay Vibhute,
 *25 Jan 2016, Modified code to work with new on-bord software to take care of bunches
 *Here I am decoding bunch information and adding them to different file.
 * */
int L1evtData::extract_data(string TCTfilename) {
    int status=0;
    int nevents=0, feb_veto_size=0, eventsize=0, vetospecstart=-1;
    int pktNo;
    int counter=0;
    int invalidForEventData=0;
    int startIndex=0, endIndex=0;
    int i=0, j=0;
    int nPkts;
    int bunch_sync=0,isbunch=0;//if event is bunch, bunch_sync will be 7 and isbunch will be 0. 
    Frame frame;
    double instrumentTime;
    double interpolatedTime; //interpolated time
    Tct TCT; //for time correlation
    Packet p1, p2;
    int iFrame=0, iPacket=0, iEvent=0, index=0, iveto=0, iSSM=0; //status variables
    long tempcount=0;
    double cztseccnt=0.0, time=0.0;
    unsigned short cztntick=0, pha=0;
    unsigned char detid=0, pixid=0, alpha=0, detx=0, dety=0, dataID, modeID, quadID;
    unsigned short veto=0;
    vector <unsigned short> tempVetospectrum, tempCztspec, tempCztspecVeto, tempCztspecAlpha;
    vector <float> tempTemperature;
    bool vetoSpecRange;

    //Initializing veto spectrum vector to be of length 256
    tempVetospectrum.resize(VETOSPEC_SIZE, 0);
    tempCztspec.resize(CZTSPEC_SIZE, 0);
    tempCztspecVeto.resize(CZTSPECVETO_SIZE, 0);
    tempCztspecAlpha.resize(CZTSPECALPHA_SIZE, 0);
    tempTemperature.resize(NUM_DET_PER_QUAD, 0.0);

    //Reading TCT file
    if (TCT.read_tct_file(TCTfilename)) {
        LOG(ERROR) << "*** Error reading TCT file***";
    }
    
    for (iFrame = 0; iFrame < totalFrames; iFrame++) {
        frame = l1evtFrames[iFrame];
        nPkts = frame.packets.size();
        
        //Check based on modeID
        switch (frame.modeID) {
            case 0:
            case 4: feb_veto_size = 240; //Detector data has 240 words for FEB header and Veto Spectrum
                eventsize = 3; //Each event consumes 3 words
                vetospecstart = 8;
                break;
            case 1:
            case 5: feb_veto_size = 8;
                eventsize = 3;
                vetospecstart = -1;
                break;
            case 2:
            case 6: eventsize = 2;
                feb_veto_size = 244;
                vetospecstart = 8;
                break;
            case 3:
            case 7: eventsize = 2;
                feb_veto_size = 8;
                vetospecstart = -1;
                break;
			// Default option added to avoid unexpected errors for other mode data
			// Mithun NPS (17/05/16) 
			default:
				eventsize = 3;
				feb_veto_size = 240;
				vetospecstart = 8;
			    break;
		} 
        
		/************** NOTE ****************************************************
		 * The values added in header file columns is correct for
		 * modeM0-M7 alone. In other cases some of the columns may be irrelevant 
		 * and can have absurd values. All hk parameters available in the 100s 
		 * frames are not decoded and stored in hdr file. 
		*************************************************************************/

        //Extracting Header Data
        hdrd.cztseccnt.push_back((double)frame.secondcount);
        hdrd.dataID.push_back(frame.dataID);
        hdrd.modeID.push_back(frame.modeID);

	int tot_words=0;
	for(i=0;i<nPkts;i++) tot_words+=frame.packets[i].wcnt;
	tot_words-=feb_veto_size;

        hdrd.wordCount.push_back(tot_words/eventsize);
        hdrd.frameNo.push_back(frame.frameno);
        hdrd.statusFlag.push_back(frame.statusFlag);
        hdrd.writePktNo.push_back(frame.writePktNo);
        hdrd.readPktNo.push_back(frame.readPktNo);
        hdrd.command.push_back(frame.command);
        hdrd.laxpcTime.push_back(frame.laxpctime);
        hdrd.dcnt.push_back(frame.dcnt);
        hdrd.errorCount.push_back(frame.errorcount);
        hdrd.errorFlag.push_back(frame.errorflag);
        hdrd.bootPageNo.push_back(frame.bootpageno);
        hdrd.cztNo.push_back(frame.packets[0].febhdr.CZTNo);
        hdrd.cztStatus.push_back(frame.packets[0].febhdr.CZTstatus);
        hdrd.evtReadMode.push_back(frame.packets[0].febhdr.EventReadMode);
        hdrd.cmdStatus.push_back(frame.packets[0].febhdr.CmdStatus); 
        hdrd.bufferNo.push_back(frame.packets[0].febhdr.BufferNo);
        hdrd.baseAdd.push_back(frame.packets[0].febhdr.BaseAdd);
        hdrd.vetoSpecRange.push_back(frame.packets[0].febhdr.VetoSpecRange);
        hdrd.channelNo.push_back(frame.packets[0].febhdr.ChannelNo);
        hdrd.ADCoutput.push_back(frame.packets[0].febhdr.ADCoutput);
        hdrd.cmd1sec.push_back(frame.packets[0].febhdr.cmd1sec);
        hdrd.vetoCount.push_back(frame.packets[0].febhdr.VetoCount);
        hdrd.alphaCount.push_back(frame.packets[0].febhdr.AlphaCount);
        hdrd.cztCount_lt_uld.push_back(frame.packets[0].febhdr.CZTcount_lt_ULD);
        hdrd.cztCount_gt_uld.push_back(frame.packets[0].febhdr.CZTcount_gt_ULD);
        hdrd.cztdataRead.push_back(frame.packets[0].febhdr.CZTdataread);

        //Interpolating time for Header Data
        instrumentTime = (double)frame.secondcount;

        if (TCT.interpolate_time(instrumentTime, interpolatedTime)) {
            LOG(ERROR) << "Error in TCT time interpolation.";
            return EXIT_FAILURE;
        }
        hdrd.time.push_back(interpolatedTime);
        //Header data extracted
        if(frame.modeID > 7 || frame.modeID <0) {
            LOG(INFO) << "Frame number: " << frame.frameno << " with Mode ID " << (int) frame.modeID << " is not valid for event data";
            invalidForEventData++;
            continue;
        }
        if ((frame.dataID >= 4 && frame.dataID <= 7) || (frame.dataID >= 12 && frame.dataID <= 16)) {
            LOG(INFO) << "Frame number: " << frame.frameno << " with Data ID " << (int) frame.dataID << " is not valid for event data";
            invalidForEventData++;
            continue;
        }
      
        if(frame.packets[0].febhdr.EventReadMode==0){
	    LOG(INFO) << "Frame number: "<< frame.frameno << " is in command mode, not valid event data";	
	    invalidForEventData++;
	    continue;
	}	
	   	bunch_sync=0;

		int nframe_bunch=0;
		 
        //Extracting event data for 
        if (frame.dataID >= 0 && frame.dataID <= 3) {
            for(iPacket=0; iPacket < nPkts; iPacket++){
                p1=frame.packets[iPacket];
                
                //Assigning start and end index for event extraction
                if(p1.pktNo==0){
                    nevents = (p1.wcnt - feb_veto_size)/eventsize;
                    startIndex=feb_veto_size;
                }
                else{
                    nevents = p1.wcnt/eventsize;
                    startIndex=0;
                }
                //start and end index for event extraction assigned
                //decoding events for normal and reduced mode
                for (index = startIndex, iEvent = 0; iEvent < nevents; index += eventsize, iEvent++) { //Each event is of three words in mode 0 & 4
                    //for three word event report
                    if (eventsize == 3) {

			if(nframe_bunch==0) bunch_sync=((p1.hbtdata[index]>>13)&0x7);

			if(bunch_sync!=7)
			{
				//Execute if event is not a bunch event
                	        cztseccnt = ((double) p1.hbtdata[index] / 50000.0)+(double) frame.secondcount; 
				//conversion of fractional time into microseconds by multipyitn it by 20 microsec
	                        cztntick = (unsigned short) p1.hbtdata[index];
                        	pha = (unsigned short) (((p1.hbtdata[index + 1] >> 4) & 0xfff) >> 2);
                	        detid = (unsigned char) (p1.hbtdata[index + 1] & 0xf);
        	                pixid = (unsigned char) ((p1.hbtdata[index + 2] >> 8) & 0xff);
        	                veto = (unsigned short) ((p1.hbtdata[index + 2] >> 1) & 0x7f);
	                        alpha = (unsigned char) (p1.hbtdata[index + 2] & 0x1);
			}
			else if(bunch_sync==7)
			{
				/*
				cztntick = (unsigned short) p1.hbtdata[index];
				LOG(INFO)<<"$$$$$$$$$$$$Found bunch sync 7...."<<iPacket<<" index "<<index <<" frmae  "<<iFrame<<" "<<cztntick;
				printf("SECCNT %lf\n",cztseccnt);
				*/
			}
			//execute if event is bunch


                    }
                    //For two word event report 
                    if (eventsize == 2) {
                        //temp = (p1.hbtdata[index] >> 7) & 0x1ff; //extracted 9 bits for time
                        cztntick = (unsigned short) ((p1.hbtdata[index] &0xff80) >> 7);
                        cztseccnt = ((double) cztntick / 50000.0)+(double) frame.secondcount;
                        pha = (unsigned short) ((((int) (p1.hbtdata[index] & 0x7f) << 2) | (int) ((p1.hbtdata[index + 1] >> 14) & 0x3)) << 1);
                        detid = (unsigned char) (p1.hbtdata[index + 1] & 0xf);
                        pixid = (unsigned char) ((p1.hbtdata[index + 1] >> 6) & 0xff);
                        veto = (unsigned short) ((p1.hbtdata[index + 1] >> 5) & 0x1);
                        alpha = (unsigned char) ((p1.hbtdata[index + 1] >> 4) & 0x1);
                    }
		    	
                   		   dataID=p1.dataID;
                   modeID=p1.modeID;
                   quadID=p1.quadID;
			    	
		    if(bunch_sync!=7)
		    {
                    
                    generate_detx_dety(detid, pixid, p1.quadID, detx, dety);
                    quadEventData.cztseccnt.push_back(cztseccnt);
                    quadEventData.cztntick.push_back(cztntick);
                    quadEventData.pha.push_back(pha);
                    quadEventData.detid.push_back(detid);
                    quadEventData.pixid.push_back(pixid);
                    quadEventData.veto.push_back(veto);
                    quadEventData.alpha.push_back(alpha);
                    quadEventData.detx.push_back(detx);
                    quadEventData.dety.push_back(dety);
                    quadEventData.dataID.push_back(dataID);
                    quadEventData.modeID.push_back(modeID);
                    quadEventData.quadID.push_back(quadID);
	
                    //Initializing vectors
                     //Interpolating itme
                    instrumentTime=cztseccnt;
					
                    if(TCT.interpolate_time(instrumentTime, interpolatedTime)){
                        LOG(ERROR) << "Error in TCT time interpolation.";
                        return EXIT_FAILURE;
                    }
					
                    quadEventData.time.push_back(interpolatedTime);

		    }
		    else if (bunch_sync==7&&frame.bootpageno==2) 
		    {

			   	if(nframe_bunch==0)
				{	
			   	//Decoding the bunch information
	 		    //quadEventData.bunch_time.push_back(interpolatedTime);
			    quadEventData.evt_row_num.push_back(total_events_extracted+(long)quadEventData.cztntick.size());
			    quadEventData.time_dfs.push_back(p1.hbtdata[index] &0x0001);
			    quadEventData.time_dsl.push_back(( p1.hbtdata[index] >>1) &0x3f);
			    quadEventData.num_bunchevents.push_back(((p1.hbtdata[index] >>7)&0x3f));
			    
			    pha = (unsigned short) (((p1.hbtdata[index + 1] >> 4) & 0xfff) >> 2);
                	    detid = (unsigned char) (p1.hbtdata[index + 1] & 0xf);
        	            pixid = (unsigned char) ((p1.hbtdata[index + 2] >> 8) & 0xff);
        	            veto = (unsigned short) ((p1.hbtdata[index + 2] >> 1) & 0x7f);
	                    alpha = (unsigned char) (p1.hbtdata[index + 2] & 0x1);
 			    generate_detx_dety(detid, pixid, p1.quadID, detx, dety);
           //                 quadEventData.cztseccnt.push_back(cztseccnt);
             //               quadEventData.cztntick.push_back(cztntick);
 	                    quadEventData.pha.push_back(pha);
          	            quadEventData.detid.push_back(detid);
                    	    quadEventData.pixid.push_back(pixid);
    	                    quadEventData.veto.push_back(veto);
                    	    quadEventData.alpha.push_back(alpha);
	                    quadEventData.detx.push_back(detx);
        	            quadEventData.dety.push_back(dety);
                    quadEventData.dataID.push_back(dataID);
                    quadEventData.modeID.push_back(modeID);
                    quadEventData.quadID.push_back(quadID);

				nframe_bunch=1;

				iEvent++;

				if(iEvent>=nevents) continue;
			    
				index += eventsize;
				}
				if (nframe_bunch==1)
				{	
			    //decoding second event in bunch	

                cztseccnt = ((double) p1.hbtdata[index] / 50000.0)+(double) frame.secondcount;
                cztntick = (unsigned short) p1.hbtdata[index];

			    pha = (unsigned short) (((p1.hbtdata[index + 1] >> 4) & 0xfff) >> 2);
                	    detid = (unsigned char) (p1.hbtdata[index + 1] & 0xf);
        	            pixid = (unsigned char) ((p1.hbtdata[index + 2] >> 8) & 0xff);
        	            veto = (unsigned short) ((p1.hbtdata[index + 2] >> 1) & 0x7f);
	                    alpha = (unsigned char) (p1.hbtdata[index + 2] & 0x1);
 			    generate_detx_dety(detid, pixid, p1.quadID, detx, dety);

			    // Writing time stamp for first event in bunch
				instrumentTime=cztseccnt-quadEventData.time_dfs.back()/50000.0;
                quadEventData.cztseccnt.push_back(cztseccnt-quadEventData.time_dfs.back()/50000.0);
                quadEventData.cztntick.push_back(cztntick-quadEventData.time_dfs.back());

                 if(TCT.interpolate_time(instrumentTime, interpolatedTime)){
                            LOG(ERROR) << "Error in TCT time interpolation.";
                            return EXIT_FAILURE;
                        }
                            quadEventData.time.push_back(interpolatedTime);
                
				//Writing time stamp for second event in bunch
			    quadEventData.cztseccnt.push_back(cztseccnt);
                quadEventData.cztntick.push_back(cztntick);
				instrumentTime=cztseccnt;

                 if(TCT.interpolate_time(instrumentTime, interpolatedTime)){
                            LOG(ERROR) << "Error in TCT time interpolation.";
                            return EXIT_FAILURE;
                        }
                            quadEventData.time.push_back(interpolatedTime);

				quadEventData.bunch_time.push_back(interpolatedTime);

  				// Writing time stamp for last event in the bunch (third event)
				instrumentTime=cztseccnt+quadEventData.time_dsl.back()/50000.0;
	             quadEventData.cztseccnt.push_back(instrumentTime);
                quadEventData.cztntick.push_back(cztntick+quadEventData.time_dsl.back());

                 if(TCT.interpolate_time(instrumentTime, interpolatedTime)){
                            LOG(ERROR) << "Error in TCT time interpolation.";
                            return EXIT_FAILURE;
                        }
                            quadEventData.time.push_back(interpolatedTime);

				// Write other details of second event of bunch
 	                    quadEventData.pha.push_back(pha);
          	            quadEventData.detid.push_back(detid);
                    	    quadEventData.pixid.push_back(pixid);
    	                    quadEventData.veto.push_back(veto);
                    	    quadEventData.alpha.push_back(alpha);
	                    quadEventData.detx.push_back(detx);
        	            quadEventData.dety.push_back(dety);
                    quadEventData.dataID.push_back(dataID);
                    quadEventData.modeID.push_back(modeID);
                    quadEventData.quadID.push_back(quadID);


				nframe_bunch=2;
                iEvent++;

				if(iEvent>=nevents) continue;

			    index += eventsize;
			   	
				}
				if (nframe_bunch==2)
			   	{
			   
			    //decoding last event in bunch
			    quadEventData.detid1.push_back((p1.hbtdata[index] &0xf));
			    quadEventData.detid2.push_back((p1.hbtdata[index]>>4 )&0xf);
			    quadEventData.detid3.push_back((p1.hbtdata[index]>>8 )&0xf);
			    quadEventData.detid4.push_back((p1.hbtdata[index]>>12)&0xf);
			    //Writing third event in bunch
			    pha = (unsigned short) (((p1.hbtdata[index + 1] >> 4) & 0xfff) >> 2);
                	    detid = (unsigned char) (p1.hbtdata[index + 1] & 0xf);
        	            pixid = (unsigned char) ((p1.hbtdata[index + 2] >> 8) & 0xff);
        	            veto = (unsigned short) ((p1.hbtdata[index + 2] >> 1) & 0x7f);
	                    alpha = (unsigned char) (p1.hbtdata[index + 2] & 0x1);
 			   
			    generate_detx_dety(detid, pixid, p1.quadID, detx, dety);
 	                    quadEventData.pha.push_back(pha);
          	            quadEventData.detid.push_back(detid);
                    	    quadEventData.pixid.push_back(pixid);
    	                    quadEventData.veto.push_back(veto);
                    	    quadEventData.alpha.push_back(alpha);
	                    quadEventData.detx.push_back(detx);
        	            quadEventData.dety.push_back(dety);
                	    quadEventData.dataID.push_back(dataID);
	                    quadEventData.modeID.push_back(modeID);
                    	    quadEventData.quadID.push_back(quadID);

					nframe_bunch=0;
				}

			}
		    	else if(bunch_sync==7&&frame.bootpageno==0)
			{
				// Don't decode the event (it is corrupted)
				LOG(WARNING)<<"FOUND A BAD EVENT";
			}
                    
                }//END LOOP SINGLE PACKET
                
            }//END LOOP PACKETS
        } //END LOGICAL IF 0<DATAID<3 
        
        //Extracting veto spectrum
        if(vetospecstart==8) {
            vsd.vetoSpectrum.push_back(tempVetospectrum);
            
            vetoSpecRange = (bool)((frame.packets[0].febhdr.VetoSpecRange>>11)&(0x1));
            vsd.vetoSpecRange.push_back(vetoSpecRange);
            if(vetoSpecRange ==0){
                for(i=vetospecstart, j=24; j< VETOSPEC_SIZE; i++, j++ ){
                    vsd.vetoSpectrum[iveto][j] = frame.packets[0].hbtdata[i];
                }
            } else if(vetoSpecRange ==1) {
                for (i = vetospecstart, j = 0; j < VETOSPEC_SIZE; i++, j++) {
                    vsd.vetoSpectrum[iveto][j] = frame.packets[0].hbtdata[i];
                }                
            }
            
            iveto++; //index of veto spectrum  
        }
        vsd.cztseccnt.push_back((double)frame.secondcount);

        //Interpolating time for veto spectrum data
        instrumentTime = (double)frame.secondcount;
       
       	if (TCT.interpolate_time(instrumentTime, interpolatedTime)) {
            LOG(ERROR) << "Error in TCT time interpolation.";
            return EXIT_FAILURE;
        }
        vsd.time.push_back(interpolatedTime);

        //Added by Mithun (08/12/15)
        vsd.quadID.push_back(frame.packets[0].quadID);


        //Veto spectrum extracted
        
        //Extracting SSM data 
        if(frame.dataID >=8 && frame.dataID<=11){
            if(nPkts!=2){
                LOG(WARNING)<< "Number of Packets in frame "<< frame.frameno << " with dataID " << frame.dataID
                        << " is "<< nPkts;
                LOG(WARNING) << "** Expected number of packets: 2 **";
                continue;
            }
            ssmd.vetoSpec.push_back(tempVetospectrum);
            ssmd.cztSpec.push_back(tempCztspec);
            ssmd.cztSpecAlpha.push_back(tempCztspecAlpha);
            ssmd.cztSpecVeto.push_back(tempCztspecVeto);
            ssmd.cztseccnt.push_back((double)frame.secondcount);
            //Added by Mithun (08/12/15)
            ssmd.quadID.push_back(frame.packets[0].quadID);

            //Extracting veto spectrum for SSM data
            for(i=0, j=184; i<VETOSPEC_SIZE; j++, i++){
                ssmd.vetoSpec[iSSM][i] = frame.packets[0].hbtdata[j];
            }
            
            //Extracting CZT spectrum for SSM data
            for(i=0, j=440; i<CZTSPEC_SIZE; j++, i++){
                ssmd.cztSpec[iSSM][i] = frame.packets[0].hbtdata[j];
            }
            
            //Extracting CZT spectrum with veto for SSM data
            for (i = 0, j = 952; i < 56; j++, i++) {
                ssmd.cztSpecVeto[iSSM][i] = frame.packets[0].hbtdata[j];
            }
            for (i = 0, j = 0; i < 456; j++, i++) {
                ssmd.cztSpecVeto[iSSM][i] = frame.packets[1].hbtdata[j];
            } 
            //Extracting CZT spectrum with alpha for SSM data
            for (i = 0, j = 456; i < CZTSPECALPHA_SIZE; j++, i++) {                
                ssmd.cztSpecAlpha[iSSM][i] = frame.packets[1].hbtdata[j];
            }
            
            //Interpolating SSM Data time
            instrumentTime = (double)frame.secondcount;

	   
            if (TCT.interpolate_time(instrumentTime, interpolatedTime)) {
                LOG(ERROR) << "Error in TCT time interpolation.";
                return EXIT_FAILURE;
            }

            ssmd.time.push_back(interpolatedTime);
            iSSM++;
        }
        
        //Extracting Temperature data
        if(frame.dataID >=8 && frame.dataID<=11){
            tempd.cztseccnt.push_back((double)frame.secondcount);
            tempd.quadID.push_back(frame.packets[0].quadID);
            
            for(i=0; i<NUM_DET_PER_QUAD; i++){
                tempTemperature[i] = (float) frame.packets[0].febhdr.temperature[i];
            }
            tempd.temperature.push_back(tempTemperature);
            
            //Interpolating time for Temperature Data
            instrumentTime = (double)frame.secondcount;
		
            if (TCT.interpolate_time(instrumentTime, interpolatedTime)) {
                LOG(ERROR) << "Error in TCT time interpolation.";
                return EXIT_FAILURE;
            }
            tempd.time.push_back(interpolatedTime);
        }
        
        if(frame.dataID>=4 && frame.dataID<=11){
            counter++;
        }

        //Temperature data extracted

        
    } //END LOOP FRAME
   
   	LOG(INFO)<<"DONE READING FRAMES**********";
	LOG(INFO)<<"seccnt "<<quadEventData.cztseccnt.size();
    LOG(INFO)<<"ntick "<<quadEventData.cztntick.size();
	LOG(INFO)<<"tim "<<quadEventData.time.size();
    LOG(INFO)<<"pha "<<quadEventData.pha.size();
    LOG(INFO)<<"quadid "<<quadEventData.quadID.size();
				
	total_events_extracted=quadEventData.cztntick.size();

    LOG(INFO) << "Number of records outside TCT: " << TCT.get_nrecOutTCT();
    
    return EXIT_SUCCESS;
}

int L1evtData::get_total_valid_packets() {
    int status=0;
    long i=0, j=0;
    int nframes=0;
    int npackets=0;
    nframes = l1evtFrames.size();
    
    if(nframes >0){
        for(i=0; i<nframes; i++){
            npackets+=l1evtFrames[i].packets.size();
        }
    }
    totalValidPackets = npackets;
    totalFrames = l1evtFrames.size();
    return npackets;
}

int L1evtData::find_missing_frames() {
    int status=0;
    int nrows=l1evtFrames.size();
    int diffFrames;
    int i=0, j=0;
    missingFrames.clear();
    for(i=0; i<(nrows-1); i++){
        diffFrames=(l1evtFrames[i+1].frameno) - (l1evtFrames[i].frameno);
        if(diffFrames>1){
            for(j=1; j<diffFrames; j++){
                missingFrames.push_back(l1evtFrames[i].frameno + j);
            }
        }
    }
    totalMissingFrames = missingFrames.size();
    
    return EXIT_SUCCESS;
    
}

//**display
int L1evtData::display_l1data_info() {
    int status=0;
    int i=0,j=0;
    int iMode=0, iData=0;
    int isModeAvailable[16];
    long nframes[16][16];
    string display;
    //Initializing nframes
    for(i=0; i<16; i++){
        for(j=0; j<16; j++){
            nframes[i][j]=0;
        }
    }
    LOG(INFO) << "TOTAL PACKETS                     : " << totalPackets;
    LOG(INFO) << "TOTAL VALID PACKETS               : " << totalValidPackets;
    LOG(INFO) << "TOTAL INVALID PACKETS             : " << totalInvalidPackets;
    if (totalInvalidPackets > 0) {
        for (i = 0; i < totalInvalidPackets; i++) {
            cout << packetsRejected[i] << "\t";
        }
        cout << endl;
    }
    LOG(INFO) << "PACKETS WITH O PKT NUMBER         : " << countP0;
    LOG(INFO) << "TOTAL VALID FRAMES                : " << totalFrames;
//    LOG(INFO) << "TOTAL MISSING FRAMES              : " << totalMissingFrames;
    
    LOG(INFO) << "  Mode  " << " DataIDs->number of frames" << endl;
    for(i=0; i<l1evtFrames.size(); i++){
        for(iMode=0; iMode<16; iMode++){
            isModeAvailable[iMode]=FALSE;
            for(iData=0; iData<16; iData++){
                if(l1evtFrames[i].modeID==iMode){
                    if(l1evtFrames[i].dataID==iData){
                        isModeAvailable[iMode]=TRUE;
                        nframes[iMode][iData]++;
                    }
                }
            }
        }
    }


    for (iMode = 0; iMode < 16; iMode++) {
        if (isModeAvailable[iMode] == TRUE) {
            display = "   " +  itoa(iMode) + "   ";
            for (iData = 0; iData < 16; iData++) {
                if (nframes[iMode][iData] > 0) {
                    display += itoa(iData) + "->" + itoa(nframes[iMode][iData]) + "   ";
                }
            }
            LOG(INFO) << display;
        }
    }
   

    return EXIT_SUCCESS;
}

int L1evtData::display_l1extracted_data() {
    int status=0;
    long i=0;
    long nrows = quadEventData.time.size();
    cout << setw(10)<< "#" << setw(15)<< " CZTSECCNT" << setw(12)<< "CZTNTICK" << setw(4)<< "DETID" << setw(6)<< "PIXID" << setw(5)<< "PHA" << setw(7) << "DATAID" << endl;
    for(i=0; i<nrows; i++){
        cout << setw(10) << i << setw(15)<< setprecision(10) << quadEventData.cztseccnt[i] << setw(12)<< quadEventData.cztntick[i] << setw(4)<< 
                (int)quadEventData.detid[i] << setw(6)<< (int)quadEventData.pixid[i] << setw(5) << quadEventData.pha[i] << setw(7) << (int)quadEventData.quadID[i]<< endl;
    }
    
    return EXIT_SUCCESS;
}
//L1EVTDATA END

//SSMDATA CLASS

int SSMdata::reset() {
    time.clear();
    cztseccnt.clear();
    vetoSpec.clear();
    cztSpec.clear();
    cztSpecVeto.clear();
    cztSpecAlpha.clear();
    quadID.clear();

    return EXIT_SUCCESS;
}

//SSMDATA END

//HEADERDATA
int HeaderData::reset() {
    time.clear();
    cztseccnt.clear();
    dataID.clear();
    modeID.clear();
    wordCount.clear();
    frameNo.clear();
    statusFlag.clear();
    writePktNo.clear();
    readPktNo.clear();
    command.clear();
    laxpcTime.clear();
    dcnt.clear();
    errorCount.clear();
    errorFlag.clear();
    bootPageNo.clear();
    cztNo.clear();
    cztStatus.clear();
    evtReadMode.clear();
    cmdStatus.clear();
    bufferNo.clear();
    baseAdd.clear();
    vetoSpecRange.clear();
    channelNo.clear();
    ADCoutput.clear();
    cmd1sec.clear();
    alphaCount.clear();
    vetoCount.clear();
    cztCount_lt_uld.clear();
    cztCount_gt_uld.clear();
    cztdataRead.clear();
    
    return EXIT_SUCCESS;
}
//HEADER END

//VETOSPECTRUMDATA CLASS
int VetoSpectrumData::reset() {
    vetoSpecRange.clear();
    vetoSpectrum.clear();
    time.clear();
    cztseccnt.clear();
    quadID.clear();

    return EXIT_SUCCESS;
}
//VETOSPECTRUMDATA END

//TEMPERATUREDATA CLASS
int TemperatureData::reset() {
    time.clear();
    cztseccnt.clear();
    temperature.clear();
    quadID.clear();
    
    return EXIT_SUCCESS;
}
//TEMPERATUREDATA END

// QUADEVTDATA CLASS
int QuadEvtData::reset() {
    time.clear();
    cztseccnt.clear();
    cztntick.clear();
    pha.clear();
    detid.clear();
    pixid.clear();
    alpha.clear();
    veto.clear();
    detx.clear();
    dety.clear();
    dataID.clear();
    modeID.clear();
    quadID.clear(); 
   //clear bunches
    bunch_time.clear();
   evt_row_num.clear(); 
    time_dfs.clear();
    time_dsl.clear();
    num_bunchevents.clear();
    detid1.clear();
    detid2.clear();
    detid3.clear();
    detid4.clear();
    return EXIT_SUCCESS;
}
// QUADEVTDATA END

//PACKET CLASS
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
    
    //Frame validity checks
    if(pktNo==0) {
    flag=(frameno>=0 && frameno<65536) ? (flag & TRUE) : (flag & FALSE);
    if(flag==FALSE){
        return flag;
    }        
    }
    
    return flag;
}

Packet::Packet(){
    isValid=-1;
}

Packet::~Packet(){
}

//PACKET END

//FRAME CLASS

Frame::Frame() {
    isvalid=-1;
}

Frame::~Frame(){

}
int Frame::validity(){
    int flag=TRUE;
    flag=(frameno>=0 && frameno<65536) ? (flag & TRUE) : (flag & FALSE);
    if(flag==FALSE){
        return flag;
    }
    
    return flag;
}
//FRAME END
//Struct FEBheader
void FEBheader::read(unsigned short *data,int size,int mode,int dataid){
    //DLOG(INFO)<<"READING FEB HEADER WITH DATAID: "<< dataid << " MODE:" << mode;
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

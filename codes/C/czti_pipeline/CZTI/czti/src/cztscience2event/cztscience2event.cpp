#include "cztscience2event.h"

using namespace std;

cztscience2event::cztscience2event() {
    strcpy(modulename, "cztscience2event_v");
    strcat(modulename, VERSION);
}

int cztscience2event::read(int argc, char **argv) {
    int status = 0;


    if (PIL_OK != (status = PILInit(argc, argv))) {
        LOG(ERROR) << "***Error Initializing PIL***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("infile", infile))) {
        LOG(ERROR) << "***Error Reading Science data file:" << infile << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("TCTfile", TCTfile))) {
        LOG(ERROR) << "***Error Reading TCT file:" << TCTfile << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("outfile", outfile))) {
        LOG(ERROR) << "***Error Reading eventfile:" << outfile << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("hdrInfoFile", hdrInfoFile))) {
        LOG(ERROR) << "***Error Reading eventfile:" << hdrInfoFile << "***";
        return status;
    }

    
    if (PIL_OK != (status = PILGetFname("bunchfile", bunchfile))) {
        LOG(ERROR) << "***Error Reading eventfile:" << hdrInfoFile << "***";
        return status;
    }
//strcpy(bunchfile,"test.bunch");
    /*
    if (PIL_OK != (status = PILGetFname("GTIfile", GTIfile))) {
        LOG(ERROR) << "***Error Reading output GTI file:" << GTIfile << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("BTIfile", BTIfile))) {
        LOG(ERROR) << "***Error Reading output BTI file:" << BTIfile << "***";
        return status;
    }
    */
    if (PIL_OK != (status = PILGetInt("nPackets", &nPackets))) {
        LOG(ERROR) << "***Error Reading number of packets " << nPackets << "***";
        return status;
    }
    
    if (PIL_OK != (status = PILGetBool("BigEndian", &BigEndian))) {
        LOG(ERROR) << "***Error Reading BigEndian parameter" << BigEndian << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("clobber", &clobber))) {
        LOG(ERROR) << "***Error Reading clobber:" << clobber << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetBool("history", &history))) {
        LOG(ERROR) << "***Error Reading history parameter" << history << "***";
        return status;
    }
    
    if (PIL_OK != (status = PILGetBool("debug", &debug))) {
        LOG(ERROR) << "***Error Reading debug parameter" << debug << "***";
        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

int cztscience2event::read(char *infile, char *TCTfile, char *outfile, char* hdrInfoFile,
        char *gtifile, char *btifile, int nPackets, int bigendian, int clobber, int history, int debug) {

    if (strcpy(this->infile, infile) == NULL) {
        LOG(ERROR) << "***Error while copying input filename***";
        return (EXIT_FAILURE);
    }

    if (strcpy(this->TCTfile, TCTfile) == NULL) {
        LOG(ERROR) << "***Error while copying TCT filename***";
        return (EXIT_FAILURE);
    }

    if (strcpy(this->outfile, outfile) == NULL) {
        LOG(ERROR) << "***Error while copying output filename***";
        return (EXIT_FAILURE);
    }

    if (strcpy(this->hdrInfoFile, hdrInfoFile) == NULL) {
        LOG(ERROR) << "***Error while copying header information filename***";
        return (EXIT_FAILURE);
    }

    /*
    if (strcpy(this->GTIfile, gtifile) == NULL) {
        LOG(ERROR) << "***Error while copying output BTI filename***";
        return (EXIT_FAILURE);
    }

    if (strcpy(this->BTIfile, btifile) == NULL) {
        LOG(ERROR) << "***Error while copying output BTI filename***";
        return (EXIT_FAILURE);
    }
    */
    this->nPackets = nPackets;
    this->clobber = clobber;
    this->history = history;
    this->BigEndian = bigendian;
    this->debug = debug;
    return (EXIT_SUCCESS);

}

void cztscience2event::display() {
    LOG(INFO) << "----------------------------------------------------------------------------";
    LOG(INFO) << "                           SCIENCEDATA2EVENT PARAMETERS                     ";
    LOG(INFO) << "----------------------------------------------------------------------------";
    LOG(INFO) << "Modulename              : " << modulename;
    LOG(INFO) << "Science data file       : " << infile; //input science data file 
    LOG(INFO) << "TCT file                : " << TCTfile;
    LOG(INFO) << "Output Event file       : " << outfile; //output event file
    LOG(INFO) << "Header Inforamtion file : " << hdrInfoFile;
//    LOG(INFO) << "Output GTI file         : " << GTIfile;
//    LOG(INFO) << "Output BTI file         : " << BTIfile;
    LOG(INFO) << "No. of Packets in buffer: " << nPackets;
    LOG(INFO) << "BigEndian               : " << BigEndian;
    LOG(INFO) << "Clobber                 : " << clobber;
    LOG(INFO) << "History                 : " << history;
    LOG(INFO) << "Debug                   : " << debug;
    LOG(INFO) << "---------------------------------------------------------------------------";

}

/* Edited by Mithun NPS (30/03/16):
 * To correct the bu while reading the level-1 file in buffer mode.
 * Now last frame read in each set is ignored and it is read in the next 
 * set again to make sure that all packets of that frame are decoded.
 * Edited relevant part in level1decode.cpp as well.
 * */

/* Edit by Mithun N P S(19/05/16):
 * To keep track of total events extracted at any point of time, events_till_now
 * is defined. It is copied to l1evtdata variable so that row number of events are 
 * correctly known. The same variable in l1evtdata (total_events_extracted) is used
 * to get events extracted in a given buffer
 * */
int cztscience2event::cztscience2eventProcess() {
    int status=0;
    int i=0, j=0;
    int height=0;
    long quadnrows=0;
    long rowsToRead=nPackets;
    long startRowno=1, endRowno=rowsToRead;
	long events_till_now;
    int iquadFrame=0;
    int lastframepackets=0;
    L1evtHandler l1evt((string) infile);
    L1evtData q0data, q1data, q2data, q3data;
    EventFileHandler l2evt;

    //FILE EXISTENCE CHECK AND UNLINKING IF IT DOES
    if (clobber == YES) {
        if (FileExists(outfile)) {
            LOG(INFO) << outfile << "  :FileExists ";
            if (unlink(outfile) != 0) {
                LOG(ERROR) << "***Error in deleting " << outfile << "***";
                return (EXIT_FAILURE);
            }
        }
        if (FileExists(bunchfile)) {
            if (unlink(bunchfile) != 0) {
                LOG(ERROR) << "***Error in deleting " << bunchfile << "***";
                return (EXIT_FAILURE);
            }
        }
        if (FileExists(hdrInfoFile)) {
            if (unlink(hdrInfoFile) != 0) {
                LOG(ERROR) << "***Error in deleting " << hdrInfoFile << "***";
                return (EXIT_FAILURE);
            }
        }
    } else {
        if (FileExists(outfile) || FileExists(bunchfile)|| FileExists(hdrInfoFile)) {
            LOG(ERROR) << "***Output file already exists***";
            LOG(ERROR) << "Use clobber=yes for overwriting the file";
            return (EXIT_FAILURE);
        }
    }
    
    //Creating empty event file
    if(l2evt.create_l2event_file((string) outfile, EVTTEMPLATE, (string)infile)){
        LOG(ERROR) << "***Error creating empty level-2 event file.***";
        return EXIT_FAILURE;
    }
    //Creating empty level2 header file
    if(l2evt.create_l2hdr_file((string) hdrInfoFile, HDRTEMPLATE, (string) infile)){
        LOG(ERROR) << "***Error in creating empty level-2 header file " << hdrInfoFile << "***";
        return EXIT_FAILURE;
    }
    
    if(l2evt.create_l2event_file((string)bunchfile, BUNCHTEMPLATE, (string)infile)){
        LOG(ERROR) << "***Error creating empty level-2 event file.***";
        return EXIT_FAILURE;
    }
    //Reading level-1 Quad0 data
    LOG(INFO) << "--------------------------------------------------";
    LOG(INFO) << "                QUADRANT 0                        ";
    LOG(INFO) << "--------------------------------------------------";
    quadnrows = l1evt.get_l1evt_ext_nrows("CZT_QUAD1");
    startRowno = 1;
    endRowno = rowsToRead;
	events_till_now=0;
    do {
        LOG(INFO) << "START ROW: " << startRowno << " END ROW: " << endRowno;
        if (l1evt.read_l1evt_file_pkts("CZT_QUAD1", startRowno, endRowno)) {
            return EXIT_FAILURE;
        }
        height = l1evt.get_pktArrayRows(); //Getting number of packets extracted from Quad 0 of level-1 data

	q0data.total_events_extracted=events_till_now;

	if(height==rowsToRead)
		lastframepackets=0;
	else if(height<rowsToRead)
		lastframepackets=-1;
        //Set Packet Array
        if (q0data.set_pktArray(l1evt.pktArray, height)) {
            LOG(ERROR) << "***Error in setting Quadrant 0 packet array into q0data (which handles packet segregation and extraction).***";
            return EXIT_FAILURE;
        }
        //Segregating level-1 data for quadrant 0
        LOG(INFO) << "Segregating data for Quad 0 of level-1 data which consists of " << height << " packets.";
        if (q0data.segregate_data(lastframepackets,debug)) {
            LOG(ERROR) << "Error in segregating data for quad 0 of Level-1 event file " << infile;
            return EXIT_FAILURE;
        }
        //Displaying level-1 data information for quadrant 0
        if (q0data.display_l1data_info()) {
            LOG(ERROR) << "Error in displaying level-1 data info for quadrant 0.";
            return EXIT_FAILURE;
        }
        //Extracting level-1 data for quadrant 0
        LOG(INFO) << "Extracting level1 data for quadrant 0.";
        if (q0data.extract_data(TCTfile)) {
            LOG(ERROR) << "Error in extracting level-1 data for quadrant 0.";
            return EXIT_FAILURE;
        }
        LOG(INFO) << "Number of events in these packets: " << q0data.quadEventData.time.size();
        LOG(INFO) << "Number of header records in these packets: " << q0data.hdrd.time.size();
        //Correcting temperature values
        LOG(INFO) << "Correcting junk temperature values.";
        try {
            q0data.tempd.remove_junk_temperature();
        } catch (ErrorHandler errHandler) {
            logError(errHandler);
            if (errHandler.severity == errERROR) {
                return EXIT_FAILURE;
            }
        }

        //Writing extracted data to level-2 file for quadrant 0;
        LOG(INFO) << "Writing extracted quadrant 0 data in level-2 event file " << outfile << ".";
        if (l2evt.write_l2event_file((string) outfile,(string)bunchfile, q0data)) {
            LOG(ERROR) << "Error in writing extracted quadrant 0 data in level-2 event file " << outfile << ".";
        }

        //Writing extracted header data to level-2 header file
        LOG(INFO) << "Writing header data...";
        if (l2evt.write_l2_hdr_file((string) hdrInfoFile, q0data)) {
            LOG(ERROR) << "Error in writing header data.";
            return EXIT_FAILURE;
        }

		events_till_now+=q0data.total_events_extracted;

		LOG(INFO)<<"Total events till now for this quadrant is: "<< events_till_now;
        //Clearing q0data & l1evt objects
        q0data.reset();
        l1evt.clear_pktArray();

/*      startRowno += rowsToRead;
        endRowno += rowsToRead;*/

        startRowno = endRowno+1-lastframepackets;
        endRowno = startRowno +rowsToRead-1;
    } while (startRowno < quadnrows);
    
    //Reading level-1 Quad1 data
    LOG(INFO) << "--------------------------------------------------";
    LOG(INFO) << "                QUADRANT 1                        ";
    LOG(INFO) << "--------------------------------------------------";
    
    quadnrows = l1evt.get_l1evt_ext_nrows("CZT_QUAD2");
    startRowno = 1;
    endRowno = rowsToRead;
    events_till_now=0;	
    do {
        
        LOG(INFO) << "START ROW: " << startRowno << " END ROW: " << endRowno;
        if (l1evt.read_l1evt_file_pkts("CZT_QUAD2", startRowno, endRowno)) {
            return EXIT_FAILURE;
        }
        height = l1evt.get_pktArrayRows(); //Getting number of packets extracted from Quad 1 of level-1 data
    
		q1data.total_events_extracted=events_till_now;


        if(height==rowsToRead)
		lastframepackets=0;
	else if(height<rowsToRead)
		lastframepackets=-1;

	//Set Packet Array
        if (q1data.set_pktArray(l1evt.pktArray, height)) {
            LOG(ERROR) << "***Error in setting Quadrant 1 packet array into q1data (which handles packet segregation and extraction).***";
            return EXIT_FAILURE;
        }
        //Segregating level-1 data for quadrant 1
        LOG(INFO) << "Segregating data for Quad 1 of level-1 data which consists of " << height << " packets.";
	if (q1data.segregate_data(lastframepackets,debug)) {            
	    LOG(ERROR) << "Error in segregating data for quad 1 of Level-1 event file " << infile;
            return EXIT_FAILURE;
        }
        //Displaying level-1 data information for quadrant 1
        if (q1data.display_l1data_info()) {
            LOG(ERROR) << "Error in displaying level-1 data info for quadrant 1.";
            return EXIT_FAILURE;
        }
        //Extracting level-1 data for quadrant 1
        LOG(INFO) << "Extracting level1 data for quadrant 1.";
        if (q1data.extract_data(TCTfile)) {
            LOG(ERROR) << "Error in extracting level-1 data for quadrant 1.";
            return EXIT_FAILURE;
        }
        
        LOG(INFO) << "Number of events in these packets: " << q1data.quadEventData.time.size();
        LOG(INFO) << "Number of header records in these packets: " << q1data.hdrd.time.size();
        //Writing extracted data to level-2 file for quadrant 1;
        LOG(INFO) << "Writing extracted quadrant 1 data in level-2 event file " << outfile << ".";
        if (l2evt.write_l2event_file((string) outfile ,(string)bunchfile , q1data)) {
            LOG(ERROR) << "Error in writing extracted quadrant 1 data in level-2 event file " << outfile << ".";
        }
        //Writing extracted header data to level-2 header file
        LOG(INFO) << "Writing header data...";
        if (l2evt.write_l2_hdr_file((string) hdrInfoFile, q1data)) {
            LOG(ERROR) << "Error in writing header data.";
            return EXIT_FAILURE;
        }
        

        events_till_now+=q1data.total_events_extracted;

        LOG(INFO)<<"Total events till now for this quadrant is: "<< events_till_now;

		//Clearing q1data & l1evt objects
        q1data.reset();
        l1evt.clear_pktArray();

    /*    startRowno += rowsToRead;
        endRowno += rowsToRead;*/
        startRowno = endRowno+1-lastframepackets;
        endRowno = startRowno +rowsToRead-1;

    } while (startRowno < quadnrows);


    //Reading level-1 Quad2 data
    LOG(INFO) << "--------------------------------------------------";
    LOG(INFO) << "                QUADRANT 2                        ";
    LOG(INFO) << "--------------------------------------------------";

    quadnrows = l1evt.get_l1evt_ext_nrows("CZT_QUAD3");
    startRowno = 1;
    endRowno = rowsToRead;
	events_till_now=0;

    do {
        LOG(INFO) << "START ROW: " << startRowno << " END ROW: " << endRowno;
        if (l1evt.read_l1evt_file_pkts("CZT_QUAD3", startRowno, endRowno)) {
            return EXIT_FAILURE;
        }
        height = l1evt.get_pktArrayRows(); //Getting number of packets extracted from Quad 2 of level-1 data

   	 	q2data.total_events_extracted=events_till_now;
		
        if(height==rowsToRead)
		lastframepackets=0;
	else if(height<rowsToRead)
		lastframepackets=-1;
		
	//Set Packet Array
        if (q2data.set_pktArray(l1evt.pktArray, height)) {
            LOG(ERROR) << "***Error in setting Quadrant 2 packet array into q2data (which handles packet segregation and extraction).***";
            return EXIT_FAILURE;
        }
        //Segregating level-1 data for quadrant 2
        LOG(INFO) << "Segregating data for Quad 2 of level-1 data which consists of " << height << " packets.";
        if (q2data.segregate_data(lastframepackets,debug)) {
            LOG(ERROR) << "Error in segregating data for quad 2 of Level-1 event file " << infile;
            return EXIT_FAILURE;
        }
        //Displaying level-1 data information for quadrant 2
        if (q2data.display_l1data_info()) {
            LOG(ERROR) << "Error in displaying level-1 data info for quadrant 2.";
            return EXIT_FAILURE;
        }
        //Extracting level-1 data for quadrant 2
        LOG(INFO) << "Extracting level1 data for quadrant 2.";
        if (q2data.extract_data(TCTfile)) {
            LOG(ERROR) << "Error in extracting level-1 data for quadrant 2.";
            return EXIT_FAILURE;
        }

        LOG(INFO) << "Number of events in these packets: " << q2data.quadEventData.time.size();
        LOG(INFO) << "Number of header records in these packets: " << q2data.hdrd.time.size();
        //Writing extracted data to level-2 file for quadrant 2;
        LOG(INFO) << "Writing extracted quadrant 2 data in level-2 event file " << outfile << ".";
        if (l2evt.write_l2event_file((string) outfile,(string)bunchfile, q2data)) {
            LOG(ERROR) << "Error in writing extracted quadrant 2 data in level-2 event file " << outfile << ".";
        }
        //Writing extracted header data to level-2 header file
        LOG(INFO) << "Writing header data...";
        if (l2evt.write_l2_hdr_file((string) hdrInfoFile, q2data)) {
            LOG(ERROR) << "Error in writing header data.";
            return EXIT_FAILURE;
        }
        
        events_till_now+=q2data.total_events_extracted;

        LOG(INFO)<<"Total events till now for this quadrant is: "<< events_till_now;

		//Clearing q2data & l1evt objects
        q2data.reset();
        l1evt.clear_pktArray();

/*	startRowno += rowsToRead;
        endRowno += rowsToRead;*/
        startRowno = endRowno+1-lastframepackets;
        endRowno = startRowno +rowsToRead-1;

    } while (startRowno < quadnrows);

    //Reading level-1 Quad3 data
    LOG(INFO) << "--------------------------------------------------";
    LOG(INFO) << "                QUADRANT 3                        ";
    LOG(INFO) << "--------------------------------------------------";
    quadnrows = l1evt.get_l1evt_ext_nrows("CZT_QUAD4");
    startRowno=1;
    endRowno=rowsToRead;
    events_till_now=0;

    do {
        LOG(INFO) << "START ROW: " << startRowno << " END ROW: " << endRowno;
        if (l1evt.read_l1evt_file_pkts("CZT_QUAD4", startRowno, endRowno)) {
            return EXIT_FAILURE;
        }
        height = l1evt.get_pktArrayRows(); //Getting number of packets extracted from Quad 3 of level-1 data

    	q3data.total_events_extracted=events_till_now;
        
        if(height==rowsToRead)
	                lastframepackets=0;
        else if(height<rowsToRead)
	                lastframepackets=-1;
	
	//Set Packet Array
        if (q3data.set_pktArray(l1evt.pktArray, height)) {
            LOG(ERROR) << "***Error in setting Quadrant 3 packet array into q3data (which handles packet segregation and extraction).***";
            return EXIT_FAILURE;
        }
        //Segregating level-1 data for quadrant 3
        LOG(INFO) << "Segregating data for Quad 3 of level-1 data which consists of " << height << " packets.";
        if (q3data.segregate_data(lastframepackets,debug)) {
            LOG(ERROR) << "Error in segregating data for quad 3 of Level-1 event file " << infile;
            return EXIT_FAILURE;
        }
        //Displaying level-1 data information for quadrant 3
        if (q3data.display_l1data_info()) {
            LOG(ERROR) << "Error in displaying level-1 data info for quadrant 3.";
            return EXIT_FAILURE;
        }
        //Extracting level-1 data for quadrant 3
        LOG(INFO) << "Extracting level1 data for quadrant 3.";
        if (q3data.extract_data(TCTfile)) {
            LOG(ERROR) << "Error in extracting level-1 data for quadrant 3.";
            return EXIT_FAILURE;
        }
        LOG(INFO) << "Number of events in these packets: " << q3data.quadEventData.time.size();
        LOG(INFO) << "Number of header records in these packets: " << q3data.hdrd.time.size();
        //Writing extracted data to level-2 file for quadrant 3;
        LOG(INFO) << "Writing extracted quadrant 3 data in level-2 event file " << outfile << ".";
        if (l2evt.write_l2event_file((string) outfile,(string)bunchfile, q3data)) {
            LOG(ERROR) << "Error in writing extracted quadrant 3 data in level-2 event file " << outfile << ".";
        }
        //Writing extracted header data to level-2 header file
        LOG(INFO) << "Writing header data...";
        if (l2evt.write_l2_hdr_file((string) hdrInfoFile, q3data)) {
            LOG(ERROR) << "Error in writing header data.";
            return EXIT_FAILURE;
        }
        
        events_till_now+=q3data.total_events_extracted;

        LOG(INFO)<<"Total events till now for this quadrant is: "<< events_till_now;

		//Clearing q3data & l1evt objects
        q3data.reset();
        l1evt.clear_pktArray();

/*	startRowno += rowsToRead;
        endRowno += rowsToRead;*/
        startRowno = endRowno+1-lastframepackets;
        endRowno = startRowno +rowsToRead-1;

    } while (startRowno < quadnrows);


    LOG(INFO) << "Copying GTI extensions ";
    if(l2evt.copy_gti_extensions(outfile, infile)){    
    LOG(ERROR) << "Error in copying GTI extensions ";
    return EXIT_FAILURE;
    }

/*    
    //Reading GTI extension & writing GTI file
    if(l1evt.read_l1evt_file_gti()){
        LOG(ERROR) << "Error in reading GTI extension of level-1 event file " << infile;
        return (EXIT_FAILURE);
    }
    if(l2evt.create_l2_gti_file((string)GTIfile, GTITEMPLATE, (string)infile)){
        LOG(ERROR) << "Error in creating GTI file " << GTIfile;
        return EXIT_FAILURE;
    }
    if(l2evt.write_l2_gti_file((string)GTIfile, l1evt.gti.tstart, l1evt.gti.tstop)){
        LOG(ERROR) << "Error in writing time data in GTI file " << GTIfile;
        return EXIT_FAILURE;
    }
    //Reading BTI extension & writing BTI file
    if(l1evt.read_l1evt_file_bti()){
        LOG(ERROR) << "Error in reading BTI extension of level-1 event file " << infile;
        return (EXIT_FAILURE);
    }
    if(l2evt.create_l2_bti_file((string)BTIfile, BTITEMPLATE, (string)infile)){
        LOG(ERROR) << "Error in creating BTI file " << BTIfile;
        return EXIT_FAILURE;
    }
    if(l2evt.write_l2_bti_file((string)BTIfile, l1evt.bti.tstart, l1evt.bti.tstop)){
        LOG(ERROR) << "Error in writing time data in BTI file " << BTIfile;
        return EXIT_FAILURE;
    }
*/

    //writing history to all HDUs of output file
    if (history == YES) {
        vector<string> vhistory;
	vhistory.clear();
        getHistory(vhistory);
        if(writeHistory(outfile, vhistory)){ //writes history to each HDU of output filE
            LOG(ERROR) << "Error in writing History to event file.";
        } 
/*        if(writeHistory(GTIfile, vhistory)){ //writes history to each HDU of output filE
            LOG(ERROR) << "Error in writing History to GTI file.";
        } 
        if(writeHistory(BTIfile, vhistory)){ //writes history to each HDU of output filE
            LOG(ERROR) << "Error in writing History to BTI file.";
        } 
*/
        if(writeHistory(hdrInfoFile, vhistory)){ //writes history to each HDU of output filE
            LOG(ERROR) << "Error in writing History to Header file.";
        } 
    }

    /*
    //UPDATING KEYWORDS
    if(updateHdrTime(outfile, "Q0", "TIME")){
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q0";
        return EXIT_FAILURE;
    }
    if(updateHdrTime(outfile, "Q1", "TIME")){
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q1";
        return EXIT_FAILURE;
    }
    if(updateHdrTime(outfile, "Q2", "TIME")){
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q2";
        return EXIT_FAILURE;
    }
    if(updateHdrTime(outfile, "Q3", "TIME")){
        LOG(ERROR) << "Error in updating TSTART, TSTOP keywords for Q3";
        return EXIT_FAILURE;
    }
*/

	LOG(INFO)<<"Updating checksum keywords";

    if(updateKeywords(outfile, modulename)){
        LOG(WARNING) << "Error in updating keywords.";
    } 


/*
    if(updateKeywords(GTIfile, modulename)){
        LOG(WARNING) << "Error in updating keywords.";
    }    
    if(updateKeywords(BTIfile, modulename)){
        LOG(WARNING) << "Error in updating keywords.";
    }    
*/
    if(updateKeywords(hdrInfoFile, modulename)){
        LOG(WARNING) << "Error in updating keywords.";
    }    
    
    return EXIT_SUCCESS;
}

int cztscience2event::getHistory(vector<string> &vhistory) {
    //char *user = getlogin();
    strcpy(modulename, "cztscience2event_v");
    strcat(modulename, VERSION);
    char *user = getenv("USER");
	vhistory.push_back("Module run by " + (string) user);
    vhistory.push_back("Parameter List START for " + (string) modulename);
    vhistory.push_back("P1 infile=" + (string) infile);
    vhistory.push_back("P2 TCTfile=" + (string) TCTfile);
    vhistory.push_back("P3 outfile=" + (string) outfile);
    vhistory.push_back("P4 hdrInfoFile=" + (string) hdrInfoFile);
    //vhistory.push_back("P5 GTIfile=" + (string) GTIfile);
    //vhistory.push_back("P6 BTIfile=" + (string) BTIfile);

    if (BigEndian == YES)
        vhistory.push_back("P7 BigEndian=yes");
    else
        vhistory.push_back("P7 BigEndian=no");

    if (clobber == YES)
        vhistory.push_back("P8 clobber=yes");
    else
        vhistory.push_back("P8 clobber=no");
    if (history == YES)
        vhistory.push_back("P9 history=yes");
    else
        vhistory.push_back("P9 history=no");
    if (debug == YES)
        vhistory.push_back("P10 debug=yes");
    else
        vhistory.push_back("P10 debug=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

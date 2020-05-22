#include "level2validation.h"
#include "ExpMap.h"
#include "maskGeometry.h"
#include "fft.h"

using namespace std;

int FitsFileInfo::get_fitsinfo(fitsfile* fptr){
    int status=0; //status variable
    int i,j=0;
    long nrows=0;
    int ncols=0;
    string errorMsg="";
    char extname[MAX_KEYWORD_SIZE];
    int hdutype=0;
    FitsInfo f;
    
    fits_get_num_hdus(fptr, &numHDU, &status);
    errorMsg = "Error in getting number of HDUs in input fits file.";
    report_error(status, errorMsg);
    
    f.hduName = "Primary";
    f.hduType = IMAGE_HDU;
    f.rows=0;
    f.cols=0;
    
    fitsFileMap[f.hduName] = f;
    
    if (numHDU>1){
        for(i=2; i<=numHDU; i++){
            fits_movabs_hdu(fptr, i, &hdutype, &status);
            errorMsg = "Error getting HDU type";
            if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
            f.hduType = hdutype;

            if(hdutype==BINARY_TBL or hdutype==ASCII_TBL) {
                fits_get_num_rows(fptr, &nrows, &status);
                errorMsg= "Error in getting number of rows.";
                if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
                f.rows = nrows;

                fits_get_num_cols(fptr, &ncols, &status);
                errorMsg= "Error in getting number of cols.";
                if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
                f.cols = ncols;
            }
            else {
                f.rows=0;
                f.cols=0;
            }
            fits_read_key(fptr, TSTRING, "EXTNAME", extname, NULL, &status);
            errorMsg= "Error reading keyname from file header";
            if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
            f.hduName = (string) extname;

            fitsFileMap.insert(pair<string, FitsInfo> (f.hduName, f));
        }
    }
    
    return status;
    
}

int FitsFileInfo::validate_fitsfile(int numHDUs, vector<string> hduNames, vector<int> hduTypes, int* status, vector<int> ncols, vector<long> nrows) {
    int i, j = 0;
    DLOG(INFO) << "Validating Input fits file:";
    
    if (hduNames.size()!=numHDUs || hduTypes.size()!=numHDUs) {
        LOG(ERROR) << "Vectors declared are not of proper size";
        LOG(ERROR) << "Claimed size: " << numHDUs << " | Size of Vector found: " << hduNames.size();
        *status = EXIT_FAILURE;
        return (EXIT_FAILURE);
    }
    if (this->numHDU != numHDUs) {
        LOG(ERROR) << "Number of HDUs do not match.";
        LOG(ERROR) << "HDUs expected: " << numHDUs << " | HDUs found: " << this->numHDU;
        *status = EXIT_FAILURE;
        return (EXIT_FAILURE);
    } else {
        for (i=0; i<numHDUs; i++) {
            if(this->fitsFileMap.find(hduNames[i]) == this->fitsFileMap.end()){
                LOG(ERROR) << "No HDU by "<< hduNames[i] << " in input file.";
                *status = EXIT_FAILURE;
                return (EXIT_FAILURE);
            }
            if (this->fitsFileMap[hduNames[i]].hduType != hduTypes[i]) {
                LOG(ERROR) << "HDU Type did not match for extension " << hduNames[i];
                LOG(ERROR) << "Expected: " << hduTypes[i] << "| Found: " << this->fitsFileMap[hduNames[i]].hduType;
                *status = EXIT_FAILURE;
                return (EXIT_FAILURE);
            }
            if (!ncols.empty()) {
                if (this->fitsFileMap[hduNames[i]].cols != ncols[i]) {
                    LOG(ERROR) << "Number of columns differ in extension: " << hduNames[i];
                    LOG(ERROR) << "Expected: " << ncols[i] << "| Found: " << this->fitsFileMap[hduNames[i]].cols;
                    *status = EXIT_FAILURE;
                    return (EXIT_FAILURE);
                }
            }
            if (!nrows.empty()) {
                if (this->fitsFileMap[hduNames[i]].rows != nrows[i]) {
                    LOG(ERROR) << "Number of rows differ in extension: " << hduNames[i];
                    LOG(ERROR) << "Expected: " << nrows[i] << "| Found: " << this->fitsFileMap[hduNames[i]].rows;
                    *status = EXIT_FAILURE;
                    return (EXIT_FAILURE);
                }
            }
        }
    }
    
    DLOG(INFO) << "Fits file is valid";
    *status=EXIT_SUCCESS;
    return (EXIT_SUCCESS);

}
   
int FitsFileInfo::get_hdu_info(string key, int &hdutype, int &ncols, long &nrows){
    if (fitsFileMap.find(key) != fitsFileMap.end()) {
        hdutype = fitsFileMap[key].hduType;
        ncols = fitsFileMap[key].cols;
        nrows = fitsFileMap[key].rows;
    }
    else {
        return (EXIT_FAILURE);
    }
    
    return EXIT_SUCCESS;
}
int FitsFileInfo::display_fitsinfo() {
    int i,j=0; //status variables
    LOG(INFO) << ">>> FITS FILE INFO <<<";
    LOG(INFO) << "Number of HDUs: " << numHDU;
    LOG(INFO) << "-------------------------------------------------------------";
    if(numHDU>1) {
        for (map<string,FitsInfo>::const_iterator p = fitsFileMap.begin(); p != fitsFileMap.end(); p++) {
            LOG(INFO) << p->first << "    ExtensionType:" << p->second.hduType << "   rows:" << p->second.rows << "   cols:" << p->second.cols;
        }
    }
    LOG(INFO) << "------------------------------------------------------------";
    LOG(INFO) << "IMAGE_HDU:0  |  ASCII_TBL:1  |  BINARY_TBL:2";
}

//CLASS EVENTFILEHANDLER

int EventFileHandler::create_l2event_file(string l2eventFilename, string eventTemplate, string l1eventFilename) {
    int status=0;
    fitsfile *fl2evt, *fl1evt;
    int i=0,j=0;
    int iquad=0;
    int hdutype=0;
    int hdunum=0;

    cztHeaderParam headKey;

    vector <string> vecKeywords;
    if(create_empty_fitsfile(l2eventFilename, eventTemplate)){
        LOG(ERROR) << "Error in creating empty level-2 event file " << l2eventFilename;
        return EXIT_FAILURE;
    }
    
    //Opening level-1 event file
    fits_open_file(&fl1evt, (char*) l1eventFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-1 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Opening level-2 event file
    fits_open_file(&fl2evt, (char*) l2eventFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-2 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Keywords to be copied
    vecKeywords.push_back("DATE-OBS");
    vecKeywords.push_back("TIME-OBS");
    vecKeywords.push_back("DATE-END");
    vecKeywords.push_back("TIME-END");
    vecKeywords.push_back("TIMESYS");
    vecKeywords.push_back("MJDREFI");
    vecKeywords.push_back("MJDREFF");
    vecKeywords.push_back("TIMEZERI");
    vecKeywords.push_back("TIMEZERF");
    vecKeywords.push_back("TIMEDEL");
    vecKeywords.push_back("TIMEUNIT");
    vecKeywords.push_back("TSTARTI");
    vecKeywords.push_back("TSTARTF");
    vecKeywords.push_back("TSTOPI");
    vecKeywords.push_back("TSTOPF");
    vecKeywords.push_back("OBJECT");
    vecKeywords.push_back("RA_PNT");
    vecKeywords.push_back("DEC_PNT");
    vecKeywords.push_back("EQUINOX");
    vecKeywords.push_back("RADECSYS");
    vecKeywords.push_back("OBSERVER");
    vecKeywords.push_back("OBS_ID");
    vecKeywords.push_back("SOURCEID");
    vecKeywords.push_back("TARGETID");
    vecKeywords.push_back("EXPOSURE");
    vecKeywords.push_back("ORB_NUM");
    
   
    if (copyUserKeyWords(fl1evt, fl2evt, "Primary", "Primary", vecKeywords)) { //copying keywords from input file
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;
    }

    //Set Deafults and Read the header Keys from level1 file
    headKey.readFromHeader(fl1evt);
    headKey.set_default_values();


    //Write to primary of level2 evt file (UPDATE)
    headKey.writeToHeader(fl2evt);

    //To copy exposure time for each quadrant extension
    vecKeywords.clear();
    vecKeywords.push_back("EXPOSURE");

    if (copyUserKeyWords(fl1evt, fl2evt, "CZT_QUAD1","Q0", vecKeywords)) { //copying keywords from input file
            LOG(ERROR) << "***Error in copying keywords***";
                return EXIT_FAILURE;
      }

    if (copyUserKeyWords(fl1evt, fl2evt, "CZT_QUAD2", "Q1", vecKeywords)) { //copying keywords from input file
            LOG(ERROR) << "***Error in copying keywords***";
                return EXIT_FAILURE;
      }

    if (copyUserKeyWords(fl1evt, fl2evt, "CZT_QUAD3", "Q2", vecKeywords)) { //copying keywords from input file
            LOG(ERROR) << "***Error in copying keywords***";
                return EXIT_FAILURE;
      }

    if (copyUserKeyWords(fl1evt, fl2evt, "CZT_QUAD4", "Q3", vecKeywords)) { //copying keywords from input file
            LOG(ERROR) << "***Error in copying keywords***";
                return EXIT_FAILURE;
      }


    // Writing other header keywords to quadrant extensions
    for(i=2;i<=5;i++)
    {       

	fits_movabs_hdu(fl1evt, i, &hdutype, &status);
	
        fits_movabs_hdu(fl2evt, i, &hdutype, &status);
        if (status) {
        LOG(ERROR) << "Error in moving to extension "<<i<<"of event file " ;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
        }
	headKey.readFromHeader(fl1evt);
        headKey.writeToHeader(fl2evt);
    }

    fits_close_file(fl1evt, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-1 event file " << l1eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fl2evt, &status);
    if (status) {
        LOG(ERROR) << "Error in closing DPI/DPH file " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


    return (EXIT_SUCCESS);    
    
    
}

int EventFileHandler::write_l2event_file(string l2eventFilename,string bunchFile, L1evtData l1data) {
    int status=0;
    fitsfile *fl2evt,*bunchptr;
    long i=0,j=0;
    long nrows=0;
    long nrows_vec=0; //number of rows in quadevent data
    double tstart=0.0, tstop=0.0;
    vector <double> quadTime;
    int colnum=0;
    int iquad=0;
    string hduname;
    vector <QuadEvtData> quadEvtData;
    QuadEvtData tempqedata;
    
    //Initializing quadEvtData
    for(iquad=0; iquad<4; iquad++){
        quadEvtData.push_back(tempqedata);
    }
    
    fits_open_file(&fl2evt, (char*) l2eventFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-2 event file: " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

 /*   
	iquad=l1data.quadEventData.quadID[0];
    
	int notsamequad=0;

	for(i=1; i<l1data.quadEventData.get_nrows(); i++){
		if(l1data.quadEventData.quadID[i] != iquad)
		{
			notsamequad=1;
			break;
		}
	}

	if(notsamequad){
*/
    //Segregatin l1data quadrant wise
    for(i=0; i<l1data.quadEventData.get_nrows(); i++){
      // DLOG(INFO) << "DEBUG:: " << (int) l1data.quadEventData.quadID[i];
        
        switch(l1data.quadEventData.quadID[i]){
            case 0:
                iquad=0;
                break;
            case 1:
                iquad=1;
                break;
            case 2:
                iquad=2;
                break;
            case 3:
                iquad=3;
                break;
            default:
                LOG(ERROR) << "Value is neither 0 or 1 or 2 or 3";
                return (EXIT_FAILURE);
        }

        quadEvtData[iquad].time.push_back(l1data.quadEventData.time[i]);
        quadEvtData[iquad].cztseccnt.push_back(l1data.quadEventData.cztseccnt[i]);
        quadEvtData[iquad].cztntick.push_back(l1data.quadEventData.cztntick[i]);
        quadEvtData[iquad].pha.push_back(l1data.quadEventData.pha[i]);
        quadEvtData[iquad].detid.push_back(l1data.quadEventData.detid[i]);
        quadEvtData[iquad].pixid.push_back(l1data.quadEventData.pixid[i]);
        quadEvtData[iquad].alpha.push_back(l1data.quadEventData.alpha[i]);
        quadEvtData[iquad].veto.push_back(l1data.quadEventData.veto[i]);
        quadEvtData[iquad].detx.push_back(l1data.quadEventData.detx[i]);
        quadEvtData[iquad].dety.push_back(l1data.quadEventData.dety[i]);
        quadEvtData[iquad].dataID.push_back(l1data.quadEventData.dataID[i]);
        quadEvtData[iquad].modeID.push_back(l1data.quadEventData.modeID[i]);
        quadEvtData[iquad].quadID.push_back(l1data.quadEventData.quadID[i]);
     
        if(i<l1data.quadEventData.bunch_time.size())
	{
	//adding information about bunches
        quadEvtData[iquad].bunch_time.push_back(l1data.quadEventData.bunch_time[i]);
        quadEvtData[iquad].evt_row_num.push_back(l1data.quadEventData.evt_row_num[i]);
        quadEvtData[iquad].time_dfs.push_back(l1data.quadEventData.time_dfs[i]);
        quadEvtData[iquad].time_dsl.push_back(l1data.quadEventData.time_dsl[i]);
        quadEvtData[iquad].num_bunchevents.push_back(l1data.quadEventData.num_bunchevents[i]);
        quadEvtData[iquad].detid1.push_back(l1data.quadEventData.detid1[i]);
        quadEvtData[iquad].detid2.push_back(l1data.quadEventData.detid2[i]);
        quadEvtData[iquad].detid3.push_back(l1data.quadEventData.detid3[i]);
        quadEvtData[iquad].detid4.push_back(l1data.quadEventData.detid4[i]);
        quadEvtData[iquad].DetId_fevt.push_back(l1data.quadEventData.DetId_fevt[i]);
        quadEvtData[iquad].DetId_sevt.push_back(l1data.quadEventData.DetId_sevt[i]);
        quadEvtData[iquad].DetId_tevt.push_back(l1data.quadEventData.DetId_tevt[i]);

        quadEvtData[iquad].PixId_fevt.push_back(l1data.quadEventData.PixId_fevt[i]);
        quadEvtData[iquad].PixId_sevt.push_back(l1data.quadEventData.PixId_sevt[i]);
        quadEvtData[iquad].PixId_tevt.push_back(l1data.quadEventData.PixId_tevt[i]);
//	LOG(INFO)<<l1data.quadEventData.bunch_time[i]<<" "<<l1data.quadEventData.evt_row_num[i]<<" "<<(int)l1data.quadEventData.time_dfs[i]<<" "<<(int)l1data.quadEventData.time_dsl[i]<<" "<<(int)l1data.quadEventData.num_bunchevents[i]<<" "<<(int)l1data.quadEventData.detid1[i];

	}
     }

/*}
	else
	{
		quadEvtData[iquad].time=l1data.quadEventData.time;
		quadEvtData[iquad].cztseccnt=l1data.quadEventData.cztseccnt;
		quadEvtData[iquad].cztntick=l1data.quadEventData.cztntick;
		quadEvtData[iquad].pha=l1data.quadEventData.pha;
		quadEvtData[iquad].detid=l1data.quadEventData.detid;
        quadEvtData[iquad].pixid=(l1data.quadEventData.pixid);
        quadEvtData[iquad].alpha=(l1data.quadEventData.alpha);
        quadEvtData[iquad].veto=(l1data.quadEventData.veto);
        quadEvtData[iquad].detx=(l1data.quadEventData.detx);
        quadEvtData[iquad].dety=(l1data.quadEventData.dety);
        quadEvtData[iquad].dataID=(l1data.quadEventData.dataID);
        quadEvtData[iquad].modeID=(l1data.quadEventData.modeID);
        quadEvtData[iquad].quadID=(l1data.quadEventData.quadID);
	}
  */

    //Writing level-2 event data for different quadrants
    for(iquad=0; iquad<4; iquad++){
        //Sorting object quadEvtData
        nrows_vec=quadEvtData[iquad].get_nrows();
        hduname = quadExtName(iquad);
      
       LOG(INFO)<<"@@@@@@@@@ Writing "<<iquad<<" rows "<<nrows_vec;

	    if(nrows_vec!=0){
	    
        fits_movnam_hdu(fl2evt, BINARY_TBL, (char*) (hduname).c_str(), NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << hduname << " of L2 event file " <<
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_get_num_rows(fl2evt, &nrows, &status);
        if (status) {
            LOG(ERROR) << "Error in getting number of rows for HDU "  << hduname << 
                    " of L2 event file " <<l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        //TIME
        fits_get_colnum(fl2evt, CASEINSEN, "Time", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for Time column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TDOUBLE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].time.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing TIME COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //CZTSECCNT
        fits_get_colnum(fl2evt, CASEINSEN, "CZTSECCNT", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for CZTSECCNT column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].cztseccnt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing CZTSECCNT COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //CZTNTICK
        fits_get_colnum(fl2evt, CASEINSEN, "CZTNTICK", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for CZTNTICK column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].cztntick.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing CZTNTICK COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //PHA
        fits_get_colnum(fl2evt, CASEINSEN, "PHA", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for PHA column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].pha.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing PHA COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //DETID
        fits_get_colnum(fl2evt, CASEINSEN, "DETID", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DETID column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].detid.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DETID COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //PIXID
        fits_get_colnum(fl2evt, CASEINSEN, "PIXID", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for PIXID column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].pixid.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing PIXID COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //DETX
        fits_get_colnum(fl2evt, CASEINSEN, "DETX", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DETX column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].detx.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DETX COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //DETY
        fits_get_colnum(fl2evt, CASEINSEN, "DETY", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DETY column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].dety.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DETY COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //VETO
        fits_get_colnum(fl2evt, CASEINSEN, "VETO", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for VETO column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].veto.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing VETO COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //ALPHA
        fits_get_colnum(fl2evt, CASEINSEN, "ALPHA", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for ALPHA column of L2 event file " << 
                    l2eventFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, quadEvtData[iquad].alpha.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing ALPHA COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        }
		   
    }//END QUADRANT EVENT DATA

    //Writing Veto Spectrum Data
    hduname="VETOSPECTRUM";
    nrows_vec = l1data.vsd.get_nrows();
    fits_movnam_hdu(fl2evt, BINARY_TBL, (char*) (hduname).c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << hduname << " of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fl2evt, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows for HDU " << hduname <<
                " of L2 event file " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //TIME
    fits_get_colnum(fl2evt, CASEINSEN, "Time", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for Time column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, l1data.vsd.time.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing TIME COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //CZTSECCNT
    fits_get_colnum(fl2evt, CASEINSEN, "CZTSECCNT", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for CZTSECCNT column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, l1data.vsd.cztseccnt.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing CZTSECCNT COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //VETO SPECTRUM
    fits_get_colnum(fl2evt, CASEINSEN, "VETOSPEC", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for VETOSPEC column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (i = 0; i < nrows_vec; i++) {
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1 + i, 1,
                VETOSPEC_SIZE, l1data.vsd.vetoSpectrum[i].data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing VETOSPECTRUM COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }

/*
 * Added by Mithun to Write QUADID to the Vetospectrum extension (08/12/15)
 * */

    fits_get_colnum(fl2evt, CASEINSEN, "QUADID", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for QUADID column of L2 event file Veto extn " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, l1data.vsd.quadID.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing QUADID COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


// End of edit

    //VETOSPECTRUM HDU DATA WRITTEN
    
    //WRITING SSM DATA
    hduname = "SSM Data";
    nrows_vec = l1data.ssmd.get_nrows();
    fits_movnam_hdu(fl2evt, BINARY_TBL, (char*) (hduname).c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << hduname << " of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fl2evt, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows for HDU " << hduname <<
                " of L2 event file " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //TIME
    fits_get_colnum(fl2evt, CASEINSEN, "Time", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for Time column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, l1data.ssmd.time.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing TIME COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //CZTSECCNT
    fits_get_colnum(fl2evt, CASEINSEN, "CZTSECCNT", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for CZTSECCNT column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, l1data.ssmd.cztseccnt.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing CZTSECCNT COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //VETO SPECTRUM  
    fits_get_colnum(fl2evt, CASEINSEN, "VETOSPEC", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for VETOSPEC column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (i = 0; i < nrows_vec; i++) {
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1 + i, 1,
                VETOSPEC_SIZE, l1data.ssmd.vetoSpec[i].data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing VETOSPEC COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    //CZTSPEC
    fits_get_colnum(fl2evt, CASEINSEN, "CZTSPEC", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for CZTSPEC column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (i = 0; i < nrows_vec; i++) {
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1 + i, 1,
                CZTSPEC_SIZE, l1data.ssmd.cztSpec[i].data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing CZTSPEC COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    //CZTSPECVETO
    fits_get_colnum(fl2evt, CASEINSEN, "CZTSPECVETO", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for CZTSPECVETO column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (i = 0; i < nrows_vec; i++) {
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1 + i, 1,
                CZTSPECVETO_SIZE, l1data.ssmd.cztSpecVeto[i].data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing CZTSPECVETO COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    //CZTSPECALPHA
    fits_get_colnum(fl2evt, CASEINSEN, "CZTSPECALPHA", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for CZTSPECALPHA column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (i = 0; i < nrows_vec; i++) {
        fits_write_col(fl2evt, TUSHORT, colnum, nrows + 1 + i, 1,
                CZTSPECALPHA_SIZE, l1data.ssmd.cztSpecAlpha[i].data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing CZTSPECALPHA COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }

/*
 * Added by Mithun to Write QUADID to the SSM data extension (08/12/15)
 * */ 

    fits_get_colnum(fl2evt, CASEINSEN, "QUADID", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for QUADID column of L2 event file SSM extn " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, l1data.ssmd.quadID.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing QUADID COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


// End of edit

    
    //SSM DATA WRITTEN
    
    //WRITING TEMPERATURE DATA
    hduname = "TEMP";
    nrows_vec = l1data.tempd.get_nrows();
    fits_movnam_hdu(fl2evt, BINARY_TBL, (char*) (hduname).c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << hduname << " of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fl2evt, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows for HDU " << hduname <<
                " of L2 event file " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //TIME
    fits_get_colnum(fl2evt, CASEINSEN, "Time", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for Time column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, l1data.tempd.time.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing TIME COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //CZTSECCNT
    fits_get_colnum(fl2evt, CASEINSEN, "CZTSECCNT", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for CZTSECCNT column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TDOUBLE, colnum, nrows + 1, 1, nrows_vec, l1data.tempd.cztseccnt.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing CZTSECCNT COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //TEMPERATURE 
    fits_get_colnum(fl2evt, CASEINSEN, "TEMPERATURE", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for TEMPERATURE column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    for (i = 0; i < nrows_vec; i++) {
        fits_write_col(fl2evt, TFLOAT, colnum, nrows + 1 + i, 1,
                NUM_DET_PER_QUAD, l1data.tempd.temperature[i].data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing TEMPERATURE COLUMN of " << l2eventFilename <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }    
    //QUADFLAG
    fits_get_colnum(fl2evt, CASEINSEN, "QUAD_FLAG", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for QUAD_FLAG column of L2 event file " <<
                l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_write_col(fl2evt, TBYTE, colnum, nrows + 1, 1, nrows_vec, l1data.tempd.quadID.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing QUAD_FLAG COLUMN of " << l2eventFilename <<
                " for HDU " << hduname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //TEMPERATURE DATA WRITTEN
    
    fits_close_file(fl2evt, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-2 event file " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
  


//writing bunch clean file
/* Modified by Ajay; Adding additional information about bunch, two columns, detid pixid
 * */




    fits_open_file(&bunchptr, (char*) bunchFile.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-2 bunch file: " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    for(iquad=0; iquad<4; iquad++)
    {
        //Sorting object quadEvtData
        nrows_vec=quadEvtData[iquad].bunch_time.size();
        hduname = quadExtName(iquad);
       
	if(nrows_vec!=0)
	{
	    
        fits_movnam_hdu(bunchptr, BINARY_TBL, (char*) (hduname).c_str(), NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << hduname << " of L2 event file " <<
                   bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_get_num_rows(bunchptr, &nrows, &status);
        if (status) {
            LOG(ERROR) << "Error in getting number of rows for HDU "  << hduname << 
                    " of L2 event file " <<bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	//TIME
        fits_get_colnum(bunchptr, CASEINSEN, "Time", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for Time column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TDOUBLE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].bunch_time.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing TIME COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	//Time_dfs	
        fits_get_colnum(bunchptr, CASEINSEN, "Time_dfs", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for Time_dfs column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].time_dfs.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing TIME_dfs COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }


	//Time_dsl	
        fits_get_colnum(bunchptr, CASEINSEN, "Time_dsl", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for Time_dsl column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].time_dsl.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing TIME_dsl COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	//NumEvent	
        fits_get_colnum(bunchptr, CASEINSEN, "NumEvent", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for NumEvent column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].num_bunchevents.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing NumEvent COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }


	//detid1	
        fits_get_colnum(bunchptr, CASEINSEN, "DetId1", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId1 column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].detid1.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing detid1 COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	//detid2	
        fits_get_colnum(bunchptr, CASEINSEN, "DetId2", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId2 column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].detid2.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing detid2 COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	//detid3	
        fits_get_colnum(bunchptr, CASEINSEN, "DetId3", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId3 column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].detid3.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing detid3 COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }


	//detid4	
        fits_get_colnum(bunchptr, CASEINSEN, "DetId4", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId4 column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].detid4.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing detid4 COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	//RowNum	
        fits_get_colnum(bunchptr, CASEINSEN, "RowNum", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for RowNum column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_write_col(bunchptr, TLONG, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].evt_row_num.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing RowNum COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
	/*Added by Ajay, 22 May 2020
	 * For new noise reduction code added new columns, detid and pixid for first three events, which will be removed in bunch cleaning
	 * */
	//RowNum	
	fits_get_colnum(bunchptr, CASEINSEN, "DetId_fevt", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId_fevt column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
	
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].DetId_fevt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DetId_fevt COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_get_colnum(bunchptr, CASEINSEN, "PixId_fevt", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for PixId_fevt column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].PixId_fevt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing PixId_fevt COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_get_colnum(bunchptr, CASEINSEN, "DetId_sevt", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId_sevt column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].DetId_sevt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DetId_sevt COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_get_colnum(bunchptr, CASEINSEN, "PixId_sevt", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for PixId_sevt column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].PixId_sevt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing PixId_sevt COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

	fits_get_colnum(bunchptr, CASEINSEN, "DetId_tevt", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for DetId_tevt column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].DetId_tevt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DetId_tevt COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_get_colnum(bunchptr, CASEINSEN, "PixId_tevt", &colnum, &status);
        if (status) {
            LOG(ERROR) << "Error in getting column number for PixId_tevt column of L2 event file " << 
                    bunchFile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_col(bunchptr, TBYTE, colnum, nrows+1,1,nrows_vec, quadEvtData[iquad].PixId_tevt.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing PixId_tevt COLUMN of " << bunchFile <<
                    " for HDU " << hduname;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }


//	detid[evt_row_num];
//	pixid[evt_row_num]

	}
    }
    fits_close_file(bunchptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-2 bunch file " << bunchFile;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
   

    quadEvtData.clear();
    return EXIT_SUCCESS;
}

int EventFileHandler::create_l2hdr_file(string l2hdrFilename, string hdrTemplate, string l1eventFilename) {
        int status=0;
    fitsfile *fl2evt, *fl1hdr;
    int i=0,j=0;
    int iquad=0;
    
    vector <string> vecKeywords;
    if(create_empty_fitsfile(l2hdrFilename, hdrTemplate)){
        LOG(ERROR) << "Error in creating empty level-2 HEADER file " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    
    //Opening level-1 event file
    fits_open_file(&fl1hdr, (char*) l1eventFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-1 header file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Opening level-2 event file
    fits_open_file(&fl2evt, (char*) l2hdrFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-2 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Keywords to be copied
    vecKeywords.push_back("DATE-OBS");
    vecKeywords.push_back("TIME-OBS");
    vecKeywords.push_back("DATE-END");
    vecKeywords.push_back("TIME-END");
    vecKeywords.push_back("TIMESYS");
    vecKeywords.push_back("MJDREFI");
    vecKeywords.push_back("MJDREFF");
    vecKeywords.push_back("TIMEZERI");
    vecKeywords.push_back("TIMEZERF");
    vecKeywords.push_back("TIMEDEL");
    vecKeywords.push_back("TIMEUNIT");
    vecKeywords.push_back("TSTARTI");
    vecKeywords.push_back("TSTARTF");
    vecKeywords.push_back("TSTOPI");
    vecKeywords.push_back("TSTOPF");
    vecKeywords.push_back("OBJECT");
    vecKeywords.push_back("RA_PNT");
    vecKeywords.push_back("DEC_PNT");
    vecKeywords.push_back("EQUINOX");
    vecKeywords.push_back("RADECSYS");
    vecKeywords.push_back("OBSERVER");
    vecKeywords.push_back("OBS_ID");
    vecKeywords.push_back("SOURCEID");
    vecKeywords.push_back("TARGETID");
    vecKeywords.push_back("EXPOSURE");
    vecKeywords.push_back("ORB_NUM");
    
   
    if (copyUserKeyWords(fl1hdr, fl2evt, "Primary", "Primary", vecKeywords)) { //copying keywords from input file
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;
    }

    fits_close_file(fl1hdr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-1 event file " << l1eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fl2evt, &status);
    if (status) {
        LOG(ERROR) << "Error in closing DPI/DPH file " << l2hdrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }



    return (EXIT_SUCCESS); 
}

int EventFileHandler::write_l2_hdr_file(string l2hdrFilename, L1evtData l1data) {
    int status=0;
    fitsfile *fhdr;
    long i=0, j=0;
    long nrows=0;
    string hduname="Header";
    
    fits_open_file(&fhdr, (char*) l2hdrFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening header file " << l2hdrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_movnam_hdu(fhdr, BINARY_TBL, (char*) hduname.c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << hduname << " of file " << l2hdrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_rows(fhdr, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows for HDU " << hduname <<
                " of L2 Header file " << l2hdrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing respective data
    if(write_fits_column(fhdr, "Time", TDOUBLE, nrows+1, 1, l1data.hdrd.time)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "Cztseccnt", TDOUBLE, nrows+1, 1, l1data.hdrd.cztseccnt)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "DATAID", TBYTE, nrows+1, 1, l1data.hdrd.dataID)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "MODEID", TBYTE, nrows+1, 1, l1data.hdrd.modeID)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "TOTALEVENTS", TUSHORT, nrows+1, 1, l1data.hdrd.wordCount)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "FRAMENO", TUSHORT, nrows+1, 1, l1data.hdrd.frameNo)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "STATUSWORD", TUSHORT, nrows+1, 1, l1data.hdrd.statusFlag)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "WRITEPKTNO", TUSHORT, nrows+1, 1, l1data.hdrd.writePktNo)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "READPKTNO", TUSHORT, nrows+1, 1, l1data.hdrd.readPktNo)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if(write_fits_column(fhdr, "LAXPCTIME", TDOUBLE, nrows+1, 1, l1data.hdrd.laxpcTime)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }    
    if(write_fits_column(fhdr, "DCNT", TUSHORT, nrows+1, 1, l1data.hdrd.dcnt)){
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }
    if (write_fits_column(fhdr, "ERRORCOUNT", TBYTE, nrows + 1, 1, l1data.hdrd.errorCount)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }    
    if (write_fits_column(fhdr, "ERRORFLAG", TBYTE, nrows + 1, 1, l1data.hdrd.errorFlag)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }    
    if (write_fits_column(fhdr, "BOOTPAGENO", TBYTE, nrows + 1, 1, l1data.hdrd.bootPageNo)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "CZTNUMBER", TBYTE, nrows + 1, 1, l1data.hdrd.cztNo)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "CZTSTATUS", TBYTE, nrows + 1, 1, l1data.hdrd.cztStatus)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "CZTEVENTREADMODE", TBYTE, nrows + 1, 1, l1data.hdrd.evtReadMode)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "COMMANDSTATUS", TBYTE, nrows + 1, 1, l1data.hdrd.cmdStatus)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "BUFFERNO", TBYTE, nrows + 1, 1, l1data.hdrd.bufferNo)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "BASEADDRESS", TBYTE, nrows + 1, 1, l1data.hdrd.baseAdd)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "VETOSPECRANGE", TBYTE, nrows + 1, 1, l1data.hdrd.vetoSpecRange)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "HKCHANNELNO", TBYTE, nrows + 1, 1, l1data.hdrd.channelNo)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "ADCOUTPUT", TUSHORT, nrows + 1, 1, l1data.hdrd.ADCoutput)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "FEBCOMMAND", TUSHORT, nrows + 1, 1, l1data.hdrd.cmd1sec)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "ALPHACOUNT", TUSHORT, nrows + 1, 1, l1data.hdrd.alphaCount)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "VETOCOUNT", TUSHORT, nrows + 1, 1, l1data.hdrd.vetoCount)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "COUNTLTULD", TUSHORT, nrows + 1, 1, l1data.hdrd.cztCount_lt_uld)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "COUNTGTULD", TUSHORT, nrows + 1, 1, l1data.hdrd.cztCount_gt_uld)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    if (write_fits_column(fhdr, "CZTDATAREAD", TUSHORT, nrows + 1, 1, l1data.hdrd.cztdataRead)) {
        LOG(ERROR) << "Error in writing header file: " << l2hdrFilename;
        return EXIT_FAILURE;
    }  
    //ALL COLUMNS WRITTEN
    
    fits_close_file(fhdr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-2 header file: " << l2hdrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
    
    
}

int EventFileHandler::create_l2_gti_file(string l2gtiFilename, string gtiTemplate, string l1eventFilename){
    int status = 0;
    fitsfile *fl2gti, *fl1hdr;
    int i = 0, j = 0;

    vector <string> vecKeywords;
    if (create_empty_fitsfile(l2gtiFilename, gtiTemplate)) {
        LOG(ERROR) << "Error in creating empty level-2 GTI file " << l2gtiFilename;
        return EXIT_FAILURE;
    }

    //Opening level-1 event file
    fits_open_file(&fl1hdr, (char*) l1eventFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-1 header file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Opening level-2 event file
    fits_open_file(&fl2gti, (char*) l2gtiFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-2 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Keywords to be copied
    vecKeywords.push_back("DATE-OBS");
    vecKeywords.push_back("TIME-OBS");
    vecKeywords.push_back("DATE-END");
    vecKeywords.push_back("TIME-END");
    vecKeywords.push_back("TIMESYS");
    vecKeywords.push_back("MJDREFI");
    vecKeywords.push_back("MJDREFF");
    vecKeywords.push_back("TIMEZERI");
    vecKeywords.push_back("TIMEZERF");
    vecKeywords.push_back("TIMEDEL");
    vecKeywords.push_back("TIMEUNIT");
    vecKeywords.push_back("TSTARTI");
    vecKeywords.push_back("TSTARTF");
    vecKeywords.push_back("TSTOPI");
    vecKeywords.push_back("TSTOPF");
    vecKeywords.push_back("OBJECT");
    vecKeywords.push_back("RA_PNT");
    vecKeywords.push_back("DEC_PNT");
    vecKeywords.push_back("EQUINOX");
    vecKeywords.push_back("RADECSYS");
    vecKeywords.push_back("OBSERVER");
    vecKeywords.push_back("OBS_ID");
    vecKeywords.push_back("SOURCEID");
    vecKeywords.push_back("TARGETID");
    vecKeywords.push_back("EXPOSURE");
    vecKeywords.push_back("ORB_NUM");


    if (copyUserKeyWords(fl1hdr, fl2gti, "Primary", "Primary", vecKeywords)) { //copying keywords from input file
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;
    }

    fits_close_file(fl1hdr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-1 event file " << l1eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fl2gti, &status);
    if (status) {
        LOG(ERROR) << "Error in closing GTI file " << l2gtiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    return (EXIT_SUCCESS);
}


int EventFileHandler::write_l2_gti_file(string l2gtiFilename, vector<double> tstart, vector<double> tstop) {
    int status=0;
    fitsfile *fgti;
    long i = 0, j = 0;
    long nrows = 0;
    string hduname = "GTI";

    fits_open_file(&fgti, (char*) l2gtiFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening header file " << l2gtiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_movnam_hdu(fgti, BINARY_TBL, (char*) hduname.c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << hduname << " of file " << l2gtiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Writing respective data
    if (write_fits_column(fgti, "START", TDOUBLE, 1, 1, tstart)) {
        LOG(ERROR) << "Error in writing gti file: " << l2gtiFilename;
        return EXIT_FAILURE;
    }
    if (write_fits_column(fgti, "STOP", TDOUBLE, 1, 1, tstop)) {
        LOG(ERROR) << "Error in writing gti file: " << l2gtiFilename;
        return EXIT_FAILURE;
    }

    fits_close_file(fgti, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-2 GTI file: " << l2gtiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
    
}


int EventFileHandler::create_l2_bti_file(string l2btiFilename, string btiTemplate, string l1eventFilename){
    int status = 0;
    fitsfile *fl2bti, *fl1hdr;
    int i = 0, j = 0;

    vector <string> vecKeywords;
    if (create_empty_fitsfile(l2btiFilename, btiTemplate)) {
        LOG(ERROR) << "Error in creating empty level-2 bti file " << l2btiFilename;
        return EXIT_FAILURE;
    }

    //Opening level-1 event file
    fits_open_file(&fl1hdr, (char*) l1eventFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-1 header file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Opening level-2 event file
    fits_open_file(&fl2bti, (char*) l2btiFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-2 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Keywords to be copied
    vecKeywords.push_back("DATE-OBS");
    vecKeywords.push_back("TIME-OBS");
    vecKeywords.push_back("DATE-END");
    vecKeywords.push_back("TIME-END");
    vecKeywords.push_back("TIMESYS");
    vecKeywords.push_back("MJDREFI");
    vecKeywords.push_back("MJDREFF");
    vecKeywords.push_back("TIMEZERI");
    vecKeywords.push_back("TIMEZERF");
    vecKeywords.push_back("TIMEDEL");
    vecKeywords.push_back("TIMEUNIT");
    vecKeywords.push_back("TSTARTI");
    vecKeywords.push_back("TSTARTF");
    vecKeywords.push_back("TSTOPI");
    vecKeywords.push_back("TSTOPF");
    vecKeywords.push_back("OBJECT");
    vecKeywords.push_back("RA_PNT");
    vecKeywords.push_back("DEC_PNT");
    vecKeywords.push_back("EQUINOX");
    vecKeywords.push_back("RADECSYS");
    vecKeywords.push_back("OBSERVER");
    vecKeywords.push_back("OBS_ID");
    vecKeywords.push_back("SOURCEID");
    vecKeywords.push_back("TARGETID");
    vecKeywords.push_back("EXPOSURE");
    vecKeywords.push_back("ORB_NUM");


    if (copyUserKeyWords(fl1hdr, fl2bti, "Primary", "Primary", vecKeywords)) { //copying keywords from input file
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;
    }

    fits_close_file(fl1hdr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-1 event file " << l1eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fl2bti, &status);
    if (status) {
        LOG(ERROR) << "Error in closing bti file " << l2btiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    return (EXIT_SUCCESS);
}

int EventFileHandler::write_l2_bti_file(string l2btiFilename, vector<double> tstart, vector<double> tstop) {
    int status=0;
    fitsfile *fbti;
    long i = 0, j = 0;
    long nrows = 0;
    string hduname = "BTI";

    fits_open_file(&fbti, (char*) l2btiFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening header file " << l2btiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_movnam_hdu(fbti, BINARY_TBL, (char*) hduname.c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << hduname << " of file " << l2btiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Writing respective data
    if (write_fits_column(fbti, "START", TDOUBLE, 1, 1, tstart)) {
        LOG(ERROR) << "Error in writing gti file: " << l2btiFilename;
        return EXIT_FAILURE;
    }
    if (write_fits_column(fbti, "STOP", TDOUBLE, 1, 1, tstop)) {
        LOG(ERROR) << "Error in writing gti file: " << l2btiFilename;
        return EXIT_FAILURE;
    }

    fits_close_file(fbti, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-2 GTI file: " << l2btiFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
    
}

int EventFileHandler::read_event_temp(fitsfile* fptr, int quad_no){
    int status=0; //status variable
    int i,j=0; //counter variables
    int colnum=0; //
    long nrows=0; //number of rows in HDU extension
    char temp_char[MAX_KEYWORD_SIZE];
    string errorMsg="";
    double *timeArray;
    double *localTimeArray;
    unsigned char* quadFlagArray;
    float* tempColArray; //NOTE: It stores value of only one record of temperature column
    string tempColName ="";
    
    
    // clearing object variables
    
    // object variables cleared
    
    fits_movnam_hdu(fptr, BINARY_TBL, "TEMP", NULL, &status);
    errorMsg = "Error in moving to HDU TEMP of event file";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in TEMP extension of Level-2 Event File.";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

   
    // Reading TIME, LOCAL_TIME & QUAD_FLAG columns from TEMP extension of Event File.
    fits_get_colnum(fptr, CASEINSEN, "TIME", &colnum, &status);
    errorMsg = "Error in getting column number for TIME column in TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    timeArray = new double[nrows];
    for(i=0; i<nrows; i++){
        timeArray[i]=0.0;
    }
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, timeArray, NULL, &status);
    errorMsg = "Error in reading TIME extension  of TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
   
    fits_get_colnum(fptr, CASEINSEN, "CZTSECCNT", &colnum, &status);
    errorMsg = "Error in getting column number for CZTSECCNT column in TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    localTimeArray = new double[nrows];
    for(i=0; i<nrows; i++){
        localTimeArray[i]=0.0;
    }
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, localTimeArray, NULL, &status);
    errorMsg = "Error in reading TIME extension of TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;} 

    fits_get_colnum(fptr, CASEINSEN, "QUAD_FLAG", &colnum, &status);
    errorMsg = "Error in getting column number for QUAD_FLAG column in TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    quadFlagArray = new unsigned char[nrows];
    for(i=0; i<nrows; i++){
        quadFlagArray[i]=0;
    }
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, quadFlagArray, NULL, &status);
    errorMsg = "Error in reading QUAD_FLAG extension  of TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    // TIME, LOCAL_TIME & QUAD_FLAG columns from TEMP extension of Event File has been read.
    
     
    fits_get_colnum(fptr, CASEINSEN, "TEMPERATURE", &colnum, &status);
    errorMsg = "Error in getting column number for " + tempColName + " column in TEMP extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    tempColArray = new float[NO_DET_PER_QUAD];
    for(i=0; i<NO_DET_PER_QUAD; i++){
       tempColArray[i]=0.0;
    }
    
 
    
    vector <float> vecTempRec(16,0.0); //temporary vector to store tempColArray. This vector will be later pushed back into the tempQuadrant object variable.

    // Reading only those records from temperature column that belongs to the requisite quadrant number (as supplied by user)
    // Pushing the pertinent information i.e. info belonging to the requisite quadrant into the object variables.
    vecTempQuadrant.clear();
    vecTimeTemp.clear();
    vecLocalTimeTemp.clear();
    for (i = 0; i < nrows; i++) {
        if ((int) quadFlagArray[i] == quad_no) {
            fits_read_col(fptr, TFLOAT, colnum, i+1, 1, 16, NULL, tempColArray, NULL, &status);
            sprintf(temp_char,"%d",i+1);
            errorMsg = "Error in reading record number " + string(temp_char) + " of column " + tempColName + " in TEMP extension of Level-2 Event File";
            if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
//            vecTempRec.insert(vecTempRec.begin(), &tempColArray[0], &tempColArray[NO_DET_PER_QUAD]);
//            vecTempRec.resize(NO_DET_PER_QUAD);
            for(j=0; j<NO_DET_PER_QUAD; j++){
                vecTempRec[j] = tempColArray[j];
            }
            vecTempQuadrant.push_back(vecTempRec);
            vecTimeTemp.push_back(timeArray[i]);
            vecLocalTimeTemp.push_back(localTimeArray[i]);
        }
    }
    // TEPERATURE column has been read.
    // tempQuadrant, time, localTime object variables have been filled by the values from the Level-2 input file.    
    
    LOG(INFO) << "Extracted " << vecTempQuadrant.size() << " Temperature values for Quadrant " << quad_no;

    // Deleting dynamically declared variables
    delete[] timeArray, localTimeArray, quadFlagArray, tempColArray;
    
    return status;
}

int EventFileHandler::read_quad_extension(fitsfile* fptr, int quad_no, string colname, long startRowno, long endRowno){
    int status=0;
    int i,j=0;
    int colnum=0;
    long nelements=0;
    long nrows=0;
    string errorMsg="";
    char filename[PIL_LINESIZE];
    string tempExtName = "";

    switch (quad_no) {
        case 0:
            tempExtName = "Q0";
            break;
        case 1:
            tempExtName = "Q1";
            break;
        case 2:
            tempExtName = "Q2";
            break;
        case 3:
            tempExtName = "Q3";
            break;
    }

    fits_file_name(fptr, filename, &status);
    if (status) {
        LOG(WARNING) << "Error in getting filename.";
    }

    fits_movnam_hdu(fptr, BINARY_TBL, (char*) tempExtName.c_str(), NULL, &status);
    errorMsg = "Error in moving to " + tempExtName + " extension of output file.";
    report_error(status, errorMsg);

    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in TEMP extension of Level-2 Event File.";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }
    if(endRowno > nrows){
        endRowno=nrows;
    }
    
    nelements = endRowno - startRowno +1;
    
    if(colname=="TIME") {
        vecTimeQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "TIME", TDOUBLE, startRowno, 1, nelements, vecTimeQuad)) {
                LOG(ERROR) << "Error in reading column Time of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    } 
    //Reading Cztseccnt
    else if(colname=="CZTSECCNT") {
        vecLocalTimeQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "CZTSECCNT", TDOUBLE, startRowno, 1, nelements, vecLocalTimeQuad)) {
                LOG(ERROR) << "Error in reading column Time of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    
    //Reading Cztnticks
    else if(colname=="CZTNTICK"){   
        vecCztntickQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "CZTNTICK", TUSHORT, startRowno, 1, nelements, vecCztntickQuad)) {
                LOG(ERROR) << "Error in reading column CZTNTICK of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    //Reading PHA
    else if(colname=="PHA"){   
        vecPHAQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "PHA", TUSHORT, startRowno, 1, nelements, vecPHAQuad)) {
                LOG(ERROR) << "Error in reading column PHA of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    //Reading DETID
    else if(colname=="DETID"){   
        vecDetIDQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "DETID", TBYTE, startRowno, 1, nelements, vecDetIDQuad)) {
                LOG(ERROR) << "Error in reading column DETID of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    //Reading PIXID
    else if(colname=="PIXID"){   
        vecPixIDQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "PIXID", TBYTE, startRowno, 1, nelements, vecPixIDQuad)) {
                LOG(ERROR) << "Error in reading column PIXID of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    //Reading DETX
    else if(colname=="DETX"){   
        vecDetXQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "DETX", TBYTE, startRowno, 1, nelements, vecDetXQuad)) {
                LOG(ERROR) << "Error in reading column DETX of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    //Reading DETY
    else if(colname=="DETY"){   
        vecDetYQuad.clear();
        try {
            if (read_fits_column_if_available(fptr, "DETY", TBYTE, startRowno, 1, nelements, vecDetYQuad)) {
                LOG(ERROR) << "Error in reading column DETY of " << tempExtName << " extension of level-2 event file.";
                return EXIT_FAILURE;
            }
        } catch (ErrorHandler errHandler) {
            LOG(INFO) << "Filename: " << eventFilename;
            logError(errHandler);
        }
    }
    
    return status;
}
int EventFileHandler::read_quad_extension(fitsfile* fptr, int quad_no){
    int status=0; //status variable
    int i,j=0; //counter variables
    int colnum=0; //
    long nrows=0; //number of rows in HDU extension
    string errorMsg="";
    char filename[PIL_LINESIZE];
    
    double *timeArray;
    double *localTimeArray;
    unsigned short* phaArray;
    unsigned char* detIDArray;
    unsigned char* pixIDArray;
    unsigned short* vetoArray;
    unsigned char* alphaArray;
    unsigned char* detXArray;
    unsigned char* detYArray;
    unsigned short* cztntickArray;
    float * energyArray;

    char tempChar[MAX_KEYWORD_SIZE];
    string tempExtName="";
    
    switch (quad_no) {
        case 0:
            tempExtName="Q0";
            break;
        case 1:
            tempExtName="Q1";
            break;
        case 2:
            tempExtName="Q2";
            break;
        case 3:
            tempExtName="Q3";
            break;   
    }
    
    fits_file_name(fptr, filename, &status);
    if (status) {
        LOG(WARNING) << "Error in getting filename.";
    }
    
    fits_movnam_hdu(fptr, BINARY_TBL, (char*) tempExtName.c_str(), NULL, &status);
    errorMsg = "Error in moving to " + tempExtName + " extension of output file.";
    report_error(status, errorMsg);
    
    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in TEMP extension of Level-2 Event File.";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;
    }

    //Initializing arrays
    timeArray = new double[nrows];
    localTimeArray = new double[nrows];
    phaArray = new unsigned short[nrows];
    detIDArray = new unsigned char[nrows];
    pixIDArray = new unsigned char[nrows];
    alphaArray = new unsigned char[nrows];
    detXArray = new unsigned char[nrows];
    vetoArray = new unsigned short[nrows];
    detYArray = new unsigned char[nrows];
    cztntickArray = new unsigned short[nrows];
    energyArray = new float [nrows];
    for(i=0; i<nrows; i++){
        timeArray[i]=0.0;
        localTimeArray[i]=0.0;
        phaArray[i]=0;
        detIDArray[i]=0;
        pixIDArray[i]=0;
        vetoArray[i]=0;
        alphaArray[i]=0;
        detXArray[i]=0;
        detYArray[i]=0;
        cztntickArray[i]=0;
        energyArray[i]=0.0;
    }
    
    //clearing class vectors
    vecTimeQuad.clear();
    vecLocalTimeQuad.clear();
    vecPHAQuad.clear();
    vecDetIDQuad.clear();
    vecPixIDQuad.clear();
    vecVetoQuad.clear();
    vecAlphaQuad.clear();
    vecDetXQuad.clear();
    vecDetYQuad.clear();
    vecCztntickQuad.clear();
    vecEnergyQuad.clear();
    //class vectors cleared
    
    // Reading TIME, CZTSECCNT, PHA, DetID, pixID, veto, alpha, DETX, DETY & CZTNTICK columns from Q_ extension of Event File.
    fits_get_colnum(fptr, CASEINSEN, "TIME", &colnum, &status);
    errorMsg = "Error in getting column number for TIME column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, timeArray, NULL, &status);
    errorMsg = "Error in reading TIME column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "CZTSECCNT", &colnum, &status);
    errorMsg = "Error in getting column number for CZTSECCNT column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, localTimeArray, NULL, &status);
    errorMsg = "Error in reading CZTSECCNT column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_colnum(fptr, CASEINSEN, "PHA", &colnum, &status);
    errorMsg = "Error in getting column number for PHA column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TUSHORT, colnum, 1, 1, nrows, NULL, phaArray, NULL, &status);
    errorMsg = "Error in reading PHA column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    

    fits_get_colnum(fptr, CASEINSEN, "DetID", &colnum, &status);
    errorMsg = "Error in getting column number for DetID column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detIDArray, NULL, &status);
    errorMsg = "Error in reading DetID column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}      

    fits_get_colnum(fptr, CASEINSEN, "pixID", &colnum, &status);
    errorMsg = "Error in getting column number for pixID column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, pixIDArray, NULL, &status);
    errorMsg = "Error in reading pixID column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "veto", &colnum, &status);
    errorMsg = "Error in getting column number for veto column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TUSHORT, colnum, 1, 1, nrows, NULL, vetoArray, NULL, &status);
    errorMsg = "Error in reading veto column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "alpha", &colnum, &status);
    errorMsg = "Error in getting column number for alpha column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, alphaArray, NULL, &status);
    errorMsg = "Error in reading alpha column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}   
    
    fits_get_colnum(fptr, CASEINSEN, "DETX", &colnum, &status);
    errorMsg = "Error in getting column number for DETX column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detXArray, NULL, &status);
    errorMsg = "Error in reading DETX column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "DETY", &colnum, &status);
    errorMsg = "Error in getting column number for DETY column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TBYTE, colnum, 1, 1, nrows, NULL, detYArray, NULL, &status);
    errorMsg = "Error in reading DETY column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_get_colnum(fptr, CASEINSEN, "CZTNTICK", &colnum, &status);
    errorMsg = "Error in getting column number for CZTNTICK column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    fits_read_col(fptr, TUSHORT, colnum, 1, 1, nrows, NULL, cztntickArray, NULL, &status);
    errorMsg = "Error in reading CZTNTICK column of " + tempExtName + " extension of Level-2 Event File";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    // TIME, CZTSECCNT, PHA, DetID, pixID, veto, alpha, DETX, DETY & CZTNTICK columns from Q_ extension have been read.
    
    //Reading Energy
    try {
        if (read_fits_column_if_available(fptr, "ENERGY", TFLOAT, 1, 1, nrows, vecEnergyQuad)) {
            LOG(ERROR) << "Error in reading column Energy of " << tempExtName << " extension of level-2 event file.";
            return EXIT_FAILURE;
        }
    } catch (ErrorHandler errHandler) {
        LOG(INFO) << "Filename: " << eventFilename;
        logError(errHandler);
    }

    //Reading PI
    try{
        if(read_fits_column_if_available(fptr, "PI", TUSHORT, 1, 1, nrows, vecPIQuad)){
            LOG(ERROR) << "Error in reading column PI of " << tempExtName << " extension of level-2 event file.";
            return EXIT_FAILURE;
        }
    } catch (ErrorHandler errHandler) {
        LOG(INFO) << "Filename: " << eventFilename;
        logError(errHandler);
    }
    
    nrowsQuad = nrows;

    for(i=0; i<nrows; i++){
        vecTimeQuad.push_back(timeArray[i]);
        vecLocalTimeQuad.push_back(localTimeArray[i]);
        vecPHAQuad.push_back(phaArray[i]);
        vecDetIDQuad.push_back(detIDArray[i]);
        vecPixIDQuad.push_back(pixIDArray[i]);
        vecVetoQuad.push_back(vetoArray[i]);
        vecAlphaQuad.push_back(alphaArray[i]);
        vecDetXQuad.push_back(detXArray[i]);
        vecDetYQuad.push_back(detYArray[i]);
        vecCztntickQuad.push_back(cztntickArray[i]);
    }
    
    //deleting dynamically declared variables;
    delete[] timeArray, localTimeArray, phaArray, detIDArray, pixIDArray, vetoArray;
    delete[] alphaArray, detXArray, detYArray, cztntickArray, energyArray;
    DLOG(INFO) << tempExtName << " extension read and data stored in class variables.";
    return status;
}

int EventFileHandler::find_tstart_tstop(fitsfile* fptr, int quadNo, 
        double *tstart, double *tstop, string l2evtFilename) {
    int status=0;
    ErrorHandler errHandler;
    char extname[PIL_LINESIZE];

    long startTimei=0; //integer start time 
    double startTimef=0.0; //fractional start time
    long stopTimei=0; //integer start time 
    double stopTimef=0.0; //fractional start time
    
    if(fptr==NULL && l2evtFilename==""){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Either provide pointer to fits file or the name of fits file.";
        
        throw errHandler;
    }
    
    if(fptr==NULL){
        fits_open_file(&fptr, l2evtFilename.c_str(), READONLY, &status);
        if (status) {
            LOG(ERROR) << "Error in opening event file: " << l2evtFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }

    quadToHDU(quadNo, extname);
    fits_movnam_hdu(fptr, BINARY_TBL, extname, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << extname << " of L2 event file ";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    try{

        readKey(fptr, TDOUBLE, "TSTART", tstart);
        readKey(fptr, TDOUBLE, "TSTOP", tstop);
       /* readKey(fptr, TLONG, "TSTARTI", &startTimei);
        readKey(fptr, TLONG, "TSTOPI", &stopTimei);
        readKey(fptr, TDOUBLE, "TSTARTF", &startTimef);
        readKey(fptr, TDOUBLE, "TSTOPF", &stopTimef);
	*/
    } catch(ErrorHandler errHandler){
        logError(errHandler);
        throw errHandler;
    }
    
//    *tstart = (double)startTimei + startTimef;
 //  *tstop = (double)stopTimei + stopTimef;
    
    return status;
}

int EventFileHandler::get_common_time(fitsfile* fptr, vector<int> quadsToProcess, double* tmin, 
        double* tmax, string l2evtFilename) {
    int status = 0;
    int i = 0, iquad = 0;
    ErrorHandler errHandler;
    char *extname;
    vector <double> tstart;
    vector <double> tstop;
    double tempstart, tempstop;

    if (fptr == NULL && l2evtFilename == "") {
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Either provide pointer to fits file or the name of fits file.";

        throw errHandler;
    }

    if (fptr == NULL) {
        fits_open_file(&fptr, l2evtFilename.c_str(), READONLY, &status);
        if (status) {
            LOG(ERROR) << "Error in opening event file: " << l2evtFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }


    //RESIZING TSTART AND TSTOP
    tstart.resize(quadsToProcess.size(), 0.0);
    tstop.resize(quadsToProcess.size(), 0.0);

    try {
        for (int i = 0; i < quadsToProcess.size(); i++) {
            iquad = quadsToProcess[i];
            find_tstart_tstop(fptr, iquad, &tempstart, &tempstop);
            tmin[i] = tempstart;
            tmax[i] = tempstop;

            /*tstart[i] = tempstart;
            tstop[i] = tempstop;*/
        }
        if (quadsToProcess.size() > 0) {
            //*tmin = *max_element(tstart.begin(), tstart.end());
            //*tmax = *min_element(tstop.begin(), tstop.end());
	   
        } else {
            errHandler.severity = errERROR;
            errHandler.errorStatus = NO_OUTPUT_POSSIBLE;
            errHandler.errorMsg = "Cannot get minimum & maximum time from event file as number of quadrants to be processed is 0.";
            throw errHandler;
        }
    } catch (ErrorHandler errHandler) {
        throw errHandler;
    }

    return status;
}

int EventFileHandler::get_min_max_time(fitsfile* fptr, vector<int> quadsToProcess, double* tmin, double* tmax, string l2evtFilename) {
    int status = 0;
    int i=0,iquad=0;
    ErrorHandler errHandler;
    char *extname;
    vector <double> tstart;
    vector <double> tstop;
    double tempstart, tempstop;

    if (fptr == NULL && l2evtFilename == "") {
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Either provide pointer to fits file or the name of fits file.";

        throw errHandler;
    }

    if (fptr == NULL) {
        fits_open_file(&fptr, l2evtFilename.c_str(), READONLY, &status);
        if (status) {
            LOG(ERROR) << "Error in opening event file: " << l2evtFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }

    
    //RESIZING TSTART AND TSTOP
    tstart.resize(quadsToProcess.size(), 0.0);
    tstop.resize(quadsToProcess.size(), 0.0);

    try{
        for(int i=0; i<quadsToProcess.size(); i++){
            iquad = quadsToProcess[i];
            find_tstart_tstop(fptr, iquad, &tempstart, &tempstop);
            tstart[i] = tempstart;
            tstop[i] = tempstop;
        }
        if(quadsToProcess.size()>0){
            *tmin = *min_element(tstart.begin(), tstart.end());
            *tmax = *max_element(tstop.begin(), tstop.end());
        } else {
            errHandler.severity = errERROR;
            errHandler.errorStatus = NO_OUTPUT_POSSIBLE;
            errHandler.errorMsg = "Cannot get minimum & maximum time from event file as number of quadrants to be processed is 0.";
            throw errHandler;
        }
    } catch(ErrorHandler errHandler){
        throw errHandler;
    }
    
    return status;
}


bool EventFileHandler::event_satisfy_criteria(int recordNo, vector <double> tstart,
        vector <double> tstop, vector <float> estart,
        vector <float> estop, eventFileQuad* evtquad,int tfilter, int efilter) {
    ErrorHandler errHandler;
    bool energyFlag = false;
    bool timeFlag = false;
    bool evtFlag=false;
    int timeLen = 0;
    int energyLen = 0;
    int itime = 0;
    int ienergy =0;
    //Verifying the input

    //ERROR Handling
    if (recordNo >= nrowsQuad) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = RECORD_DOES_NOT_EXIST;
        errHandler.errorMsg = "Record number " + itoa(recordNo) + " is out of range.";
        throw errHandler;
    }

    //WARNING
    if (tstart.size() != tstop.size()) {
        LOG(WARNING) << "Vectors tstart and tstop are of different size.";
        LOG(WARNING) << "Considering the shortest of two vectors.";
    }
    if (estart.size() != estop.size()) {
        LOG(WARNING) << "Vectors estart and estop are of different size.";
        LOG(WARNING) << "Considering the shortest of two vectors.";
    }

    //Main code
    timeLen = (tstart.size() <= tstop.size()) ? tstart.size() : tstop.size();
    energyLen = (estart.size() <= estop.size()) ? estart.size() : estop.size();

    if (tfilter==0)
        timeFlag=true;
    else if(timeLen>0){
        for (itime = 0; itime < timeLen; itime++) {
            if (vecTimeQuad[recordNo] >= tstart[itime] && vecTimeQuad[recordNo] <= tstop[itime]) {
                timeFlag = true;
                break;
            }
        }
    }
    else
    {
        LOG(ERROR)<<"Time selection failed";
        return(EXIT_FAILURE);
    } 
   
    if (efilter==0)
        energyFlag=true;
    else if(energyLen>0){
        for(ienergy = 0; ienergy<energyLen; ienergy++){
            if (vecEnergyQuad[recordNo] >= estart[ienergy] && vecEnergyQuad[recordNo] <= estop[ienergy]){
             energyFlag= true;
             break;
            }
        }
    }
    else
    {
        LOG(ERROR)<<"Energy selection failed";
        return(EXIT_FAILURE);    
    } 
    

    evtFlag = (timeFlag && energyFlag);
    if (evtFlag==true) {
        evtquad->time = vecTimeQuad[recordNo];
        evtquad->detID = vecDetIDQuad[recordNo];
        evtquad->pixID = vecPixIDQuad[recordNo];
        evtquad->veto = vecVetoQuad[recordNo];
        evtquad->alpha = vecAlphaQuad[recordNo];
        evtquad->detX = vecDetXQuad[recordNo];
        evtquad->detY = vecDetYQuad[recordNo];
        evtquad->energy = vecEnergyQuad[recordNo];
        evtquad->PHA = vecPHAQuad[recordNo];
        evtquad->PI = vecPIQuad[recordNo];
    }

    return evtFlag;

}

int EventFileHandler::is_column_empty(string colname) {
    int emptyFlag=1;
    
    if(colname=="PI"){
        emptyFlag = vecPIQuad.empty();
    } else if (colname=="Energy"){
        emptyFlag = vecEnergyQuad.empty();
    } else {
        emptyFlag = -1; //column does not exist.
    }
    
    return emptyFlag;
}

vector<float> EventFileHandler::get_quadrant_temperature(double eventTime, int& status) {
    long nrows = vecTimeTemp.size();
    int j=0; //counter variable
    long i=0; //counter variables
    float slope_temp=0.0;
    float offset_temp=0.0;
    float temperature=0.0;
    vector<float> quadrant_temperature(16,1.0);
    
    long indexMinTime;
    long indexMaxTime;

    quadrant_temperature.clear();     
    //event time is less than earliest temperature record
    if(eventTime<= vecTimeTemp[0]){
        quadrant_temperature = vecTempQuadrant[0];
        //LOG(WARNING)<<"Event Time " << eventTime << " is less than earliest temperature record.";
        return quadrant_temperature;
    }
    //event time is greater then latest temperature record
    else if(eventTime > vecTimeTemp[nrows-1]){
        quadrant_temperature = vecTempQuadrant[nrows-1];
        //LOG(WARNING) << "Event Time " << eventTime << " is greater than latest temperature record.";
        return quadrant_temperature;
    }


    for(i=0; i<nrows-1; i++){
        if(eventTime> vecTimeTemp[i] && eventTime <=vecTimeTemp[i+1]){
            indexMinTime = i;
            indexMaxTime = i+1;
            //LOG(INFO)<< "Event Time " << eventTime << " lies between " << indexMinTime << " and " << indexMaxTime;
            break;
        }
    }
   
    for(j=0; j<NO_DET_PER_QUAD; j++){

        slope_temp = (vecTempQuadrant[indexMaxTime][j] - vecTempQuadrant[i][j])/((float)(vecTimeTemp[indexMaxTime] - vecTimeTemp[indexMinTime]));
        offset_temp = vecTempQuadrant[indexMaxTime][j] - (slope_temp*vecTimeTemp[indexMaxTime]);
        temperature = (slope_temp*eventTime) + offset_temp;
        quadrant_temperature.push_back((float)temperature);
    }


    return quadrant_temperature;    
       
}

// Function to copy 5 GTI extensions from level 1 file to event file

int EventFileHandler::copy_gti_extensions(string l2eventFilename, string l1scienceFilename)
{
    fitsfile *fl1sci,*fptr;
    int status=0,hdutype=0;
    int num_hdus;

    fits_open_file(&fptr, (char*) l2eventFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in re-opening level-2 event file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


    fits_open_file(&fl1sci, (char*) l1scienceFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening level-1 science file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_get_num_hdus(fl1sci,&num_hdus,&status);
    if (status) {
        LOG(ERROR) << "Error in getting num_hdus in level-1 science file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    

    if(num_hdus>7)
    {
        
    fits_movnam_hdu(fl1sci, BINARY_TBL, "GTI", 0, &status);
    if (status) {
        LOG(INFO) << "GTI extension not present in level-1 science file.";
        return (EXIT_SUCCESS);
    }    

    fits_copy_hdu(fl1sci, fptr, 0,&status);    
    if (status) {
        LOG(ERROR) << "Error in copying GTI extension to evt file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    

    fits_movnam_hdu(fl1sci, BINARY_TBL, "Q0_GTI", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in finding Q0_GTI extension in level-1 science file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_copy_hdu(fl1sci, fptr, 0,&status);  
    if (status) {
        LOG(ERROR) << "Error in copying Q0_GTI extension to evt file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }



    fits_movnam_hdu(fl1sci, BINARY_TBL, "Q1_GTI", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in finding Q1_GTI extension in level-1 science file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_copy_hdu(fl1sci, fptr, 0,&status);  
    if (status) {
        LOG(ERROR) << "Error in copying Q1_GTI extension to evt file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_movnam_hdu(fl1sci, BINARY_TBL, "Q2_GTI", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in finding Q2_GTI extension in level-1 science file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_copy_hdu(fl1sci, fptr, 0,&status);  
    if (status) {
        LOG(ERROR) << "Error in copying Q2_GTI extension to evt file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


    fits_movnam_hdu(fl1sci, BINARY_TBL, "Q3_GTI", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in finding Q3_GTI extension in level-1 science file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_copy_hdu(fl1sci, fptr, 0,&status);  
    if (status) {
        LOG(ERROR) << "Error in copying Q3_GTI extension to evt file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    }
    else
    {
        LOG(INFO) << "All GTI extensions are not present..Skipping copying GTIs";
    }


    fits_close_file(fl1sci, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-1 science file " << l1scienceFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing level-2 event file " << l2eventFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
   
    return(EXIT_SUCCESS); 
}


//EVENTFILEHANDLER END

//ASPECT CLASS

AspectFileHandler::AspectFileHandler() {
    nrows = 0;
    RA = 0.0;
    DEC = 0.0;
    TWIST = 0.0;
    teldefFilename = "";
}
int AspectFileHandler::read_aspect_file(string aspectFilename){
    int status = 0; // status variable
    int i, j = 0; // counter variables
    int colnum = 0;
    string errorMsg = "";
    char temp_char[MAX_KEYWORD_SIZE];
    fitsfile *fptr; // Pointer to level-2 ASPECT FILE.
    double *timeArray;
    float *nXArray;
    float *nYArray;
    float *nZArray;
    float *nXtArray;
    float *nYtArray;
    float *nZtArray;
    fitsfile *fptrAspect;
    

    //reading aspect file;
    fits_open_file(&fptrAspect, aspectFilename.c_str(), READONLY, &status );
    errorMsg = "*** Error in opening ASPECT file: " + aspectFilename + " ***";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    fits_movnam_hdu(fptrAspect, BINARY_TBL, "ASPECT", 0, &status);
    errorMsg = "Error in reading ASPECT extension in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    fits_get_num_rows(fptrAspect, &nrows, &status);
    errorMsg = "Error in getting number of rows in ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    //reading pertinent keys from Aspect file
    //RA
    try{
        readFitsKey(fptrAspect, TFLOAT, "RA", &RA, NULL);
        readFitsKey(fptrAspect, TFLOAT, "DEC", &DEC, NULL);
        readFitsKey(fptrAspect, TFLOAT, "TWIST", &TWIST, NULL);
        readFitsKeyStr(fptrAspect, "TELDEF", &teldefFilename, NULL);
    }
    catch (ErrorHandler errHandler){
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    // pertinent keys read.
    
    //reading Time
    fits_get_colnum(fptrAspect, CASEINSEN, "TIME", &colnum, &status);
    errorMsg = "Error in getting column number of TIME column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    timeArray = new double[nrows];
    for(i=0; i<nrows; i++){
        timeArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TDOUBLE, colnum, 1, 1, nrows, NULL, timeArray, NULL, &status);
    errorMsg = "Error in reading TIME column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    //reading nX
    fits_get_colnum(fptrAspect, CASEINSEN, "Nx", &colnum, &status);
    errorMsg = "Error in getting column number of Nx column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    nXArray = new float[nrows];
    for(i=0; i<nrows; i++){
        nXArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TFLOAT, colnum, 1, 1, nrows, NULL, nXArray, NULL, &status);
    errorMsg = "Error in reading Nx column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

    //reading nY
    fits_get_colnum(fptrAspect, CASEINSEN, "Ny", &colnum, &status);
    errorMsg = "Error in getting column number of Ny column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    nYArray = new float[nrows];
    for(i=0; i<nrows; i++){
        nYArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TFLOAT, colnum, 1, 1, nrows, NULL, nYArray, NULL, &status);
    errorMsg = "Error in reading Ny column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    
    //reading nZ
    fits_get_colnum(fptrAspect, CASEINSEN, "Nz", &colnum, &status);
    errorMsg = "Error in getting column number of Nz column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    nZArray = new float[nrows];
    for(i=0; i<nrows; i++){
        nZArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TFLOAT, colnum, 1, 1, nrows, NULL, nZArray, NULL, &status);
    errorMsg = "Error in reading Nz column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    
    //reading nXt
    fits_get_colnum(fptrAspect, CASEINSEN, "Nxt", &colnum, &status);
    errorMsg = "Error in getting column number of Nxt column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    nXtArray = new float[nrows];
    for(i=0; i<nrows; i++){
        nXtArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TFLOAT, colnum, 1, 1, nrows, NULL, nXtArray, NULL, &status);
    errorMsg = "Error in reading Nxt column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    
    //reading nYt
    fits_get_colnum(fptrAspect, CASEINSEN, "Nyt", &colnum, &status);
    errorMsg = "Error in getting column number of Nyt column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    nYtArray = new float[nrows];
    for(i=0; i<nrows; i++){
        nYtArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TFLOAT, colnum, 1, 1, nrows, NULL, nYtArray, NULL, &status);
    errorMsg = "Error in reading Nyt column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    
    //reading nZt
    fits_get_colnum(fptrAspect, CASEINSEN, "Nxt", &colnum, &status);
    errorMsg = "Error in getting column number of Nzt column in ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}    
    nZtArray = new float[nrows];
    for(i=0; i<nrows; i++){
        nZtArray[i] =0.0;
    }
    fits_read_col(fptrAspect, TFLOAT, colnum, 1, 1, nrows, NULL, nZtArray, NULL, &status);
    errorMsg = "Error in reading Nzt column of ASPECT extension of ASPECT file: " + aspectFilename;
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    //aspect file read and relevant information stored in object variables.
    
    fits_close_file(fptrAspect, &status);
    if (status) {
        LOG(ERROR) << "Error in closing Aspect file: " << aspectFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //clearing class vectors
    vecTime.clear();
    vecNx.clear();
    vecNy.clear();
    vecNz.clear();
    vecNxt.clear();
    vecNyt.clear();
    vecNzt.clear();
    // Assigning data read from aspect file to object variables
    for(i=0; i<nrows; i++){
        vecTime.push_back(timeArray[i]);
        vecNx.push_back(nXArray[i]);
        vecNy.push_back(nYArray[i]);
        vecNz.push_back(nZArray[i]);
        vecNxt.push_back(nXtArray[i]);
        vecNyt.push_back(nYtArray[i]);
        vecNzt.push_back(nZtArray[i]);
    }
    // Data assigned to object variables
    
    return status;
}

int AspectFileHandler::read_multiple_aspect_files(vector <string> aspectFilenames) {
    int status=0;
    long i=0,j=0;
    long nfiles = aspectFilenames.size(); //number of aspect files
    vector <double> temptime;
    vector <float> tempnx; //temporary vectors to store nx, ny, nz, nxt, nyt, nzt data
    vector <float> tempny; //from all aspect files.
    vector <float> tempnz; 
    vector <float> tempnxt; 
    vector <float> tempnyt; 
    vector <float> tempnzt; 
    
    for(i=0; i<nfiles; i++){
        if(read_aspect_file(aspectFilenames[i])){
            LOG(ERROR) << "Error in reading aspect file " << i << " : " << aspectFilenames[i];
            return EXIT_FAILURE;
        }
        DLOG(INFO) << "Aspect file " << aspectFilenames[i] << " with nrows="<< nrows << " read.";
        temptime.insert(temptime.end(), vecTime.begin(), vecTime.end());
        tempnx.insert(tempnx.end(), vecNx.begin(), vecNx.end());
        tempny.insert(tempny.end(), vecNy.begin(), vecNy.end());
        tempnz.insert(tempnz.end(), vecNz.begin(), vecNz.end());
        tempnxt.insert(tempnxt.end(), vecNxt.begin(), vecNxt.end());
        tempnyt.insert(tempnyt.end(), vecNyt.begin(), vecNyt.end());
        tempnzt.insert(tempnzt.end(), vecNzt.begin(), vecNzt.end());   

    }
    vecTime.clear();
    vecNx.clear();
    vecNy.clear();
    vecNz.clear();
    vecNxt.clear();
    vecNyt.clear();
    vecNzt.clear();
    vecTime.insert(vecTime.end(), temptime.begin(), temptime.end());
    vecNx.insert(vecNx.begin(), tempnx.begin(), tempnx.end());
    vecNy.insert(vecNy.begin(), tempny.begin(), tempny.end());
    vecNz.insert(vecNz.begin(), tempnz.begin(), tempnz.end());
    vecNxt.insert(vecNxt.begin(), tempnxt.begin(), tempnxt.end());
    vecNyt.insert(vecNyt.begin(), tempnyt.begin(), tempnyt.end());
    vecNzt.insert(vecNzt.begin(), tempnzt.begin(), tempnzt.end());

    return EXIT_SUCCESS;
}


int AspectFileHandler::get_avg_pointing_and_twist_vector(float& avg_nX, float& avg_nY, float& avg_nZ, float& avg_nXt, 
                                                            float& avg_nYt, float& avg_nZt, int& status){
    
    if(vecNx.size()==0 || vecNy.size()==0 || vecNz.size()==0 || vecNxt.size()==0 || vecNyt.size()==0 || vecNzt.size()==0){
        LOG(ERROR)<< "Vectors containing Nx, Ny, Nz, Nxt, Nyt & Nzt empty().";
        status = EXIT_FAILURE;
    }
    else{
        avg_nX = (accumulate(vecNx.begin(), vecNx.end(), 0.0)) / vecNx.size();
        avg_nY = (accumulate(vecNy.begin(), vecNy.end(), 0.0)) / vecNy.size();
        avg_nZ = (accumulate(vecNz.begin(), vecNz.end(), 0.0)) / vecNz.size();
        avg_nXt = (accumulate(vecNxt.begin(), vecNxt.end(), 0.0)) / vecNxt.size();
        avg_nYt = (accumulate(vecNyt.begin(), vecNyt.end(), 0.0)) / vecNyt.size();
        avg_nZt = (accumulate(vecNzt.begin(), vecNzt.end(), 0.0)) / vecNzt.size();
    }
    return status;
}

//EXPOSUREMAP CLASS
ExposuremapFileHandler::ExposuremapFileHandler() {

}

int ExposuremapFileHandler::read_exposure_map_file(string expMapFilename) {
    int status=0;
    int colnum=0;
    long nrows=0;
    fitsfile *fexpmap;
    
    fits_open_file(&fexpmap, (char*) expMapFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_movnam_hdu(fexpmap, BINARY_TBL, "CZTEXPOSURE", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to CZTEXPOSURE hdu in file : " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    } 
    
    //Resetting exposure table
    expTab.reset();
    
    if(read_fits_column(fexpmap, "DetX", TBYTE, 1, 1,-1, expTab.detX)){
        LOG(ERROR) << "Error in reading column DetX of " << expMapFilename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fexpmap, "DetY", TBYTE, 1, 1,-1, expTab.detY)){
        LOG(ERROR) << "Error in reading column DetY of " << expMapFilename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fexpmap, "PixX", TBYTE, 1, 1,-1, expTab.pixX)){
        LOG(ERROR) << "Error in reading column PixX of " << expMapFilename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fexpmap, "PixY", TBYTE, 1, 1,-1, expTab.pixY)){
        LOG(ERROR) << "Error in reading column PixY of " << expMapFilename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fexpmap, "X", TBYTE, 1, 1,-1, expTab.locX)){
        LOG(ERROR) << "Error in reading column X of " << expMapFilename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fexpmap, "Y", TBYTE, 1, 1,-1, expTab.locY)){
        LOG(ERROR) << "Error in reading column Y of " << expMapFilename;
        return EXIT_FAILURE;
    }
    if(read_fits_column(fexpmap, "Exposure", TFLOAT, 1, 1,-1, expTab.openfrac)){
        LOG(ERROR) << "Error in reading column Exposure of " << expMapFilename;
        return EXIT_FAILURE;
    }
    
    fits_close_file(fexpmap, &status);
    if(status){
        LOG(ERROR) << "Error in closing exposure map file :" << expMapFilename;
        return EXIT_FAILURE;
    }
    
    return status;
}

int ExposuremapFileHandler::write_exposure_map_file(string expMapFilename, string expMapTemplate){
    int status=0;
    int colnum=0;
    long nrows=0;
    fitsfile *fexpmap;
    
    nrows = expTab.get_nrows();
    
    if(create_empty_fitsfile(expMapFilename, expMapTemplate)){
        LOG(ERROR)<< "Error in creating exposure map file from corresponding template";
        return EXIT_FAILURE;
    }
    
    fits_open_file(&fexpmap, (char*) expMapFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) <<"Error in opening file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(fexpmap, BINARY_TBL, "CZTEXPOSURE", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to CZTEXPOSURE hdu in file : " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing DetX
    fits_get_colnum(fexpmap, CASEINSEN, "DetX", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for DetX column in CZTEXPOSURE extension of"
                << " file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_write_col(fexpmap,TBYTE,colnum,1,1,nrows,expTab.detX.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing data in DetX column in file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing DetY
    fits_get_colnum(fexpmap, CASEINSEN, "DetY", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for DetY column in CZTEXPOSURE extension of"
                << " file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_write_col(fexpmap,TBYTE,colnum,1,1,nrows,expTab.detY.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing data in DetY column in file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }   
    //Writing PixX
    fits_get_colnum(fexpmap, CASEINSEN, "PixX", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for PixX column in CZTEXPOSURE extension of"
                << " file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_write_col(fexpmap,TBYTE,colnum,1,1,nrows,expTab.pixX.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing data in PixX column in file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Writing PixY
    fits_get_colnum(fexpmap, CASEINSEN, "PixY", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for PixY column in CZTEXPOSURE extension of"
                << " file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_write_col(fexpmap, TBYTE, colnum, 1, 1, nrows, expTab.pixY.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing data in PixY column in file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Writing Exposure
    fits_get_colnum(fexpmap, CASEINSEN, "Exposure", &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for Exposure column in CZTEXPOSURE extension of"
                << " file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_write_col(fexpmap, TFLOAT, colnum, 1, 1, nrows, expTab.openfrac.data(), &status);
    if (status) {
        LOG(ERROR) << "Error in writing data in Exposure column in file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    if(write_fits_column(fexpmap, "x", TBYTE, 1, 1, expTab.locX)){
        LOG(ERROR) << "Error in writing column x of exposure file " << expMapFilename;
        return EXIT_FAILURE;
    }
    
    if(write_fits_column(fexpmap, "y", TBYTE, 1, 1, expTab.locY)) {
        LOG(ERROR) << "Error in writing column x of exposure file " << expMapFilename;
        return EXIT_FAILURE;
    }
    fits_close_file(fexpmap, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file: " << expMapFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}

int ExposuremapFileHandler::read_exposure_array_file(string expArrFilename) {
    int status = 0;
    int colnum = 0;
    long nrows = 0;
    char keyValue[PIL_LINESIZE];
    fitsfile *fexpmap;

    fits_open_file(&fexpmap, (char*) expArrFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file: " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Resetting exposure table
    expTab.reset();
    
    //READING EXPOSURE ARRAY EXTENSION
    fits_movnam_hdu(fexpmap, BINARY_TBL, "EXPOSURE_ARRAY", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to CZTEXPOSURE hdu in file : " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    if (read_fits_column(fexpmap, "DetX", TBYTE, 1, 1, -1, expTab.detX)) {
        LOG(ERROR) << "Error in reading column DetX of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "DetY", TBYTE, 1, 1, -1, expTab.detY)) {
        LOG(ERROR) << "Error in reading column DetY of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "PixX", TBYTE, 1, 1, -1, expTab.pixX)) {
        LOG(ERROR) << "Error in reading column PixX of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "PixY", TBYTE, 1, 1, -1, expTab.pixY)) {
        LOG(ERROR) << "Error in reading column PixY of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "X", TBYTE, 1, 1, -1, expTab.locX)) {
        LOG(ERROR) << "Error in reading column X of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "Y", TBYTE, 1, 1, -1, expTab.locY)) {
        LOG(ERROR) << "Error in reading column Y of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_array_column(fexpmap, "EXPOSURE", TFLOAT, 1, 1, -1, expTab.openfracArray)) {
        LOG(ERROR) << "Error in reading column Exposure of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_array_column(fexpmap, "WEIGHTS", TFLOAT, 1, 1, -1, expTab.weightsArray)) {
        LOG(ERROR) << "Error in reading column Weights of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "BADPIXFLAG", TBYTE, 1, 1, -1, expTab.badpixFlag)) {
        LOG(ERROR) << "Error in reading column BADPIXFLAG of " << expArrFilename;
        return EXIT_FAILURE;
    }    

    //READING RENORM_OFFSETS EXTENSION
    fits_movnam_hdu(fexpmap, BINARY_TBL, "RENORM_OFFSETS", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to RENORM_OFFSETS hdu in file : " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    if (read_fits_column(fexpmap, "Energy", TFLOAT, 1, 1, -1, expTab.energies)) {
        LOG(ERROR) << "Error in reading column Energy of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "Channel", TUSHORT, 1, 1, -1, expTab.PIs)) {
        LOG(ERROR) << "Error in reading column Channel of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "Renormalization_offset", TFLOAT, 1, 1, -1, expTab.renormOffsetArray)) {
        LOG(ERROR) << "Error in reading column Renormalization_offset of " << expArrFilename;
        return EXIT_FAILURE;
    }
    if (read_fits_column(fexpmap, "Area", TFLOAT, 1, 1, -1, expTab.areaArray)) {
        LOG(ERROR) << "Error in reading column Area of " << expArrFilename;
        return EXIT_FAILURE;
    }
    
    //READING NECESSARY KEYWORDS
    fits_movnam_hdu(fexpmap, BINARY_TBL, "EXPOSURE_ARRAY", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to CZTEXPOSURE hdu in file : " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    try{
        
        readKey(fexpmap, TFLOAT, "THETAX", &expTab.thetaxd);
        readKey(fexpmap, TFLOAT, "THETAY", &expTab.thetayd);
        readKey(fexpmap, TINT, "NCHANNEL", &expTab.nchannels);
        readKey(fexpmap, TINT, "EFLAG", &expTab.eboundFlag);
        readKey(fexpmap, TSTRING, "ERANGE", &keyValue);
        expTab.energyRange = (string) keyValue;
    } catch (ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }
    
    fits_close_file(fexpmap, &status);
    if (status) {
        LOG(ERROR) << "Error in closing exposure map file :" << expArrFilename;
        return EXIT_FAILURE;
    }

    return status;
}

int ExposuremapFileHandler::write_exposure_array_file(string expArrFilename, string expArrTemplate){
    int status = 0;
    int colnum = 0;
    long nrows = 0;
    fitsfile *fexpmap;
    vector <float> tempofArray;
    vector <string> keys;
    vector <float> valuesF;
    vector <string> valuesS;

    nrows = expTab.get_nrows();

    if (create_empty_fitsfile(expArrFilename, expArrTemplate)) {
        LOG(ERROR) << "Error in creating exposure array file from corresponding template";
        return EXIT_FAILURE;
    }

    fits_open_file(&fexpmap, (char*) expArrFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file: " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_movnam_hdu(fexpmap, BINARY_TBL, "EXPOSURE_ARRAY", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to EXPOSURE_ARRAY in file : " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Writing DetX
    if(write_fits_column(fexpmap, "DetX", TBYTE, 1, 1, expTab.detX)) {
        LOG(ERROR) << "Error in writing column DetX of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing DetY
    if(write_fits_column(fexpmap, "DetY", TBYTE, 1, 1, expTab.detY)) {
        LOG(ERROR) << "Error in writing column DetY of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing PixX
    if(write_fits_column(fexpmap, "PixX", TBYTE, 1, 1, expTab.pixX)) {
        LOG(ERROR) << "Error in writing column PixX of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing PixY
    if(write_fits_column(fexpmap, "PixY", TBYTE, 1, 1, expTab.pixY)) {
        LOG(ERROR) << "Error in writing column PixY of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing x
    if (write_fits_column(fexpmap, "x", TBYTE, 1, 1, expTab.locX)) {
        LOG(ERROR) << "Error in writing column x of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing y
    if (write_fits_column(fexpmap, "y", TBYTE, 1, 1, expTab.locY)) {
        LOG(ERROR) << "Error in writing column x of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing Exposure
    if(write_fits_array_column(fexpmap, "EXPOSURE", TFLOAT, 1, 1, expTab.openfracArray)){
        LOG(ERROR) << "Error in writing column EXPOSURE of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Writing Weights
    if(write_fits_array_column(fexpmap, "WEIGHTS", TFLOAT, 1, 1, expTab.weightsArray)){
        LOG(ERROR) << "Error in writing column EXPOSURE of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    
    //Writing Flags
    if(!(expTab.badpixFlag.empty())){
        if(write_fits_column(fexpmap, "badpixFlag", TBYTE, 1, 1, expTab.badpixFlag)){
            LOG(ERROR) << "Error in writing column badpixFlag of exposure file " << expArrFilename;
            return EXIT_FAILURE;
        }
    }
    //Updating keys
    string keysArray[] = {"THETAX", "THETAY", "NCHANNEL"};
    float valuesFArray[] = {expTab.thetaxd, expTab.thetayd, expTab.nchannels};
    keys.assign(keysArray, keysArray + 3);
    valuesF.assign(valuesFArray, valuesFArray + 3);
    if(updateKeys(fexpmap, "EXPOSURE_ARRAY", keys, valuesF, TFLOAT)){
        LOG(ERROR) << "Error in updating keys.";
        return EXIT_FAILURE;
    }
    
    fits_update_key(fexpmap, TINT, "EFLAG", &(expTab.eboundFlag), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in updating EFLAG keyword";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_update_key_str(fexpmap, "ERANGE", (char*)(expTab.energyRange).c_str(), NULL, &status);
    if(status){
        LOG(WARNING) << "Unable to update key ERANGE";
    }
    //keys updated

    //Moving to HDU Renorm_offsets
    fits_movnam_hdu(fexpmap, BINARY_TBL, "RENORM_OFFSETS", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to RENORM_OFFSETS in file : " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    //Energy 
    if(write_fits_column(fexpmap, "Energy", TFLOAT, 1, 1, expTab.energies)){
        LOG(ERROR) << "Error in writing Energy column of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Channel
    if (write_fits_column(fexpmap, "Channel", TUSHORT, 1, 1, expTab.PIs)) {
        LOG(ERROR) << "Error in writing Channel column of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Renormalization offsets
    if(write_fits_column(fexpmap, "Renormalization_offset", TFLOAT, 1, 1, expTab.renormOffsetArray)){
        LOG(ERROR) << "Error in writing Renormalization_offset column of exposure file " << expArrFilename;
        return EXIT_FAILURE;
    }
    //Areas
    if(write_fits_column(fexpmap, "Area", TFLOAT, 1, 1, expTab.areaArray)) {
        LOG(ERROR) << "Error in writing Area column of exposure file " << expArrFilename;
        return EXIT_FAILURE;        
    }
    fits_close_file(fexpmap, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file: " << expArrFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    return EXIT_SUCCESS;    
}
int ExposuremapFileHandler::set_exposuretable(ExposureTable expTable) {
    this->expTab = expTable;
    return EXIT_SUCCESS;
}


//EXPOSURE TABLE STRUCT

ExposureTable::ExposureTable() {
    eboundFlag = false;
}

int ExposureTable::set_exposure(unsigned char detx, unsigned char dety, unsigned char pixx, 
        unsigned char pixy, float openfraction, float weight){
    unsigned char x,y;
    detX.push_back(detx);
    detY.push_back(dety);
    pixX.push_back(pixx);
    pixY.push_back(pixy);
    openfrac.push_back(openfraction);
    
    x = detx*NO_PIX_X_PER_DET + pixx;
    y = dety*NO_PIX_Y_PER_DET + pixy;
    
    locX.push_back(x);
    locY.push_back(y);
    weights.push_back(weight);
    return EXIT_SUCCESS;
}

int ExposureTable::generate_full_exposure_image(){
    int status=0;
    long i=0, ix=0, iy=0;
    long nrows=get_nrows();
    vector <float> tempExposure;
    //Initializing vectors
    fullExposure.clear();
    for(ix=0; ix<XSIZE; ix++){
        tempExposure.push_back(0.0);
    }
    for(iy=0; iy<YSIZE; iy++){
        fullExposure.push_back(tempExposure);
    }
    //vectors initialized
    for (i = 0; i < nrows; i++) {
        //Calculating pixel co-ordinates for each open fraction table record.
        ix = detX[i] * PIXELS_PER_ROW_DET + pixX[i];
        iy = detY[i] * PIXELS_PER_COL_DET + pixY[i];

        //Full Open fraction 2D vector 
        fullExposure[iy][ix] = openfrac[i];
    }  
    
    return EXIT_SUCCESS;
}

int ExposureTable::generate_index_table(bool regenerateFlag){
    int status=0;
    ErrorHandler errHandler;
    long irecord=0;
    vector <long> tempIndexX;
    
    //Checking if indexTable already available
    if (regenerateFlag == false) {
        if (indexTable.size() == YSIZE) {
            if (indexTable[0].size() == XSIZE) {
                LOG(INFO) << "Index Table is already available.";
                return EXIT_SUCCESS;
            }
        }
    }
    
    //Checking for errors in exposure table
    if(locX.size() != locY.size()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Size of locX and locY in exposure table is not equal";
        throw errHandler;
    } else if(locX.size()!=NO_PIX_ALL_QUADS){
        errHandler.severity = errERROR;
        errHandler.errorStatus = VECTOR_2D_GENERATION_ERROR;
        errHandler.errorMsg = "Cannot generate index table as number of records is less than " + itoa(NO_PIX_ALL_QUADS);
        throw errHandler;
    }
    
    indexTable.clear();
    //Resizing 2D vector to store 128x128 index values
    tempIndexX.resize(XSIZE, 0);
    indexTable.resize(YSIZE, tempIndexX);
    for(irecord=0; irecord<NO_PIX_ALL_QUADS; irecord++){
        indexTable[locY[irecord]][locX[irecord]]=irecord;
    }
    
    return EXIT_SUCCESS;
}

int ExposureTable::reset() {
    thetaxd=-99.9999;
    thetayd=-99.9999;
    detX.clear();
    detY.clear();
    pixX.clear();
    pixY.clear();
    locX.clear();
    locY.clear();
    openfrac.clear();
    fullExposure.clear();
    weights.clear();
    openfracArray.clear();
    weightsArray.clear();
    PIs.clear();
    energies.clear();
    badpixFlag.clear();
    renormOffsetArray.clear();
    areaArray.clear();
    indexTable.clear();
    eboundFlag=false;
}

float ExposureTable::get_weight(unsigned char locx, unsigned char locy, unsigned short PIchannel){
    float weight=0.0;
    int numrows = XSIZE*YSIZE;
    
}

// END EXPOSURE TABLE

//DPI CLASS

DPI::DPI() {
    //Initializing certain DPI parameters
    eMin=0.0;
    eMax=0.0;
    tStart=0.0;
    tStop=0.0;
    totalDPHcount = 0;
    totalDPIcount =0.0;
    badpixFile = "";
    badpixFlag = false;
    badpixThreshold = GOODPIX;
    //Initializing full dph and dpi 2d
    vector <long> tempFullDPH;
    vector <float> tempFullDPI;
    tempFullDPH.resize(XSIZE, 0.0);
    tempFullDPI.resize(XSIZE, 0.0);
    fullDPH2D.resize(YSIZE, tempFullDPH);
    fullDPI2D.resize(YSIZE, tempFullDPI);
}

vector<vector<float> > DPI::get_dpi_image(string dpiFilename, string extname) {
    int status = 0;
    fitsfile *fdpih;
    vector <vector <float> > dpi;
    ErrorHandler errHandler;
    int dataType = 0;
    fits_open_file(&fdpih, (char*) dpiFilename.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening DPI/DPH file.";
        throw errHandler;
    }
    
    if(read_img(dpi, fdpih, extname, TFLOAT)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_READING_FITS_IMAGE ;
        errHandler.errorMsg = "Error in reading dpi for hduname " + extname;
        throw errHandler;
    }
    
    fits_close_file(fdpih, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in closing fits file " + dpiFilename;
        throw errHandler;
    }
    
    return dpi;
    
}

vector<vector<long> > DPI::get_dph_image(string dphFilename, string extname) {
    int status = 0;
    fitsfile *fdpih;
    vector <vector <long> > dph;
    ErrorHandler errHandler;
    int dataType = 0;
    fits_open_file(&fdpih, (char*) dphFilename.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening DPI/DPH file.";
        throw errHandler;
    }
    
    if(read_img(dph, fdpih, extname, TLONG)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_READING_FITS_IMAGE ;
        errHandler.errorMsg = "Error in reading dpi for hduname " + extname;
        throw errHandler;
    }
    
    fits_close_file(fdpih, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in closing fits file " + dphFilename;
        throw errHandler;
    }
    
    return dph;
    
}


int DPI::add_quad_to_fulldpih(int quadNo, string fileType){
    int status=0;
    if (fileType == "DPI") {
        if (rearrange_quads(&fullDPI2D, dpi2D, quadNo)) {
            LOG(ERROR) << "Error in properly placing quadDPI into fullDPI";
            return EXIT_FAILURE;
        }
    } else if (fileType == "DPH") {
        if (rearrange_quads(&fullDPH2D, dph2D, quadNo)) {
            LOG(ERROR) << "Error in properly placing quadDPH into fullDPH";
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}

int DPI::create_dpi2d(string effAreaFilename, int quad_no){
    int status=0;
    int ix, iy=0;
    long i=0;
    long nrows=0;
    float totcount=0.0;
    EffArea effectiveArea;
    char* HDUname= new char[3];
    //string effAreaFilepath="";
    vector <vector <float> > normEffArea;
    vector <float> tempDPI;
    normEffArea.clear();
    //effAreaFilepath = caldb_full_path_generator(effAreaFilename);
    quadToHDU(quad_no, HDUname);

    LOG(INFO) << "Calculating normalized effective area for Quadrant " << quad_no;
    if(effectiveArea.read_effarea_file(effAreaFilename, (string) HDUname)){
        LOG(ERROR)<<"Error in reading Effective Area File : " << effAreaFilename;
    }
    if(effectiveArea.calculate_normalized_effarea(eMin, eMax)){
        LOG(ERROR)<<"Error in calculating normalized effective area.";
    }

    dpi2D.clear();
    dpi2D.clear();
    tempDPI.resize(XPIX_QUAD, 0);
    dpi2D.resize(YPIX_QUAD, tempDPI);
    
    normEffArea = effectiveArea.get_normalized_effarea();
    for(ix=0; ix<XPIX_QUAD; ix++){
        for(iy=0; iy<YPIX_QUAD; iy++){
            dpi2D[ix][iy]=dph2D[ix][iy]/normEffArea[ix][iy];
            totcount+=dpi2D[ix][iy];
        }
    }
    LOG(INFO) << setprecision(20) << totcount;
    totalDPIcount = totcount;
    
    return EXIT_SUCCESS;
}

int DPI::create_dph2d(string eventFilename, int quad_no, float emin, float emax, 
        double tstart, double tstop,int timefilter,int energyfilter){
    int status=0;
    int ix, iy=0;
    long i=0;
    long nrows=0;
    EventFileHandler evt;
    ErrorHandler errHandler;
    Badpix badpix; //to read badpix file
    BadpixTable bptable; //badpix table for quadrant defined by quad_no
    vector < vector <unsigned char> > quadBadpixMap;
    long totcount=0;
    vector <long> tempDPH;
    fitsfile *fevt; //event file pointer
    
    eMin = emin;
    eMax = emax;
    tStart = tstart;
    tStop = tstop;
    
    //Validating input
    if(quad_no<0 || quad_no>3){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Quadrant ID should be in the range 0-3";
        throw errHandler;
    }
    //Opening event file
    fits_open_file(&fevt, (char*) eventFilename.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening input event file: "+ eventFilename;
        throw errHandler;
    }
    
    if(evt.read_quad_extension(fevt, quad_no)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading event data for Quadrant " + itoa(quad_no) + " from event file " + eventFilename;
        throw errHandler;
        return EXIT_FAILURE;
    }
    
    //Reading badpixel file if File exists 
    if(badpixFile!=""){
        if(!FileExists((char*) badpixFile.c_str())){
			LOG(ERROR) << "Error in opening badpix file: " << badpixFile;
        	return (EXIT_FAILURE);

        }
    }
    if(badpixFile!=""){
        if(badpix.read_badpix_file(badpixFile)){
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in reading badpixel file " + badpixFile;
            throw errHandler;
        }
        bptable = badpix.get_badpix_table(quad_no);   
        bptable.create_quadrant_badPixMap();
    }
    
    dph2D.clear();
    tempDPH.resize(XPIX_QUAD, 0);
    dph2D.resize(YPIX_QUAD, tempDPH);

    
    vecDetX = evt.get_vecDetXQuad();
    vecDetY = evt.get_vecDetYQuad();
    vecEnergy = evt.get_vecEnergyQuad();
    vecTime = evt.get_vecTimeQuad();
    vecDetID = evt.get_vecDetIDQuad();
    vecPixID = evt.get_vecPixIDQuad();
    
    nrows = vecTime.size();

    //Counting events to make DPH
    for(i=0; i<nrows; i++){
        if(vecTime[i]>=tStart && vecTime[i]<=tStop || timefilter==0){
            if(vecEnergy[i]>=eMin && vecEnergy[i]<=eMax ||energyfilter==0){
                dph2D[vecDetY[i]][vecDetX[i]] +=1;
                totcount=totcount+1;
            }
        }
    }
    
    // Setting bad pixels to 0 [a pixel is defined as bad if its flag is 
    // greater than user defined threshold value.]
    for(iy=0; iy<YPIX_QUAD; iy++){
        for(ix=0; ix<XPIX_QUAD; ix++){
            if(bptable.vecBadpixMap[iy][ix] > badpixThreshold){
                dph2D[iy][ix]=0;
            }
        }
    }
    
    totalDPHcount=totcount;
    LOG(INFO)<<"DPH evaluated for Quadrant " << quad_no << " with following parameters : ";
    LOG(INFO) << "MINIMUM TIME:" << setprecision(20) << tstart;
    LOG(INFO) << "MAXIMUM TIME:" << setprecision(20) << tstop;
    LOG(INFO) << "MINIMUM ENERGY:" << emin;
    LOG(INFO) << "MAXIMUM ENERGY:" << emax;    
    LOG(INFO)<<"TOTAL COUNTS:"<<totcount;

    fits_close_file(fevt, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in closing event file "+ eventFilename;
        throw errHandler;
    }
    
    return EXIT_SUCCESS;
    
}

int DPI::create_dpi_file(string outDPIfilename, string DPItemplate, string inEvtFilename){
    int status=0;
    fitsfile *fevt, *fdpi;
    
    vector <string> vecKeywords;
    if(create_empty_fitsfile(outDPIfilename, DPItemplate)){
        LOG(ERROR)<<"Error in creating empty DPI file " << outDPIfilename;
        return (EXIT_FAILURE);
    }
    
    fits_open_file(&fevt, (char*) inEvtFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening input event file";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_open_file(&fdpi, (char*) outDPIfilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening output DPI file";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Keywords to be copied
    vecKeywords.push_back("DATE-OBS");
    vecKeywords.push_back("TIME-OBS");
    vecKeywords.push_back("DATE-END");
    vecKeywords.push_back("TIME-END");
    vecKeywords.push_back("TIMESYS");
    vecKeywords.push_back("MJDREFI");
    vecKeywords.push_back("MJDREFF");
    vecKeywords.push_back("TIMEZERI");
    vecKeywords.push_back("TIMEZERF");
    vecKeywords.push_back("TIMEDEL");
    vecKeywords.push_back("TIMEUNIT");
    vecKeywords.push_back("TSTARTI");
    vecKeywords.push_back("TSTARTF");
    vecKeywords.push_back("TSTOPI");
    vecKeywords.push_back("TSTOPF");
    vecKeywords.push_back("TSTART");
    vecKeywords.push_back("TSTOP");
    vecKeywords.push_back("TELAPSE");
				
    //add GTI_TYPE 
    vecKeywords.push_back("GTITYPE");
    vecKeywords.push_back("OBJECT");
    vecKeywords.push_back("RA_PNT");
    vecKeywords.push_back("DEC_PNT");
    vecKeywords.push_back("EQUINOX");
    vecKeywords.push_back("RADECSYS");
    vecKeywords.push_back("OBSERVER");
    vecKeywords.push_back("OBS_ID");
    vecKeywords.push_back("SOURCEID");
    vecKeywords.push_back("TARGETID");
    vecKeywords.push_back("EXPOSURE");
    vecKeywords.push_back("ORB_NUM");
 
    
    if (copyUserKeyWords(fevt, fdpi, "Primary", "Primary", vecKeywords)) { //copying keywords from input file
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;
    }

    int i;
    //Copy GTITYPE keywords to all quad extensions
    vecKeywords.clear();
    vecKeywords.push_back("GTITYPE");

    LOG(INFO)<<"Copying heder keyowrds";
    for(i=0;i<4;i++)
    {
        char dpiextname[20];
        sprintf(dpiextname,"Q%d",i);
        if (copyUserKeyWords(fevt, fdpi, dpiextname, dpiextname, vecKeywords)) {
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;        
        }
    }    
    //Copy keywords to full_dph extension
    if (copyUserKeyWords(fevt, fdpi, "Primary", "FULL_DPH", vecKeywords)) { 
        LOG(ERROR) << "***Error in copying keywords***";
        return EXIT_FAILURE;
    }

    
    fits_close_file(fevt, &status);
    if (status) {
        LOG(ERROR) << "Error in closing event file "<< inEvtFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_close_file(fdpi, &status);
    if(status){
        LOG(ERROR) << "Error in closing DPI/DPH file " << outDPIfilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);        
    }
    


    return (EXIT_SUCCESS);
}

int DPI::write_dpih_quad(string DPIfilename, string filetype, string extname) {
    int status=0;
    fitsfile *fdpih;
    vector<long> dph;
    vector<float> dpi;
    long fpixel[2];
    int bitpix=LONG_IMG;
    long naxis1=64;
    long naxis2=64;

    
    fpixel[0]=fpixel[1]=1;
    
    fits_open_file(&fdpih, (char*) DPIfilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening output dph file: " << DPIfilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }   

    fits_movnam_hdu(fdpih, IMAGE_HDU, (char *) extname.c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to " << extname << " HDU of file " << DPIfilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_update_key(fdpih, TFLOAT, "EMIN", &eMin, NULL, &status);
    fits_update_key(fdpih, TFLOAT, "EMAX", &eMax, NULL, &status);
    fits_update_key(fdpih, TDOUBLE, "TSTART", &tStart, NULL, &status);
    fits_update_key(fdpih, TDOUBLE, "TSTOP", &tStop, NULL, &status);
    headerParam.writeToHeader(fdpih); 
    if(filetype=="DPH"){
        fits_update_key(fdpih, TLONG, "DPHCOUNT", &totalDPHcount, NULL, &status);
    }
    else if(filetype=="DPI"){
        fits_update_key(fdpih, TFLOAT, "DPICOUNT", &totalDPIcount, NULL, &status);
    }
   headerParam.writeToHeader(fdpih); 
    if(filetype=="DPH"){
        if(linearize_2Dvector(dph2D, &dph)){
            LOG(ERROR)<<"Error in linearizing 2D vector for DPH";
            return EXIT_FAILURE;
        }
        fits_write_pix(fdpih, TLONG, fpixel, dph.size(), dph.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DPH to file : " << DPIfilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    else if(filetype=="DPI") {
        if (linearize_2Dvector(dpi2D, &dpi)) {
            LOG(ERROR) << "Error in linearizing 2D vector for DPI";
            return EXIT_FAILURE;
        }
        fits_write_pix(fdpih, TFLOAT, fpixel, dpi.size(), dpi.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing DPI to file : " << DPIfilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    

    
    fits_close_file(fdpih, &status);
    if (status) {
        LOG(ERROR) << "Error in closing DPI";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    
    
    return EXIT_SUCCESS;
}

int DPI::write_dpih_full(string DPIfilename, string filetype, string extname) {
    int status = 0;
    fitsfile* fdpih;
    vector<long> fulldph;
    vector<float> fulldpi;
    
    long fpixel[2];
    int bitpix = LONG_IMG;

    fpixel[0] = fpixel[1] = 1;

    fits_open_file(&fdpih, (char*) DPIfilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening output dph file: " << DPIfilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    } 

    fits_movnam_hdu(fdpih, IMAGE_HDU, (char *) extname.c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to " << extname << " HDU of file " << DPIfilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_update_key(fdpih, TFLOAT, "EMIN", &eMin, NULL, &status);
    fits_update_key(fdpih, TFLOAT, "EMAX", &eMax, NULL, &status);
    fits_update_key(fdpih, TDOUBLE, "TSTART", &tStart, NULL, &status);
    fits_update_key(fdpih, TDOUBLE, "TSTOP", &tStop, NULL, &status);
    
    if (filetype == "DPH") {
        if (linearize_2Dvector(fullDPH2D, &fulldph)) {
            LOG(ERROR) << "Error in linearizing 2D vector for DPH";
            return EXIT_FAILURE;
        }
        fits_write_pix(fdpih, TLONG, fpixel, fulldph.size(), fulldph.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing FULL_DPH to file : " << DPIfilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    } else if (filetype == "DPI") {
        if (linearize_2Dvector(fullDPI2D, &fulldpi)) {
            LOG(ERROR) << "Error in linearizing 2D vector for DPI";
            return EXIT_FAILURE;
        }
        fits_write_pix(fdpih, TFLOAT, fpixel, fulldpi.size(), fulldpi.data(), &status);
        if (status) {
            LOG(ERROR) << "Error in writing FULL_DPI to file : " << DPIfilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }



    fits_close_file(fdpih, &status);
    if (status) {
        LOG(ERROR) << "Error in closing DPI";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    } 

    return EXIT_SUCCESS;
}

//END DPI CLASS

//SHADOW TABLE STRUCT
int ShadowTable::reset(){
    detX.clear();
    detY.clear();
    pixX.clear();
    pixY.clear();
    locx.clear();
    locy.clear();
    shadow.clear();
    shadow2D.clear();
    thetaXr=0.0;
    thetaYr=0.0;
    energyMin=0.0;
    energyMax=0.0;
    RA=0.0;
    DEC=0.0;
    nBins=1;
}

int ShadowTable::set_shadow(int detx, int dety, int pixx, int pixy, int x, int y, float shadow){
    detX.push_back(detx);
    detY.push_back(dety);
    pixX.push_back(pixx);
    pixY.push_back(pixy);
    locx.push_back(x);
    locy.push_back(y);
    this->shadow.push_back(shadow);
    
    return EXIT_SUCCESS;
}

int ShadowTable::is_empty() {
    int flag=FALSE;
    if(detX.size()==0){
        LOG(WARNING) << "DETX IS EMPTY.";
        flag &= TRUE;
    }
    else if(detY.size()==0){
        LOG(WARNING) << "DETY IS EMPTY.";
        flag &= TRUE;
    }
    else if(pixX.size()==0){
        LOG(WARNING) << "PIXX IS EMPTY.";
        flag &= TRUE;
    }
    else if(pixY.size()==0){
        LOG(WARNING) << "PIXY IS EMPTY.";
        flag &= TRUE;
    }
    else if(locx.size()==0){
        LOG(WARNING) << "LOCX IS EMPTY.";
        flag &= TRUE;
    }
    else if(locy.size()==0){
        LOG(WARNING) << "LOCY IS EMPTY.";
        flag &= TRUE;
    }
    else if(shadow.size()==0){
        LOG(WARNING) << "SHADOW VECTOR IS EMPTY.";
        flag &= TRUE;
                
    }

    return flag;
}
//SHADOW TABLE STRUCT END
//SHADOW FILE HANDLER

int ShadowFileHandler::create_shadow_file(string shadowFilename, string shadowTemplate, vector<string> extNames) {
    int status = 0;
    fitsfile *fshadow, *fshadowcpy;
    int nhdus=0; //number of hdus.
    int hdunum=0;
    int i=0;

    //number of hdus to be created in empty shadow file
    nhdus=extNames.size();
    
    if (create_empty_fitsfile(shadowFilename, shadowTemplate)) {
        LOG(ERROR) << "Error in creating empty Shadow file " << shadowFilename;
        return (EXIT_FAILURE);
    }

    if (nhdus > 1) {
        fits_open_file(&fshadow, (char*) shadowFilename.c_str(), READWRITE, &status);
        if (status) {
            LOG(ERROR) << "Error in opening shadow file" << shadowFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_open_file(&fshadowcpy, (char*) shadowFilename.c_str(), READWRITE, &status);
        if (status) {
            LOG(ERROR) << "Error in opening shadow file" << shadowFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_movabs_hdu(fshadow, 2, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU number " << 2 << " of output shadow file: " << shadowFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        
        //copying shadow HDU
        for (i = 0; i < (nhdus-1); i++) {
            fits_copy_hdu(fshadow, fshadowcpy, NULL, &status);
            if (status) {
                LOG(ERROR) << "Error in copying SHADOW HDU to create multiple shadow extensions" <<
                        " in output shadow file: " << shadowFilename;
                fits_report_error(stderr, status);
                return (EXIT_FAILURE);
            }
        }
        
        //renaming shadow extension
        for(i=0; i<nhdus; i++) {
            hdunum=2+i;
            fits_movabs_hdu(fshadowcpy, hdunum, NULL, &status);
            if (status) {
                LOG(ERROR) << "Error in moving to HDU number " << 2 << " of output shadow file: " << shadowFilename;
                fits_report_error(stderr, status);
                return (EXIT_FAILURE);
            }
            //Updating Extension names for SHADOW extensions
            fits_update_key(fshadowcpy, TSTRING, "EXTNAME", (char*) extNames[i].c_str(), NULL, &status);
            if (status) {
                LOG(ERROR) << "Error in updating extension name for shadow extension " << hdunum;
                fits_report_error(stderr, status);
                return (EXIT_FAILURE);
            }
        }

        fits_close_file(fshadow, &status);
        if(status) {
            LOG(ERROR) << "Error in closing output shadow file" << shadowFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_close_file(fshadowcpy, &status);
        if(status) {
            LOG(ERROR) << "Error in closing output shadow file" << shadowFilename;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

    }

    return (EXIT_SUCCESS);
}


int ShadowFileHandler::write_shadow_file(string shadowFilename, string extName) {
    int status=0;
    fitsfile *fshadow;
    
    fits_open_file(&fshadow, (char*) shadowFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening shadow file " << shadowFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(fshadow, BINARY_TBL, (char*)extName.c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << extName << " in shadow file " << shadowFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    if(shadowtab.is_empty()){
        LOG(ERROR) << "Improper shadow data.";
        return EXIT_FAILURE;
    } 
    else {
        if(write_fits_column(fshadow, "DETX", TINT, 1, 1, shadowtab.detX)){
            LOG(ERROR) << "Error in writing DetX column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
        if(write_fits_column(fshadow, "DETY", TINT, 1, 1, shadowtab.detY)){
            LOG(ERROR) << "Error in writing DETY column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
        if(write_fits_column(fshadow, "PIXX", TINT, 1, 1, shadowtab.pixX)){
            LOG(ERROR) << "Error in writing PIXX column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
        if(write_fits_column(fshadow, "PIXY", TINT, 1, 1, shadowtab.detX)){
            LOG(ERROR) << "Error in writing PIXY column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
        if(write_fits_column(fshadow, "x", TINT, 1, 1, shadowtab.locx)){
            LOG(ERROR) << "Error in writing x column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
        if(write_fits_column(fshadow, "y", TINT, 1, 1, shadowtab.locy)){
            LOG(ERROR) << "Error in writing y column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
        if(write_fits_column(fshadow, "SHADOW", TFLOAT, 1, 1, shadowtab.shadow)){
            LOG(ERROR) << "Error in writing DetX column of HDU " << extName << " in shadow file " << shadowFilename;
            return EXIT_FAILURE;
        }
    }
    
    fits_close_file(fshadow, &status);
    if (status) {
        LOG(ERROR) << "Error in closing shadow file: " << shadowFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    
    
    return EXIT_SUCCESS;
}

int ShadowFileHandler::read_shadow_file(string shadowFilename, string extName) {
    int status=0;
    fitsfile *fshadow;
    ErrorHandler errHandler;
    fits_open_file(&fshadow, (char*) shadowFilename.c_str(), READWRITE, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening shadow file " + shadowFilename;
        throw errHandler;
    }

    fits_movnam_hdu(fshadow, BINARY_TBL, (char*) extName.c_str(), NULL, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg ="Error in moving to HDU " + extName + " in shadow file " + shadowFilename;
        throw errHandler;
    }

    try{
        read_fits_columnN(fshadow, "DETX", TINT, 1, 1, -1, shadowtab.detX);
        read_fits_columnN(fshadow, "DETY", TINT, 1, 1, -1, shadowtab.detY);
        read_fits_columnN(fshadow, "PIXX", TINT, 1, 1, -1, shadowtab.pixX);
        read_fits_columnN(fshadow, "PIXY", TINT, 1, 1, -1, shadowtab.pixY);
        read_fits_columnN(fshadow, "x", TINT, 1, 1, -1, shadowtab.locx);
        read_fits_columnN(fshadow, "y", TINT, 1, 1, -1, shadowtab.locy);
        read_fits_columnN(fshadow, "SHADOW", TFLOAT, 1, 1, -1, shadowtab.shadow);
    } catch(ErrorHandler errHandler){
        throw errHandler;
    }
    
    //Reading header keywords 
    try {
        readFitsKey(fshadow, TFLOAT, "THETAX", &shadowtab.thetaXr, NULL);
        readFitsKey(fshadow, TFLOAT, "THETAY", &shadowtab.thetaYr, NULL);
        readFitsKey(fshadow, TFLOAT, "FLUX", &shadowtab.sourceFlux, NULL);
        readFitsKey(fshadow, TFLOAT, "ESTART", &shadowtab.energyMin, NULL);
        readFitsKey(fshadow, TFLOAT, "ESTOP", &shadowtab.energyMax, NULL);
        readFitsKey(fshadow, TINT, "NBINS", &shadowtab.nBins, NULL);
        readFitsKey(fshadow, TFLOAT, "RA", &shadowtab.RA, NULL);
        readFitsKey(fshadow, TFLOAT, "DEC", &shadowtab.DEC, NULL);
        headerReadFlag = true;
    } catch (ErrorHandler errHandler) {
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    
    fits_close_file(fshadow, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in clsoing shadow file: " + shadowFilename;
        throw errHandler;
    }
    return status;
}

vector < vector <float> > ShadowFileHandler::get_shadow(string quadname){
    int status=0;
    ErrorHandler errHandler;
    vector <float> tempshadow;
    long nrows=0;
    long i=0;
    int quadid=0;
    int pixelX=0;
    int pixelY=0;
    if(shadowtab.is_empty()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Shadow data is incomplete.";
        throw errHandler;
    }
    nrows = shadowtab.get_nrows();
    if (quadname == "Q0") {
        shadowtab.shadow2D.clear();
        tempshadow.resize(XPIX_QUAD, 0.0);
        shadowtab.shadow2D.resize(YPIX_QUAD, tempshadow);
        for (i = 0; i < nrows; i++) {
            if (shadowtab.locx[i] < XPIX_QUAD && shadowtab.locy[i] >= YPIX_QUAD) {
                pixelX = shadowtab.locx[i];
                pixelY = shadowtab.locy[i]-YPIX_QUAD;
                shadowtab.shadow2D[pixelY][pixelX] = shadowtab.shadow[i];
    
            }
        }
    } else if (quadname == "Q1") {
        shadowtab.shadow2D.clear();
        tempshadow.resize(XPIX_QUAD, 0.0);
        shadowtab.shadow2D.resize(YPIX_QUAD, tempshadow);
        for (i = 0; i < nrows; i++) {
            if (shadowtab.locx[i] >= XPIX_QUAD && shadowtab.locy[i] >= YPIX_QUAD) {
                pixelX = shadowtab.locx[i] - XPIX_QUAD;
                pixelY = shadowtab.locy[i] - YPIX_QUAD;
                shadowtab.shadow2D[pixelY][pixelX] = shadowtab.shadow[i];
            }
        }
    } else if (quadname == "Q2") {
        shadowtab.shadow2D.clear();
        tempshadow.resize(XPIX_QUAD, 0.0);
        shadowtab.shadow2D.resize(YPIX_QUAD, tempshadow);
        for (i = 0; i < nrows; i++) {
            if (shadowtab.locx[i] >= XPIX_QUAD && shadowtab.locy[i] < YPIX_QUAD) {
                pixelX = shadowtab.locx[i] - XPIX_QUAD;
                pixelY = shadowtab.locy[i];
                shadowtab.shadow2D[pixelY][pixelX] = shadowtab.shadow[i];
            }
        }
    } else if (quadname == "Q3") {
        shadowtab.shadow2D.clear();
        tempshadow.resize(XPIX_QUAD, 0.0);
        shadowtab.shadow2D.resize(YPIX_QUAD, tempshadow);
        for (i = 0; i < nrows; i++) {
            if (shadowtab.locx[i] < XPIX_QUAD && shadowtab.locy[i] < YPIX_QUAD) {
                pixelX = shadowtab.locx[i];
                pixelY = shadowtab.locy[i];
                shadowtab.shadow2D[pixelY][pixelX] = shadowtab.shadow[i];
            }
        }
    } else if (quadname == "QALL") {
        shadowtab.shadow2D.clear();
        tempshadow.resize(XSIZE, 0.0);
        shadowtab.shadow2D.resize(YSIZE, tempshadow);
        for (i = 0; i < nrows; i++) {
            pixelX = shadowtab.locx[i];
            pixelY = shadowtab.locy[i];
            shadowtab.shadow2D[pixelY][pixelX] = shadowtab.shadow[i];
        }
    }

    return shadowtab.shadow2D;
}
int ShadowFileHandler::reset() {
    shadowtab.reset();
    shadowFilename="";
    extname="";
    
    return EXIT_SUCCESS;
}


//SHADOW FILE HANDLER END

//SPECTRUMLC CLASS

int SpectrumLc::bin_spectrum(string eventFileName, vector<double> tstart, vector<double> tstop, 
                     vector <float> estart, vector <float> estop, vector <int> quadsToProcess,
                        ExposureTable &expTab,int tfilter,int efilter,int maskWeight){
    int status = 0;
    EventFileHandler evt; //to read and store event file.
    ErrorHandler errh; //error handler for this function.
    fitsfile *fevt; //pointer to event file
    int i = 0, iquad = 0, ievent = 0;
    string errorMsg = "";
    eventFileQuad evtQuad;
    unsigned char locx, locy; //actual location x and y of an event pixel.
    long indexPix; //index record of pixel with locx,locy in exposure table weights array.
    vector <float> sumwc; //Flux calculation
    vector <float> sumw2c; //Poisson standard deviation calculation
    Badpix badpix;

    //Checking if exposure array is generated for ebounds table or not.
    if(maskWeight==true&&expTab.eboundFlag==false){
        errh.severity = errERROR;
        errh.errorStatus = IMPROPER_INPUT;
        errh.errorMsg = "Spectrum cannot be generated for channels other than those provided in CALDB Ebounds file.";
        throw errh;
    }
    
    //Opening energy added event file.
    fits_open_file(&fevt, (char*) eventFileName.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening event file " << eventFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    // Read the badpix file
    if(badpix.read_badpix_file(badpixfile)){
        LOG(ERROR) << "Error in reading bad pixel file " << badpixfile;
        return EXIT_FAILURE;
    }

    vector < vector <unsigned char> > badpixMap;

    badpixMap = badpix.get_badpix_map();
    if(badpixMap.empty()){
        LOG(INFO) << "Creating badpix map.";
        try {badpix.create_badPixMap();}
        catch(string s){
            LOG(ERROR)<< s;
            return EXIT_FAILURE;
        }
        badpixMap = badpix.get_badpix_map();
    }

    //Get the avg exposure time 
    float exp_time;
    getAvgExposureTime(fevt,quadsToProcess,&exp_time);
    exposure_time=exp_time;

//    printf("EXPOSURE TIME %f\n",exp_time);

    //Initializing sumwc and sumw2c
    sumwc.resize(expTab.PIs.size(),0.0);
    sumw2c.resize(expTab.PIs.size(),0.0);
   
	 
    //Generating exposure index table for faster processing
    if(maskWeight==true)
	{
	try{
        expTab.generate_index_table(false);
    } catch (ErrorHandler errHandler){
        throw errHandler;
    }

	}
    // exposure index table generated
    
    //Reading event file for quadrants to be processed
    for (i = 0; i < quadsToProcess.size(); i++) {
        iquad = quadsToProcess[i];

        //Read quadrant extension of event file
        if (evt.read_quad_extension(fevt, iquad)) {
            LOG(ERROR) << "Error in reading Quadrant " << iquad << " data from file: " << eventFileName;
            return EXIT_FAILURE;
        }
        
        //Handling error if event file does not have column Energy or PI.
        if(evt.is_column_empty("Energy")==1 || evt.is_column_empty("PI")){
            errh.severity = errERROR;
            errh.errorStatus = IMPROPER_INPUT;
            errh.errorMsg = "Energy/PI field values for quadrant " + itoa(iquad) + " data is empty.";
            throw errh;
        }
        
		if(maskWeight==true)
		{
		for (ievent = 0; ievent < evt.get_nrowsQuad(); ievent++) {
            try {
                if (evt.event_satisfy_criteria(ievent, tstart, tstop, estart, estop, &evtQuad,tfilter,efilter)) {
                    
                    //Generating location x and location y of this quadrant pixel on 128x128 array
                    generate_locx_locy(evtQuad.detX, (unsigned char)iquad, evtQuad.detY, locx, locy);
                    indexPix = expTab.indexTable[locy][locx];
                    sumwc[evtQuad.PI]+=expTab.weightsArray[indexPix][evtQuad.PI];
                    sumw2c[evtQuad.PI]+=(expTab.weightsArray[indexPix][evtQuad.PI])*
                            (expTab.weightsArray[indexPix][evtQuad.PI]);
                }
            } catch(ErrorHandler errHandler){
                throw errHandler;
            }
        }
		flux = sumwc; //assigning counts 
        transform(sumw2c.begin(), sumw2c.end(), sumw2c.begin(), (float(*)(float))sqrt );
        error = sumw2c;

		}
		else
		{
		for (ievent = 0; ievent < evt.get_nrowsQuad(); ievent++) {
            try {
                if (evt.event_satisfy_criteria(ievent, tstart, tstop, estart, estop, &evtQuad,tfilter,efilter)) {

                    generate_locx_locy(evtQuad.detX, (unsigned char)iquad, evtQuad.detY, locx, locy);
                    
					if(badpixMap[locy][locx] <= badpixthresh) sumwc[evtQuad.PI]+=1;
					
                }
            } catch(ErrorHandler errHandler){
                throw errHandler;
            }
        }
				
    	flux = sumwc; //assigning counts 
    	transform(sumwc.begin(), sumwc.end(), sumwc.begin(), (float(*)(float))sqrt );
		
    	error = sumwc;

			
		}
    }
  
    PI.assign(expTab.PIs.begin(), expTab.PIs.end());
    energy = expTab.energies;
    thetaxd = expTab.thetaxd;
    thetayd = expTab.thetayd;
    energyRange = expTab.energyRange;
    nchannels = expTab.nchannels;

    fits_close_file(fevt,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing event file";
        return (EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}


int SpectrumLc::bin_lc(string eventFileName, vector<double> tstart, vector<double> tstop, 
        vector<float> estart, vector<float> estop, vector<int> quadsToProcess, 
        ExposureTable& expTab, double binsize,int tfilter,int efilter,char *livetimefile,int maskWeight) {
    int status = 0;
    EventFileHandler evt; //to read and store event file.
    ErrorHandler errh; //error handler for this function.
    fitsfile *fevt; //pointer to event file
    fitsfile *flive;
    int i = 0, iquad = 0, ievent = 0;
    string errorMsg = "";
    double timerange=0.0; //maximum(tstop) - minimum(tstart) time specified by user
    double mintstart, mintstop, maxtstart, maxtstop;
    long nbins=0;
    long binindex=0;
    eventFileQuad evtQuad;
    unsigned char locx, locy; //actual location x and y of an event pixel.
    long indexPix; //index record of pixel with locx,locy in exposure table weights array.
    vector <float> sumwc; //Flux calculation
    vector <float> sumw2c; //Poisson standard deviation calculation
    double binsizeInv=1.0/binsize;
    long Tstarti,Tstopi;
    int hdutype=0;
    Badpix badpix;

    this->timedel=binsize;

    //Opening energy added event file.
    fits_open_file(&fevt, (char*) eventFileName.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening event file " << eventFileName;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //getting timerange;

/*
    if(tstart.size()>0 && tstop.size()>0){
        compute_min_and_max(tstart.begin(), tstart.end(), mintstart, maxtstart);
        compute_min_and_max(tstop.begin(), tstop.end(), mintstop, maxtstop);
    } else {
        try{
            evt.get_min_max_time(NULL, quadsToProcess, &mintstart, &maxtstop, eventFileName);
        } catch (ErrorHandler errHandler){
            throw errHandler;
        }
    }
*/
   
    fits_movabs_hdu(fevt, iquad+2, &hdutype, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in moving to quadrant HDU";
        return (EXIT_FAILURE);
    }

		    
    fits_read_key(fevt,TLONG,"TSTARTI",&Tstarti,NULL,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in reading TSTARTI in event file header";
        return (EXIT_FAILURE);
    }

    fits_read_key(fevt,TLONG,"TSTOPI",&Tstopi,NULL,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in reading TSTOPI in event file header";
        return (EXIT_FAILURE);
    }
    
    timerange = (double)Tstopi+1.0 - (double)Tstarti;
    nbins = ceil(timerange/binsize);
    UT.resize(nbins, 0.0);
    fracexp.resize(nbins,1.0);

	LOG(INFO)<<"Numbins "<<nbins;

    for(i=0; i<nbins; i++){
        UT[i] = (double)Tstarti + (i+0.5)*binsize;
    }

    fits_open_file(&flive, livetimefile, READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening livetime file " << livetimefile;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    long ntbins,k;
    double *timearray;
    double *livetime;
    int l;

    if(floor(binsize)!=binsize)
    {
        LOG(INFO)<<"--------------------------------------------------------------------------------------------------------";
        LOG(INFO)<<"   The binsize provided is not an integer. Fractional exposure is incorrect by few percentage ";
        LOG(INFO)<<"--------------------------------------------------------------------------------------------------------";
    }
    
    // Read the badpix file
    if(badpix.read_badpix_file(badpixfile)){
        LOG(ERROR) << "Error in reading bad pixel file " << badpixfile;
        return EXIT_FAILURE;
    }

    vector < vector <unsigned char> > badpixMap;
    badpixMap = badpix.get_badpix_map();
    if(badpixMap.empty()){
        LOG(INFO) << "Creating badpix map.";
        try {badpix.create_badPixMap();}
        catch(string s){
            LOG(ERROR)<< s;
            return EXIT_FAILURE;
        }
        badpixMap = badpix.get_badpix_map();
    }

    
    
    //Initializing sumwc and sumw2c

    sumwc.resize(nbins, 0.0);
    sumw2c.resize(nbins, 0.0);

    //Generating exposure index table for faster processing
    if(maskWeight==true)
    {
        try {
            expTab.generate_index_table(false);
        } catch (ErrorHandler errHandler) {
            throw errHandler;
        }
    }

    //Reading event file for quadrants to be processed
    for (i = 0; i < quadsToProcess.size(); i++) {
        iquad = quadsToProcess[i];
        if (evt.read_quad_extension(fevt, iquad)) {
            LOG(ERROR) << "Error in reading Quadrant " << iquad << " data from file: " << eventFileName;
            return EXIT_FAILURE;
        }

        if(floor(binsize)==binsize)
        {
        //Read livetime
        try{
        fits_movabs_hdu(flive, iquad+2, &hdutype, &status);

        fits_get_num_rows(flive, &ntbins, &status);
        }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }

        livetime=(double*)malloc(sizeof(double)*ntbins);
        timearray=(double*)malloc(sizeof(double)*ntbins);

        try{
        fits_read_col(flive, TDOUBLE, 1, 1, 1, ntbins, NULL, timearray,NULL, &status);
        fits_read_col(flive, TDOUBLE, 2, 1, 1, ntbins, NULL, livetime,NULL, &status);
        }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }

        for(k=0;k<nbins;k++)
        {
            fracexp[k]=0;
            for(l=Tstarti+k*binsize-timearray[0]-0.5;l<Tstarti+k*binsize-timearray[0]-0.5+binsize;l++)
                fracexp[k]+=livetime[l];

            fracexp[k]/=binsize;
			if(fracexp[k]==0) fracexp[k]=1.0;
        }

		}

        //Handling error if event file does not have column Energy or PI.
        if (evt.is_column_empty("Energy") == 1 || evt.is_column_empty("PI")) {
            errh.severity = errERROR;
            errh.errorStatus = IMPROPER_INPUT;
            errh.errorMsg = "Energy/PI field values for quadrant " + itoa(iquad) + " data is empty.";
            throw errh;
        }

        LOG(INFO) << "Binning events for quadrant " << iquad;

        for (ievent = 0; ievent < evt.get_nrowsQuad(); ievent++) {
            try {
                if (evt.event_satisfy_criteria(ievent, tstart, tstop, estart, estop, &evtQuad,tfilter,efilter)) {
                    binindex = (evtQuad.time - Tstarti)/binsize;
                   
				   	//LOG(INFO)<<"event "<<ievent<<" binindex "<<binindex; 
                    if(binindex<0||binindex>=nbins) {

						LOG(INFO)<<"Event before tstart or after tstop";
						LOG(INFO)<<"Evt time "<<setprecision(20)<<evtQuad.time<<" tstarti "<<Tstarti<<" tstopi "<<Tstopi;

						continue;}
                    
                    if(maskWeight==true)
                    {
                    generate_locx_locy(evtQuad.detX, (unsigned char) iquad, evtQuad.detY, locx, locy);
                    indexPix = expTab.indexTable[locy][locx];
                    sumwc[binindex] += expTab.weightsArray[indexPix][evtQuad.PI];
                    sumw2c[binindex] += (expTab.weightsArray[indexPix][evtQuad.PI])*
                            (expTab.weightsArray[indexPix][evtQuad.PI]);
                    }
                    else
                    {
                    generate_locx_locy(evtQuad.detX, (unsigned char)iquad, evtQuad.detY, locx, locy);

                    if(badpixMap[locy][locx] <= badpixthresh){

                        sumwc[binindex] +=1;
                        sumw2c[binindex]+=1;
                        }
                    }

                }
            } catch (ErrorHandler errHandler) {
                throw errHandler;
            }
        }

    }

    
    //Convert to rate and compute error
    transform(sumwc.begin(), sumwc.end(), sumwc.begin(), bind1st(multiplies<double>(), binsizeInv) );
    flux = sumwc;
    transform(sumw2c.begin(), sumw2c.end(), sumw2c.begin(), (float(*)(float))sqrt );
    transform(sumw2c.begin(), sumw2c.end(), sumw2c.begin(), bind1st(multiplies<double>(), binsizeInv) );
    error = sumw2c;
 
    thetaxd = expTab.thetaxd;
    thetayd = expTab.thetayd;
    energyRange = expTab.energyRange;
    nchannels = expTab.nchannels;

    if (floor(binsize)==binsize)
    {
        free(timearray);
        free(livetime);
    }

    fits_close_file(flive,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing event file";
        return (EXIT_FAILURE);
    }

    fits_close_file(fevt,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing event file";
        return (EXIT_FAILURE);
    }

    return EXIT_SUCCESS;    
}

int SpectrumLc::verify_spectrum() {
    int status=0;
    ErrorHandler errHandler;
    
    if(flux.size()!=PI.size()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Improper spectrum input as number of records in flux and PI column are not equal.";
        throw errHandler;
    }
    if(error.size()!=PI.size()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Improper spectrum input as number of records in error and PI column are not equal.";
        throw errHandler;
    }
    
    return EXIT_SUCCESS;
}
int SpectrumLc::verify_lc() {
    int status=0;
    ErrorHandler errHandler;
    
    if(flux.size()!=UT.size()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Improper lc input as number of records in flux and UT column are not equal.";
        throw errHandler;
    }
    if(error.size()!=UT.size()){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Improper lc input as number of records in error and UT column are not equal.";
        throw errHandler;
    }
    
    return EXIT_SUCCESS;
}

int SpectrumLc::clear()
{

    flux.clear();
    energy.clear();
    PI.clear();
    UT.clear();
    error.clear();

    return(EXIT_SUCCESS);
}
//SPECTRUMLC CLASS END


//SPECTRUMFILEHANDLER CLASS

int SpectrumFileHandler::write_spectrum_file(string specFilename, string specTemplate,string infile) {
    int status = 0;
    int colnum = 0;
    long nrows = 0;
    fitsfile *fspec;
    vector <string> keys;
    vector <float> valuesF;
    vector <string> valuesS;

    int hdunum=0;
    int hdutype=0;
    fitsfile *fevt;

	spec.PI=spec.channels;

    hdunum=spec.quadrantid+2;

    try{
	           fits_open_file(&fevt, (char*)infile.c_str(), READONLY, &status);
	                  fits_movabs_hdu(fevt,hdunum,&hdutype,&status);

	    } catch (ErrorHandler errHandler) {
                      logError(errHandler);
                              return EXIT_FAILURE;
	                                }


    try{
        spec.verify_spectrum();
    } catch (ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }

    if (create_empty_fitsfile(specFilename, specTemplate)) {
        LOG(ERROR) << "Error in creating spectrum file from corresponding template";
        return EXIT_FAILURE;
    }

    fits_open_file(&fspec, (char*) specFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file: " << specFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_movnam_hdu(fspec, BINARY_TBL, "SPECTRUM", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to SPECTRUM in file : " << specFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
   
  
	 
    cztHeaderParam headkey;
    headkey.readFromHeader(fevt);

    headkey.writeToHeader(fspec);

	

    //Writing PI channel
    if (write_fits_column(fspec, "CHANNEL", TINT, 1, 1, spec.PI)) {
        LOG(ERROR) << "Error in writing column CHANNEL of spectrum file " << specFilename;
        return EXIT_FAILURE;
    }
    //Writing flux
    if (write_fits_column(fspec, "COUNTS", TFLOAT, 1, 1, spec.flux)) {
        LOG(ERROR) << "Error in writing column COUNTS of spectrum file " << specFilename;
        return EXIT_FAILURE;
    }
    //Writing Error
    if (write_fits_column(fspec, "STAT_ERR", TFLOAT, 1, 1, spec.error)) {
        LOG(ERROR) << "Error in writing column ERROR of spectrum file " << specFilename;
        return EXIT_FAILURE;
    }


    //Updating keys
    string keysArray[] = {"THETAX", "THETAY", "NCHANNEL","EXPOSURE","MASKWT","PIXTHRES","QUADID"};
    float valuesFArray[] = {spec.thetaxd, spec.thetayd, spec.nchannels,spec.exposure_time,spec.maskweighted,spec.badpixthresh,spec.quadrantid};
    keys.assign(keysArray, keysArray + 6);
    valuesF.assign(valuesFArray, valuesFArray + 6);
    if (updateKeys(fspec, "SPECTRUM", keys, valuesF, TFLOAT)) {
        LOG(ERROR) << "Error in updating keys.";
        return EXIT_FAILURE;
    }
    fits_update_key_str(fspec, "ERANGE", (char*) (spec.energyRange).c_str(), NULL, &status);
    if (status) {
        LOG(WARNING) << "Unable to update key ERANGE";
    }

    fits_update_key_str(fspec, "GTITYPE", (char*) (spec.gtitype.c_str()), NULL, &status);
    if (status) {
        LOG(WARNING) << "Unable to update key ERANGE";
    }
	
        if(spec.allquad==0)
            updateKey(fspec,TINT,"QUADID",&spec.quadrantid);
	
    //keys updated
    
    fits_close_file(fspec, &status);
    if (status) {
        LOG(ERROR) << "Error in closing spectrum file " << specFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fevt,&status);
        if (status) {
		                    LOG(ERROR) << "Error in closing event file " ;
				                                fits_report_error(stderr, status);
								                                    return (EXIT_FAILURE);
												                                            }


    return EXIT_SUCCESS;
}

//SPECTRUMFILEHANDLER CLASS END

//LIGHTCURVEFILEHANDLER CLASS

int LightCurveFileHandler::write_lc_file(string lcFilename, string lcTemplate,string infile) {
    int status = 0;
    int colnum = 0;
    long nrows = 0;
    fitsfile *flc;
    char* keyValue;
    int hdunum=0;
    int hdutype=0;
    fitsfile *fevt;
	
    hdunum=lc.quadrantid+2;

    try{
	   fits_open_file(&fevt, (char*)infile.c_str(), READONLY, &status);
           fits_movabs_hdu(fevt,hdunum,&hdutype,&status); 

    } catch (ErrorHandler errHandler) {
	            logError(errHandler);
		            return EXIT_FAILURE;
			        }


    try {
        lc.verify_lc();
    } catch (ErrorHandler errHandler) {
        logError(errHandler);
        return EXIT_FAILURE;
    }

    if (create_empty_fitsfile(lcFilename, lcTemplate)) {
        LOG(ERROR) << "Error in creating light curve file from corresponding template";
        return EXIT_FAILURE;
    }

    fits_open_file(&flc, (char*) lcFilename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file: " << lcFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_movnam_hdu(flc, BINARY_TBL, "LIGHTCURVE", NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to SPECTRUM in file : " << lcFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Read and write keywords
    cztHeaderParam headkey;
    headkey.readFromHeader(fevt);
    headkey.writeToHeader(flc);
    
    double telapse;
    headkey.getTelapse(telapse);

    fits_update_key(flc,TDOUBLE,"TELAPSE",&telapse,"Elapsed time",&status);


    //Writing PI channel
    if (write_fits_column(flc, "TIME", TDOUBLE, 1, 1, lc.UT)) {
        LOG(ERROR) << "Error in writing column TIME of light curve file " << lcFilename;
        return EXIT_FAILURE;
    }
    //Writing flux
    if (write_fits_column(flc, "RATE", TFLOAT, 1, 1, lc.flux)) {
        LOG(ERROR) << "Error in writing column RATE of light curve file " << lcFilename;
        return EXIT_FAILURE;
    }
    //Writing Error
    if (write_fits_column(flc, "ERROR", TFLOAT, 1, 1, lc.error)) {
        LOG(ERROR) << "Error in writing column ERROR of light curve file " << lcFilename;
        return EXIT_FAILURE;
    }

    //Write fracexp
    if (write_fits_column(flc, "FRACEXP", TFLOAT, 1, 1, lc.fracexp)) {
        LOG(ERROR) << "Error in writing column ERROR of light curve file " << lcFilename;
        return EXIT_FAILURE;
    }
    

    LOG(INFO)<<"Timedel is "<<lc.timedel;

    //Updating keys
    try{
        updateKey(flc, TFLOAT, "THEATAX", &lc.thetaxd);
        updateKey(flc, TFLOAT, "THEATAY", &lc.thetaxd);
        keyValue = (char*)(lc.energyRange).c_str();
        updateKey(flc, TSTRING, "ERANGE", &keyValue);
        updateKey(flc,TDOUBLE,"TIMEDEL",&lc.timedel);
		updateKey(flc,TINT,"MASKWT",&lc.maskweighted);
		updateKey(flc,TINT,"PIXTHRES",&lc.badpixthresh);
		if(lc.allquad==0)
	        updateKey(flc,TINT,"QUADID",&lc.quadrantid);
				
    }
    catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }
    //keys updated

    fits_close_file(flc, &status);
    if (status) {
        LOG(ERROR) << "Error in closing LIGHT CURVE file " << lcFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_close_file(fevt,&status);
    if (status) {
	            LOG(ERROR) << "Error in closing event file " ;
		            fits_report_error(stderr, status);
			            return (EXIT_FAILURE);
				        }


    return EXIT_SUCCESS;
}
//LIGHTCURVEFILEHANDLER CLASS END

//IMAGE
Image::Image(){
    teldefFile = "";
    aspectFile = "";
    oversamplingfactor=1;
    ysize=0;
    xsize=0;
    //DPI KEYWORDS
    eMin = 0.0;
    eMax = 0.0;
    tStart=0.0;
    tStop=0.0;
    badpixThreshold = GOODPIX;
    //WCS KEYWORDS
    //RA DEC
    crpix1=0.0;
    crpix2=0.0;
    crval1=0.0;
    crval2=0.0;
    cdelt1=0.0;
    cdelt2=0.0;
    crota1=0.0;
    crota2=0.0;
    cd1_1=0.0;
    cd1_2=0.0;
    cd2_1=0.0;
    cd2_2=0.0;
    //Camera coordinates
    crpix1b=0.0;
    crpix2b=0.0;
    crval1b=0.0;
    crval2b=0.0;
    cdelt1b=0.0;
    cdelt2b=0.0;
    crota1b=0.0;
    crota2b=0.0;
    cd1_1b=0.0;
    cd1_2b=0.0;
    cd2_1b=0.0;
    cd2_2b=0.0;
}
void Image::evaluate_wcs_coordinates(int oversamplingFactor, string aspectFile ) {
    int status=0;
    AspectFileHandler aspect;
    ErrorHandler errHandler;
    float RA=0.0, DEC=0.0, TWIST=0.0;
    float scale=0.0;
    
    this->aspectFile = aspectFile;
    this->oversamplingfactor = oversamplingFactor;
    //reading aspect file
    if(!FileExists((char*) aspectFile.c_str())){
        aspectFile = "";
    }
    if (aspectFile != "") {
        if (aspect.read_aspect_file(aspectFile)) {
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error reading aspect file " + aspectFile;
            throw errHandler;
        }
        RA = aspect.get_RA();
        DEC = aspect.get_DEC();
        TWIST = aspect.get_TWIST();
        this->teldefFile = aspect.get_teldefFilename();
    }
    scale = atan(XLDET/(MASKHEIGHT*16*oversamplingfactor)) * TODEG; //scale in degrees
    //evaluate wcs coordinates for ICRS
    cdelt1 = scale;
    cdelt2 = scale;
    crota1 = -TWIST;
    crota2 = -TWIST;
    crpix1 = 15*oversamplingfactor +1;
    crpix2 = 15*oversamplingfactor +1;
    crval1 = RA;
    crval2 = DEC;
    cd1_1 = scale * cos(crota1*TORAD);
    cd1_2 = scale * sin(crota1*TORAD);
    cd2_1 = -cd1_2;
    cd2_2 = cd1_1;
    //evaluate wcs coordinates for LOCALE CAMERA COORDINATES
    cdelt1b = scale;
    cdelt2b = scale;
    crota1b = 0.0;
    crota2b = 0.0;
    crpix1b = 15 * oversamplingfactor + 1;
    crpix2b = 15 * oversamplingfactor + 1;
    crval1b = 0.0;
    crval2b = 0.0;
    cd1_1b = scale;
    cd1_2b = 0;
    cd2_1b = -cd1_2b;
    cd2_2b = cd1_1b;
}

int Image::is_empty(){
    ErrorHandler errHandler;
    int ysize=crossCorrelationImage.size();
    if(ysize==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_STRUCTURE_NOT_COMPLETE;
        errHandler.errorMsg = "Image 2d vector is empty.";
        throw errHandler;
    }
}
//IMAGE
void Image::generate_cross_correlation_image(vector<vector<float> >* dpi, vector<vector<float> >* mask, int oversampfactor, float xshift,float yshift) {
    vector <float> dpi1d;
    vector <float> mask1d;
    vector <float> correlationMatrix1d;
    ErrorHandler errHandler;
    int ydpisize = 0;
    int xdpisize = 0;
    int ymasksize = 0;
    int xmasksize = 0;
    
    ydpisize = (*dpi).size();
    ymasksize = (*mask).size();
    if(ydpisize==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "DPI provided by user to generate image is empty.";
        throw errHandler;
    }
    if(ymasksize==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Mask provided by user to generate image is empty.";
        throw errHandler;
    }
    xdpisize = (*dpi)[0].size();
    xmasksize = (*mask)[0].size();
    if(ydpisize!=ymasksize || xdpisize!=xmasksize){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Mask and DPI dimension provide by user do not match.";
        throw errHandler;
    }
    
    try {
        dpi1d.clear();
        linearize_2Dvector((*dpi), &dpi1d);
        mask1d.clear();
        linearize_2Dvector((*mask), &mask1d);
        DLOG(INFO) << "Oversampling DPI.";
        oversample((&dpi1d), ydpisize, xdpisize, oversampfactor);
        DLOG(INFO) << "Oversampling mask.";

        oversample((&mask1d), ymasksize, xmasksize, oversampfactor);
        crossCorrelationImage = generate_correlation_matrix2D(dpi1d, mask1d,
      	&correlationMatrix1d, oversampfactor,xshift,yshift);
        ysize = crossCorrelationImage.size();
        if(ysize!=0){
            xsize = crossCorrelationImage[0].size();
        }
    } catch(ErrorHandler errHandler){
        throw errHandler;
    }

}


//IMAGE END

//Image File Handler

ImageFileHandler::ImageFileHandler() {
    this->imageFilename = "";
}

void ImageFileHandler::create_image_file(string outImageFilename, string imageTemplate){
    int status=0;
    ErrorHandler errHandler;
    if(create_empty_fitsfile(outImageFilename, imageTemplate)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = UNABLE_TO_WRITE;
        errHandler.errorMsg = "Error in creating empty image file " + outImageFilename;
        throw errHandler;
    }
    
}

int ImageFileHandler::write_image_file(string imgFilename, string extname, Image *img){
    int status=0;
    fitsfile *fimg;
    ErrorHandler errHandler;
    vector <float> img1d;
    char *keyValue;
    int keyIntValue;
    long fpixel[2];
    long naxes[2];
    long naxis1=0;
    long naxis2=0;
    fpixel[0]=1;
    fpixel[1]=1;
    
    fits_open_file(&fimg, (char*) imgFilename.c_str(), READWRITE, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error opening image file " + imgFilename;
        throw errHandler;
    }
    
    //Checking whether image data is proper or not.
    try{
        (*img).is_empty();
    } catch(ErrorHandler errHandler){
        throw(errHandler);
    }

    //Writing image data
    naxis1 = (*img).get_imgysize();
    if(naxis1==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Error: Image data empty.";
        throw errHandler;
    }
    naxis2 = (*img).get_imgxsize();
    naxes[0]=naxis1;
    naxes[1]=naxis2;
    //linearizing image data
    try{
        linearize_2Dvector((*img).crossCorrelationImage, &img1d);
    } catch(ErrorHandler errHandler){
        logError(errHandler);
        throw(errHandler);
    }
    
//    }
    fits_create_img(fimg, FLOAT_IMG, 2, naxes, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in creating image extension";
        throw errHandler;
    }
    DLOG(INFO) << "[DEBUG] Writing image with dimensions " << naxis1 << " x " << naxis2;

    fits_write_pix(fimg, TFLOAT, fpixel, img1d.size(), img1d.data(), &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in writing image data to extension " + extname + " of image file " + imgFilename;
        throw errHandler;
    }
    
    //Updating basic header keywords
    try{
        updateKey(fimg, TINT, "OVERSAMP", &(*img).oversamplingfactor);
    } catch(ErrorHandler errHandler){
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    //Updating header keywords (WCS)
    try {
        //RA DEC
        keyValue = (char*) extname.c_str();
        updateKey(fimg, TSTRING, "EXTNAME", keyValue);
        keyValue = "ICRS";
        updateKey(fimg, TSTRING, "WCSNAME", keyValue);
        keyIntValue = 2;
        updateKey(fimg, TINT, "WCSAXES", &keyIntValue);
        keyValue = "RA---TAN";
        updateKey(fimg, TSTRING, "CTYPE1", keyValue);
        keyValue = "DEC--TAN";
        updateKey(fimg, TSTRING, "CTYPE2", keyValue);
        updateKey(fimg, TFLOAT, "CRPIX1", &(*img).crpix1);
        updateKey(fimg, TFLOAT, "CRPIX2", &(*img).crpix2);
        updateKey(fimg, TFLOAT, "CRVAL1", &(*img).crval1);
        updateKey(fimg, TFLOAT, "CRVAL2", &(*img).crval2);
        updateKey(fimg, TFLOAT, "CDELT1", &(*img).cdelt1);
        updateKey(fimg, TFLOAT, "CDELT2", &(*img).cdelt2);
        updateKey(fimg, TFLOAT, "CROTA1", &(*img).crota1);
        updateKey(fimg, TFLOAT, "CROTA2", &(*img).crota2);
        updateKey(fimg, TFLOAT, "CD1_1", &(*img).cd1_1);
        updateKey(fimg, TFLOAT, "CD1_2", &(*img).cd1_2);
        updateKey(fimg, TFLOAT, "CD2_1", &(*img).cd2_1);
        updateKey(fimg, TFLOAT, "CD2_2", &(*img).cd2_2);
        //Camera coordinates
        keyValue = "camera";
        updateKey(fimg, TSTRING, "WCSNAMEb", keyValue);
        keyIntValue = 2;
        updateKey(fimg, TINT, "WCSAXESb", &keyIntValue);
        keyValue = "X----TAN";
        updateKey(fimg, TSTRING, "CTYPE1b", keyValue);
        keyValue = "Y----TAN";
        updateKey(fimg, TSTRING, "CTYPE2b", keyValue);
        updateKey(fimg, TFLOAT, "CRPIX1b", &(*img).crpix1b);
        updateKey(fimg, TFLOAT, "CRPIX2b", &(*img).crpix2b);
        updateKey(fimg, TFLOAT, "CRVAL1b", &(*img).crval1b);
        updateKey(fimg, TFLOAT, "CRVAL2b", &(*img).crval2b);
        updateKey(fimg, TFLOAT, "CDELT1b", &(*img).cdelt1b);
        updateKey(fimg, TFLOAT, "CDELT2b", &(*img).cdelt2b);
        updateKey(fimg, TFLOAT, "CROTA1b", &(*img).crota1b);
        updateKey(fimg, TFLOAT, "CROTA2b", &(*img).crota2b);
        updateKey(fimg, TFLOAT, "CD1_1b", &(*img).cd1_1b);
        updateKey(fimg, TFLOAT, "CD1_2b", &(*img).cd1_2b);
        updateKey(fimg, TFLOAT, "CD2_1b", &(*img).cd2_1b);
        updateKey(fimg, TFLOAT, "CD2_2b", &(*img).cd2_2b);
        keyValue = (char*) ((*img).teldefFile).c_str();
        updateKey(fimg, TSTRING, "TELDEF", keyValue);
        updateKey(fimg, TINT, "OVERSAMP", &(*img).oversamplingfactor);
    	headerParam.writeToHeader(fimg); 
    }    catch (ErrorHandler errHandler) {
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    //keys updated
    
    fits_close_file(fimg, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in closing image file " + imgFilename;
        throw errHandler;
    }

    return status;
}

//Image File Handler End


//INDEPENDENT FUNCTIONS
int verify_num_hdus(fitsfile* fptr, int expectedNumHDUs, int status){
    int numhdu;
    string errorMsg="";
    
    fits_get_num_hdus(fptr, &numhdu, &status);
    errorMsg = "Error in getting number of HDUs in input file";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    
    if (numhdu != expectedNumHDUs) {
        LOG(ERROR) << "Number of HDUs found in input file: " << numhdu;
        LOG(ERROR) << "Number of HDUs expected in input file: " << expectedNumHDUs;
        status = EXIT_FAILURE;
        return status;
    }
    
    return status;
}


// Functions added by Mithun

// Compute the avg exposure time
int getAvgExposureTime(fitsfile *fevt,vector <int>quadsToProcess,float *exp_time)
{

    int status=0;
    int i,qid,j;
    float qexp_time;
    double exposure_fraction[4096];
    int exp_col;

    double tot_exp=0.;
    long npix=0;


    for(i=0;i<quadsToProcess.size();i++)
    {

    qid=quadsToProcess[i];

    char quadextname[20];
    sprintf(quadextname,"Q%d",qid);

    fits_movnam_hdu(fevt,BINARY_TBL,quadextname,0,&status);
    if (status) {
        LOG(ERROR) << "Error in opening quad extension of event file";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_read_key(fevt,TFLOAT,"EXPOSURE",&qexp_time,NULL,&status);
    if (status) {
        LOG(ERROR) << "Error in reading exp_time from event file";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    char colname[20];
    sprintf(colname,"EXPOSURE_Q%d",qid);

    fits_movnam_hdu(fevt,BINARY_TBL,"EXPOSURE",0,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in moving to EXPOSURE extension of fits file";
        return (EXIT_FAILURE);
    }

    fits_get_colnum(fevt,CASEINSEN,colname,&exp_col,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Not finding exposure fraction in the event file";
        return (EXIT_FAILURE);
    }

    fits_read_col(fevt, TDOUBLE, exp_col, 1, 1, 4096, NULL, &exposure_fraction,NULL, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in moving to EXPOSURE extension of fits file";
        return (EXIT_FAILURE);
    }

    for(j=0;j<4096;j++)
    {
        tot_exp+=exposure_fraction[j];

        if(exposure_fraction[j]!=0)
            npix+=1;  
    }
    
    }

    //printf("IN EXP %f\t%f\t%d\n",qexp_time,tot_exp,npix);

    //Compute the avg exposure time
    tot_exp/=(float)npix;

    *exp_time=qexp_time*tot_exp;
    
}


int writelivetime(double *timearray,double*livetime,long num_timebins,int qid,char*livetimefile)
{
    int status=0;

    char quadrant_ext[4];
    fitsfile *fout;

    sprintf(quadrant_ext,"Q%d",qid);

    fits_open_file(&fout, livetimefile, READWRITE, &status);
    if (status) {
        LOG(ERROR) <<"Error in opening livetime file";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


    fits_movnam_hdu(fout, BINARY_TBL, quadrant_ext, 0, &status);
    if (status) {
        LOG(ERROR) <<"Error in moving to hdu "<<quadrant_ext;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }


    fits_write_col(fout, TDOUBLE, 1, 1,1, num_timebins, timearray, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_col(fout, TDOUBLE, 2, 1,1, num_timebins, livetime, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_close_file(fout,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    return (EXIT_SUCCESS);
}

//int create_livetime_file(char *livetimefile)//commented by Mayuri,16th Dec,2017
int create_livetime_file(char *livetimefile,double livetime_binsize)
{

    int status=0,bitpix=8,naxis=0;
    int tfields=2,qid=0;
    char telescope[10]="ASTROSAT";
    char instrument[10]="CZTI";
    char comment[100],quadrant_ext[4];
    fitsfile *fout;


    // Create file
    fits_create_file(&fout,livetimefile,&status);
    if(status) {fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    // Create Primary extension  
    fits_create_img(fout,bitpix,naxis,NULL,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    strcpy(comment,"Name of mission/satellite");

    fits_write_key(fout,TSTRING,"TELESCOP",telescope,comment,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_key(fout,TDOUBLE,"LV_BINSIZE",&livetime_binsize,"Livetime binsize",&status);//added by Mayuri,16th Dec 2017
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }//added by Mayuri,16th Dec 2017

    strcpy(comment,"Name of Instrument/detector");
    fits_write_key(fout,TSTRING,"INSTRUME",instrument,comment,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_write_chksum(fout, &status);
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }

     char *ttype[] = {"TIME","FRACEXP"};
     char *tform[] = {"D","E"};

       for(qid=0;qid<4;qid++)
       {

       sprintf(quadrant_ext,"Q%d",qid);

       fits_create_tbl(fout, BINARY_TBL, 0,  tfields,ttype,tform,NULL, quadrant_ext, &status);
       if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

       fits_write_chksum(fout, &status);
       if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }

	fits_write_key(fout,TDOUBLE,"LV_BINSIZE",&livetime_binsize,"Livetime binsize",&status);//added by Mayuri,16th Dec 2017
    	if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }//added by Mayuri,16th Dec 2017

       }

       fits_close_file(fout,&status);
       if(status){fits_report_error(stderr,status); return(EXIT_FAILURE);}


    return(0);
}


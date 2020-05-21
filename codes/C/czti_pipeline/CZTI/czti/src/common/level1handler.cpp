#include "level1handler.h"

using namespace std;

//TCT
Tct::Tct(){
    nrecOutTCT=0;
	last_index=0;
}

//Function to compute average UT/CZT time fraction
int Tct::get_clockCorrectionFactor(double &corFac){
    
    int i;
    long nrows = vec_instrumentTime.size();
    double instrDiff,spsDiff;

    corFac=0;

        for (i = 1; i < nrows; i++)
        {
            instrDiff = vec_instrumentTime[i] - vec_instrumentTime[i-1];
            spsDiff = vec_spsTime[i] - vec_spsTime[i-1];

            if(spsDiff==0||instrDiff==0)
            {
                LOG(ERROR) << "Wrong UT/Instrument time";            
                return EXIT_FAILURE;
            }
            else
            corFac+=(spsDiff/instrDiff);
        }
        corFac = (corFac/(nrows-1));

        return EXIT_SUCCESS;
}


int Tct::clean_tct_data() {
    int status=0;
    long i=0;
    int j=0;
    long nrows = vec_instrumentTime.size();
    double delta = 16.0;
    double difference = 0.0;
    vector <double> delta1;
    vector <unsigned int> spuriousRecords;
    delta1.resize(nrows - 1);
    do {
        spuriousRecords.clear();
        for (i = 1; i < nrows; i++) {
            difference = vec_instrumentTime[i] - vec_instrumentTime[i-1];
            if (abs(difference - (1.0 * delta)) > DELTA_INSTRUMENT_TIME) {
                LOG(INFO) << i << " " << difference << " " << abs(difference - (1.0 * delta));
                spuriousRecords.push_back(i-1);
            }
        }
        LOG(INFO) << spuriousRecords.size() << "--------------";
        //print_vector(spuriousRecords);
        for (j = 0; j < spuriousRecords.size(); j++) {
            vec_instrumentTime[spuriousRecords[j]] = vec_instrumentTime[spuriousRecords[j]+1] - 16.0;
            vec_spsTime[spuriousRecords[j]] = vec_spsTime[spuriousRecords[j]+1] - 16.0;
        }
    } while (spuriousRecords.size() != 0);


    return status;
}

int Tct::read_tct_file(string tctFilename){
    int status=0; //status variable
    long i,j=0; //counter variable
    long nrows=0;
    fitsfile *ftct;
    
    //reading TCT file
    fits_open_file(&ftct, tctFilename.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening TCT file : " << tctFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_movnam_hdu(ftct, BINARY_TBL, "TCT", 0, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to TCT hdu in TCT file : " << tctFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    fits_get_num_rows(ftct, &nrows, &status);
    if (status) {
        LOG(ERROR) << "Error in getting number of rows from TCT file : " << tctFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    //Reading SPS_TIME
    if(read_fits_column(ftct, "SPS_TIME", TDOUBLE, 1, 1, nrows, vec_spsTime)) {
        LOG(ERROR) << "Error in reading SPS_TIME column of TCT file: " << tctFilename;
        return (EXIT_FAILURE);
    }
    //Reading SPS_OBT  
    if(read_fits_column(ftct, "SPS_OBT",TDOUBLE, 1, 1, nrows, vec_spsObt)) {
        LOG(ERROR) << "Error in reading SPS_OBT column of TCT file: " << tctFilename;
        return (EXIT_FAILURE);
    }
    //Reading INSTRUMENT_TIME  
    if(read_fits_column(ftct, "INSTRUMENT_TIME", TDOUBLE, 1, 1, nrows, vec_instrumentTime)) {
        LOG(ERROR) << "Error in reading INSTRUMENT_TIME column of TCT file: " << tctFilename;
        return (EXIT_FAILURE);
    }
    
    //Calculating slope and offset for each 100 second interval of tct
    
    fits_close_file(ftct, &status);
    if (status) {
        LOG(ERROR) << "Error in closing TCT file: " << tctFilename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }    
    
   
   	this->is_reset();
	 
    return status;
}

int Tct::is_reset(){
    long nrows;
    long i,j=0;
    int resetSize=0;
    nrows=vec_instrumentTime.size();
    resetFlag=NO;
    resetIndex.clear();
    for(i=1; i<nrows; i++){
        if(vec_instrumentTime[i] < vec_instrumentTime[i-1]){
            resetFlag = YES;
            resetIndex.push_back(i);
            LOG(INFO)<<"Reset found at index : " << i;
        }
    }
    resetSize=resetIndex.size();
    if(resetSize>1){
        LOG(ERROR)<<"Error for TCT interpolation. Reset is found twice in TCT file.";
		LOG(ERROR)<<"I can't handle more than one reset in TCT file. Exiting....";
        LOG(ERROR)<<setw(10)<<"Index"<< setw(10)<<"Instrument Time";
        for(i=0; i<resetSize; i++){
            LOG(ERROR)<<setw(10)<<resetIndex[i]<< setw(10)<<vec_instrumentTime[resetIndex[i]];
        }
        return (EXIT_FAILURE);
    }
	else if (resetSize==1)
	{
		if(vec_instrumentTime[nrows-1]>=vec_instrumentTime[0])
		{
			LOG(ERROR)<<"Confusion exists in TCT interpolation.";
			LOG(ERROR)<<"Unable to perforam TCT interpolation.";
			return(EXIT_FAILURE);
		}
		
	}
    
    return EXIT_SUCCESS;
}

int Tct::interpolate_time(double instrumentTime, double &intTime, bool isReset){
    long nrows=0;
    long i,j=0;
    int interpolateFlag=0;
    long startIndex, endIndex;
    nrows = vec_instrumentTime.size();
    
    if(nrows==0){
        LOG(ERROR) << "No TCT data.";
        return EXIT_FAILURE;
    }

/*	if(isReset==YES && resetFlag==NO){
        LOG(ERROR)<< "Reset observed in event file but no reset in TCT file";
        return (EXIT_FAILURE);
    }
    
    if(isReset == YES) {
        startIndex = resetIndex[1];
        endIndex = nrows-1;
    }
    else if(isReset == NO){
        startIndex =0;
        endIndex = nrows -1;
    }

*/

    startIndex =0;
	endIndex = nrows -1;


    /* Modified by ajay and mithun to take care of
     * reset of intrument time present in tct and 
     * the 2s delay of time tagging of frames from CZTI FEB
     * Dec 7 2015, 09:50 PM.
     */
    instrumentTime-=2;
    instrumentTime-=((int)(instrumentTime/1048576)*1048576.0);
    //End of fix

  
	/* Modified by Mithun to improve the performance
	 * Jan 7 2016
	 *  
	 * Modified by Mithun to take care of reset within a TCT file 
	 * Jan 15 2016
	 *
	 * Bug fix to account for interpolation of time stamps between the 
	 * roll over. Added one missing else if statement below.
	 *
	 * - Mithun N P S (16/06/16)
 	 */

  if(resetFlag==NO)
  {
    // Check if the instrument time is out of the tct range and interpolate with the first two/last two
	if (instrumentTime < vec_instrumentTime[0]) {
		nrecOutTCT++;
        if (find_interpolated_value(instrumentTime, vec_instrumentTime[0], vec_instrumentTime[1],
                vec_spsTime[0], vec_spsTime[1], intTime)) {
            LOG(ERROR) << "***Error in TCT time interpolation.***";
            return EXIT_FAILURE;
        }
    } else if (instrumentTime >= vec_instrumentTime[endIndex]) {
		nrecOutTCT++;
        if (find_interpolated_value(instrumentTime, vec_instrumentTime[endIndex-1], vec_instrumentTime[endIndex],
                vec_spsTime[endIndex-1], vec_spsTime[endIndex], intTime)) {
            LOG(ERROR) << "***Error in TCT time interpolation.***";
            return EXIT_FAILURE;
        }

    }
	//Else find the nearest records in tct file
	else
	{
	startIndex=this->last_index;

    for (i = startIndex; i < endIndex; i++) {
        if (instrumentTime >= vec_instrumentTime[i] && instrumentTime < vec_instrumentTime[i + 1]) {
            if (find_interpolated_value(instrumentTime, vec_instrumentTime[i], vec_instrumentTime[i + 1],
                    vec_spsTime[i], vec_spsTime[i + 1], intTime)) {
                LOG(ERROR) << "***Error in TCT time interpolation.***";
                return EXIT_FAILURE;
            }
            interpolateFlag = 1;
		    this->last_index=i;

            break;
        }
    }

	if(interpolateFlag!=1) 
	{
		startIndex=0;

    	for (i = startIndex; i < endIndex; i++) {
        	if (instrumentTime >= vec_instrumentTime[i] && instrumentTime < vec_instrumentTime[i + 1]) {
            	if (find_interpolated_value(instrumentTime, vec_instrumentTime[i], vec_instrumentTime[i + 1],
                	    vec_spsTime[i], vec_spsTime[i + 1], intTime)) {
                	LOG(ERROR) << "***Error in TCT time interpolation.***";
                	return EXIT_FAILURE;
            	}
            	interpolateFlag = 1;
		    	this->last_index=i;

            	break;
        	}
    		}
	
	}

	}

   }
   else //THERE IS RESET IN TCT
   {
	   if(instrumentTime < vec_instrumentTime[0]&& instrumentTime > vec_instrumentTime[endIndex])
	   {
		   nrecOutTCT++;
		   if(this->last_index==0)
		   {
		        if (find_interpolated_value(instrumentTime, vec_instrumentTime[0], vec_instrumentTime[1],
        	        vec_spsTime[0], vec_spsTime[1], intTime)) {
            		LOG(ERROR) << "***Error in TCT time interpolation.***";
 		           	return EXIT_FAILURE;
        		}
			   
		   		this->last_index=0;
		   }
		   else //if (this->lastindex==endIndex)
		   {
        		if (find_interpolated_value(instrumentTime, vec_instrumentTime[endIndex-1], vec_instrumentTime[endIndex],
                	vec_spsTime[endIndex-1], vec_spsTime[endIndex], intTime)) {
            		LOG(ERROR) << "***Error in TCT time interpolation.***";
           	 		return EXIT_FAILURE;
        		}

				this->last_index=endIndex;
		   }
			/*
		   else
		   {
			   LOG(ERROR)<<"This is upexpected ERROR. ";
			   LOG(ERROR)<<"CHECK WHAT IS WRONG";
			   RETURN(EXIT_FAILURE);
		   }
			*/
	   }
	   else if (instrumentTime < vec_instrumentTime[resetIndex[0]]) 
	   {

                if (find_interpolated_value(instrumentTime, vec_instrumentTime[resetIndex[0]-1]-1048576.0, vec_instrumentTime[resetIndex[0]],
                        vec_spsTime[resetIndex[0]-1], vec_spsTime[resetIndex[0]], intTime)) {
                    LOG(ERROR) << "***Error in TCT time interpolation.***";
                    return EXIT_FAILURE;
                }
                interpolateFlag = 1;
                this->last_index=resetIndex[0]-1;			

	   }
	   else if (instrumentTime > vec_instrumentTime[resetIndex[0]-1])
	   {
		   if (find_interpolated_value(instrumentTime, vec_instrumentTime[resetIndex[0]-1], vec_instrumentTime[resetIndex[0]]+1048576.0,
					   vec_spsTime[resetIndex[0]-1], vec_spsTime[resetIndex[0]], intTime)) {
                    LOG(ERROR) << "***Error in TCT time interpolation.***";
	                     return EXIT_FAILURE;
		                 }
                interpolateFlag = 1;
	                 this->last_index=resetIndex[0]-2;	   
	   }
	   else
	   {		
	    startIndex=this->last_index;

    	for (i = startIndex; i < endIndex; i++) {
        	if (instrumentTime >= vec_instrumentTime[i] && instrumentTime < vec_instrumentTime[i + 1]) {
            	if (find_interpolated_value(instrumentTime, vec_instrumentTime[i], vec_instrumentTime[i + 1],
                	    vec_spsTime[i], vec_spsTime[i + 1], intTime)) {
                	LOG(ERROR) << "***Error in TCT time interpolation.***";
                	return EXIT_FAILURE;
            	}
            	interpolateFlag = 1;
            	this->last_index=i;

            	break;
        	}
    	}

	    if(interpolateFlag!=1)
    	{
        	startIndex=0;

	        for (i = startIndex; i < endIndex; i++) {
    	        if (instrumentTime >= vec_instrumentTime[i] && instrumentTime < vec_instrumentTime[i + 1]) {
        	        if (find_interpolated_value(instrumentTime, vec_instrumentTime[i], vec_instrumentTime[i + 1],
            	            vec_spsTime[i], vec_spsTime[i + 1], intTime)) {
                	    LOG(ERROR) << "***Error in TCT time interpolation.***";
                    	return EXIT_FAILURE;
                	}
                	interpolateFlag = 1;
                	this->last_index=i;

                	break;
            	}
            	}
   
    	}

	   }
   }

		
    return EXIT_SUCCESS;
}
//Attitude
//Comparison function to check sorting
//These template function will be used to check whether attitude data is sorted or not.

bool whether_ascending_att(attStruct a, attStruct b) {
    return a.time < b.time;
}

bool whether_descending_att(attStruct a, attStruct b) {
    return a.time > b.time;
}

attStruct::attStruct(){
    time=0.0;
    rollRA=0.0;
    rollDEC=0.0;
    rollRot=0.0;
    qSat.resize(4,0.0);
}
Attitude::Attitude() {
    attFilename = "";
    attFileRead = false;
    nrows=0;
    lastRecord=0;
    fileSorted=false;
    headerReadFlag=false;
    tstarti=0;
    tstartf=0.0;
    tstopi=0;
    tstopf=0.0;
    RAPnt=0.0;
    DECPnt=0.0;
    obsID="";
    minSamplingInterval=0.0;
    maxSamplingInterval=0.0;
    dataMinTime=0.0;
    dataMaxTime=0.0;
}

int Attitude::read_attitude_file(string attFilename){
    int status=0;
    fitsfile *fatt;
    ErrorHandler errHandler;
    long i=0, nelements=0;
    vector <double> vecTime;
    vector < vector <double> > vecQsat;
    vector <double> vecRollRa;
    vector <double> vecRollDEC;
    vector <double> vecRollRot;
    
    fits_open_file(&fatt, (char*) attFilename.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening fits file: " + attFilename;
        throw errHandler;
    }

    //Reading header keywords [PRIMARY HEADER]
    try {
        readFitsKey(fatt, TLONG, "TSTARTI", &tstarti, NULL);
        readFitsKey(fatt, TDOUBLE, "TSTARTF", &tstartf, NULL);
        readFitsKey(fatt, TLONG, "TSTOPI", &tstopi, NULL);
        readFitsKey(fatt, TDOUBLE, "TSTOPF", &tstopf, NULL);
        readFitsKey(fatt, TDOUBLE, "RA_PNT", &RAPnt, NULL);
        readFitsKey(fatt, TDOUBLE, "DEC_PNT", &DECPnt, NULL);
        readFitsKeyStr(fatt, "OBS_ID", &obsID, NULL);
        headerReadFlag = true;
    } catch (ErrorHandler errHandler) {
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    
    fits_movnam_hdu(fatt, BINARY_TBL, "ATTITUDE", 0, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in moving to HDU ATTITUDE in fits file: " + attFilename;
        throw errHandler;
    }
    
    if(read_fits_column(fatt, "TIME", TDOUBLE, 1, 1, -1, vecTime)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column TIME of attitude file: " + attFilename;
        throw errHandler;
    }
    if(read_fits_column(fatt, "ROLL_RA", TDOUBLE, 1, 1, -1, vecRollRa)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column ROLL_RA of attitude file: " + attFilename;
        throw errHandler;
    }
    if(read_fits_column(fatt, "ROLL_DEC", TDOUBLE, 1, 1, -1, vecRollDEC)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column ROLL_DEC of attitude file: " + attFilename;
        throw errHandler;
    }
    if(read_fits_column(fatt, "ROLL_ROT", TDOUBLE, 1, 1, -1, vecRollRot)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column ROLL_ROT of attitude file: " + attFilename;
        throw errHandler;
    }
    if(read_fits_array_column(fatt,"Q_SAT",TDOUBLE, 1, 1, -1, vecQsat)){
        errHandler.severity=errERROR;
        errHandler.errorMsg = "Error in reading column Q_SAT of attitude file: " + attFilename;
        throw errHandler;       
    }
    
    //Assigning data to attitude structure
    nelements = vecTime.size();
    this->att.resize(nelements);
    for(i=0; i<nelements; i++){
        (this->att[i]).time = vecTime[i];
        (this->att[i]).qSat = vecQsat[i];
        (this->att[i]).rollRA = vecRollRa[i];
        (this->att[i]).rollDEC = vecRollDEC[i];
        (this->att[i]).rollRot = vecRollRot[i];
    }
    //Checking whether the data is sorted by time or not
    if (adjacent_find(att.begin(), att.end(), whether_descending_att) == att.end()) {
        fileSorted = true;
    }
    
    this->attFilename = attFilename;
    attFileRead = true;
    this->nrows = att.size();
    try{
        get_min_max_sampling_interval();
        get_min_max_dataTime();
    } catch(ErrorHandler errHandler){
        logError(errHandler);
    }
    
    return EXIT_SUCCESS;
    
}

attStruct Attitude::get_interpolated_attitude(double time, bool extrapolateFlag) {
    long i = 0;
    long indexMin = 0;
    long indexMax = 0;
    ErrorHandler errHandler;
    bool attTimeFlag = false;
    attStruct attOut;
    double q0, q1, q2, q3, rollRA, rollDEC, rollRot;

    if (fileSorted == false) {
        //sort orbit data
        sort(att.begin(), att.end());
        fileSorted = true;
    } else if (fileSorted == true) {
        for (i = lastRecord; i < nrows - 1; i++) {
            if (att[i + 1].time > time && att[i].time <= time) {
                indexMin = i;
                indexMax = i + 1;
                this->lastRecord = indexMin;
                //DLOG(INFO) << "Orb Time " << time << " lies between " << indexMin << " and " << indexMax;
                attTimeFlag = true;
                break;
            }
        }
    }

    if (attTimeFlag == true) {
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].qSat[0],
                (double) att[indexMax].qSat[0],q0 );
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].qSat[1],
                (double) att[indexMax].qSat[1], q1);
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].qSat[2],
                (double) att[indexMax].qSat[2], q2);
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].qSat[3],
                (double) att[indexMax].qSat[3], q3);
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].rollRA,
                (double) att[indexMax].rollRA, rollRA);
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].rollDEC,
                (double) att[indexMax].rollDEC, rollDEC);
        find_interpolated_value(time, (double) att[indexMin].time,
                (double) att[indexMax].time, (double) att[indexMin].rollRot,
                (double) att[indexMax].rollRot, rollRot);


        attOut.qSat[0] = q0; 
        attOut.qSat[1] = q1; 
        attOut.qSat[2] = q2; 
        attOut.qSat[3] = q3; 
        attOut.rollRA = rollRA;
        attOut.rollDEC = rollDEC;
        attOut.rollRot = rollRot;
    }

    if (attTimeFlag == false) {
        LOG(WARNING) << " Attitude EXTRAPOLATION REQUIRED " << setprecision(20) << time;
        errHandler.severity = errERROR;
        errHandler.errorStatus = EXTRAPOLATION_REQUIRED;
        errHandler.errorMsg = "Time " + itoa(time, 20) + " is outside attitude file time range. Extrapolation required.";
        throw errHandler;
    }


    return attOut;
}

void Attitude::display_header_keywords() {
    if (attFileRead) {
        LOG(INFO) << "-------------------------------------------";
        LOG(INFO) << "Attitude File           : " << attFilename;
        LOG(INFO) << "All header keywords read: " << headerReadFlag;
        LOG(INFO) << "TSTARTI                 : " << tstarti;
        LOG(INFO) << "TSTARTF                 : " << tstartf;
        LOG(INFO) << "TSTOPI                  : " << tstopi;
        LOG(INFO) << "TSTOPF                  : " << tstopf;
        LOG(INFO) << "RA Pointing (degrees)   : " << RAPnt;
        LOG(INFO) << "DEC Pointing (degrees)  : " << DECPnt;
        LOG(INFO) << "Observation ID          : " << obsID;
        LOG(INFO) << "Minimum sampling        : " << minSamplingInterval;
        LOG(INFO) << "Maximum sampling        : " << maxSamplingInterval;
        LOG(INFO) << setprecision(20) << "Minimum data time       : " << dataMinTime;
        LOG(INFO) << setprecision(20) << "Maximum data time       : " << dataMaxTime;
    }
}

void Attitude::get_min_max_sampling_interval(){
    vector <double> samplingInterval;
    ErrorHandler errHandler;
    long i=0;
    
    samplingInterval.resize(nrows-1, 0.0);
    if(nrows>0){
        for(i=1; i<nrows; i++){
            samplingInterval[i-1]=att[i].time - att[i-1].time;
        }
        compute_min_and_max(samplingInterval.begin(), samplingInterval.end(),
            minSamplingInterval, maxSamplingInterval);
    }
    else{
        errHandler.severity = errWARNING;
        errHandler.errorMsg = "No attitude data to find sampling interval.";
        throw errHandler;
    }  
}

void Attitude::get_min_max_dataTime() {
    if (attFileRead && fileSorted) {
        dataMinTime = att[0].time;
        dataMaxTime = att[nrows - 1].time;
    } else {
        LOG(INFO) << "Sorting attitude table...";
        sort(att.begin(), att.end());
        fileSorted=true;
        dataMinTime = att[0].time;
        dataMaxTime = att[nrows - 1].time;
    }
}
//ORBIT FILE
//Comparison function to check sorting
//These template function will be used to check whether orbit data is sorted or not.

bool whether_ascending_orb(orbStruct a, orbStruct b) {
    return a.time < b.time;
}

bool whether_descending_orb(orbStruct a, orbStruct b) {
    return a.time > b.time;
}

orbStruct::orbStruct(){
    time=0.0;
    x=0.0;
    y=0.0;
    vX=0.0;
    vY=0.0;
    vZ=0.0;
    lat = 0.0;
    lon = 0.0;
    altitude = 0.0;
    semimajorAxis = 0.0;
    ecentricity = 0.0;
    inclination = 0.0;
    aop=0.0;
    raan=0.0;
    trueAnomaly =0.0;
}
Orbit::Orbit() {
    orbFilename = "";
    orbFileRead=false;
    nrows=0;
    lastRecord=0;
    fileSorted=false;
    headerReadFlag=false;
    tstarti=0;
    tstartf=0.0;
    tstopi=0;
    tstopf=0.0;
    RAPnt = 0.0;
    DECPnt = 0.0;
    obsID = "";
    minSamplingInterval = 0.0;
    maxSamplingInterval = 0.0;
    dataMinTime=0.0;
    dataMaxTime=0.0;
}

int Orbit::read_orbit_file(string orbFilename) {
    int status = 0;
    fitsfile *forb;
    ErrorHandler errHandler;
    long i = 0, nelements = 0;
    vector <double> vecTime;
    vector <float>  vecX;
    vector <float>  vecY;
    vector <float>  vecZ;
    vector <float>  vecvX;
    vector <float>  vecvY;
    vector <float>  vecvZ;
    vector <float>  vecLat;
    vector <float>  vecLon;
    vector <float>  vecAlt;
    vector <float>  vecSemiMajorAxis;
    vector <float>  vecEccentricity;
    vector <float>  vecInclination;
    vector <float>  vecAOP;
    vector <float>  vecRAAN;
    vector <float>  vecTrueAnomaly;

    fits_open_file(&forb, (char*) orbFilename.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening fits file: " + orbFilename;
        throw errHandler;
    }
    
    //Reading header keywords [PRIMARY HEADER]
    try{
        readFitsKey(forb, TLONG, "TSTARTI", &tstarti, NULL);
        readFitsKey(forb, TDOUBLE, "TSTARTF", &tstartf, NULL);
        readFitsKey(forb, TLONG, "TSTOPI", &tstopi, NULL);
        readFitsKey(forb, TDOUBLE, "TSTOPF", &tstopf, NULL);
        readFitsKey(forb, TDOUBLE, "RA_PNT", &RAPnt, NULL);
        readFitsKey(forb, TDOUBLE, "DEC_PNT", &DECPnt, NULL);
        readFitsKeyStr(forb, "OBS_ID", &obsID, NULL);
        headerReadFlag=true;
    } catch(ErrorHandler errHandler){
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    
    fits_movnam_hdu(forb, BINARY_TBL, "ORBIT", 0, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in moving to HDU ORBIT in fits file: " + orbFilename;
        throw errHandler;
    }
    
    //TIME
    if (read_fits_column(forb, "TIME", TDOUBLE, 1, 1, -1, vecTime)) {
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column TIME of orbit file: " + orbFilename;
        throw errHandler;
    }
    //X
    if (read_fits_column(forb, "x", TFLOAT, 1, 1, -1, vecX)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column x of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Y
    if (read_fits_column(forb, "Y", TFLOAT, 1, 1, -1, vecY)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column Y of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Z
    if (read_fits_column(forb, "Z", TFLOAT, 1, 1, -1, vecZ)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column Z of orbit file: " + orbFilename;
        throw errHandler;
    }
    //vX
    if (read_fits_column(forb, "vx", TFLOAT, 1, 1, -1, vecvX)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column vx of orbit file: " + orbFilename;
        throw errHandler;
    }
    //vY
    if (read_fits_column(forb, "vy", TFLOAT, 1, 1, -1, vecvY)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column vy of orbit file: " + orbFilename;
        throw errHandler;
    }
    //vz
    if (read_fits_column(forb, "vz", TFLOAT, 1, 1, -1, vecvZ)){
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading column vz of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Latitude
    if(read_fits_column(forb, "Latitude", TFLOAT, 1, 1, -1, vecLat)){
        errHandler.severity=errERROR;
        errHandler.errorMsg = "Error in reading column Latitude of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Longitude
    if(read_fits_column(forb, "Longitude", TFLOAT, 1, 1, -1, vecLon)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Longitude of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Longitude
    if(read_fits_column(forb, "Longitude", TFLOAT, 1, 1, -1, vecLon)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Longitude of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Altitude
    if(read_fits_column(forb, "Altitude", TFLOAT, 1, 1, -1, vecAlt)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Altitude of orbit file: " + orbFilename;
        throw errHandler;
    }
    //SemiMajor_Axis
    if(read_fits_column(forb, "SemiMajor_Axis", TFLOAT, 1, 1, -1, vecSemiMajorAxis)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column SemiMajor_Axis of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Eccentricity
    if(read_fits_column(forb, "Ecentricity", TFLOAT, 1, 1, -1, vecEccentricity)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Ecentricity of orbit file: " + orbFilename;
        throw errHandler;
    }
    //Inclination
    if(read_fits_column(forb, "Inclination", TFLOAT, 1, 1, -1, vecInclination)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Inclination of orbit file: " + orbFilename;
        throw errHandler;
    }
    //ArgOfPerigee
    if(read_fits_column(forb, "ArgOfPerigee", TFLOAT, 1, 1, -1, vecAOP)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column ArgOfPerigee of orbit file: " + orbFilename;
        throw errHandler;
    }
    //RAAN
    if(read_fits_column(forb, "RAAN", TFLOAT, 1, 1, -1, vecRAAN)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column RAAN of orbit file: " + orbFilename;
        throw errHandler;
    }
    //TrueAnomaly
    if(read_fits_column(forb, "TrueAnomaly", TFLOAT, 1, 1, -1, vecTrueAnomaly)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column TrueAnomaly of orbit file: " + orbFilename;
        throw errHandler;
    }

    //Assigning data to attitude structure
    nelements = vecTime.size();
    this->orb.resize(nelements);
    for (i = 0; i < nelements; i++) {
        (this->orb[i]).time = vecTime[i];
        (this->orb[i]).x = vecX[i];
        (this->orb[i]).y = vecY[i];
        (this->orb[i]).z = vecZ[i];
        (this->orb[i]).vX = vecvX[i];
        (this->orb[i]).vY = vecvY[i];
        (this->orb[i]).vZ = vecvZ[i];
        (this->orb[i]).lat = vecLat[i];
        (this->orb[i]).lon = vecLon[i];
        (this->orb[i]).altitude = vecAlt[i];
        (this->orb[i]).semimajorAxis = vecSemiMajorAxis[i];
        (this->orb[i]).ecentricity = vecEccentricity[i];
        (this->orb[i]).inclination = vecInclination[i];
        (this->orb[i]).aop = vecAOP[i];
        (this->orb[i]).raan = vecRAAN[i];
        (this->orb[i]).trueAnomaly = vecTrueAnomaly[i];
    }

    //Checking whether the data is sorted by time or not
    if (adjacent_find(orb.begin(), orb.end(), whether_descending_orb) == orb.end()) {
        fileSorted = true;
    }
    
    this->orbFilename = orbFilename;
    orbFileRead = true; //attitude file flag
    this->nrows = orb.size();
    try {
        get_min_max_sampling_interval();
        get_min_max_dataTime();
    } catch (ErrorHandler errHandler) {
        logError(errHandler);
    }
    
    return EXIT_SUCCESS;
    
}

void Orbit::sort_orbit_data() {
    if(orbFileRead==true){
        sort(orb.begin(), orb.end());
        fileSorted=true;
    }
}

void Orbit::display_header_keywords() {
    if (orbFileRead) {
        LOG(INFO) << "-------------------------------------------";
        LOG(INFO) << "Orbit File              : " << orbFilename;
        LOG(INFO) << "All header keywords read: " << headerReadFlag;
        LOG(INFO) << "TSTARTI                 : " << tstarti;
        LOG(INFO) << "TSTARTF                 : " << tstartf;
        LOG(INFO) << "TSTOPI                  : " << tstopi;
        LOG(INFO) << "TSTOPF                  : " << tstopf;
        LOG(INFO) << "RA Pointing (degrees)   : " << RAPnt;
        LOG(INFO) << "DEC Pointing (degrees)  : " << DECPnt;
        LOG(INFO) << "Observation ID          : " << obsID;
        LOG(INFO) << "Minimum sampling        : " << minSamplingInterval;
        LOG(INFO) << "Maximum sampling        : " << maxSamplingInterval;
        LOG(INFO) << setprecision(20) << "Minimum data time       : " << dataMinTime;
        LOG(INFO) << setprecision(20) << "Maximum data time       : " << dataMaxTime;
    }
}

void Orbit::get_min_max_sampling_interval() {
    vector <double> samplingInterval;
    ErrorHandler errHandler;
    long i = 0;

    samplingInterval.resize(nrows - 1, 0.0);
    if (nrows > 0) {
        for (i = 1; i < nrows; i++) {
            samplingInterval[i - 1] = orb[i].time - orb[i - 1].time;
        }
        compute_min_and_max(samplingInterval.begin(), samplingInterval.end(),
                minSamplingInterval, maxSamplingInterval);
    } else {
        errHandler.severity = errWARNING;
        errHandler.errorMsg = "No orbit data to find sampling interval.";
        throw errHandler;
    }
}

void Orbit::get_min_max_dataTime() {
    if (orbFileRead && fileSorted) {
        dataMinTime = orb[0].time;
        dataMaxTime = orb[nrows - 1].time;
    } else {
        LOG(INFO) << "Sorting lbt table...";
        sort(orb.begin(), orb.end());
        fileSorted = true;
        dataMinTime = orb[0].time;
        dataMaxTime = orb[nrows - 1].time;
    }
}

orbStruct Orbit::get_interpolated_orbit(double time, bool extrapolateFlag) {
    long i=0;
    long indexMin=0;
    long indexMax=0;
    ErrorHandler errHandler;
    bool orbitTimeFlag=false;
    orbStruct orbOut;
    double x, y, z, vX, vY, vZ, lat, lon, altitude;

    if (fileSorted == false) {
        //sort orbit data
        sort(orb.begin(), orb.end());
        fileSorted = true;
    } else if (fileSorted == true) {
        for (i = lastRecord; i < nrows - 1; i++) {
            if (orb[i + 1].time > time && orb[i].time <= time) {
                indexMin = i;
                indexMax = i + 1;
                this->lastRecord = indexMin;
                //DLOG(INFO) << "Orb Time " << time << " lies between " << indexMin << " and " << indexMax;
                orbitTimeFlag = true;
                break;
            }
        }
    }
    
    if(orbitTimeFlag==true) {
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].x,
                (double) orb[indexMax].x, x);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].y,
                (double) orb[indexMax].y, y);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].z,
                (double) orb[indexMax].z, z);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].vX,
                (double) orb[indexMax].vX, vX);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].vY,
                (double) orb[indexMax].vY, vY);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].vZ,
                (double) orb[indexMax].vZ, vZ);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].lat,
                (double) orb[indexMax].lat, lat);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].lon,
                (double) orb[indexMax].lon, lon);
        find_interpolated_value(time, (double) orb[indexMin].time,
                (double) orb[indexMax].time, (double) orb[indexMin].altitude,
                (double) orb[indexMax].altitude, altitude);

        orbOut.x = (float) x;
        orbOut.y = (float) y;
        orbOut.z = (float) z;
        orbOut.vX = (float) vX;
        orbOut.vY = (float) vY;
        orbOut.vZ = (float) vZ;
        orbOut.lat = (float) lat;
        orbOut.lon = (float) lon;
        orbOut.altitude = (float) altitude;
    }
    
    if (orbitTimeFlag==false) {
        LOG(WARNING) << " ORB EXTRAPOLATION REQUIRED " << setprecision(20) << time;
        errHandler.severity = errERROR;
        errHandler.errorStatus = EXTRAPOLATION_REQUIRED;
        errHandler.errorMsg = "Time " + itoa(time, 20)+ " is outside orbit file time range. Extrapolation required.";
        throw errHandler;
    }
    

    return orbOut;
}


//LBT FILE
//Comparison function to check sorting
//These template function will be used to check whether lbt data is sorted or not.

bool whether_ascending_lbt(lbtStruct a, lbtStruct b) {
    return a.time < b.time;
}

bool whether_descending_lbt(lbtStruct a, lbtStruct b) {
    return a.time > b.time;
}

lbtStruct::lbtStruct() {
    time = 0.0;
    //Quadrant A
    q1modeID = 5;
    q1Pos5VMonitor = 0.0;
    q1CZTCounter = 0;
    q1Temperature1 = 0.0;
    q1VetoHVMonitor = 0.0;
    q1Pos2dot5VMonitor = 0.0;
    q1CZTHVMonitor = 0.0;
    q1VetoCounter = 0;
    q1vetoLLD = 0.0;
    q1DVDD = 0.0;
    q1AlphaCounter = 0;
    //Quadrant B
    q2modeID = 5;
    q2Pos5VMonitor = 0.0;
    q2CZTCounter = 0;
    q2Temperature1 = 0.0;
    q2VetoHVMonitor = 0.0;
    q2Pos2dot5VMonitor = 0.0;
    q2CZTHVMonitor = 0.0;
    q2VetoCounter = 0;
    q2vetoLLD = 0.0;
    q2DVDD = 0.0;
    q2AlphaCounter = 0;
    //Quadrant C
    q3modeID = 5;
    q3Pos5VMonitor = 0.0;
    q3CZTCounter = 0;
    q3Temperature1 = 0.0;
    q3VetoHVMonitor = 0.0;
    q3Pos2dot5VMonitor = 0.0;
    q3CZTHVMonitor = 0.0;
    q3VetoCounter = 0;
    q3vetoLLD = 0.0;
    q3DVDD = 0.0;
    q3AlphaCounter = 0;
    //Quadrant D
    q4modeID = 5;
    q4Pos5VMonitor = 0.0;
    q4CZTCounter = 0;
    q4Temperature1 = 0.0;
    q4VetoHVMonitor = 0.0;
    q4Pos2dot5VMonitor = 0.0;
    q4CZTHVMonitor = 0.0;
    q4VetoCounter = 0;
    q4vetoLLD = 0.0;
    q4DVDD = 0.0;
    q4AlphaCounter = 0;

    cztMemLvl=0;
    cztBootPage=0;
    cztPEErrCnt=0;
    cztLastCmd=0;
    cpmRate=0;

}

bool is_event_modeID(unsigned char modeID){
    bool eventModeFlag=false;
    if(modeID>=0 || modeID<=7){
        eventModeFlag = true;
    }
    return eventModeFlag;
}
LBT::LBT(){
    lbtFileRead=false;
    lbtFilename="";
    nrows=0;
    lastRecord=0;
    fileSorted=false;
    headerReadFlag=false;
    tstarti=0;
    tstartf=0.0;
    tstopi=0;
    tstopf=0.0;
    RAPnt=0.0;
    DECPnt=0.0;
    obsID="";
    minSamplingInterval=0.0;
    maxSamplingInterval=0.0;
}

int LBT::read_lbt_file(string lbtFilename) {
    int status = 0;
    fitsfile *flbt;
    ErrorHandler errHandler;
    long i = 0, nelements = 0;
    vector <double> vecTime;
    //Quadrant 1
    vector<unsigned char> vecq1modeID;
    vector<float> vecq1Pos5VMonitor;
    vector<unsigned long> vecq1CZTCounter;
    vector<float> vecq1Temperature1;
    vector<float> vecq1VetoHVMonitor;
    vector<float> vecq1Pos2dot5VMonitor;
    vector<float> vecq1CZTHVMonitor;
    vector<unsigned long> vecq1VetoCounter;
    vector<float> vecq1vetoLLD;
    vector<float> vecq1DVDD;
    vector<unsigned long> vecq1AlphaCounter;
    //Quadrant 2
    vector<unsigned char> vecq2modeID;
    vector<float> vecq2Pos5VMonitor;
    vector<unsigned long> vecq2CZTCounter;
    vector<float> vecq2Temperature1;
    vector<float> vecq2VetoHVMonitor;
    vector<float> vecq2Pos2dot5VMonitor;
    vector<float> vecq2CZTHVMonitor;
    vector<unsigned long> vecq2VetoCounter;
    vector<float> vecq2vetoLLD;
    vector<float> vecq2DVDD;
    vector<unsigned long> vecq2AlphaCounter;
    //Quadrant 3
    vector<unsigned char> vecq3modeID;
    vector<float> vecq3Pos5VMonitor;
    vector<unsigned long> vecq3CZTCounter;
    vector<float> vecq3Temperature1;
    vector<float> vecq3VetoHVMonitor;
    vector<float> vecq3Pos2dot5VMonitor;
    vector<float> vecq3CZTHVMonitor;
    vector<unsigned long> vecq3VetoCounter;
    vector<float> vecq3vetoLLD;
    vector<float> vecq3DVDD;
    vector<unsigned long> vecq3AlphaCounter;
    //Quadrant 4
    vector<unsigned char> vecq4modeID;
    vector<float> vecq4Pos5VMonitor;
    vector<unsigned long> vecq4CZTCounter;
    vector<float> vecq4Temperature1;
    vector<float> vecq4VetoHVMonitor;
    vector<float> vecq4Pos2dot5VMonitor;
    vector<float> vecq4CZTHVMonitor;
    vector<unsigned long> vecq4VetoCounter;
    vector<float> vecq4vetoLLD;
    vector<float> vecq4DVDD;
    vector<unsigned long> vecq4AlphaCounter;

    vector<unsigned char> veccztMemLvl;
    vector<unsigned char> veccztBootPage;
    vector<unsigned char> veccztPEErrCnt;
    vector<unsigned long> veccztLastCmd;
    vector<unsigned long> veccpmRate;

    fits_open_file(&flbt, (char*) lbtFilename.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in opening fits file: " + lbtFilename;
        throw errHandler;
    }
    
    //Reading header keywords [PRIMARY HEADER]
    try {
        readFitsKey(flbt, TLONG, "TSTARTI", &tstarti, NULL);
        readFitsKey(flbt, TDOUBLE, "TSTARTF", &tstartf, NULL);
        readFitsKey(flbt, TLONG, "TSTOPI", &tstopi, NULL);
        readFitsKey(flbt, TDOUBLE, "TSTOPF", &tstopf, NULL);
        readFitsKey(flbt, TDOUBLE, "RA_PNT", &RAPnt, NULL);
        readFitsKey(flbt, TDOUBLE, "DEC_PNT", &DECPnt, NULL);
        readFitsKeyStr(flbt, "OBS_ID", &obsID, NULL);
        headerReadFlag = true;
    } catch (ErrorHandler errHandler) {
        errHandler.severity = errWARNING;
        logError(errHandler);
    }
    fits_movnam_hdu(flbt, BINARY_TBL, "CZT-LBTHK", 0, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in moving to HDU CZT-LBTHK in fits file: " + lbtFilename;
        throw errHandler;
    }
    
    //TIME
    if (read_fits_column(flbt, "TIME", TDOUBLE, 1, 1, -1, vecTime)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column TIME of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1_MODEID
    if (read_fits_column(flbt, "Q1_MODEID", TBYTE, 1, 1, -1, vecq1modeID)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1_MODEID of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2_MODEID
    if (read_fits_column(flbt, "Q2_MODEID", TBYTE, 1, 1, -1, vecq2modeID)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2_MODEID of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3_MODEID
    if (read_fits_column(flbt, "Q3_MODEID", TBYTE, 1, 1, -1, vecq3modeID)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3_MODEID of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4_MODEID
    if (read_fits_column(flbt, "Q4_MODEID", TBYTE, 1, 1, -1, vecq4modeID)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4_MODEID of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1_POS_5V_MONITOR
    if (read_fits_column(flbt, "Q1_POS_5V_MONITOR", TFLOAT, 1, 1, -1, vecq1Pos5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1_POS_5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2_POS_5V_MONITOR
    if (read_fits_column(flbt, "Q2_POS_5V_MONITOR", TFLOAT, 1, 1, -1, vecq2Pos5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2_POS_5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3_POS_5V_MONITOR
    if (read_fits_column(flbt, "Q3_POS_5V_MONITOR", TFLOAT, 1, 1, -1, vecq3Pos5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3_POS_5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4_POS_5V_MONITOR
    if (read_fits_column(flbt, "Q4_POS_5V_MONITOR", TFLOAT, 1, 1, -1, vecq4Pos5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4_POS_5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1_POS_2DOT5V_MONITOR
    if (read_fits_column(flbt, "Q1_POS_2DOT5V_MONITOR", TFLOAT, 1, 1, -1, vecq1Pos2dot5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1_POS_2DOT5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2_POS_2DOT5V_MONITOR
    if (read_fits_column(flbt, "Q2_POS_2DOT5V_MONITOR", TFLOAT, 1, 1, -1, vecq2Pos2dot5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2_POS_2DOT5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3_POS_2DOT5V_MONITOR
    if (read_fits_column(flbt, "Q3_POS_2DOT5V_MONITOR", TFLOAT, 1, 1, -1, vecq3Pos2dot5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3_POS_2DOT5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4_POS_2DOT5V_MONITOR
    if (read_fits_column(flbt, "Q4_POS_2DOT5V_MONITOR", TFLOAT, 1, 1, -1, vecq4Pos2dot5VMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4_POS_2DOT5V_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1_CZT_COUNTER
    if (read_fits_column(flbt, "Q1_CZT_COUNTER", TULONG, 1, 1, -1, vecq1CZTCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1_CZT_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2_CZT_COUNTER
    if (read_fits_column(flbt, "Q2_CZT_COUNTER", TULONG, 1, 1, -1, vecq2CZTCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2_CZT_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3_CZT_COUNTER
    if (read_fits_column(flbt, "Q3_CZT_COUNTER", TULONG, 1, 1, -1, vecq3CZTCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3_CZT_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4_CZT_COUNTER
    if (read_fits_column(flbt, "Q4_CZT_COUNTER", TULONG, 1, 1, -1, vecq4CZTCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4_CZT_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1TEMPERATURE1
    if (read_fits_column(flbt, "Q1TEMPERATURE1", TFLOAT, 1, 1, -1, vecq1Temperature1)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1TEMPERATURE1 of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2TEMPERATURE1
    if (read_fits_column(flbt, "Q2TEMPERATURE1", TFLOAT, 1, 1, -1, vecq2Temperature1)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2TEMPERATURE1 of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3TEMPERATURE1
    if (read_fits_column(flbt, "Q3TEMPERATURE1", TFLOAT, 1, 1, -1, vecq3Temperature1)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3TEMPERATURE1 of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4TEMPERATURE1
    if (read_fits_column(flbt, "Q4TEMPERATURE1", TFLOAT, 1, 1, -1, vecq4Temperature1)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4TEMPERATURE1 of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1VETOHV_MONITOR
    if (read_fits_column(flbt, "Q1VETOHV_MONITOR", TFLOAT, 1, 1, -1, vecq1VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1VETOHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2VETOHV_MONITOR
    if (read_fits_column(flbt, "Q2VETOHV_MONITOR", TFLOAT, 1, 1, -1, vecq2VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2VETOHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3VETOHV_MONITOR
    if (read_fits_column(flbt, "Q3VETOHV_MONITOR", TFLOAT, 1, 1, -1, vecq3VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3VETOHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4VETOHV_MONITOR
    if (read_fits_column(flbt, "Q4VETOHV_MONITOR", TFLOAT, 1, 1, -1, vecq4VetoHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4VETOHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1CZTHV_MONITOR
    if (read_fits_column(flbt, "Q1CZTHV_MONITOR", TFLOAT, 1, 1, -1, vecq1CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1CZTHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2CZTHV_MONITOR
    if (read_fits_column(flbt, "Q2CZTHV_MONITOR", TFLOAT, 1, 1, -1, vecq2CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2CZTHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3CZTHV_MONITOR
    if (read_fits_column(flbt, "Q3CZTHV_MONITOR", TFLOAT, 1, 1, -1, vecq3CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3CZTHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4CZTHV_MONITOR
    if (read_fits_column(flbt, "Q4CZTHV_MONITOR", TFLOAT, 1, 1, -1, vecq4CZTHVMonitor)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4CZTHV_MONITOR of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1_VETOCOUNTER
    if (read_fits_column(flbt, "Q1_VETOCOUNTER", TULONG, 1, 1, -1, vecq1VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1VETO_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2_VETOCOUNTER
    if (read_fits_column(flbt, "Q2_VETOCOUNTER", TULONG, 1, 1, -1, vecq2VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2VETO_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3_VETOCOUNTER
    if (read_fits_column(flbt, "Q3_VETOCOUNTER", TULONG, 1, 1, -1, vecq3VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3VETO_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4_VETOCOUNTER
    if (read_fits_column(flbt, "Q4_VETOCOUNTER", TULONG, 1, 1, -1, vecq4VetoCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4VETO_COUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1VETOLLD
    if (read_fits_column(flbt, "Q1VETOLLD", TFLOAT, 1, 1, -1, vecq1vetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1VETOLLD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2VETOLLD
    if (read_fits_column(flbt, "Q2VETOLLD", TFLOAT, 1, 1, -1, vecq2vetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2VETOLLD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3VETOLLD
    if (read_fits_column(flbt, "Q3VETOLLD", TFLOAT, 1, 1, -1, vecq3vetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3VETOLLD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4VETOLLD
    if (read_fits_column(flbt, "Q4VETOLLD", TFLOAT, 1, 1, -1, vecq4vetoLLD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4VETOLLD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1DVDD
    if (read_fits_column(flbt, "Q1DVDD", TFLOAT, 1, 1, -1, vecq1DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1DVDD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2DVDD
    if (read_fits_column(flbt, "Q2DVDD", TFLOAT, 1, 1, -1, vecq2DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2DVDD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3DVDD
    if (read_fits_column(flbt, "Q3DVDD", TFLOAT, 1, 1, -1, vecq3DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3DVDD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4DVDD
    if (read_fits_column(flbt, "Q4DVDD", TFLOAT, 1, 1, -1, vecq4DVDD)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4DVDD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q1_ALPHACOUNTER
    if (read_fits_column(flbt, "Q1_ALPHACOUNTER", TULONG, 1, 1, -1, vecq1AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q1_ALPHACOUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q2_ALPHACOUNTER
    if (read_fits_column(flbt, "Q2_ALPHACOUNTER", TULONG, 1, 1, -1, vecq2AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q2_ALPHACOUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q3_ALPHACOUNTER
    if (read_fits_column(flbt, "Q3_ALPHACOUNTER", TULONG, 1, 1, -1, vecq3AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q3_ALPHACOUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //Q4_ALPHACOUNTER
    if (read_fits_column(flbt, "Q4_ALPHACOUNTER", TULONG, 1, 1, -1, vecq4AlphaCounter)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column Q4_ALPHACOUNTER of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //CZT_MEM_LVL
    if(read_fits_column(flbt, "CZT_MEM_LVL", TBYTE, 1, 1, -1, veccztMemLvl)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column CZT_MEM_LVL of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //CZT_BOOT_PAGE
    if(read_fits_column(flbt, "CZT_BOOT_PAGE", TBYTE, 1, 1, -1, veccztBootPage)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column CZT_BOOT_PAGE of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //CZT_PE_ERR_CNT
    if(read_fits_column(flbt, "CZT_PE_ERR_CNT", TBYTE, 1, 1, -1, veccztPEErrCnt)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column CZT_PE_ERR_CNT of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //CZT_LAST_CMD
    if(read_fits_column(flbt, "CZT_LAST_CMD", TULONG, 1, 1, -1, veccztLastCmd)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column CZT_LAST_CMD of lbt file: " + lbtFilename;
        throw errHandler;
    }
    //CPM_RATE
    if(read_fits_column(flbt, "CPM_RATE", TULONG, 1, 1, -1, veccpmRate)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error in reading column CPM_RATE of lbt file: " + lbtFilename;
        throw errHandler;
    }

    //Assigning data to attitude structure
    nelements = vecTime.size();
    this->lbt.resize(nelements);
    for (i = 0; i < nelements; i++) {
        (this->lbt[i]).time = vecTime[i];
        (this->lbt[i]).q1modeID = vecq1modeID[i];
        (this->lbt[i]).q2modeID = vecq2modeID[i];
        (this->lbt[i]).q3modeID = vecq3modeID[i];
        (this->lbt[i]).q4modeID = vecq4modeID[i];
        (this->lbt[i]).q1Pos5VMonitor = vecq1Pos5VMonitor[i];
        (this->lbt[i]).q1Pos2dot5VMonitor = vecq1Pos2dot5VMonitor[i];
        (this->lbt[i]).q1CZTCounter=vecq1CZTCounter[i];
        (this->lbt[i]).q1Temperature1 = vecq1Temperature1[i];
        (this->lbt[i]).q1VetoHVMonitor=vecq1VetoHVMonitor[i];
        (this->lbt[i]).q1CZTHVMonitor=vecq1CZTHVMonitor[i];
        (this->lbt[i]).q1VetoCounter=vecq1VetoCounter[i];
        (this->lbt[i]).q1vetoLLD=vecq1vetoLLD[i];
        (this->lbt[i]).q1DVDD=vecq1DVDD[i];
        (this->lbt[i]).q1AlphaCounter=vecq1AlphaCounter[i];
        (this->lbt[i]).q2Pos5VMonitor = vecq2Pos5VMonitor[i];
        (this->lbt[i]).q2Pos2dot5VMonitor = vecq2Pos2dot5VMonitor[i];
        (this->lbt[i]).q2CZTCounter=vecq2CZTCounter[i];
        (this->lbt[i]).q2Temperature1 = vecq2Temperature1[i];
        (this->lbt[i]).q2VetoHVMonitor=vecq2VetoHVMonitor[i];
        (this->lbt[i]).q2CZTHVMonitor=vecq2CZTHVMonitor[i];
        (this->lbt[i]).q2VetoCounter=vecq2VetoCounter[i];
        (this->lbt[i]).q2vetoLLD=vecq2vetoLLD[i];
        (this->lbt[i]).q2DVDD=vecq2DVDD[i];
        (this->lbt[i]).q2AlphaCounter=vecq2AlphaCounter[i];
        (this->lbt[i]).q3Pos5VMonitor = vecq3Pos5VMonitor[i];
        (this->lbt[i]).q3Pos2dot5VMonitor = vecq3Pos2dot5VMonitor[i];
        (this->lbt[i]).q3CZTCounter=vecq3CZTCounter[i];
        (this->lbt[i]).q3Temperature1 = vecq3Temperature1[i];
        (this->lbt[i]).q3VetoHVMonitor=vecq3VetoHVMonitor[i];
        (this->lbt[i]).q3CZTHVMonitor=vecq3CZTHVMonitor[i];
        (this->lbt[i]).q3VetoCounter=vecq3VetoCounter[i];
        (this->lbt[i]).q3vetoLLD=vecq3vetoLLD[i];
        (this->lbt[i]).q3DVDD=vecq3DVDD[i];
        (this->lbt[i]).q3AlphaCounter=vecq3AlphaCounter[i];
        (this->lbt[i]).q4Pos5VMonitor = vecq4Pos5VMonitor[i];
        (this->lbt[i]).q4Pos2dot5VMonitor = vecq4Pos2dot5VMonitor[i];
        (this->lbt[i]).q4CZTCounter=vecq4CZTCounter[i];
        (this->lbt[i]).q4Temperature1 = vecq4Temperature1[i];
        (this->lbt[i]).q4VetoHVMonitor=vecq4VetoHVMonitor[i];
        (this->lbt[i]).q4CZTHVMonitor=vecq4CZTHVMonitor[i];
        (this->lbt[i]).q4VetoCounter=vecq4VetoCounter[i];
        (this->lbt[i]).q4vetoLLD=vecq4vetoLLD[i];
        (this->lbt[i]).q4DVDD=vecq4DVDD[i];
        (this->lbt[i]).q4AlphaCounter=vecq4AlphaCounter[i];
        (this->lbt[i]).cztMemLvl=veccztMemLvl[i];
        (this->lbt[i]).cztBootPage = veccztBootPage[i];
        (this->lbt[i]).cztPEErrCnt= veccztPEErrCnt[i];
        (this->lbt[i]).cztLastCmd=veccztLastCmd[i];
        (this->lbt[i]).cpmRate=veccpmRate[i];
    }

    //Checking whether the data is sorted by time or not
    if (adjacent_find(lbt.begin(), lbt.end(), whether_descending_lbt) == lbt.end()) {
        fileSorted = true;
    }
    this->lbtFilename = lbtFilename;
    lbtFileRead = true; //attitude file flag
    nrows = lbt.size();
    try {
        get_min_max_sampling_interval();
        get_min_max_dataTime();
    } catch (ErrorHandler errHandler) {
        logError(errHandler);
    }
    return EXIT_SUCCESS;
}

void LBT::display_header_keywords() {
    if (lbtFileRead == true) {
        LOG(INFO) << "---------------------------------------------";
        LOG(INFO) << "LBT File                : " << lbtFilename;
        LOG(INFO) << "All header keywords read: " << headerReadFlag;
        LOG(INFO) << "TSTARTI                 : " << tstarti;
        LOG(INFO) << "TSTARTF                 : " << tstartf;
        LOG(INFO) << "TSTOPI                  : " << tstopi;
        LOG(INFO) << "TSTOPF                  : " << tstopf;
        LOG(INFO) << "RA Pointing (degrees)   : " << RAPnt;
        LOG(INFO) << "DEC Pointing (degrees)  : " << DECPnt;
        LOG(INFO) << "Observation ID          : " << obsID;
        LOG(INFO) << "Minimum sampling        : " << minSamplingInterval;
        LOG(INFO) << "Maximum sampling        : " << maxSamplingInterval;
        LOG(INFO) << setprecision(20) << "Minimum data time       : " << dataMinTime;
        LOG(INFO) << setprecision(20) << "Maximum data time       : " << dataMaxTime;
    }
}

void LBT::get_min_max_sampling_interval() {
    vector <double> samplingInterval;
    ErrorHandler errHandler;
    long i = 0;

    samplingInterval.resize(nrows - 1, 0.0);
    if (nrows > 0) {
        for (i = 1; i < nrows; i++) {
            samplingInterval[i - 1] = lbt[i].time - lbt[i - 1].time;
        }
        compute_min_and_max(samplingInterval.begin(), samplingInterval.end(),
                minSamplingInterval, maxSamplingInterval);
    } else {
        errHandler.severity = errWARNING;
        errHandler.errorMsg = "No orbit data to find sampling interval.";
        throw errHandler;
    }
}

void LBT::get_min_max_dataTime(){
    if(lbtFileRead && fileSorted){
        dataMinTime = lbt[0].time;
        dataMaxTime = lbt[nrows-1].time;
    } else {
        LOG(INFO) << "Sorting lbt table...";
        sort(lbt.begin(), lbt.end());
        dataMinTime = lbt[0].time;
        dataMaxTime = lbt[nrows - 1].time;
    }
}
lbtStruct LBT::get_lbt(double time, bool interpolateFlag){
    long i=0;
    long indexMin=0;
    long indexMax=0;
    ErrorHandler errHandler;
    bool lbtTimeFlag=false;
    lbtStruct lbtOut;
    unsigned char modeIDfirst;
    unsigned char modeIDsecond;
    
    if(fileSorted==false){
        //sort lbt data
        sort(lbt.begin(), lbt.end());
        fileSorted=true;
    } 
    if(fileSorted==true){
        for(i=lastRecord; i<nrows-1; i++){
		//cout<<"LBT TIme:"<<lbt[i].time <<" Time:"<<time;
            if(lbt[i+1].time > time && lbt[i].time <= time){
                indexMin=i;
                indexMax= i+1;
                this->lastRecord=indexMin;
                //DLOG(INFO) << "LBT time " << time << " lies between " << indexMin << " and " << indexMax;
                lbtTimeFlag = true;
                break;
            }
        }
    }
    
    if(lbtTimeFlag==true){
        //Getting parameters for quadrant A
        modeIDfirst = lbt[indexMin].q1modeID;
        modeIDsecond = lbt[indexMax].q1modeID;
        if(is_event_modeID(modeIDfirst) && !is_event_modeID(modeIDsecond)){
            lbtOut = lbt[indexMax];
        } else {
            lbtOut = lbt[indexMin];
        }
        //Getting parameters for quadrant B
        modeIDfirst = lbt[indexMin].q2modeID;
        modeIDsecond = lbt[indexMax].q2modeID;
        if(is_event_modeID(modeIDfirst) && !is_event_modeID(modeIDsecond)){
            lbtOut = lbt[indexMax];
        } else {
            lbtOut = lbt[indexMin];
        }
        //Getting parameters for quadrant C
        modeIDfirst = lbt[indexMin].q3modeID;
        modeIDsecond = lbt[indexMax].q3modeID;
        if(is_event_modeID(modeIDfirst) && !is_event_modeID(modeIDsecond)){
            lbtOut = lbt[indexMax];
        } else {
            lbtOut = lbt[indexMin];
        }
        //Getting parameters for quadrant D
        modeIDfirst = lbt[indexMin].q4modeID;
        modeIDsecond = lbt[indexMax].q4modeID;
        if(is_event_modeID(modeIDfirst) && !is_event_modeID(modeIDsecond)){
            lbtOut = lbt[indexMax];
        } else {
            lbtOut = lbt[indexMin];
        }
    } else if (lbtTimeFlag==false) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = EXTRAPOLATION_REQUIRED;
        errHandler.errorMsg = "Time is outside lbt file time range. Extrapolation required.";
        throw errHandler;
    }
    
    return lbtOut;
}


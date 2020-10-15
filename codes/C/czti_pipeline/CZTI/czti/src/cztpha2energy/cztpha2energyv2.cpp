#include "cztpha2energyv2.h"

Cztpha2energy::Cztpha2energy() {
	  
	    strcpy(modulename, "cztpha2energy_v");
	    strcat(modulename, VERSION);
}

void Cztpha2energy::display(){
	   
	    LOG(INFO)  << "----------------------------------------------------------------------------";
	    LOG(INFO)  << "                           CZTPHA2ENERGY PARAMETERS                         ";
	    LOG(INFO)  << "----------------------------------------------------------------------------";
		LOG(INFO)  << "Modulename                   : "<<modulename;
	    LOG(INFO)  << "Input Event File             : " << inputEvtFile;
	    //LOG(INFO)  << "Ebounds File                 : " << eboundsFile;	  
	    LOG(INFO)  << "Output Event File            : " << outputEvtFile;
	    LOG(INFO)  << "Temperature Extension file   : " << tempExtFile;
		LOG(INFO)  << "Buffer                       : " << buffer;
	    LOG(INFO)  << "Clobber                      : " << clobber;
	    LOG(INFO)  << "History                      : " << history;
	    LOG(INFO)  << "----------------------------------------------------------------------------";   
}

int Cztpha2energy::read(int argc, char **argv) {
	    int status=0; //status variable
	    if (PIL_OK != (status= PILInit(argc, argv))) {
		LOG(ERROR) << "Error Initializing PIL";
		return (EXIT_FAILURE);
	    }
	    if (PIL_OK != (status = PILGetFname("par_infile", inputEvtFile))) {
		LOG(ERROR) << "Error Reading Input Event file";
		return (EXIT_FAILURE);
	    }
		
	    /*if (PIL_OK != (status = PILGetFname("par_eboundsfile", eboundsFile))) {
		LOG(ERROR) << "Error in reading Reading Ebounds File";
		return (EXIT_FAILURE);
	    }*/ //already commented
	    
	    /*if (PIL_OK != (status = PILGetFname("par_gofile", goFile))) {
		LOG(ERROR) << "Error in reading Reading Gain File";
		return (EXIT_FAILURE);
	    }*/ //already commented
	    
	    if (PIL_OK != (status = PILGetFname("par_outfile", outputEvtFile))) {
		LOG(ERROR) << "Error Reading Output Event File";
		return (EXIT_FAILURE);
	    }
	    
	    if (PIL_OK != (status = PILGetFname("par_tempextfile", tempExtFile))) {
		LOG(ERROR) << "Error Reading Temperature Extension File";
		return (EXIT_FAILURE);
	    }

		if (PIL_OK != (status = PILGetInt("par_buffer", &buffer))) {
		LOG(ERROR) << "Error Reading buffer size";
		return (EXIT_FAILURE);
	    }

	    if (PIL_OK != (status = PILGetBool("par_clobber", &clobber))) {
		return (EXIT_FAILURE);
	    }
	    if (PIL_OK != (status = PILGetBool("par_history", &history))) {
		return (EXIT_FAILURE);
	    }
	    PILClose(status);
	    return (status);
}

int Cztpha2energy::read(char* inputEvtFile, char* outputEvtFile, char* tempextfile,int buffer, int clobber, int history) {
	
	    int status=0;
	    strcpy(this->inputEvtFile, inputEvtFile);
	    //strcpy(this-> eboundsFile, eboundsFile);
	    //strcpy(this->goFile, goFile);
	    strcpy(this->outputEvtFile, outputEvtFile);
	    strcpy(this->tempExtFile, tempextfile);
		this->buffer = buffer;
	    this->clobber = clobber;
	    this->history = history;
	    return (status);
}

int Cztpha2energy::cztpha2energy_process(){


	    fitsfile *fptrIn, *fptrtempExt, *fptrOut; //Input and Output event file pointers
	    int status=0; // status variable
	    int i,j,k=0; // counter variables
	    int hdunum=0; //stores current hdu number
	    int colnum=0; //stores column number of fits extension
		int energycolnum=0,picolnum=0,timecolnum,cztseccntcolnum,cztntickcolnum,phacolnum,detidcolnum,pixidcolnum,detxcolnum,detycolnum,vetocolnum,alphacolnum,intnull;
	    int hdutype=0; //stores hdutype
	    int ncols=0; //store number of columns in current hdu
	    long nrows=0; //store number of rows in current hdu
	    string errorMsg="";
	    char tempChar[MAX_KEYWORD_SIZE];
	    double doublenull;

	    FitsFileInfo fInInfo; //to store info of any opened file
	    
	    // Variables for validating input files
	    int numhdus;
	    vector<string> v_hdunames;
	    vector<int> v_hdutypes;
	    vector<int> v_ncols;
	    vector<long> v_nrows;
	    
	    // To read and store gain, offset and ebounds information
	    Ebounds ebounds; //to read ebounds file
	    GainOffset gainOffset; //to read CALDB GAIN file
	    EventFileHandler eventFile; //to read Event file and its various extensions.

	    vector <vector <float> > vGainT1_2D;
	    vector <vector <float> > vGainT2_2D;
	    vector <vector <float> > vGainT3_2D;
	    vector <vector <float> > vGainT4_2D;
	    vector <vector <float> > vGainT5_2D;
	    vector <vector <float> > vOffsetT1_2D;
	    vector <vector <float> > vOffsetT2_2D;
	    vector <vector <float> > vOffsetT3_2D;
	    vector <vector <float> > vOffsetT4_2D;
	    vector <vector <float> > vOffsetT5_2D;
	    
	    char *extnames[]={"Q0", "Q1", "Q2", "Q3"}; //can't be set from outside.
	    int eboundsExtnum,frow,felem;
	    double tstart,tstop;
		//already commented
	    //commented for 536
	    // Checking whether output event file exists or not
	    // If yes then deletes it (for clobber=yes)
	    // Otherwise raises an error
	    if (FileExists(outputEvtFile)){
		if (clobber==YES) {
		    if (unlink(outputEvtFile)!=0){
		        LOG(ERROR) << "Error in deleting Output Event File: " << outputEvtFile;
		    }
		}
		else {
		    LOG(INFO) << outputEvtFile << " already exists.";
		    LOG(INFO) << "Use clobber=yes for overwriting the file.";
		    return (EXIT_FAILURE);
		}
	    }
	   // Output event File existence check finished.
	    
	       
	    //Reading input event file
	    //fits_open_file(&fptrIn, inputEvtFile, READWRITE, &status);     //commented by mayuri
	    fits_open_file(&fptrIn, inputEvtFile, READONLY, &status);
	    
	    errorMsg = "Error in opening INPUT EVENT FILE: " + (string) inputEvtFile;
	    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

	    fits_read_key(fptrIn, TDOUBLE, "TSTART", &tstart, NULL, &status);
	    report_error(status, (string) "Error in reading keyword TSTART");

	    fits_read_key(fptrIn, TDOUBLE, "TSTOP", &tstop, NULL, &status);
	    report_error(status, (string) "Error in reading keyword TSTOP");


	    fInInfo.get_fitsinfo(fptrIn); //storing info of input event file in variable f.
	    //commented for 536
	    // Creating output file and copying data from input event file into it.
	
        fits_create_file(&fptrOut, outputEvtFile, &status);
	    errorMsg="Unable to create output fits file: " + (string)outputEvtFile;
	    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}


		 //added by mayuri,5 Oct 2017
	    fits_movabs_hdu(fptrIn,1,&hdutype,&status);
		//fits_movnam_hdu(fptrIn,BINARY_TBL ,"Primary", 0, &status);
		errorMsg="Error in moving to Primary HDU in input file: " + (string)inputEvtFile;
		if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

		/*fits_movnam_hdu(fptrOut, hdutype,1, 0, &status);
		errorMsg= "Error in moving to Primary HDU in output file: " + outputEvtFile;
		if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
 		*/

		
	    fits_copy_hdu(fptrIn, fptrOut, NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in copying Primary extension to output event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }


	    fits_movnam_hdu(fptrIn, BINARY_TBL, "SSM Data", NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in moving to SSM Data extensio of input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }
	    // Adding PI & Energy columns to extensions Q0, Q1, Q2, and Q3 of output file.

	    char *ttype[]={"PI", "ENERGY"};
	    char *tform[]={"I", "E"};

	 	

	    LOG(INFO) << "Adding Columns  PI & ENERGY to output file: "<< outputEvtFile;
	    for(i=0; i<4; i++)
		{
			if (fInInfo.get_hdu_info(extnames[i], hdutype, ncols, nrows)) {
			    LOG(ERROR) << "There is no key/hduname by the name " << extnames[i] << " in the file.";
			    status = EXIT_FAILURE;
			    return status;
			};
			fits_movnam_hdu(fptrIn, hdutype, extnames[i], 0, &status);
			errorMsg= "Error in moving to HDU "+ string(extnames[i]) + " in output file: " + outputEvtFile;
			if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

			//mayuri
			 fits_copy_hdu(fptrIn, fptrOut, NULL, &status);
			    if (status) {
				LOG(ERROR) << "Error in copying Q%d extension to input event file.";;
				fits_report_error(stderr, status);
				return (EXIT_FAILURE);
			    }  

			fits_insert_cols(fptrOut, ncols+1, 2, ttype, tform, &status);
			errorMsg = "Error in inserting new columns viz. PI & ENERGY in " + string(extnames[i]) + "of Output Event File " + outputEvtFile;
			if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

		
	      }
	    
        //mayuri
        // EXTENSIONS VETOSPECTRUM IS COPIED HERE.
        fits_movnam_hdu(fptrIn, BINARY_TBL, "VETOSPECTRUM", NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in moving to VETOSPECTRUM extension of input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }    
	    //mayuri
	    fits_copy_hdu(fptrIn, fptrOut, NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in copying VETOSPECTRUM extension to input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }
 
	    LOG(INFO) << "Columns  PI & ENERGY added to output file: "<< outputEvtFile;    
	    // PI & Energy columns added to extension Q0, Q1, Q2 and Q3 of output file.

	    
	    //Copying SSM and TEMP extensions from SSM event data into other event files.

	    LOG(INFO) << "Copying SSM Data and TEMP extensions from SSM file ";

	    fits_open_file(&fptrtempExt, tempExtFile, READONLY, &status);
	    if (status) {
		LOG(ERROR) <<"Error in opening ModeSS data with TEMP extension file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }


	    fits_movnam_hdu(fptrtempExt, BINARY_TBL, "SSM Data", NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in moving to SSM Data extension of modeSS file";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }   


	    fits_copy_hdu(fptrtempExt, fptrOut, NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in copying SSM Data extension to input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }   


	    fits_movnam_hdu(fptrtempExt, BINARY_TBL, "TEMP", NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in moving to TEMP extension of modeSS file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }


	    fits_copy_hdu(fptrtempExt, fptrOut, NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in copying Temperature extension to input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }

	   fits_close_file(fptrtempExt, &status);
	    if (status) {
		LOG(ERROR) << "Error in closing Temperature Extension file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }

	    //Copying extensions from ssm data is done
	    //Copying GTI extensions from the input modeM0 event file

	    LOG(INFO)<<"Copying GTI extensions ";

	    fits_open_file(&fptrIn, inputEvtFile, READONLY, &status);
	    errorMsg = "Error in opening INPUT EVENT FILE: " + (string) inputEvtFile;
	    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

	    fits_movnam_hdu(fptrIn, BINARY_TBL, "GTI", NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in moving to GTI  extension of input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }    

	    //mayuri
	    fits_copy_hdu(fptrIn, fptrOut, NULL, &status);
	    if (status) {
		LOG(ERROR) << "Error in copying GTI extension to input event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }

	    fits_copy_file(fptrIn, fptrOut, 0, 0 , 1, &status);
	    if (status) {
		LOG(ERROR) << "Error in copyin GTI extensions to output event file.";
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }

	    //fits_close_file(fptrIn, &status);
	    errorMsg = "Error in closing input event file.";
	    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
		//commented for 536
	    //fits_close_file(fptrOut, &status);
	    errorMsg = "Error in closing output fits file.";
	    if (report_error(status, errorMsg)) {
		return EXIT_FAILURE;
	    }





	    //Copying GTI extensions is done

	  
	    // Reading EBOUNDS file.
	    LOG(INFO) << "Reading Ebounds File.";
	    
	    if(QueryCaldb("ASTROSAT","CZTI","-","EBOUNDS",tstart,tstop,eboundsFile,eboundsExtnum))
	    {
		LOG(ERROR) << "Not able to get CALDB EBOUNDS file";
		return (EXIT_FAILURE);
	    }
	    
	    if (ebounds.read_ebounds_file(eboundsFile)) {
		LOG(ERROR) << "Unable to read CALDB EBOUNDS file: " << eboundsFile;
		return (EXIT_FAILURE);
	    }
	    LOG(INFO) << "EBOUNDS File read.";
	    // EBOUNDS file read.
		    
	    // Initializing 2D vectors for gain and offset values
	    vector <float> temp(256, 0.0);
	    vGainT1_2D.resize(NO_DET_PER_QUAD, temp);
	    vGainT2_2D.resize(NO_DET_PER_QUAD, temp);
	    vGainT3_2D.resize(NO_DET_PER_QUAD, temp);
	    vGainT4_2D.resize(NO_DET_PER_QUAD, temp);
	    vGainT5_2D.resize(NO_DET_PER_QUAD, temp);
	    vOffsetT1_2D.resize(NO_DET_PER_QUAD, temp);
	    vOffsetT2_2D.resize(NO_DET_PER_QUAD, temp);
	    vOffsetT3_2D.resize(NO_DET_PER_QUAD, temp);
	    vOffsetT4_2D.resize(NO_DET_PER_QUAD, temp);
	    vOffsetT5_2D.resize(NO_DET_PER_QUAD, temp);



		//Reopening energy added event file
	    fits_open_file(&fptrOut, outputEvtFile, READWRITE, &status);
	    if(status){
		LOG(ERROR) << "Error in opening output event file: " << outputEvtFile;
		fits_report_error(stderr, status);
		return (EXIT_FAILURE);
	    }

		//mayuri
		int bufsize=0,anynull,*arr_pi;
		double *arr_time;
		unsigned short *arr_pha,PI=0;
		unsigned char *arr_detid,*arr_pixid;
		float *arr_energy,gain=0.0,offset=0.0,energy=0.0;
		char quadname[20];
		vector <float> quadTemperature(16,0.0); //to store quadrant temperature evaluated after fitting with TEMP extension.

		bufsize=buffer;
	   	arr_time  = (double*)malloc(bufsize * sizeof(double));
		arr_pha  = (unsigned short*)malloc(bufsize * sizeof(unsigned short));
		arr_detid  = (unsigned char*)malloc(bufsize * sizeof(unsigned char));
	 	arr_pixid  =(unsigned char*) malloc(bufsize * sizeof(unsigned char));

		arr_energy=(float*)malloc(sizeof(float)*bufsize);
		arr_pi=(int*)malloc(sizeof(int)*bufsize);

		if(arr_pi==NULL || arr_energy==NULL)
		{
			LOG(ERROR)<<"Unable to allocate memory";
		}

	    // Running loop for all 4 quadrants
	    for(i=0; i<4; i++)
		{
       


			bufsize=buffer;

			//Reading Gain and Offset
			LOG(INFO) << "Reading Gain and offset files for Quadrant " << i;
		 

			sprintf(quadname,"QUADRANT%d",i);

		    	if(QueryCaldb("ASTROSAT","CZTI",(string)quadname,"GAIN",tstart,tstop,goFile,eboundsExtnum))
		    	{
				LOG(ERROR) << "Not able to get CALDB GAIN file";
				return (EXIT_FAILURE);
		    	}

			if (gainOffset.read_gainoffset_file(goFile, extnames[i])) {
			    LOG(ERROR) << "Unable to read CALDB GAINS file: " << goFile;
			    return (EXIT_FAILURE);
			}


			// Rearranging values in gain and offset 2D vectors
			DLOG(INFO) << "Rearranging gains AND offsets to created 2D vectors";
		
			// Storing gain and offset for TEMPERATURE 1
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_gain(1, status), &vGainT1_2D);
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_offset(1, status), &vOffsetT1_2D);
			//TEMPERATURE 2
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_gain(2, status), &vGainT2_2D);
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_offset(2, status), &vOffsetT2_2D);
			//TEMPERATURE 3
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_gain(3, status), &vGainT3_2D);
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_offset(3, status), &vOffsetT3_2D);
			//TEMPERATURE 4
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_gain(4, status), &vGainT4_2D);
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_offset(4, status), &vOffsetT4_2D);
			//TEMPERATURE 5
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_gain(5, status), &vGainT5_2D);
			create_2d_vector(gainOffset.get_detID(), gainOffset.get_PixID(), gainOffset.get_offset(5, status), &vOffsetT5_2D);        
		
			DLOG(INFO) << "2D gains and offsets created.";
			// Rearranged gain and offset 2d vectors created.
			// Gain and Offset read.
		       
		
			// Reading and storing event file and its various extensions in object eventFile of Event File Handler.
			if(eventFile.read_event_temp(fptrOut, i)){
			    LOG(ERROR) << "Error in reading temperature extension for quadrant " << i << "." ;
			    return EXIT_FAILURE;
			}


			gain=0.0;
			offset=0.0;
			energy=0.0;
			PI=0;
			
		   	status=0;    //status changed by mayuri
			fits_movabs_hdu(fptrOut, i+2, NULL, &status);
			sprintf(tempChar, "%d", i+2);
			errorMsg = "Error in moving to HDU number " + string(tempChar) + " of Output Event File.";
			if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
		
			fits_get_colnum(fptrOut, CASEINSEN, "ENERGY", &energycolnum, &status);
			sprintf(tempChar, "%d", i + 2);
			errorMsg = "Error getting column number of ENERGY column in HDU number " + string(tempChar) + " of Output Event File.";
			if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

			fits_get_colnum(fptrOut, CASEINSEN, "PI", &picolnum, &status);
			sprintf(tempChar, "%d", i + 2);
			errorMsg = "Error getting column number of PI column in HDU number " + string(tempChar) + " of Output Event File.";
			if(report_error(status, errorMsg)) {return EXIT_FAILURE;}

			fits_get_num_rows(fptrOut,&nrows,&status);
			
			fits_get_colnum(fptrOut, CASEINSEN, "Time",&timecolnum, &status);
			fits_get_colnum(fptrOut, CASEINSEN, "PHA",&phacolnum, &status);
			fits_get_colnum(fptrOut, CASEINSEN, "DetID", &detidcolnum, &status);
			fits_get_colnum(fptrOut, CASEINSEN, "pixID",&pixidcolnum, &status);
			
			felem     = 1;

			for(k=0;k<nrows;k+=bufsize)
			{
				if(nrows-k<bufsize)	
						bufsize=nrows-k;
				
			  	frow      = k+1;

				fits_read_col(fptrOut, TDOUBLE, timecolnum, frow, felem, bufsize, NULL, arr_time,&anynull, &status); 
				fits_read_col(fptrOut, TUSHORT, phacolnum, frow, felem, bufsize, NULL, arr_pha, NULL, &status);
				fits_read_col(fptrOut, TBYTE, detidcolnum, frow, felem, bufsize, NULL, arr_detid,&anynull, &status);  
				fits_read_col(fptrOut, TBYTE, pixidcolnum, frow, felem, bufsize, NULL, arr_pixid,&anynull, &status);        

				for(j=0; j<bufsize;j++)                 
				{
					    //Getting temperature of all 16 detectors of quadrant i at time quadTime[j].
						status=0;
					  	quadTemperature = eventFile.get_quadrant_temperature(arr_time[j], status);

					    if(status){
							LOG(ERROR) << "Error in evaluating Quadrant Temperature. ";
							return (status);
					    }
					    
					    if (quadTemperature[arr_detid[j]] < 2.5) {
							gain = vGainT1_2D[arr_detid[j]][arr_pixid[j]];
							offset = vOffsetT1_2D[arr_detid[j]][arr_pixid[j]];
					    }
					    else if (quadTemperature[arr_detid[j]] >= 2.5 && quadTemperature[arr_detid[j]] < 7.5){
							gain = vGainT2_2D[arr_detid[j]][arr_pixid[j]];
							offset = vOffsetT2_2D[arr_detid[j]][arr_pixid[j]];                
					    } 
					    else if (quadTemperature[arr_detid[j]] >= 7.5 && quadTemperature[arr_detid[j]] < 12.5) {
							gain = vGainT3_2D[arr_detid[j]][arr_pixid[j]];
							offset = vOffsetT3_2D[arr_detid[j]][arr_pixid[j]];
					    } 
					    else if (quadTemperature[arr_detid[j]] >= 12.5 && quadTemperature[arr_detid[j]] < 17.5) {
							gain = vGainT4_2D[arr_detid[j]][arr_pixid[j]];
							offset = vOffsetT4_2D[arr_detid[j]][arr_pixid[j]];
					    }
					    else if (quadTemperature[arr_detid[j]] >= 17.5) {
							gain = vGainT5_2D[arr_detid[j]][arr_pixid[j]];
							offset = vOffsetT5_2D[arr_detid[j]][arr_pixid[j]];
					    }

						
					    if(arr_pha[j]==0){
							energy=0.0;
					    }
					    else {
					    	energy = gain*arr_pha[j] + offset;
					    }


					   
					    PI = ebounds.find_channel_from_energy(energy, status);
						if(PI<0)
						{
							printf("%d\t%d\t%f\n",PI,j,energy);
						}
					    if(status){
							LOG(INFO) << "PHA: " << arr_pha[j] << " Energy: " << energy;
							LOG(ERROR)<< "Cannot calculate PI as energy " << energy << " is out of Ebounds";
							return EXIT_FAILURE;
					    }


						arr_energy[j]=energy;
						arr_pi[j]=PI;
		
				}//end loop for j, used for buffering

				fits_write_col(fptrOut, TFLOAT, energycolnum, k+1, 1, bufsize, arr_energy, &status);
				sprintf(tempChar, "%d", i + 2);
				errorMsg = "Error in writing ENERGY values in HDU number " + string(tempChar) + " of Output Event File.";
				if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
		
				fits_write_col(fptrOut, TINT, picolnum, k+1, 1, bufsize, arr_pi, &status);
				sprintf(tempChar, "%d", i + 2);
				errorMsg = "Error in writing PI values in HDU number " + string(tempChar) + " of Output Event File.";
				if(report_error(status, errorMsg)) {return EXIT_FAILURE;}



			}//end loop for k, iterates over nrows

			//mayuri, 23 Aug 2017
         }//end of loop i, iterates over quadrants


    // Read and stored gain, offset and ebounds in their respective class objects.


    status=0;         //mayuri changed the status
    // Updating keywords and writing history if required
    for (i = 1; i <= 8; i++) 
    {
        
        fits_movabs_hdu(fptrOut, i, NULL, &status);
        if (status) {
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_write_date(fptrOut, &status);
        if (status) {
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        fits_update_key(fptrOut, TSTRING, "CREATOR", modulename, "Module that created this file", &status);
        if (status) {
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        if (history == YES) {
            vector<string> vhistory;
            get_history(vhistory);
            for (int j = 0; j < vhistory.size(); j++) {
                fits_write_history(fptrOut, vhistory[j].c_str(), &status);
                if (status) {
                    fits_report_error(stderr, status);
                    return (EXIT_FAILURE);
                }
            }
        }
        fits_write_chksum(fptrOut, &status);
        if (status) {
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }


	free(arr_time);
	free(arr_pha);
	free(arr_detid);
	free(arr_pixid);
	free(arr_energy);
    // printf("status:%d\n",status); 
    status=0;
    fits_close_file(fptrOut, &status);
    errorMsg = "Error in closing output fits file.";
    if(report_error(status, errorMsg)) {return EXIT_FAILURE;}
    // printf("status:%d\n"); 
    return status;

}

int Cztpha2energy::get_history(vector<string> &vhistory) {

    //char *user = getlogin();
    strcpy(modulename, "cztpha2energy_v");
    strcat(modulename, VERSION);
    char *user = getenv("USER");
    string str = "Module run by " + (string) user;
    vhistory.push_back(str);
    vhistory.push_back("Parameter List START for " + (string) modulename);
    vhistory.push_back("P1 infile=" + (string) inputEvtFile);
    vhistory.push_back("P2 Ebounds file="+(string)eboundsFile);
    vhistory.push_back("P3 Gain file="+(string)goFile);
    vhistory.push_back("P5 outfile file="+(string)outputEvtFile);

    if (clobber == YES)
        vhistory.push_back("P5 clobber=yes");
    else
        vhistory.push_back("P5 clobber=no");
    if (history == YES)
        vhistory.push_back("P6 history=yes");
    else
        vhistory.push_back("P6 history=no");
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

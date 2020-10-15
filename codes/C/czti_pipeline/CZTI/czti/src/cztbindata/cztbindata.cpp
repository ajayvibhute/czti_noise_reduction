/* cztbindata.cpp
 *
 * This routine generates mask-weighted/non-maskweighted spectrum and light curves
 * for Cadmium Zinc Telluride Imager on-board AstroSat. Current version is rewritten 
 * based on legacy code (by Tanul Gupta+SAC team). The salient features are:
 *
 * # Maskweight is computed based on open fraction of pixel considering 
 *   the mask alone. 
 *
 * # Three options are expected to be present for usage of aspect information for computing
 *   the mask-weights:
 *
 * 	  $ Average aspect in the Good Time interval used for data selection
 *	  $ Mask-weights are recomputed for each time bin, where time binsize is a parameter
 *	  $ Mask-weights are computed when either tx or ty changes beyond a threhsold value which
 *		is also a parameter
 *
 * 	  At present only the last one is incorporated as it is found to be most suitable	
 *
 * # Mask-weights are renormalized at each module level. This takes into account the 
 * 	 effective area of pixels, exposure fraction of pixels and the background non-uniformity
 *   in the detector plane.
 * 
 * # Also corrections for partial coding is applied such that counts are rescaled to 
 *   "counts per fully illuminated unit area in the detector plane". Further it is 
 *   multiplied by area to rescale it to mimic total counts. 
 *
 * # Offsets between boresight of quadrants with the satellite axis is measured from 
 *   mask-weighted flux profiles. These measured offset corrections are applied here.
 *
 * This code still uses some structures and classes from earlier version and in some 
 * cases there is a change in conventions used from the earlier versions. 
 *
 * Mithun N P S (August 2016)
 *
 *
 * Edits:
 *
 * Added Input/output for inclusion in pipeline (07/09/16)
 * 
 * Ajay Vibhute (Sept 25, 2017)
 * Changed livetime correction
 * Sept 26, 2017
 * Added livetime calculation of fractional bins
 * Oct 3, 2017
 * If fractional exposures are not available then changed the binsize to closed binsize for which
 * Fractional exposures are available
 * */


#include <cztbindata.h>

cztbindata::cztbindata(){
    strcpy(modulename,"cztbindata_v");
    strcat(modulename,VERSION);
    strcpy(quadsToProcess,"-");
    binsize=1.0;
    maskWeight=YES;
}


cztbindata::~cztbindata(){
}

/* *************************  VERY IMPORTANT NOTE  *****************************************
 * 
 * In actual definition of ExposureTable structure (level2validation.h), detx dety meant 
 * det no from 0 to 7 in detector plane which is inconsistent with convention used every 
 * other place. Also there is pixx, and pixy which given pixel location within detector.
 *
 * In this new implementation of cztbindata, instead of this detx vector in ExposureTable 
 * is actually detx (0-63) and similarly dety. pixx, pixy, locx and locy 
 * are of no importance now and can have arbitary values.
 * Also ExposureTable structure is now ONLY FOR ONE QUADRANT. Many variables in the 
 * structure like weightsArray, openfracArray etc are irrelevant now.
 *
 * For each quadrant exposure fractions and weights are stored in vector which 
 * is indexed as detx+dety*64
 *
 * Functions defined in Exposure.cpp are no longer used in cztbindata. All the dependant
 * shadow computation in now present in this folder of cztbindata only.
 *
 * At some point of time efforts should be made to clean up this mess.
 *
 * Mithun N P S (11/08/16)
 *
 * ******************************************************************************************/

int cztbindata::read(int argc,char **argv){
    int status=0;
    vector <string> inStrings;
    
    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error Initializing PIL***";
        return status;
    }
    
    //Input event file
    if(PIL_OK!=(status=PILGetFname("inevtfile",infile))){
        LOG(ERROR)<<"***Error reading input energy added event file"<<infile<<"***";
        return status;
 	}
    
    //MKF file
	if(PIL_OK!=(status=PILGetFname("mkffile",mkffile))){
        LOG(ERROR)<<"***Error getting input MKF file name"<<mkffile<<"***";
        return status;
    }

    //PIXEL QUALITY FILE
    if(PIL_OK!=(status=PILGetFname("badpixfile",pixelqualityfile))){
        LOG(ERROR)<<"***Error getting badpix file name***";
        return status;
    }
   
	//Livetime file
    if (PIL_OK != (status = PILGetFname("livetimefile", livetimefile))) {
        LOG(ERROR) << "***Error reading livetime file name:" << livetimefile<< "***";
        return status;
    }

    //OUTFILE PREFIX: out (spectrum name: out.spec; light curve name: lc.spec)
    if (PIL_OK != (status = PILGetFname("outfile", outfile))) {
        LOG(ERROR) << "***Error reading output LC/SPEC file name:" << outfile << "***";
        return status;
    }

	//output event file name
    if (PIL_OK != (status = PILGetFname("outevtfile", outevtfile))) {
        LOG(ERROR) << "***Error reading output weight added evt file name:" << outevtfile << "***";
        return status;
    }

	//Whether mask-weighted or not
    if (PIL_OK != (status = PILGetBool("maskWeight", &maskWeight))) {
        LOG(ERROR) << "***Error reading maskWeight parameter***";

        return status;
    }

	// RA and DEC of src
    if (PIL_OK != (status = PILGetReal4("rasrc", &RAsrc))) {
        LOG(ERROR) << "***Error reading RA of source***";
        return status;
    }

    if (PIL_OK != (status = PILGetReal4("decsrc", &DECsrc))) {
        LOG(ERROR) << "***Error reading DEC of source***";
        return status;
    }


    //BADPIX THRESHOLD
    if (PIL_OK != (status=PILGetInt("badpixThreshold", &badpixThreshold))){
        LOG(ERROR) << "***Error reading badpix threshold.***";
        return  EXIT_FAILURE;
    }


    //OUTTYPE: SPEC, LC, BOTH
    inStrings.clear();
    inStrings.push_back("lc");
    inStrings.push_back("spec");
    inStrings.push_back("both");
    do {
        if (PIL_OK != (status = PILGetString("outputtype", outputtype))) {
            LOG(ERROR) << "***Error reading output type***";
            return status;
        }
    } while (isStringValid((string) outputtype, inStrings));

    if (strcasecmp(outputtype, "lc") == 0) {
        outtypeflag = LC;
    } else if (strcasecmp(outputtype, "spec") == 0) {
        outtypeflag = SPEC;
    } else if (strcasecmp(outputtype, "both") == 0) {
        outtypeflag = SPECLC;
    }


    //Get time binsize and energy range for lc 
    if(outtypeflag==LC || outtypeflag==SPECLC) {
        if (PIL_OK != (status = PILGetReal4("timebinsize", &binsize))) {
            LOG(ERROR) << "Error reading time binsize***";
            return status;
        }
        do {

            if (PIL_OK != (status = PILGetString("energyrange", energyrange))) {
                LOG(ERROR) << "***Error Reading energyrange:" << energyrange << "***";
                return status;
            }
        } while (isBinsValid(energyrange));

    }


    //QUADS TO PROCESS
    do {
        if (PIL_OK != (status = PILGetString("quadsToProcess", quadsToProcess))) {
            LOG(ERROR) << "***Error reading quadrant:" << quadsToProcess << "***";
            return status;
        }
    } while (isQuadrantNumbersValid(quadsToProcess));



/*
 CALDB FILES (for querying from CIF give input as CALDB)
*/

    //Input compressed mask file
    if(PIL_OK!=(status=PILGetFname("maskfile",compmaskfile))){
        LOG(ERROR)<<"***Error reading mask file"<<compmaskfile<<"***";
        return status;
    }

    //EBOUNDS FILE
    if(PIL_OK!=(status=PILGetFname("eboundsfile",eboundsfile))){
        LOG(ERROR)<<"***Error reading energy bounds CALDB file***";
        return status;
    }
    
    //CAMERA GEOMETRY FILE
    if(PIL_OK!=(status=PILGetFname("camgeomfile", cameraGeomFile))) {
        LOG(ERROR) << "***Error reading camera geometry CALDB file***";
        return status;
    }
    //Effective Area file
    if (PIL_OK != (status = PILGetFname("effareafile", effectiveAreafile))) {
        LOG(ERROR) << "***Error reading effective area CALDB file***";
        return status;
    }

    //LLD file
    if (PIL_OK != (status = PILGetFname("lldfile", LLDfile))) {
        LOG(ERROR) << "***Error reading LLD CALDB file***";
        return status;
    }

/*
    END OF CALDB FILES
*/   

    if(PIL_OK!=(status=PILGetBool("applyLLD",&applyLLD))){
        LOG(ERROR)<<"***Error reading generate_eventfile parameter***";
        return status;
    }

    if(PIL_OK!=(status=PILGetBool("generate_eventfile",&generate_eventfile))){
        LOG(ERROR)<<"***Error reading generate_eventfile parameter***";
        return status;
    }
 
    if(PIL_OK!=(status=PILGetBool("history",&history))){
        LOG(ERROR)<<"***Error reading history parameter***";
        return status;
    }
    if(PIL_OK!=(status=PILGetBool("clobber",&clobber))){
        LOG(ERROR)<<"***Error Reading clobber***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("debug", &debug))) {
        LOG(ERROR) << "***Error reading debugmode parameter***";

        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

void cztbindata::clear(){
    strcpy(infile,"");
    //strcpy(compmaskfile,"");
    strcpy(aspectfileQ0,"");
    strcpy(aspectfileQ1,"");
    strcpy(aspectfileQ2,"");
    strcpy(aspectfileQ3,"");
    //strcpy(cameraGeomFile,"");
    //strcpy(effectiveAreafile,"");
    strcpy(pixelqualityfile, "");
    //strcpy(eboundsFile,"");
    strcpy(outfile,"");
    strcpy(outexparrayfile,"");
//    strcpy(timerange,"");
    strcpy(energyrange,"");
//    strcpy(gtifile,"");
    strcpy(quadsToProcess, "");
    binsize=1.0; 
    strcpy(outputtype,"");
}



int cztbindata::read(char *infile, char* mkffile,char *pixelqualityfile, char *livetimefile, char * outfile, 
		char *outevtfile,int maskWeight,float RAsrc,float DECsrc,int badpixThreshold,char *outputtype,float binsize,
		char *energyrange,char *quadsToProcess,char *compmaskfile,char *eboundsfile, char *cameraGeomFile,
        char *effectiveAreafile,char *LLDfile,int applyLLD,int generate_eventfile, int clobber, int history, int debug)

{
    strcpy(this->infile,infile);
    strcpy(this->mkffile,mkffile);
    strcpy(this->pixelqualityfile,pixelqualityfile);
    strcpy(this->livetimefile,livetimefile);
    strcpy(this-> outfile, outfile);
    strcpy(this->outevtfile,outevtfile);

    strcpy(this->outputtype,outputtype);
    strcpy(this->quadsToProcess,quadsToProcess);
    strcpy(this->compmaskfile,compmaskfile);
    strcpy(this->eboundsfile,eboundsfile);
    strcpy(this->cameraGeomFile,cameraGeomFile);
    strcpy(this->effectiveAreafile,effectiveAreafile);
    strcpy(this->LLDfile,LLDfile);
    //strcpy(this->);
	
    if (strcasecmp(outputtype, "lc") == 0) {
        outtypeflag = LC;
    } else if (strcasecmp(outputtype, "spec") == 0) {
        outtypeflag = SPEC;
    } else if (strcasecmp(outputtype, "both") == 0) {
        outtypeflag = SPECLC;
    }

    if(outtypeflag==LC || outtypeflag==SPECLC){
        strcpy(this->energyrange, energyrange);
        this->binsize=binsize;
    }

	this->maskWeight=maskWeight;
	this->RAsrc=RAsrc;	
    this->DECsrc=DECsrc;
    this->badpixThreshold=badpixThreshold;
    //this->=;

	this->applyLLD=applyLLD;
 	this->generate_eventfile=generate_eventfile;

    this->debug=debug;
    this->clobber=clobber;
    this->history=history;
    this->maskWeight=maskWeight;

    return (EXIT_SUCCESS);

}



void cztbindata::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                          CZTBINDATA PARAMETERS                            ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
	LOG(INFO)<<  "Modulename                        : " <<modulename;
    LOG(INFO) << "Input Event file                  : " << infile;
    LOG(INFO) << "Attitude file                     : " << mkffile;
    LOG(INFO) << "Input badpix file                 : " << pixelqualityfile;
    LOG(INFO) << "Input livetime file               : " << livetimefile;
    LOG(INFO) << "Output file                       : " << outfile;
    LOG(INFO) << "Output Event file                 : " << outevtfile;	
	if(maskWeight==YES)
    LOG(INFO)<<  "Maskweight                        : YES";
    else
    LOG(INFO)<<  "Maskweight                        : NO";
    LOG(INFO) << "Output Type                       : " << outputtype;
    LOG(INFO) << "RA of source (degrees)            : " << RAsrc;
    LOG(INFO) << "Dec of source (degrees)           : " << DECsrc;
    LOG(INFO) << "Badpixel threshold                : " << badpixThreshold;
    if(outtypeflag==LC || outtypeflag==SPECLC) {
    LOG(INFO) << "Energy Range (for lc generation)  : " << energyrange;
    LOG(INFO) << "Time Bin Size (for lc generation) : " << binsize;
    }
	LOG(INFO) << "Ebounds File                      : " << eboundsfile;
    LOG(INFO) << "Camera Geom File                  : " << cameraGeomFile;
    LOG(INFO) << "Effective Area File               : " << effectiveAreafile;
    LOG(INFO) << "Mask pattern file                 : " << compmaskfile;
    LOG(INFO) << "LLD file                          : " << LLDfile;

    if(applyLLD==YES)
    LOG(INFO)<< "Apply LLD cut                      : YES";
    else
    LOG(INFO)<< "Apply LLD cut                      : NO";

    if(generate_eventfile==YES)
    LOG(INFO)<< "Generate event file with weight    : YES";
    else
    LOG(INFO)<< "Generate event file with weight    : NO";

    if(debug==YES)
    LOG(INFO)<< "Debug                              : YES";
    else
    LOG(INFO)<< "Debug                              : NO";

    if(clobber==YES)
    LOG(INFO)<< "Clobber                            : YES";
    else
    LOG(INFO)<< "Clobber                            : NO";
    if(history==YES)
    LOG(INFO)<< "History                            : YES";
    else
    LOG(INFO)<< "History                            : NO";

    LOG(INFO)<<"----------------------------------------------------------------------------";

    //LOG(INFO) << "Quadrants to be processed        :" << quadsToProcess;
    
}


int cztbindata::cztbindataProcess(){


	int status=0,hdutype=0;
	int qid=0;
	long i,j,k,l;
	
    string compMaskFile;
    string camGeomFile;
    string eboundsFile;
    string effAreaFile;
	string lldFile;
    char gtitype[20];
    fitsfile *fevt,*flive,*fmkf;
   
    Attitude attfilehandler;
    vector <attStruct> attvec;	
	double *att_thx,*att_thy;
	long nrows_att;
	double att_startT,att_stopT,att_dur;
	
	double avgthetax=0,avgthetay=0;
	vector <int> quadVec;
   
   	ExposureTable exptable;
    int *maskElements;
    //ajay, 25 Sept
    int pistart,piend;	
    //Getting quadsToProcess
    if(get_quadsToProcessVec((string)quadsToProcess, &quadVec)){
        LOG(ERROR) << "Error getting quadrants to be processed.";
        return EXIT_FAILURE;
    }


	// ********* Measured offsets in boresight of quadrants **************
	// *******************************************************************
    float tx_shift[4]={0.015,-0.125,-0.08,+0.035};
    float ty_shift[4]={-0.025,-0.025,0.205,0.215};
	
	//float tx_shift[4]={0.,-0.12,-0.07,0.03};
	//float ty_shift[4]={0.0,-0.05,0.185,0.17};
	// *******************************************************************
	// *******************************************************************


    //Querying CALDB CIF file to get the file names


    fits_open_file(&fevt, infile, READONLY, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in opening event file";
        return (EXIT_FAILURE);
    }
	
    double head_tstart,head_tstop;

    int  effareaExtnum,eboundsExtnum,geometryExtnum,compmaskExtnum,lldExtnum;

    fits_read_key(fevt, TDOUBLE, "TSTART", &head_tstart, NULL, &status);
    report_error(status, (string) "Error in reading keyword TSTART");

    fits_read_key(fevt, TDOUBLE, "TSTOP", &head_tstop, NULL, &status);
    report_error(status, (string) "Error in reading keyword TSTOP");


	if (strcasecmp(eboundsfile, "CALDB") == 0) 
	{
		if(QueryCaldb("ASTROSAT","CZTI","-","EBOUNDS",head_tstart,head_tstop,eboundsFile,eboundsExtnum))
		{
			LOG(ERROR) << "Not able to get CALDB EBOUNDS file";
			return (EXIT_FAILURE);
    	}
	}
	else
	{
		eboundsFile=(string)eboundsfile;
	}	

    if (strcasecmp(effectiveAreafile, "CALDB") == 0)
    {
		if(QueryCaldb("ASTROSAT","CZTI","QUADRANT0","EFF_AREA",head_tstart,head_tstop,effAreaFile,effareaExtnum))
		{
        LOG(ERROR) << "Not able to get CALDB EFF_AREA file";
        return (EXIT_FAILURE);
    	}
	}
    else
    {
        effAreaFile=(string)effectiveAreafile;
    }
	

    if (strcasecmp(cameraGeomFile, "CALDB") == 0)
    {
	    if(QueryCaldb("ASTROSAT","CZTI","-","GEOMETRY",head_tstart,head_tstop,camGeomFile,geometryExtnum))
    	{
        LOG(ERROR) << "Not able to get CALDB GEOMETRY file";
        return (EXIT_FAILURE);
    	}
	}
    else
    {   
        camGeomFile=(string)cameraGeomFile;
    }     

    if (strcasecmp(compmaskfile, "CALDB") == 0)
    {
    	if(QueryCaldb("ASTROSAT","CZTI","-","MASK_OVERSAMPLED",head_tstart,head_tstop,compMaskFile,compmaskExtnum))
    	{
       	LOG(ERROR) << "Not able to get CALDB MASK_OVERSAMPLED file";
        return (EXIT_FAILURE);
    	}
	}
    else
    {
        compMaskFile=(string)compmaskfile;
    }

    if (strcasecmp(LLDfile, "CALDB") == 0)
	{
        if(QueryCaldb("ASTROSAT","CZTI","QUADRANT0","LLD",head_tstart,head_tstop,lldFile,lldExtnum))
        {
            LOG(ERROR) << "Not able to get CALDB LLD file";
            return (EXIT_FAILURE);
        }
     
	}
	else
	{
		lldFile=(string)LLDfile;
	}		

	LOG(INFO)<<"CALDB LLD file "<<lldFile<<" is used";

    
    fits_read_key(fevt,TSTRING,"GTITYPE",gtitype,NULL,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in reading GTITYPE in event file header";
        return (EXIT_FAILURE);
    }


    fits_close_file(fevt,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing event file";
        return (EXIT_FAILURE);
    }

	// END OF CALDB Query


	//Computing camera cooridnates of the source in case of mask-weighted spec/lc
	
	if(maskWeight==true) 
	{

	// Getting attitude information

/*		
	status=attfilehandler.read_attitude_file((string)/attitudefile);
	if(status) {LOG(ERROR)<<"Error in reading attitude file"; return EXIT_FAILURE;}

	attvec=attfilehandler.get_att();
	nrows_att=attvec.size();
    
	att_thx=(double*)malloc(sizeof(double)*nrows_att);
	att_thy=(double*)malloc(sizeof(double)*nrows_att);

	att_startT=attvec[0].time;
	att_stopT=attvec[nrows_att-1].time;
	att_dur=att_stopT-att_startT;

*/
	double* mkf_time;
	float *rollRA,*rollDEC,*rollRot;	
    fits_open_file(&fmkf, mkffile, READONLY, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in opening MKF file";
        return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fmkf, 2, &hdutype, &status);

    fits_get_num_rows(fmkf, &nrows_att, &status);

	mkf_time=(double*)malloc(sizeof(double)*nrows_att);
	rollRA=(float*)malloc(sizeof(float)*nrows_att);
	rollDEC=(float*)malloc(sizeof(float)*nrows_att);
	rollRot=(float*)malloc(sizeof(float)*nrows_att);

	att_thx=(double*)malloc(sizeof(double)*nrows_att);
    att_thy=(double*)malloc(sizeof(double)*nrows_att);

	int mkf_time_col,mkf_ra_col,mkf_dec_col,mkf_rot_col;


    fits_get_colnum(fmkf,CASEINSEN,"TIME",&mkf_time_col,&status);
    fits_get_colnum(fmkf,CASEINSEN,"Roll_RA",&mkf_ra_col,&status);
    fits_get_colnum(fmkf,CASEINSEN,"Roll_DEC",&mkf_dec_col,&status);
    fits_get_colnum(fmkf,CASEINSEN,"ROLL_ROT",&mkf_rot_col,&status);

    fits_read_col(fmkf, TDOUBLE, mkf_time_col, 1, 1, nrows_att, NULL, mkf_time,NULL, &status);
    fits_read_col(fmkf, TFLOAT, mkf_ra_col, 1, 1, nrows_att, NULL, rollRA,NULL, &status);
    fits_read_col(fmkf, TFLOAT, mkf_dec_col, 1, 1, nrows_att, NULL, rollDEC,NULL, &status);
    fits_read_col(fmkf, TFLOAT, mkf_rot_col, 1, 1, nrows_att, NULL, rollRot,NULL, &status);

    fits_close_file(fmkf,&status);
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }

    att_startT=mkf_time[0];
    att_stopT=mkf_time[nrows_att-1];
    att_dur=att_stopT-att_startT;
	

	LOG(INFO)<<"Attitude information for duration of "<<att_dur<<" seconds is read from the MKF file";
    
	// Completed reading attitude information
	

	// Compute camera coordinates for the source as function of time.
	// Angle is with respect to defined CZTI coordinates (actual boresight shifts are taken care of separately)
	for (i=0;i<nrows_att;i++)
	{
		double RA=RAsrc, DEC=DECsrc;

		//fcompute_tx_ty_(&attvec[i].rollRA,&attvec[i].rollDEC,&attvec[i].rollRot,&RA,&DEC,&att_thx[i],&att_thy[i]);


		double rapnt=rollRA[i];
		double decpnt=rollDEC[i];
		double rotpnt=rollRot[i];

		fcompute_tx_ty_(&rapnt,&decpnt,&rotpnt,&RA,&DEC,&att_thx[i],&att_thy[i]);

//		LOG(INFO)<<"TX TY "<<att_thx[i]<<"  "<<att_thy[i];

		avgthetax+=att_thx[i];
		avgthetay+=att_thy[i];

   /*     if(abs(att_thx[i]) > 5.8||abs(att_thy[i]) > 5.8)
        {
            LOG(ERROR)<<"*****Thetax or Thetay greater than 5.8. Source outside FOV or incorrect attitude*******";
            LOG(ERROR)<<"Thetax = "<<att_thx[i];
            LOG(ERROR)<<"Thetay = "<<att_thy[i];
            LOG(ERROR)<<"Exiting without generating spectrum/lc";
            return EXIT_FAILURE;
        }
	*/	
	}
	
	avgthetax/=(double)nrows_att;
        avgthetay/=(double)nrows_att;

	LOG(INFO)<<"*********AVERAGE TX TY ********* "<<avgthetax<<"   "<<avgthetay;

	}

    //*******Generation of spectrum/lc with/without mask-weighting******


    SpectrumLc spec;
    SpectrumLc lc;
    SpectrumFileHandler specfh;
    LightCurveFileHandler lcfh;
    spec.quadrantid=qid;
    lc.quadrantid=qid;
    char outfileLC[PIL_LINESIZE];
    char outfileSPEC[PIL_LINESIZE];
    vector <string> outfilenames;

	int nebins=0;
	vector <float> userEstart, userEstop;
	long evtnrows;
	double *evttime;
	unsigned char *evtdetx,*evtdety,*evtdetid,*evtpixid;
	int *evtPI;
	float *evtfrac,*evtweight;
	int LLD[4096];
	int evttime_col,evtdetx_col,evtdety_col,evtPI_col,evtdetid_col,evtpixid_col;
	double evtstartT,evtstopT,lvbinsize=0.0,*tmplv_timearray;
	long lcstartT,lcstopT;
	long current_att_tbin=-1,att_tbin,ievt;
	float weight;
	float exp_time,tmpbinsize=0.0;
	int index,mulfact=1.0;
	float shadow_pixels[4096],tmpbinrem=0.0;
	double r2d=180.0/3.14159265;
	double d2r=3.14159265/180.0;
	double *livetime,*lv_timearray;
	long nlvtbins;

    //Opening energy added event file.
    fits_open_file(&fevt, infile, READONLY, &status);
    if (status) {
        LOG(ERROR) << "Error in opening event file " << infile;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

	//Open livetime file
    fits_open_file(&flive, livetimefile, READONLY, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in opening livetime file";
        return (EXIT_FAILURE);
    }

    //Copying input event file to output event file in case maskweighting is applied

    if(maskWeight==true&&generate_eventfile==true)
    {

	  outfilenames.clear();	
      outfilenames.push_back((string)outevtfile);
      if (deleteFile(outfilenames, clobber)) {
        LOG(ERROR) << "Error in deleting output files.";
            return EXIT_FAILURE;
      }

     fitsfile *foutevt;
    
     fits_create_file(&foutevt,outevtfile,&status);
      if(status){
         fits_report_error(stderr,status);
         LOG(ERROR)<<"***Error in creating file "<<outevtfile<<"***";
         return (EXIT_FAILURE);
     }
     LOG(INFO) << "Output evt file created" ;
     fits_copy_file(fevt,foutevt,1,1,1,&status);
     if(status){
         fits_report_error(stderr,status);
         LOG(ERROR)<<"***Error in copying input to output file***";
         return (EXIT_FAILURE);
     }
     LOG(INFO)<<"File copied";

    fits_close_file(foutevt,&status);  
    if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }
	
	}

    //Get energy range for lc generation
    if(outtypeflag==LC || outtypeflag==SPECLC){
        nebins=get_nbins(energyrange);
        if (getVecBins((string) energyrange, &userEstart, &userEstop, nebins)) {
            LOG(ERROR) << "Error extracting Energy bins from energy range provided by user.";
            return EXIT_FAILURE;
        }
	//Ajay 25 Sept 2017
	if(userEstart.size()!=0)
	{
		pistart=(int)((userEstart[0]-5)*2);
		piend=(int)( (userEstop[0]-5)*2);
		if(piend==0)
			piend=512;
	}
	else
	{
		pistart=0;
		piend=512;
	}


	//Ajay, 3 Oct 2017
	//Code to modify binsize on the fly
	tmplv_timearray=(double*)malloc(sizeof(double)*4);
	
	fits_movabs_hdu(flive, 2, &hdutype, &status);
        try
	{
        	fits_read_col(flive, TDOUBLE, 1, 1, 1, 2, NULL, tmplv_timearray,NULL, &status);
        }catch(ErrorHandler errHandler)
	{
        	logError(errHandler);
        	return EXIT_FAILURE;
	}
	//lvbinsize=(tmplv_timearray[1]-tmplv_timearray[0]);
	//read lvbinsize from header
	fits_read_key(flive,TDOUBLE,"LV_BINSIZE",&lvbinsize,NULL,&status);//added by Mayuri
    	if(status) 
	{ 
		fits_report_error(stderr,status);
		LOG(WARNING)<<"Unable to read lvbinsize keyword from livetime header file";
		LOG(WARNING)<<"Calculating the livetime binsize, there may be round up errors in that";
		lvbinsize=(double)((double)tmplv_timearray[1]*1000-(double)tmplv_timearray[0]*1000);
		lvbinsize/=1000;
		//printf("%lf\t%lf\n",tmplv_timearray[0],tmplv_timearray[1]);
	}//added by Mayuri

	tmpbinsize=binsize;
	if(binsize>=1)
	{
		tmpbinrem=fmodf(binsize,lvbinsize);
		tmpbinrem=floor(tmpbinrem);
	//	printf("%.10f\t%.10f\t%.10f\n",tmpbinrem,binsize,lvbinsize);
		//scanf("%d");
		
/*
		if(tmpbinrem<0.0)
			tmpbinrem*=-1;
		if(tmpbinrem>0.000001)
		{
			binsize=floor(binsize/lvbinsize);
			LOG(WARNING)<<"Unable to calculate fractional exposures for binsize "<<tmpbinsize<<", changing binsize to "<<binsize;
			//printf("change binsize to: %f\n",floor(binsize/lvbinsize) );
		}


				
*/
	}
	else
	{
		
		float tmpdiv=lvbinsize/binsize;	
		tmpbinrem=fmodf(tmpdiv,1.0);
		tmpbinrem=floor(tmpbinrem);
	//	printf("%.10f\t%.10f\t%.10f\n",tmpbinrem,binsize,lvbinsize);
		//scanf("%d");

		if(tmpbinrem!=0.0)
		{
			while(tmpbinsize<lvbinsize)
			{	
				tmpbinsize*=10.0;
				mulfact*=10;	
			}
			binsize=((long)(binsize*mulfact))*1.0/mulfact;
			LOG(WARNING)<<"Unable to calculate fractional exposures for binsize "<<tmpbinsize<<", changing binsize to "<<binsize;
		}
	}


	
    }

    spec.clear();
    lc.clear();

    Badpix badpix;

    vector < vector <unsigned char> > badpixMap; //badpixel map for all 128x128 pixels

    if(badpix.read_badpix_file((string)pixelqualityfile)){
        LOG(ERROR) << "Error in reading bad pixel file " << pixelqualityfile;
        return EXIT_FAILURE;
    }

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

	spec.clear();
	lc.clear();


    if(!(strcasecmp(gtitype,"COMMON")==0||strcasecmp(gtitype,"QUAD")==0))
	{

      	LOG(ERROR)<<"Unknown GTITYPE in event file header";
      	return(EXIT_FAILURE);		
	} 

    if(strcasecmp(gtitype,"COMMON")==0) // Generate for all quadrants added
    {
		 LOG(INFO)<<"Event file has COMMON GTI. Spectrum/lc of quadrants will be added together";
         sprintf(outfileLC,"%s.lc",outfile);
         sprintf(outfileSPEC,"%s.pha",outfile);
    }
	else if (strcasecmp(gtitype,"QUAD")==0)
	{

		//added by Ajay Vibhute, Sept 17 2016
		spec.allquad=0;
		lc.allquad=0;
        LOG(INFO)<<"Event file has QUAD GTI. Spectrum/lc of quadrants will be generated separately";
	}

	int qseq,qno;

	for (qseq=0;qseq<quadVec.size();qseq++)
	{
		qid=quadVec[qseq];

        vector <int> quadsToProcess;
        quadsToProcess.clear();
        if(strcasecmp(gtitype,"COMMON")==0) quadsToProcess=quadVec;
        else quadsToProcess.push_back(qid);

		// Initialize spectrum and light curve
		spec.clear();
		lc.clear();
        //lcstartT=(long)head_tstart;//evttime[0];
        //lcstopT=(long)head_tstop+1;//evttime[evtnrows-1];

	//Ajay, 27 Sept 2017
        long nlcbins=0;
	vector <float> tot_lc;
	vector <float> tot_lc_err;
	vector <float> fracexp;





        int nPIch=512;
        vector <float> sumwc,sumw2c;
        sumwc.resize(nPIch);
        sumw2c.resize(nPIch);

        for (i=0;i<nPIch;i++) {sumwc[i]=0;sumw2c[i]=0;}

		for (qno=0;qno<quadsToProcess.size();qno++) // Loop for adding quadrant data in case of common GTI
		{

		qid=quadsToProcess[qno];
		LOG(INFO)<<"******QUAD****** IS "<<qid	;

		//READ mask pattern for this quadrant
	    	maskElements=(int*)malloc(sizeof(int)*TOTALROWS*COMP_COLS);
    		getMaskPattern((char *)compMaskFile.c_str(),maskElements,qid+2);
	
		//READ LLD file for this quadrant
		if(read_lldfile(lldFile,qid,LLD))
		{
			LOG(ERROR)<<"Error in reading LLD file";
			return EXIT_FAILURE;
		}
	
		try{

    	fits_movabs_hdu(fevt, qid+2, &hdutype, &status);
		fits_get_num_rows(fevt, &evtnrows, &status);
        if(evtnrows==0)
        {
            LOG(WARNING)<<"No events in quadrant "<<qid;
            continue;
        }
        fits_get_colnum(fevt,CASEINSEN,"TIME",&evttime_col,&status);
        fits_get_colnum(fevt,CASEINSEN,"DETID",&evtdetid_col,&status);		
        fits_get_colnum(fevt,CASEINSEN,"PIXID",&evtpixid_col,&status);
		fits_get_colnum(fevt,CASEINSEN,"DETX",&evtdetx_col,&status);
        fits_get_colnum(fevt,CASEINSEN,"DETY",&evtdety_col,&status);
        fits_get_colnum(fevt,CASEINSEN,"PI",&evtPI_col,&status);

		}catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }
	
		evttime=(double*)malloc(sizeof(double)*evtnrows);
		evtdetid=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
        evtpixid=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
        evtdetx=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
        evtdety=(unsigned char*)malloc(sizeof(unsigned char)*evtnrows);
        evtPI=(int*)malloc(sizeof(int)*evtnrows);
		evtfrac=(float*)malloc(sizeof(float)*evtnrows);
		evtweight=(float*)malloc(sizeof(float)*evtnrows);
		
		try{
        fits_read_col(fevt, TDOUBLE, evttime_col, 1, 1, evtnrows, NULL, evttime,NULL, &status);
        fits_read_col(fevt, TBYTE, evtdetid_col, 1, 1, evtnrows, NULL, evtdetid,NULL, &status);
        fits_read_col(fevt, TBYTE, evtpixid_col, 1, 1, evtnrows, NULL, evtpixid,NULL, &status);
		fits_read_col(fevt, TBYTE, evtdetx_col, 1, 1, evtnrows, NULL, evtdetx,NULL, &status);
        fits_read_col(fevt, TBYTE, evtdety_col, 1, 1, evtnrows, NULL, evtdety,NULL, &status);
        fits_read_col(fevt, TINT, evtPI_col, 1, 1, evtnrows, NULL, evtPI,NULL, &status);
		}catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }

	//Ajay, Sept 27, 2017
        lcstartT=(long)evttime[0];
        lcstopT=(long)evttime[evtnrows-1]+1;

	nlcbins=ceil((double)(lcstopT-lcstartT)/(double)binsize);
	tot_lc.resize(nlcbins,0.0);
        tot_lc_err.resize(nlcbins,0.0);
	fracexp.resize(nlcbins,0.0);



		LOG(INFO)<<"Read event file for quadrant "<<qid;	
				
		//LOG(INFO)<<"DIFF "<<evtstartT-att_startT<<" "<<evtstopT-att_stopT;

        float prev_thx=1000.0,prev_thy=1000.0;
        long ntotshadowbins=0;
        current_att_tbin=-1;
        float prev_theta=100,theta;
		int detpix_index;

		for (ievt=0;ievt<evtnrows;ievt++)
		{
			detpix_index=evtdetid[ievt]*256+evtpixid[ievt];

			if(maskWeight)
			{
			att_tbin=(long)(evttime[ievt]-att_startT);

			if(att_tbin<0||att_tbin>=nrows_att) 
			{
				LOG(ERROR)<<"*********Attitude file doesn't span the event file********";
				LOG(ERROR)<<"******Unable to compute maskweights.. Exiting.... ******";
				return EXIT_FAILURE;
			}

			theta=atan(sqrt((tan(att_thx[att_tbin]*d2r))*(tan(att_thx[att_tbin]*d2r))+(tan(att_thy[att_tbin]*d2r))*(tan(att_thy[att_tbin]*d2r))))*r2d;

			if(abs(att_thx[att_tbin]) > 5.8||abs(att_thy[att_tbin]) > 5.8) 
			{
				LOG(ERROR)<<"*****Thetax or Thetay greater than 5.8. Source outside FOV*******";
				LOG(ERROR)<<"Tbin = "<<att_tbin<<" out of "<<nrows_att;
				LOG(ERROR)<<"Thetax = "<<att_thx[att_tbin];
				LOG(ERROR)<<"Thetay = "<<att_thy[att_tbin];
				LOG(ERROR)<<"Exiting without generating spectrum/lc";
				return EXIT_FAILURE;
			}

			//if((att_tbin-current_att_tbin)>=50000 || current_att_tbin==-1)
			
			if(abs(att_thx[att_tbin]-prev_thx)>.03||abs(att_thy[att_tbin]-prev_thy)>.03)
			{
				//compute maskweights for this time duration
				
				exptable.reset();			
				if(debug) LOG(INFO)<<"THETX AND THETAY "<<(float)att_thx[att_tbin]+tx_shift[qid]<<"   "<<(float)att_thy[att_tbin]+ty_shift[qid];	
				getShadow((float)att_thx[att_tbin]+tx_shift[qid],(float)att_thy[att_tbin]+ty_shift[qid],qid,shadow_pixels,maskElements,exptable);
				calculate_renormalized_weights(exptable,badpix,badpixThreshold,qid,effAreaFile,(string)infile);
				ntotshadowbins++;
				if(debug) LOG(INFO)<<"Shadow computation: "<<ntotshadowbins<<" txdiff "<<abs(att_thx[att_tbin]-prev_thx)<<" tydiff "<<abs(att_thy[att_tbin]-prev_thy);
				current_att_tbin=att_tbin;
	            prev_theta=atan(sqrt((tan(att_thx[att_tbin]*d2r))*(tan(att_thx[att_tbin]*d2r))+(tan(att_thy[att_tbin]*d2r))*(tan(att_thy[att_tbin]*d2r))))*r2d;
    	        prev_thx=att_thx[att_tbin];
        	    prev_thy=att_thy[att_tbin];
			}

			index=evtdetx[ievt]+evtdety[ievt]*64;
			weight=exptable.weights[index];
			evtfrac[ievt]=exptable.openfrac[index];
			evtweight[ievt]=weight;			
           		

			//Ajay, adding filter on energy, 25 Sept 2017
			//if(  ((evtPI[ievt]>LLD[detpix_index])||(applyLLD==0)) && (evtPI[ievt]> pistart && evtPI[ievt]< piend ) )//Commented by Mayuri,9th Dec 2017 ....energy condition is applied to both lc and spec
			if((evtPI[ievt]>LLD[detpix_index])||(applyLLD==0))
			{
				sumwc[evtPI[ievt]]+=weight;
				sumw2c[evtPI[ievt]]+=weight*weight;
		
			//	if((long)((evttime[ievt]-lcstartT)/binsize)<nlcbins)//Commented by Mayuri,9th Dec 2017
				if((long)((evttime[ievt]-lcstartT)/binsize)<nlcbins && (evtPI[ievt]> pistart && evtPI[ievt]< piend ))//Added by Mayuri,9th Dec 2017...energy condition is applicable only for lc				
				{	
					tot_lc[(long)((evttime[ievt]-lcstartT)/binsize)]+=weight;
					tot_lc_err[(long)((evttime[ievt]-lcstartT)/binsize)]+=weight*weight;
				}

				//module_lc[evtdetid[ievt]][(long)(evttime[ievt]-evtstartT)]+=weight;
			}

			}
			else
			{
				//Ajay 25, Sept 2017
            			if((evtPI[ievt]>LLD[detpix_index])||(applyLLD==0))
            			//if( ( (evtPI[ievt]>LLD[detpix_index])||(applyLLD==0) ) &&  (evtPI[ievt]> pistart && evtPI[ievt]< piend ) )//Commented by Mayuri,9th Dec 2017 ....energy condition is applied to both lc and spec
                   		{				
					//Add this event to spec/lc
					sumwc[evtPI[ievt]]+=1.0;
					sumw2c[evtPI[ievt]]+=1.0;
					//if((long)((evttime[ievt]-lcstartT)/binsize)<nlcbins)//Commented by Mayuri,9th Dec 2017
					if((long)((evttime[ievt]-lcstartT)/binsize)<nlcbins &&  (evtPI[ievt]> pistart && evtPI[ievt]< piend ))//Added by Mayuri,9th Dec 2017...energy condition is applicable only for lc 
					{
						tot_lc[(long)((evttime[ievt]-lcstartT)/binsize)]+=1.0;
						tot_lc_err[(long)((evttime[ievt]-lcstartT)/binsize)]+=1.0;
					}
			
				}

			}

			
		}
        LOG(INFO)<<"Total shadow computations "<<ntotshadowbins;

		//Write weight evt file 

		if(maskWeight&&generate_eventfile)
		{
			if(write_weight_evtfile(outevtfile,qid,evtfrac,evtweight,evtnrows))
			{
				LOG(ERROR)<<"Error in writing event file";
				return EXIT_FAILURE;
			}
		}


		} // End of loop on adding over quadsToProcess in case of common GTI 


        if (strcasecmp(gtitype,"QUAD")==0) // Quadrant wise products
        {
            spec.clear();
			lc.clear();
			sprintf(outfileLC,"%s_Q%d.lc",outfile,qid);
            sprintf(outfileSPEC,"%s_Q%d.pha",outfile,qid);
		}
						
        //Checking existence of outfiles and deletes them if clobber is yes
        outfilenames.clear();
        outfilenames.push_back((string)outfileLC);
        outfilenames.push_back((string)outfileSPEC);
        if (deleteFile(outfilenames, clobber)) {
        LOG(ERROR) << "Error in deleting output files.";
            return EXIT_FAILURE;
        }

		// Write spectrum and light curve as required

        if(outtypeflag==SPEC||outtypeflag==SPECLC)
		{

			spec.flux=sumwc;	
        	transform(sumw2c.begin(), sumw2c.end(), sumw2c.begin(), (float(*)(float))sqrt);
			spec.error = sumw2c;
			spec.maskweighted=maskWeight;
			spec.badpixthresh=badpixThreshold;	
	    	getAvgExposureTime(fevt,quadsToProcess,&exp_time);
    	    spec.exposure_time=exp_time;
			spec.nchannels=nPIch;
			spec.thetaxd=avgthetax;
    	    spec.thetayd=avgthetay;
			spec.gtitype=(string)gtitype;

			spec.channels.resize(nPIch,0);
			for (i=0;i<nPIch;i++) spec.channels[i]=i;
	
        	spec.quadrantid=qid;

	        specfh.set_spectrum(spec);
    	    if(specfh.write_spectrum_file((string)outfileSPEC, "spectrumTemplate",(string)infile)){
        	    LOG(ERROR) << "Error writing spectrum file.";
            	return EXIT_FAILURE;
	        	}
        
		}



        if(outtypeflag==LC||outtypeflag==SPECLC)
        {


        	//Read livetime

			try{
				fits_movabs_hdu(flive, qid+2, &hdutype, &status);
				fits_get_num_rows(flive, &nlvtbins, &status);
			}catch(ErrorHandler errHandler){
				logError(errHandler);
				return EXIT_FAILURE;
			}

        	livetime=(double*)malloc(sizeof(double)*nlvtbins);
        	lv_timearray=(double*)malloc(sizeof(double)*nlvtbins);

        	try{
        		fits_read_col(flive, TDOUBLE, 1, 1, 1, nlvtbins, NULL, lv_timearray,NULL, &status);
        		fits_read_col(flive, TDOUBLE, 2, 1, 1, nlvtbins, NULL, livetime,NULL, &status);
        	}catch(ErrorHandler errHandler){
        		logError(errHandler);
        		return EXIT_FAILURE;
			}
			//getting binsize for the entered livetime file
			double frac_binsize=lvbinsize;//(lv_timearray[1]-lv_timearray[0]);
			//checking we can divide bin or not	
			float binrem=0;//fmod((float)binsize,(float)frac_binsize);
			//calculating no of livetime bins per light-curve bin
			long binratio=0;
			int ii=0,mulfact=1;
		


			//AJay 26 Sept 2017
			if(binsize>=frac_binsize)
			{
				//no of livetime bins per lc
				binratio=(long)(binsize/frac_binsize);
				//binrem=fmodf((float)binsize,(float)frac_binsize);
			//	printf("Bin ratio:%ld\n",binratio);
				binrem=(binratio*frac_binsize)-binsize;
			}
			else
			{
				//no of lc bins per livetime bin
				binratio=(long)ceil(frac_binsize/binsize);
				binrem=fmod((float)frac_binsize/(float)binsize,1.0);
			}
			if(binrem<0.0)
				binrem*=-1;
			//ajay 26 Sept 2017
			if(binrem<=0.000001)
			{
				//getting first bin of the livetime
				long fr_startbin=0;
				//ajay 26 Sept 2017
				if(binsize>=1)
				{
					//getting first livetime bin
					//fr_startbin=(long)((((double)(lcstartT-lv_timearray[0])-(frac_binsize/2.0))));///frac_binsize));
		 			fr_startbin=(long)round((((double)(lcstartT-lv_timearray[0])-(frac_binsize/2.0))/frac_binsize)+1);
				}
				else
					fr_startbin=(long)round((((double)lcstartT-lv_timearray[0]-frac_binsize/2.0)/frac_binsize)+1);
	
			//	printf("fr_startbin:%d\n",fr_startbin);
			//	scanf("%d");
				//Ajay Sept 27,2017 
				if(fr_startbin<=-1) 
				{
					LOG(INFO)<<"Unable to get fractional exposures";
					LOG(INFO)<<"Fractional exposures are set to 1.0";
	                		for(k=0;k<nlcbins;k++)
    	            			{
        	        		    fracexp[k]=1.0;
            	    			} 
				}
				
				//ajay 26 Sept 2017
				if(binsize>=lvbinsize)
				{
					for(k=0;k<nlcbins;k++)
					{
						fracexp[k]=0.0;
						for(l=fr_startbin+k*binratio;l<fr_startbin+(k+1)*binratio;l++)
							fracexp[k]+=livetime[l];
	
					
						/*if(binsize<=1.0)
							fracexp[k]/=binsize;
						else*/
						fracexp[k]/=binsize;
						if(fracexp[k]==0) fracexp[k]=1.0;
				
					}
					//printf("Binratio:%d\t%f\n",binratio,binsize);
				}
				else
				{
				/*	//Added by ajay 26 Sept 2017
					//Livetime calculation for bins < 0 and multiple of binsize
					l=fr_startbin;
					for(k=0;k<nlcbins;l++)
					{
					
					fracexp[k]=0.0;
						for(ii=0;ii<binratio;ii++,k++)
						{

							fracexp[k]=(float)(livetime[l]/binratio);
					
						//	printf("%d\t%d\t%d\t%f\n",ii,k,l,fracexp[k]);
						}
						
					}

					if(fracexp[k]==0) fracexp[k]=1.0;
				*/
	

					LOG(INFO)<<"Unable to get fractional exposures";
					LOG(INFO)<<"Fractional exposures are set to 1.0";
					for(k=0;k<nlcbins;k++)
					{
						fracexp[k]=1.0;
					}
				}
				

			}
		//ajay
			else
			{
				LOG(INFO)<<"Fractional exposures are not available as lc bin size in not multiple of livetime binsize";
				LOG(INFO)<<"Fractional exposures are set to 1.0";
				for(k=0;k<nlcbins;k++)
				{
					fracexp[k]=1.0;
				}	
			}	
    			//Convert to rate and compute error	
			//condition added by ajay.
			//If binsize is < 0 write counts per bin instead of count rate
			//if(binsize>=0)

			for(k=0;k<nlcbins;k++)
			{
				if(fracexp[k]!=0)
				{
					tot_lc[k]/=fracexp[k];
					tot_lc_err[k]/=fracexp[k];
				}
				
			}
		    	transform(tot_lc.begin(), tot_lc.end(), tot_lc.begin(), bind1st(multiplies<double>(), 1.0/binsize) );
		
			lc.flux = tot_lc;
			transform(tot_lc_err.begin(), tot_lc_err.end(), tot_lc_err.begin(), (float(*)(float))sqrt );
			
			transform(tot_lc_err.begin(), tot_lc_err.end(), tot_lc_err.begin(), bind1st(multiplies<double>(), 1.0/binsize) );
			lc.error = tot_lc_err;
				

			lc.UT.clear();
			lc.UT.resize(nlcbins);
			for(i=0;i<nlcbins;i++) lc.UT[i]=(double)lcstartT+(double)i*binsize+binsize/2.0;
			lc.fracexp=fracexp;
            
			lc.timedel=binsize;

			lcfh.set_lc(lc);
        	if(lcfh.write_lc_file((string)outfileLC, "lcTemplate",(string)infile)){
            	LOG(ERROR) << "Error writing light curve file.";
            	return EXIT_FAILURE;
        	}


		}


		if(strcasecmp(gtitype,"COMMON")==0) break;

			/*
			 *************** Include pixel exposure fractions in renorm calculation ***************
			 * */

	}


    fits_close_file(fevt,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing event file";
        return (EXIT_FAILURE);
    }


    fits_close_file(flive,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in closing livetime file";
        return (EXIT_FAILURE);
    }

	

	return EXIT_SUCCESS;

	// FOR MODULE LC 

        /*  
        float **module_lc;  
        module_lc=(float**)malloc(sizeof(float*)*16);
        if(module_lc==NULL) 
        {
            LOG(ERROR)<<"Unable to allocate memory!!!!!";
            LOG(ERROR)<<"Exiting ...";
            return EXIT_FAILURE;
        }
        
        for(i=0;i<16;i++)
        {
             module_lc[i]=(float*)malloc(sizeof(float)*nlcbins);
             if(module_lc[i]==NULL) 
             {
                 LOG(ERROR)<<"Unable to allocate memory!!!!!";
                 LOG(ERROR)<<"Exiting ...";
                 return EXIT_FAILURE;
             }
             for(j=0;j<nlcbins;j++)
                 module_lc[i][j]=0.0;
        }
        */



}



int cztbindata::write_weight_evtfile(char *outevtfile,int qid,float *evtfrac, float *evtweight,long nevt)
{

	fitsfile *foutevt;
	int status=0;	

    char *ttype[]={"OPEN_FRACTION", "WEIGHT"};
    char *tform[]={"E", "E"};

    long nrows;
    int ncols;

    fits_open_file(&foutevt, outevtfile, READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening event file " << outevtfile;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

	char quadnam[10];
    sprintf(quadnam,"Q%d",qid);

    try{
    fits_movnam_hdu(foutevt, BINARY_TBL, quadnam, 0, &status);
    fits_get_num_cols(foutevt,&ncols,&status);
    fits_get_num_rows(foutevt, &nrows, &status);
    fits_insert_cols(foutevt,ncols+1,2,ttype,tform,&status);

	if(nevt!=nrows) {LOG(ERROR)<<"Incorrect rows in event file, exiting..."; return EXIT_FAILURE;}

    fits_write_col(foutevt, TFLOAT, ncols+1, 1,1, nrows, evtfrac, &status);
    fits_write_col(foutevt, TFLOAT, ncols+2, 1,1, nrows, evtweight, &status);

    }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }


    fits_close_file(foutevt,&status); 
    if(status) { 
		LOG(ERROR) << "Error in closing event file " << outevtfile;
		fits_report_error(stderr,status); return (EXIT_FAILURE); 
	}

	LOG(INFO)<<"Writing weight column for quadrant "<<qid<<" is complete.";

	return EXIT_SUCCESS;
}



int cztbindata::read_lldfile(string caldb_lld,int qid,int *lld)
{
    fitsfile *fptr;
    int status=0,hdutype=0;
    int lld_col,i;
    long nrows;

    fits_open_file(&fptr, (char *)caldb_lld.c_str(), READONLY, &status);
    if (status)
    {
    printf("Error (%s:%d): Error while opening %s file\n",__FILE__,__LINE__,(char*)caldb_lld.c_str());
    fits_report_error(stderr, status);
    return (EXIT_FAILURE);
    }

    fits_movabs_hdu(fptr, qid+2, &hdutype, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_num_rows(fptr, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_get_colnum(fptr,CASEINSEN,"LLD",&lld_col,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    for(i=0;i<4096;i++)
    {

        fits_read_col(fptr, TINT, lld_col, i+1, 1, 1, NULL, &lld[i],NULL, &status);
        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    }

    fits_close_file(fptr,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    return(EXIT_SUCCESS);

}

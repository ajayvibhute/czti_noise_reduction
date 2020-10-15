#include "cztimage_v2.h"
#include "maskGeometry.h"

Cztimage::Cztimage() {
    strcpy(modulename, "cztimage_v");
    strcat(modulename, VERSION);
}

int Cztimage::read(int argc, char** argv) {
    int status = 0;
    vector <string> inStrings;
    if (PIL_OK != (status = PILInit(argc, argv))) {
        LOG(ERROR) << "***Error Initializing PIL***";
        return status;
    }
    
    //INPUT TYPE: DPI/DPH/SHADOW
    inStrings.clear();
    inStrings.push_back("dpi");
    inStrings.push_back("dph");
    inStrings.push_back("shadow");
    do{
    if (PIL_OK != (status = PILGetString("par_intype", inputType))) {
        LOG(ERROR) << "***Error reading input type:" << inputType << "***";
        return status;
    }
    } while (isStringValid((string) inputType, inStrings));
    if(strcasecmp(inputType, "dpi")==0){
        inFlag = DPIFILE;
    } 
    else if(strcasecmp(inputType, "dph")==0){
        inFlag = DPHFILE;
    } 
    else if(strcasecmp(inputType, "shadow")==0){
        inFlag = SHADOWFILE;
    } 
    //INPUT TYPE SET
    
    if (PIL_OK != (status = PILGetFname("par_infile", infile))) {
        LOG(ERROR) << "***Error reading input input file:" << infile << "***";
        return status;
    }
    
    /*if (PIL_OK != (status = PILGetFname("par_maskfile64x64", maskfile64x64))) {
        LOG(ERROR) << "***Error reading mask file with 64x64 pixels:" << maskfile64x64 << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetFname("par_maskfile8100", maskfile8100))) {
        LOG(ERROR) << "***Error reading mask file:" << maskfile8100 << "***";
        return status;
    }*/

    if (PIL_OK != (status = PILGetFname("par_shadowfile", shadowfile))) {
        LOG(ERROR) << "***Error reading shadow file:" << shadowfile << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_aspectfileQ0", aspectfileQ0))) {
        LOG(ERROR) << "***Error reading aspect file:" << aspectfileQ0 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_aspectfileQ1", aspectfileQ1))) {
        LOG(ERROR) << "***Error reading aspect file:" << aspectfileQ1 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_aspectfileQ2", aspectfileQ2))) {
        LOG(ERROR) << "***Error reading aspect file:" << aspectfileQ2 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_aspectfileQ3", aspectfileQ3))) {
        LOG(ERROR) << "***Error reading aspect file:" << aspectfileQ3 << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_catalogfile", catalogfile))) {
        LOG(ERROR) << "***Error reading catalog file:" << catalogfile << "***";
        return status;
    }
    if (PIL_OK != (status = PILGetInt("par_extnum", &catalogExtnum))) {
        LOG(ERROR) << "***Error Reading extension number of catalog file:" << catalogExtnum << "***";
        return status;
    }

    if (PIL_OK != (status = PILGetString("par_RAcolname", RA_colname))) {
        LOG(ERROR) << "***Error reading RA column name:***";
        return status;
    }

    if (PIL_OK != (status = PILGetString("par_DECcolname", DEC_colname))) {
        LOG(ERROR) << "***Error reading DEC column name:***";
        return status;
    }

    if (PIL_OK != (status = PILGetString("par_FLUXcolname", FLUX_colname))) {
        LOG(ERROR) << "***Error reading FLUX column name:***";
        return status;
    }
    if (PIL_OK != (status = PILGetFname("par_outImgFile", outImgFile))) {
        LOG(ERROR) << "***Error reading output file:" << outImgFile << "***";
        return status;
    }


    do {
        if (PIL_OK != (status = PILGetString("par_quadsToProcess", quadsToProcess))) {
            LOG(ERROR) << "***Error reading quadrant:" << quadsToProcess << "***";
            return status;
        }
    } while (isQuadrantNumbersValid(quadsToProcess));

    do {
        if (PIL_OK != (status = PILGetReal4("par_resolutionLimit", &resolutionLimit))) {
            LOG(ERROR) << "***Error reading resolution_limit***";
            return status;
        }
    } while (isResolutionLimitValid(resolutionLimit));

    if (PIL_OK != (status = PILGetReal4("par_threshold", &threshold))) {
        LOG(ERROR) << "***Error reading threshold***";
        return status;
    }



    do {
        if (PIL_OK != (status = PILGetInt("par_oversamplingfactor", &oversamplingfactor))) {
            LOG(ERROR) << "***Error reading oversampling factor***";
            return status;
        }
    } while (isOversamplingfactorValid(oversamplingfactor));


    if (PIL_OK != (status = PILGetFname("par_sourcelist", sourcelist))) {
        LOG(ERROR) << "***Error reading sourcelist file name***";
        return status;
    }
    if (PIL_OK != (status = PILGetReal4("par_sourceStrengthThr", &sourceStrengthThr))) {
        LOG(ERROR) << "***Error reading threshold for source strengths***";
        return status;
    }

    do {
        if (PIL_OK != (status = PILGetInt("par_nBkgdCoef", &nBkgdCoef))) {
            LOG(ERROR) << "***Error reading number of background coefficients --Allowed range is 1-6***";
            return status;
        }
    } while (isNBkgdCoefValid(nBkgdCoef));



    if (PIL_OK != (status = PILGetBool("par_debugmode", &debugmode))) {
        LOG(ERROR) << "***Error reading debugmode parameter***";

        return status;
    }
    if (PIL_OK != (status = PILGetBool("par_history", &history))) {
        LOG(ERROR) << "***Error reading history parameter***";
        return status;
    }
    if (PIL_OK != (status = PILGetBool("par_clobber", &clobber))) {
        LOG(ERROR) << "***Error Reading clobber***";
        return status;
    }
    PILClose(status);
    return (EXIT_SUCCESS);
}

void Cztimage::display() {
    LOG(INFO) << "----------------------------------------------------------------------------";
    LOG(INFO) << "                          CZTIMAGE PARAMETERS                            ";
    LOG(INFO) << "----------------------------------------------------------------------------";
	LOG(INFO)<<  "Modulename                   : " <<modulename;
    LOG(INFO) << "Input type                   : " << inputType;
    LOG(INFO) << "Input file                   : " << infile;
    //LOG(INFO) << "Input Mask file            : " << maskfile64x64;
    //LOG(INFO) << "Compressed mask file       : " << maskfile8100;
    LOG(INFO) << "Shadow file             	   : " << shadowfile;
    LOG(INFO) << "Input Aspect file Q0         : " << aspectfileQ0;
    LOG(INFO) << "Input Aspect file Q1         : " << aspectfileQ1;
    LOG(INFO) << "Input Aspect file Q2         : " << aspectfileQ2;
    LOG(INFO) << "Input Aspect file Q3         : " << aspectfileQ3;
    LOG(INFO) << "Catalog file                 : " << catalogfile;
    LOG(INFO) << "Catalog Extension number     : " << catalogExtnum;
    LOG(INFO) << "RA column name               : " << RA_colname;
    LOG(INFO) << "DEC column name              : " << DEC_colname;
    LOG(INFO) << "FLUX column name             : " << FLUX_colname;
    LOG(INFO) << "Resolution limit             : " << resolutionLimit;
    LOG(INFO) << "Threshold                    : " << threshold;
    LOG(INFO) << "Oversampling factor          : " << oversamplingfactor;
    LOG(INFO) << "Output file                  : " << outImgFile;
    LOG(INFO) << "Quadrants to be processed    : " << quadsToProcess;
    LOG(INFO) << "No of background coef        : " << nBkgdCoef;
    LOG(INFO) << "Source list filename         : " << sourcelist;
    LOG(INFO) << "Source strength threshold    : " << sourceStrengthThr;
    LOG(INFO) << "Debug Mode                   : " << debugmode;
    LOG(INFO) << "Clobber                      : " << clobber;
    LOG(INFO) << "History                      : " << history;
    LOG(INFO) << "----------------------------------------------------------------------------";
}

int Cztimage::cztimageProcess(){
    int status=0,quadwise=0;
    long i=0, j=0,ii=0,jj=0;
    DPI dpi;
    ErrorHandler errHandler;
    ShadowFileHandler shdow;
    TanMaskGeometry mask;
    Image img;
    ImageFileHandler imgfh;

	vector < vector <float> > allQuadImage;
    char gtitype[100];
    string key="GTITYPE";
    string extname ="";
    string infilePath = (string) infile;
	float xshift=0,yshift=0;
    //string maskFilePath = caldb_full_path_generator((string) maskfile64x64);
    //string compMaskFilePath = caldb_full_path_generator((string)maskfile8100);
    string shadowfilePath = (string) shadowfile;
    string aspectfilePath = "";
//commented by ajay 
/*
    string aspectfileQ0Path = caldb_full_path_generator((string) aspectfileQ0);
    string aspectfileQ1Path = caldb_full_path_generator((string) aspectfileQ1);
    string aspectfileQ2Path = caldb_full_path_generator((string) aspectfileQ2);
    string aspectfileQ3Path = caldb_full_path_generator((string) aspectfileQ3); 
*/
    string tmp;
    string aspectfileQ0Path = aspectfileQ0;
    string aspectfileQ1Path = aspectfileQ1;
    string aspectfileQ2Path = aspectfileQ2;
    string aspectfileQ3Path = aspectfileQ3; 
    //string catalogfilePath = caldb_full_path_generator((string) catalogfile);
    string outfilePath = (string) outImgFile;
    vector <string> inputfilenames;
    vector <string> outputfilenames;
    vector <string> nonExistentFiles;
    vector < vector <float> > dpiImage;
    vector < vector <long> > dphImage;
    vector < vector <unsigned char> > charMask;
    vector < vector <float> > fltMask;
    stringstream stringbuff;

    //Checking existence of input files for generating cross-correlation image
    inputfilenames.clear();
    inputfilenames.push_back(infilePath);
    //inputfilenames.push_back(maskFilePath);
    files_exist(inputfilenames, &nonExistentFiles);
    if (nonExistentFiles.size() > 0) {
        LOG(ERROR) << "Following input files do not exist:";
        for (i = 0; i < nonExistentFiles.size(); i++) {
            LOG(ERROR) << i + 1 << ". " << nonExistentFiles[i];
        }
        return EXIT_FAILURE;
    }

    //Checking existence of cross-correlation output file and deleting it if clobber
    //is set to yes.
    if(clobber==1){
   		string temp_str;
    	temp_str=outfilePath;
		remove(temp_str.c_str());
	    temp_str.append("_Q0");
    	remove(temp_str.c_str());
	    temp_str=outfilePath;
    	temp_str.append("_Q1");
	    remove(temp_str.c_str());
    	temp_str=outfilePath;
	    temp_str.append("_Q2");
    	remove(temp_str.c_str());
	    temp_str=outfilePath;
    	temp_str.append("_Q3");
	    remove(temp_str.c_str());
 	 }

/*    outputfilenames.push_back(outfilePath);
    if (deleteFile(outputfilenames, clobber)) {
        LOG(ERROR) << "Error in deleting output files.";
        return EXIT_FAILURE;
    }
*/   
   fitsfile *fptrIn; 
   double head_tstart,head_tstop;   
   int maskExtnum;
   string maskFile;
    
    //Opening input event file and checking its validity.
    fits_open_file(&fptrIn,infile,READONLY,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"Error in opening input dpi/dph file:"<<infile;
        return (EXIT_FAILURE);
    }    
    fits_movnam_hdu(fptrIn, IMAGE_HDU,"Q0",0, &status);
    /**
     *Added by Ajay Vibhute, 16 Dec 2015, 10:36PM
     *Added to read the header paramter
     */
    imgfh.headerParam.readFromHeader(fptrIn);	
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"Error in moving hdu file:"<<infile;
        return (EXIT_FAILURE);
    }    

    fits_read_key(fptrIn, TDOUBLE, "TSTART", &head_tstart, NULL, &status);
    report_error(status, (string) "Error in reading keyword TSTART");

    fits_read_key(fptrIn, TDOUBLE, "TSTOP", &head_tstop, NULL, &status);
    report_error(status, (string) "Error in reading keyword TSTOP");
    
    fits_close_file(fptrIn,&status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"Error while closing the  file:"<<infile;
        return (EXIT_FAILURE);
    }    
    if(QueryCaldb("ASTROSAT","CZTI","QUADRANT0","MASK_PATTERN",head_tstart,head_tstop,maskFile,maskExtnum))
    {
         LOG(ERROR) << "Not able to get CALDB MASK_PATTERN file";
         return (EXIT_FAILURE);
    }

    
    
    
    try {
        //Reading mask file
        if (mask.read_64x64_mask(maskFile)) {
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in reading mask file: " + maskFile;
            throw errHandler;
        }

        if(inFlag==SHADOWFILE){
        shdow.read_shadow_file(infilePath, "SHADOW");
        }
        
	//Get Type of GTI
        getHeaderKey(infilePath.c_str(),1,key.c_str(),gtitype);


	if (strcasecmp (gtitype,"common") == 0) 
	{
		quadwise=0;
        	//Creating empty image file
	        imgfh.create_image_file(outfilePath, IMAGETEMPLATE);
	        //Empty image file created.
	}
	else if(strcasecmp (gtitype,"quad") == 0) 
	{
		quadwise=1;	

		for(i=0;i<=NUMQUAD-quadwise;i++)
		{
			stringbuff <<i;
			tmp=outfilePath+"_Q"+stringbuff.str();
			LOG(INFO)<<"Create outfile:"<< tmp;
			stringbuff.str("");
	        	imgfh.create_image_file(tmp, IMAGETEMPLATE);
		}
	}
	else
	{
        	LOG(ERROR) << "ERROR::: UNKNOWN GTI TYPE";
	        return EXIT_FAILURE;
	}
      

	for(i=0; i<=NUMQUAD-quadwise; i++) {
            extname = quadExtName(i);

        if (inFlag == DPIFILE) 
	    {
                if(i==NUMQUAD && extname=="")
				{
		       	extname = "FULL_DPI";
				}
                dpiImage.clear();
                dpiImage = dpi.get_dpi_image(infilePath, extname);

        } 
	    else if(inFlag==DPHFILE)
	    {
                if(i==NUMQUAD && extname=="")
				{
		       		extname = "FULL_DPH";
				}
                dphImage.clear();
                dpiImage.clear();
                dphImage = dpi.get_dph_image(infilePath, extname);
                assign_vector(dphImage, &dpiImage);
                //convert this vector into float here.
        } 
	    else if(inFlag==SHADOWFILE)
	    {
                if(i==NUMQUAD && extname=="")
				{ 
					extname = "QALL";
				}
                dpiImage.clear();
                dpiImage = shdow.get_shadow(extname);
        }
            
            //Getting mask for respective quadrants and assigning aspect file path
            charMask.clear();
            if(i==0){ //Quadrant 0

		/*Ajay Vibhute, 8 Dec 2015 12:05 
		 * Create a single function which will accept qid as an input 
		 * and return corrosponding mask
		 * */    
                charMask = mask.get_q0Mask();
                aspectfilePath = aspectfileQ0Path;
				xshift=0.0;
				yshift=0.0;
            } 
            else if(i==1){ //Quadrant 1
                charMask = mask.get_q1Mask();
                aspectfilePath = aspectfileQ1Path;

				xshift=0.0;
				yshift=0.0;
            } 
            else if(i==2){ //Quadrant 2
                charMask = mask.get_q2Mask();
                aspectfilePath = aspectfileQ2Path;
				xshift=0.0;
				yshift=1.0;

            } 
            else if(i==3){ //Quadrant 3
                charMask = mask.get_q3Mask();
                aspectfilePath = aspectfileQ3Path;
				xshift=0.0;
				yshift=1.0;
            } 
            else if(i==4){ //All Quadrants
                mask.get_full64x64_mask();
                charMask = mask.get_fullMask();
                aspectfilePath = aspectfileQ0Path;
            }
            //Mask extractd and aspect file path set.
            
            //converting mask vector into float here.
            fltMask.clear();
            assign_vector(charMask, &fltMask);
			if(i!=4)
	            img.generate_cross_correlation_image(&dpiImage, &fltMask, oversamplingfactor,xshift,yshift);
			

           	if(i==0)
			{
				allQuadImage.resize(img.crossCorrelationImage.size());
				for(j=0;j<img.crossCorrelationImage.size();j++)
				{
					allQuadImage[j].resize(img.crossCorrelationImage[j].size(),0.0);
	//				allQuadImage[i][j]=0;
				}

			}	   
			if(i!=4)
			{
			for(ii=0;ii<img.crossCorrelationImage.size();ii++)
			{
					for(jj=0;jj<img.crossCorrelationImage[ii].size();jj++)
					{
							allQuadImage[ii][jj]=allQuadImage[ii][jj]+img.crossCorrelationImage[ii][jj];
					}		
			}
			}
			if(i==4)
			{
				for(ii=0;ii<img.crossCorrelationImage.size();ii++)
				{
					for(jj=0;jj<img.crossCorrelationImage[ii].size();jj++)
					{
							img.crossCorrelationImage[ii][jj]=allQuadImage[ii][jj];
					}		
				}

			}
            img.evaluate_wcs_coordinates(oversamplingfactor, aspectfilePath);
            if(i==NUMQUAD && extname==""){ extname = "IMAGE";}




			if(quadwise==0)
		    tmp=outfilePath;
	    else
	    {

			stringbuff.str("");
			stringbuff <<i;
			tmp=outfilePath+"_Q"+stringbuff.str();
	    }

            imgfh.write_image_file(tmp, extname, &img);
        }
    } catch (ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }
    
    
    return status;
}

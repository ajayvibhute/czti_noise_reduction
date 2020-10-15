/* cztrspgen:
 *
 * Module for generating response file for CZTI spectrum
 * Writes response in .rsp (fits) format

See the documents
* A note on cztrspgen module
* Spectral extraction and response matrix generation for CZTI

Mithun N P S(01/07/2015)

Edit:10/08/2015
LLD file is read and respose below LLD is set to zero

Edit: Made changes to include in the pipeline
15/01/2016

Edit: Made changes to conform with updates in cztbindata.
08/09/2016

*/


#include<cztrspgen.h>


cztrspgen::cztrspgen(){
    strcpy(modulename, "cztrspgen_v");
	strcat(modulename,VERSION);
}

cztrspgen::~cztrspgen(){}

int cztrspgen::read(int argc,char **argv)
{

    int status=0;

    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error Initializing PIL***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("phafile",phafile))){
        LOG(ERROR)<<"***Error reading input Spectrum file"<<phafile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("rspfile",rspfile))){
        LOG(ERROR)<<"***Error reading response file"<<rspfile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("evtfile",evtfile))){
        LOG(ERROR)<<"***Error reading Event file"<<evtfile<<"***";
        return status;
    }

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

    //Resp_par file
    if (PIL_OK != (status = PILGetFname("respparfile", respparFile))) {
        LOG(ERROR) << "***Error reading LLD CALDB file***";
        return status;
    }

    //pixresp file
    if (PIL_OK != (status = PILGetFname("pixrespfile", pixrespFile))) {
        LOG(ERROR) << "***Error reading LLD CALDB file***";
        return status;
    }


/*
    END OF CALDB FILES
*/





/*	
    if(PIL_OK!=(status=PILGetFname("expmapfile",expmapfile))){
        LOG(ERROR)<<"***Error reading exposure map file"<<expmapfile<<"***";
        return status;
    }
*/
    if(PIL_OK!=(status=PILGetFname("badpixfile",badpixfile))){
        LOG(ERROR)<<"***Error reading badpix file"<<badpixfile<<"***";
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

    PILClose(status);

    return(EXIT_SUCCESS);
}

int cztrspgen::read(char *phafile,char *rspfile,char *evtfile, char *compmaskfile,char *eboundsfile, char *cameraGeomFile,
        char *effectiveAreafile,char *LLDfile, char *respparFile,char*pixrespFile,char *badpixfile,int clobber,int history)
{
    strcpy(this->phafile,phafile);
    strcpy(this->rspfile,rspfile);
    strcpy(this->evtfile,evtfile);
    strcpy(this->badpixfile,badpixfile);
    strcpy(this->compmaskfile,compmaskfile);
    strcpy(this->eboundsfile,eboundsfile);
    strcpy(this->cameraGeomFile,cameraGeomFile);
    strcpy(this->effectiveAreafile,effectiveAreafile);
    strcpy(this->LLDfile,LLDfile);
	strcpy(this->respparFile,respparFile);
    strcpy(this->pixrespFile,pixrespFile);
    this->clobber = clobber;
    this->history = history;
    return (EXIT_SUCCESS);
}

void cztrspgen::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                          CZTRSPGEN PARAMETERS                            ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
	LOG(INFO) << "Modulename              : " << modulename;
    LOG(INFO) << "Input spectrum file     : " << phafile;
    LOG(INFO) << "Response matrix file    : " << rspfile;
    LOG(INFO) << "Event file              : " << evtfile;
    LOG(INFO) << "Input badpix file       : " << badpixfile;

    LOG(INFO) << "Ebounds File            : " << eboundsfile;
    LOG(INFO) << "Camera Geom File        : " << cameraGeomFile;
    LOG(INFO) << "Effective Area File     : " << effectiveAreafile;
    LOG(INFO) << "Mask pattern file       : " << compmaskfile;
    LOG(INFO) << "LLD file                : " << LLDfile;
    LOG(INFO) << "Response parameter file : " << respparFile;
    LOG(INFO) << "Pixel response file     : " << pixrespFile;

    if(clobber==YES)
    LOG(INFO) << "Clobber                 : YES";
    else
    LOG(INFO) << "Clobber                 : NO";
    if(history==YES)
    LOG(INFO) << "History                 : YES";
    else
    LOG(INFO) << "History                 : NO";


}

int cztrspgen::cztrspgenProcess()
{

    int status=0;
    int i,j,k,T;

    // Files required by this routine
    fitsfile *fpha;
    fitsfile *fexpmap;
	fitsfile *fbadpix;
    fitsfile *fpixresp;
    fitsfile *fresp_par;
    fitsfile *feff_area;
    fitsfile *febounds;
    fitsfile *flld;

	string pixresp_file,resp_par_file,lld_file;
    string eff_area_file,febounds_name;
	string compMaskFile;
	string camGeomFile;

    double thetax,thetay,theta=0.;
	float quad_tx,quad_ty;
    long nresp_group;
    char Erange_string[50];
	char gtitype[20]; 
    float Emin, Emax;
    int n_Ebin,n_PIbin;
    int matrix_n;
    int det,pix,quad,row;
    int detx,dety,locx,locy;
	int indexPix;
	int QUAD_ID=0;
    double Tfrac[5];
    int pixflag_n,eff_area_n;
    int pix_flag;
    int resp_id_n,lld_n;
	int pixthreshold,maskwt,quadrantid;
    double Time_frac[16][5];
    int Ebin,PIbin;
    int hdutype=0;
    double head_tstart;
    double head_tstop;
	TanMaskGeometry tanmask;
	vector <vector <unsigned char> > fullfinermask;
	ExposureTable RespExpTable;
	int quadstart,quadstop;
    int *maskElements;

    // ********* Measured offsets in boresight of quadrants **************
    // *******************************************************************
    float tx_shift[4]={0.015,-0.125,-0.08,+0.035};
    float ty_shift[4]={-0.025,-0.025,0.205,0.215};
    // *******************************************************************
    // *******************************************************************


    //FILE EXISTENCE CHECK AND UNLINKING IF IT DOES
    if (clobber == YES) {
        if (FileExists(rspfile)) {
            LOG(INFO) << rspfile << "  :FileExists.. Replacing the old file";
            if (unlink(rspfile) != 0) {
                LOG(ERROR) << "***Error in deleting " << rspfile << "***";
                return (EXIT_FAILURE);
            }
        }
    }
    else {
        if (FileExists(rspfile)) {
            LOG(ERROR) << "***Output file already exists***";
            LOG(ERROR) << "Use clobber=yes for overwriting the file";
            return (EXIT_FAILURE);
        }
    }


    //Read spectrum header
    try
    {
    fits_open_file(&fpha,phafile,READONLY,&status);
    fits_movabs_hdu(fpha, 2, &hdutype, &status);
    fits_read_key(fpha, TINT, "MASKWT", &maskwt, NULL, &status);
    fits_read_key(fpha, TINT, "PIXTHRES", &pixthreshold, NULL, &status);
    fits_read_key(fpha, TINT, "QUADID", &quadrantid, NULL, &status);
    fits_read_key(fpha, TSTRING, "GTITYPE", gtitype, NULL, &status);
	fits_read_key(fpha, TDOUBLE, "THETAX", &thetax, NULL, &status);
    fits_read_key(fpha, TDOUBLE, "THETAY", &thetay, NULL, &status);
    fits_read_key(fpha, TDOUBLE, "TSTART", &head_tstart, NULL, &status);
    fits_read_key(fpha, TDOUBLE, "TSTOP", &head_tstop, NULL, &status);
    fits_close_file(fpha,&status);
    }catch (ErrorHandler errHandler){
        throw errHandler;
    }

	//Get the CALDB file names from CALDB index file
    int effareaExtnum,eboundsExtnum,pixrespExtnum,resp_parExtnum,lldExtnum,geometryExtnum,compmaskExtnum;

    if (strcasecmp(eboundsfile, "CALDB") == 0)
    {
    	if(QueryCaldb("ASTROSAT","CZTI","-","EBOUNDS",head_tstart,head_tstop,febounds_name,eboundsExtnum))
    	{
        	LOG(ERROR) << "Not able to get CALDB EBOUNDS file";
        	return (EXIT_FAILURE);
    	}
	}
	else
    {
        febounds_name=(string)eboundsfile;
    }


    if (strcasecmp(effectiveAreafile, "CALDB") == 0)
    {
    	if(QueryCaldb("ASTROSAT","CZTI","QUADRANT0","EFF_AREA",head_tstart,head_tstop,eff_area_file,effareaExtnum))
    	{
       	 	LOG(ERROR) << "Not able to get CALDB EFF_AREA file";
        	return (EXIT_FAILURE);
    	}
	}
    else
    {
        eff_area_file=(string)effectiveAreafile;
    }

    if (strcasecmp(respparFile, "CALDB") == 0)
    {
    	if(QueryCaldb("ASTROSAT","CZTI","QUADRANT0","RESP_PAR",head_tstart,head_tstop,resp_par_file,resp_parExtnum))
    	{
       	 	LOG(ERROR) << "Not able to get CALDB RESP_PAR file";
        	return (EXIT_FAILURE);
    	}
	}
	else
	{
		resp_par_file=(string)respparFile;
	}


    if (strcasecmp(pixrespFile, "CALDB") == 0)
    {
	    if(QueryCaldb("ASTROSAT","CZTI","-","PIXRESP",head_tstart,head_tstop,pixresp_file,pixrespExtnum))
    	{
        	LOG(ERROR) << "Not able to get CALDB PIXRESP file";
        	return (EXIT_FAILURE);
    	}
	}
    else
    {
        pixresp_file=(string)pixrespFile;
    }


    if (strcasecmp(LLDfile, "CALDB") == 0)
    {

	    if(QueryCaldb("ASTROSAT","CZTI","QUADRANT0","LLD",head_tstart,head_tstop,lld_file,lldExtnum))
    	{
        	LOG(ERROR) << "Not able to get CALDB LLD file";
        	return (EXIT_FAILURE);
	    }
	}
    else
    {
        lld_file=(string)LLDfile;
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
	

	//Open pixresp file
    try{
    fits_open_file(&fpixresp,(char *)pixresp_file.c_str(),READONLY,&status);

    fits_movabs_hdu(fpixresp, 2, &hdutype, &status);

    fits_read_key(fpixresp, TINT, "N_EBIN", &n_Ebin, NULL, &status);

    fits_read_key(fpixresp, TINT, "N_PIBIN", &n_PIbin, NULL, &status);

    fits_read_key(fpixresp, TSTRING, "CBD20001", Erange_string, NULL, &status);

    fits_get_num_rows(fpixresp, &nresp_group, &status);

    fits_get_colnum(fpixresp,CASEINSEN,"MATRIX",&matrix_n,&status);

    }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }

    // Get the Emin value frm the string
    Emin=atof(Erange_string+5);

    char tmp=*" ";

    i=0;

    while(tmp!=*"-")
    {
        tmp=*(Erange_string+i);
        i+=1;
    }

    // Get the Emax value after the '-'
    Emax=atof(Erange_string+i);

	float Ebinsize=(Emax-Emin)/(n_Ebin);


	//Initialize total_response and eff_area
    float **total_response,*eff_area;
    total_response=(float**)malloc(sizeof(float*)*n_Ebin);
	eff_area=(float*)malloc(sizeof(float)*n_Ebin);
    for(i=0;i<n_Ebin;i++)
    {
		eff_area[i]=0;
        total_response[i]=(float*)malloc(sizeof(float)*n_PIbin);
        for(j=0;j<n_PIbin;j++)
            total_response[i][j]=0;
    }

	//Read the mask pattern 

/*	
	LOG(INFO) << "Reading compressed mask file to create 18200x18200 full uncompressed mask.";
	if (tanmask.read_compressed_mask(compMaskFile)) {
    	return EXIT_FAILURE;
	}
	if (tanmask.get_full_finer_mask()) {
    	LOG(ERROR) << "Error in generating full finer mask from individual quadrant masks.";
	    return EXIT_FAILURE;
	}

	//storing full finer mask
	fullfinermask = tanmask.get_fullUnMask();
	//clearing tanmask stored variables to save memory
	tanmask.reset_full_uncomressed_mask();

	LOG(INFO)<<"Computing open fractions of pixels at response energy bins";
	string noebounds="";
	if(create_exposure_array(thetax, thetay,camGeomFile,fullfinermask, RespExpTable, 1,10 ,noebounds, Emin, Emax , n_Ebin))
	{
		LOG(ERROR)<<"Unable to compute exposure fractions for Response energy bins";
		return (EXIT_FAILURE);
	}

    //Generating exposure index table for faster processing
    try{
        RespExpTable.generate_index_table(false);
    } catch (ErrorHandler errHandler){
        throw errHandler;
    }
*/
	
    Badpix badpix;

    if(badpix.read_badpix_file((string)badpixfile)){
        LOG(ERROR) << "Error in reading bad pixel file " << badpixfile;
        return EXIT_FAILURE;
    }

	float weight,openfrac;
    float shadow_pixels[4096];
    ExposureTable exptable;

/*	
	FILE *ff;
	ff=fopen("mask_weights_406.txt","r");
	
	for(i=0;i<64;i++)
	{
		for(j=0;j<64;j++)
		{
			fscanf(ff,"%f\n",&wt);
			weight[j][i]=wt;
		}
	}
	fclose(ff);
*/

	if(strcasecmp(gtitype,"COMMON")==0) 
	{
		quadstart=0;
		quadstop=4;
	}
	else
	{
		quadstart=quadrantid;
		quadstop=quadrantid+1;
	}

    for (quad=quadstart;quad<quadstop;quad++)    // Loop on quadrant 
    {


        //READ mask pattern for this quadrant
        maskElements=(int*)malloc(sizeof(int)*TOTALROWS*COMP_COLS);
        getMaskPattern((char *)compMaskFile.c_str(),maskElements,quad+2);

        exptable.reset();
		quad_tx=thetax+tx_shift[quad];
		quad_ty=thetay+ty_shift[quad];
        getShadow(quad_tx,quad_ty,quad,shadow_pixels,maskElements,exptable);
		calculate_renormalized_weights(exptable,badpix,pixthreshold,quad,eff_area_file,(string)evtfile);	

        LOG(INFO)<<"Processing quadrant "<<quad;

		//Open CALDB files and badpix file
        int resp_par_hdu=quad+2,lld_hdu=quad+2,eff_area_hdu=quad+2,badpix_hdu=quad+2;
        try{
        fits_open_file(&fresp_par,(char *)resp_par_file.c_str(),READONLY,&status);
        fits_movabs_hdu(fresp_par, resp_par_hdu, &hdutype, &status);
        fits_get_colnum(fresp_par,CASEINSEN,"RESP_ID",&resp_id_n,&status);

        fits_open_file(&flld,(char *)lld_file.c_str(),READONLY,&status);
        fits_movabs_hdu(flld, lld_hdu, &hdutype, &status);
        fits_get_colnum(flld,CASEINSEN,"LLD",&lld_n,&status);

        fits_open_file(&feff_area,(char *)eff_area_file.c_str(),READONLY,&status);
        fits_movabs_hdu(feff_area, eff_area_hdu, &hdutype, &status);
        fits_get_colnum(feff_area,CASEINSEN,"AREA",&eff_area_n,&status);
        
		fits_open_file(&fbadpix,badpixfile,READONLY,&status);
        fits_movabs_hdu(fbadpix, badpix_hdu, &hdutype, &status);
        fits_get_colnum(fbadpix,CASEINSEN,"PIX_FLAG",&pixflag_n,&status);
        }catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
        }

		LOG(INFO)<<"Read CALDB information";
        //Compute the time fraction at each temperature for all modules in the quadrant
        //compute_exposure_frac(Time_frac,quad);

        for(det=0;det<16;det++)  // For each detector in the quadrant
        {

            for (i=0;i<5;i++) Time_frac[det][i]=0.;
            Time_frac[det][2]=1.0;                  //All except 10 deg set to zero for testing

            for (pix=0;pix<256;pix++)   // For each pixel in the module
            {
                row=det*256+pix;        // Row of this pixel in caldb file tables (0-4095)
                

				fits_read_col(fbadpix,TINT,pixflag_n,row+1,1,1,NULL,&pix_flag,NULL,&status);
				if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

				if(pix_flag<=pixthreshold)
			    {	
				generate_locx_locy(det, pix,quad, detx, dety, locx, locy);
				//indexPix = RespExpTable.indexTable[locy][locx];
                
				weight=exptable.weights[detx+dety*64];

				openfrac=exptable.openfrac[detx+dety*64];
 
                fits_read_col(feff_area,TFLOAT,eff_area_n,row+1,1,n_Ebin,NULL,eff_area,NULL,&status);
                if(status){ LOG(ERROR)<<"Unable to read effective area"; return(EXIT_FAILURE); }

                // Read the response ID for the pixel at five temperatures (Group number of response)
                int resp_id[5],LLD[5];

                fits_read_col(fresp_par,TINT,resp_id_n,row+1,1,5,NULL,&resp_id,NULL,&status);
                if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }
                
				fits_read_col(flld,TINT,lld_n,row+1,1,5,NULL,&LLD,NULL,&status);
                if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }      				

                for (T=0;T<5;T++) // For each temperature (index 0-4 corresponds to 0,5,10,15,20 deg C temp)
                {
                    if(Time_frac[det][T] > 0) // If the fraction of total time at this temperature is non zero
                    {

					//fprintf(ff,"%d\t%d\t%f\n",det,pix,RespExpTable.openfracArray[indexPix][40]);	

                    for (Ebin=0;Ebin<n_Ebin;Ebin++) // for each incident photon energy bin
                        {
                        float resp[n_PIbin];

                        // Read the redistribution probability for current energy bin
                        fits_read_col(fpixresp,TFLOAT,matrix_n,resp_id[T]+1,(Ebin*n_PIbin+1),n_PIbin,NULL,resp,NULL,&status);
                        if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

                        //Set response to zero below the LLD
                        for(i=0;i<LLD[T];i++)
                            resp[i]=0;

					    float alphaTa=absco_Ta(Ebin*Ebinsize+Emin+Ebinsize/2.0);
						double trTam=exp(-1*0.05*alphaTa);

                        // Add the response multiplied with the corresponding weight to the total response
                        	if(maskwt==1)
							{
								for(PIbin=0;PIbin<n_PIbin;PIbin++)
                            		total_response[Ebin][PIbin]+=resp[PIbin]*eff_area[Ebin]*Time_frac[det][T]*weight*(openfrac+trTam*(1-openfrac));
							}
							else
							{
								for(PIbin=0;PIbin<n_PIbin;PIbin++)
									total_response[Ebin][PIbin]+=resp[PIbin]*eff_area[Ebin]*Time_frac[det][T]*(openfrac+trTam*(1-openfrac));
							}
                        }

                    }
               	}
			
			  	} //pixthreshold logic if end	
			}

		} //det loop end


		//Close CALDB files
	    try{
    	fits_close_file(fresp_par,&status);
    	fits_close_file(flld,&status);
    	fits_close_file(feff_area,&status);
    	}catch(ErrorHandler errHandler){
    		logError(errHandler);
    		return EXIT_FAILURE;
    	}

	} 

    LOG(INFO)<<"Response calculation completed. Writing rsp file";

	//fclose(ff);
    // Write the response matrix to text file 
/*
    FILE *ftxt;
    ftxt=fopen("response.txt","w");

    for (i=0;i<n_Ebin;i++)
    {
        for (j=0;j<n_PIbin;j++)
            fprintf(ftxt,"%15.12f\t",total_response[i][j]);
        fprintf(ftxt,"\n");
    }

    fclose(ftxt);
*/
    
	float LowEnergy[n_Ebin],HighEnergy[n_Ebin];
    int NumberGroups[n_Ebin];
    int NumberChannelsGroup[n_Ebin];
    int FirstChannelGroup[n_Ebin];

    for(i=0;i<n_Ebin;i++)
    {
        // Compute the lower and higher energies of bins (incident photon energy)
        LowEnergy[i]=Emin+i*Ebinsize;
        HighEnergy[i]=Emin+(i+1)*Ebinsize;
        // Set the variables defining the rmf compression so that full 2-d matrix 
        NumberGroups[i]=1;
        FirstChannelGroup[i]=0;
        NumberChannelsGroup[i]=n_PIbin;
    }

    write_rsp(LowEnergy,HighEnergy,NumberGroups,FirstChannelGroup,NumberChannelsGroup,total_response,n_Ebin,n_PIbin,(char *)febounds_name.c_str());

	//Write response file name in pha header

	
    // Close the caldb pixresp file 
    fits_close_file(fpixresp,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

	return(EXIT_SUCCESS);
}


int cztrspgen::write_rsp(float *Elo,float *Ehi,int *ngrp,int *fchan,int *nchan,float **rspmatrix,int n_Ebin,int n_PIbin,char *febounds_name)
{
    int status=0,i;

    LOG(INFO)<<"Reached writing routine";

    if (create_empty_fitsfile(rspfile, "rspTemplate")) {
        LOG(ERROR) << "Error in creating response file from corresponding template";
        return EXIT_FAILURE;
    }

    fitsfile *frsp,*febounds;

    fits_open_file(&frsp, rspfile, READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening response file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    // Open the caldb ebounds file
    fits_open_file(&febounds,febounds_name,READONLY,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_movnam_hdu(febounds, BINARY_TBL, "EBOUNDS", 0, &status);

    // Read the ebounds file into the vectors defining the ChannelLowEnergy and ChannelHighEnergy

    float PImin[n_PIbin],PImax[n_PIbin];

    int channel[n_PIbin];

    for (i=0;i<n_PIbin;i++) channel[i]=i;

//   fits_write_col(fout,TDOUBLE,1,i+1,1,1,&(gti_final[i].tstart),&status); 

    try{
    fits_read_col(febounds,TFLOAT,2,1,1,n_PIbin,NULL,&PImin,NULL,&status);

    fits_read_col(febounds,TFLOAT,3,1,1,n_PIbin,NULL,&PImax,NULL,&status);

    fits_movnam_hdu(frsp, BINARY_TBL, "EBOUNDS", 0, &status);

    fits_write_col(frsp,TINT,1,1,1,n_PIbin,channel,&status);
    fits_write_col(frsp,TFLOAT,2,1,1,n_PIbin,PImin,&status);
    fits_write_col(frsp,TFLOAT,3,1,1,n_PIbin,PImax,&status);

    fits_movnam_hdu(frsp, BINARY_TBL, "MATRIX", 0, &status);


    for(i=0;i<n_Ebin;i++)
    {
    fits_write_col(frsp, TFLOAT, 1, i+1,1, 1,&Elo[i], &status);
    fits_write_col(frsp, TFLOAT, 2, i+1,1, 1,&Ehi[i], &status);
    fits_write_col(frsp, TINT, 3, i+1,1, 1,&ngrp[i], &status);
    fits_write_col(frsp, TINT, 4, i+1,1,1,&fchan[i], &status);
    fits_write_col(frsp, TINT, 5, i+1,1,1,&nchan[i], &status);
    fits_write_col(frsp, TINT, 5, i+1,1,1,&nchan[i], &status);
    fits_write_col(frsp, TFLOAT, 6, i+1,1, n_PIbin,rspmatrix[i], &status);
    }


    } catch(ErrorHandler errHandler){
        logError(errHandler);
        return EXIT_FAILURE;
    }

    // Close the caldb ebounds
    fits_close_file(febounds,&status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }


    fits_close_file(frsp,&status);
    if (status) {
        LOG(ERROR) << "Error in closing response file.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    return(EXIT_SUCCESS);
}


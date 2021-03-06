#include "gtiHandler.h"

GTIhandler::GTIhandler() {

}

/* *
 * Edited by Mithun to include reading gtiextname
 * and column names changed to start and stop 
 * */
int GTIhandler::read_gti_file(string gtiFilename,string gtiextname){
    int status=0;
    ErrorHandler errHandler;
    fitsfile *fgti;
    string gtiFilepath="";
    vector <double> tstart;
    vector <double> tstop;
    long nrows=0;
    if(gtiFilename=="" && this->gtiFilename==""){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "No GTI filename provided by user.";
        throw errHandler;
    }
    gtiFilepath=(gtiFilename=="")?this->gtiFilename:gtiFilename;
   

    fits_open_file(&fgti, (char*)gtiFilepath.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error opening gti file: " + gtiFilepath ;
        throw errHandler;
    }
   
    fits_movnam_hdu(fgti, BINARY_TBL, (char *)gtiextname.c_str(), 0, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in moving to HDU GTI.";
        throw errHandler;
    }
   
    if(read_fits_column(fgti, "START", TDOUBLE, 1, 1 , -1, tstart)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error reading start column of GTI file " + gtiFilepath;
        throw errHandler;
    }
   

    if(read_fits_column(fgti, "STOP", TDOUBLE, 1, 1 , -1, tstop)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error reading stop column of GTI file " + gtiFilepath;
        throw errHandler;
    }
   
    fits_close_file(fgti, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error closing fits file: " + gtiFilepath;
        throw errHandler;
    }
    
    set_GTI(tstart, tstop);
   
    return EXIT_SUCCESS;
}

//SETTERS
void GTIhandler::set_GTI(vector<double> tstart, vector<double> tstop){
    int nstart=0, nstop=0; //number of start and stop records
    int i=0;
    ErrorHandler errHandler;
    
    nstart=tstart.size();
    nstop=tstop.size();
    
    if(nstart!=nstop){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Size of tstart and tstop vectors do not match.";
        throw errHandler;
    }
    
    GTI.clear();
    GTI.resize(nstart);
    for(i=0; i<nstart; i++){
        GTI[i].tstart = tstart[i];
        GTI[i].tstop = tstop[i];
    }
}


//Function to remove reduntant entries from the GTI vector
//Mithun(08/12/15)

int GTIhandler::cleanup_GTI()
{
    long i;

    
    vector <double> new_tstart;
    vector <double> new_tstop;
    for(i=0;i<GTI.size();i++)
    {
	    if(GTI[i].tstart>=GTI[i].tstop) // THIS SHOULD ACTUALLY BE ==. For bypassing errors it is set to >=
	    {
		    GTI.erase(GTI.begin()+i);
		    i-=1;
	    }
    }

    long nrows=GTI.size();


    new_tstart.push_back(GTI[0].tstart);

    for(i=0;i<nrows-1;i++)
    {
        if(GTI[i].tstop!=GTI[i+1].tstart) 
        {
            new_tstop.push_back(GTI[i].tstop);
            new_tstart.push_back(GTI[i+1].tstart);
        }

    }

    new_tstop.push_back(GTI[nrows-1].tstop);

    set_GTI(new_tstart,new_tstop);

    return(EXIT_SUCCESS);    
}

// Function to write to GTI file

int GTIhandler::write_gti_file(fitsfile *fgti,char *extname)
{
    int status=0;
    long nrows=GTI.size(),i;

    fits_movnam_hdu(fgti, BINARY_TBL, extname, 0, &status);
    
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in moving to extensions in GTI file";
        return (EXIT_FAILURE);
    }


    double exposure_time=0;
   
    for(i=0;i<nrows;i++)
    {
        fits_write_col(fgti,TDOUBLE,1,i+1,1,1,&GTI[i].tstart,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in writing start column";
            return (EXIT_FAILURE);
        }

        fits_write_col(fgti,TDOUBLE,2,i+1,1,1,&GTI[i].tstop,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in writing stop column***";
            return (EXIT_FAILURE);
        }

    	exposure_time+=GTI[i].tstop-GTI[i].tstart;
	}

    fits_update_key(fgti, TLONG, "NAXIS2",&nrows, NULL, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in updating naxis2 for GTI***";
        return (EXIT_FAILURE);
    }

	//Update exposure time

    fits_update_key(fgti,TDOUBLE,"EXPOSURE",&exposure_time, NULL, &status);
    if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

    cztHeaderParam headKey;
    headKey.writeTimekey(GTI[0].tstart,GTI[GTI.size()-1].tstop,fgti);

    return(EXIT_SUCCESS); 
}


//Independent functions

/* Function re-written by Rakesh (13/12/15)
   based on the earlier version of the function
 * */
int find_intersecting_range(vector<GTIrecord> gti1, vector<GTIrecord> gti2, vector<GTIrecord>* gtiout)
{
    long nsize1=0,nsize2=0;
    
    nsize1=gti1.size();
    nsize2=gti2.size();

    bool gtiValidityFlag=true;
    vector<GTIrecord> mergedVec;

    GTIrecord intersec;
    int i,j,k,l,i1,i2;
    double index[4],subgti[4];
    

    if(nsize1>0){
        for(i1=0; i1<nsize1; i1++){
            if(gti1[i1].tstart >gti1[i1].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            LOG(ERROR)<<"GTI1 tstart should be less than tstop";
            return(EXIT_FAILURE);
        }
    }
    if(nsize2>0){
        for(i2=0; i2<nsize2; i2++){
            if(gti2[i2].tstart >gti2[i2].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            LOG(ERROR)<<"GTI2 tstart should be less than tstop";
            return(EXIT_FAILURE);
        }
    }


    if(nsize1==0 && nsize2==0){
        LOG(ERROR)<<"Empty GTIs.... exiting..";
        return(EXIT_FAILURE);
    }
    else if(nsize1==0)
    {
        *gtiout = gti2;
        return(EXIT_SUCCESS);
    } else if(nsize2==0)
    {
        *gtiout=gti1;
        return(EXIT_SUCCESS);
    } 

    for(i=0;i<gti1.size();i++)
       {
            for(j=0;j<gti2.size();j++)
            {
                subgti[0]=gti1[i].tstart;
                subgti[1]=gti1[i].tstop;
                subgti[2]=gti2[j].tstart;
                subgti[3]=gti2[j].tstop;
                
                index[0]=0;
                index[1]=1;
                index[2]=0;
                index[3]=1;

                for(k=0;k<4;k++)
                {
                    double gtiTemp=0;
                    int indTemp=0;

                    for(l=k+1;l<4;l++)
                    {
                        if(subgti[k]>subgti[l])
                        {
                        gtiTemp=subgti[k];
                        indTemp=index[k];
                        subgti[k]=subgti[l];
                        index[k]=index[l];
                        subgti[l]=gtiTemp;
                        index[l]=indTemp;
                        }
                    }
                }

                    if(index[1]==0 && index[2]==1)
                    {
                        intersec.tstart=subgti[1];
                        intersec.tstop=subgti[2];
                        mergedVec.push_back(intersec);
                    }                        

            }
        
       }

        
    //FIX to Remove the reduntant entries with tstart==tstop
        for(i=0;i<mergedVec.size();i++)
        {
            if(mergedVec[i].tstart==mergedVec[i].tstop)
            {
                mergedVec.erase(mergedVec.begin()+i);
                i-=1;
            }
        }


    *gtiout=mergedVec;

    return(EXIT_SUCCESS);
}    


/*

void find_intersecting_range(vector<GTIrecord> gti1, vector<GTIrecord> gti2, vector<GTIrecord>* gtiout){
    long nsize1=0; //size of gti record vector 1
    long nsize2=0; //size of gti record vector 2
    long i1=0, i2=0, i=0;
    ErrorHandler errHandler;
    GTIrecord intersection;
    GTIrecord merged;
    vector<GTIrecord> intersectionVec;
    vector<GTIrecord> mergedVec;
    bool stopFlag=true;
    bool gtiValidityFlag=true;
    nsize1=gti1.size();
    nsize2=gti2.size();

    //Checking whether tstart and tstop of gti1 are in proper format i.e. tstart <=tstop
    if(nsize1>0){
        for(i1=0; i1<nsize1; i1++){
            if(gti1[i1].tstart >gti1[i1].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            errHandler.severity = errERROR;
            errHandler.errorStatus =IMPROPER_GTI;
            errHandler.errorMsg = "GTI1tstart value is greater than tstop.";
            throw errHandler;
        }
    }
    if(nsize2>0){
        for(i2=0; i2<nsize2; i2++){
            if(gti2[i2].tstart >gti2[i2].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            errHandler.severity = errERROR;
            errHandler.errorStatus =IMPROPER_GTI;
            errHandler.errorMsg = "GTI2 tstart value is greater than tstop.";
            throw errHandler;
        }
    }
    if(nsize1==0 && nsize2==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = NO_OUTPUT_POSSIBLE;
        errHandler.errorMsg = "Error in generating intersecting GTI range as both GTI inputs are empty";
        throw errHandler;
    } else if(nsize1==0){
        *gtiout = gti2;
    } else if(nsize2==0){
        *gtiout=gti1;
    } else{
        //Finding GTI intersection range
    
        intersectionVec.clear();
        for(i1=0; i1<nsize1; i1++){
            for(i2=0; i2<nsize2; i2++){
                intersection.tstart = max(gti1[i1].tstart, gti2[i2].tstart);
                intersection.tstop = min(gti1[i1].tstop, gti2[i2].tstop);
                if(intersection.tstop<intersection.tstart){
                    continue;
                }
                intersectionVec.push_back(intersection);
            }
        }
        sort(intersectionVec.begin(), intersectionVec.end(), min_tstart_value);

        for(i=0; i<(intersectionVec.size()-1); i++){
            if(stopFlag==true){
                merged.tstart=intersectionVec[i].tstart;
                stopFlag=false;
            }
            if(intersectionVec[i].tstop<intersectionVec[i+1].tstart){
                merged.tstop=intersectionVec[i].tstop;
                mergedVec.push_back(merged);
                stopFlag=true;
            }
        }
        
        merged.tstop = intersectionVec[intersectionVec.size()-1].tstop;
        mergedVec.push_back(merged);

       //FIX to Remove the reduntant entries with tstart==tstop
        for(i=0;i<mergedVec.size();i++)
        {
            if(mergedVec[i].tstart==mergedVec[i].tstop)
            {
                mergedVec.erase(mergedVec.begin()+i);
                i-=1;
            }
        }
        // End of fix

        *gtiout = mergedVec;
    }
    
}

//compare function to sort gti
bool min_tstart_value(GTIrecord r1, GTIrecord r2){
    return r1.tstart < r2.tstart;
}
*/



//original code 


//commented by mayuri
/*#include "gtiHandler.h"

GTIhandler::GTIhandler() {

}*/

/* *
 * Edited by Mithun to include reading gtiextname
 * and column names changed to start and stop 
 * */
 
//commented by mayuri 
 /*
int GTIhandler::read_gti_file(string gtiFilename,string gtiextname){
    int status=0;
    ErrorHandler errHandler;
    fitsfile *fgti;
    string gtiFilepath="";
    vector <double> tstart;
    vector <double> tstop;
    long nrows=0;
    if(gtiFilename=="" && this->gtiFilename==""){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "No GTI filename provided by user.";
        throw errHandler;
    }
    gtiFilepath=(gtiFilename=="")?this->gtiFilename:gtiFilename;
   

    fits_open_file(&fgti, (char*)gtiFilepath.c_str(), READONLY, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error opening gti file: " + gtiFilepath ;
        throw errHandler;
    }
   
    fits_movnam_hdu(fgti, BINARY_TBL, (char *)gtiextname.c_str(), 0, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in moving to HDU GTI.";
        throw errHandler;
    }
   
    if(read_fits_column(fgti, "START", TDOUBLE, 1, 1 , -1, tstart)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error reading start column of GTI file " + gtiFilepath;
        throw errHandler;
    }
   

    if(read_fits_column(fgti, "STOP", TDOUBLE, 1, 1 , -1, tstop)){
        errHandler.severity = errERROR;
        errHandler.errorStatus = TABLE_READ_ERROR;
        errHandler.errorMsg = "Error reading stop column of GTI file " + gtiFilepath;
        throw errHandler;
    }
   
    fits_close_file(fgti, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error closing fits file: " + gtiFilepath;
        throw errHandler;
    }
    
    set_GTI(tstart, tstop);
   
    return EXIT_SUCCESS;
}

//SETTERS
void GTIhandler::set_GTI(vector<double> tstart, vector<double> tstop){
    int nstart=0, nstop=0; //number of start and stop records
    int i=0;
    ErrorHandler errHandler;
    
    nstart=tstart.size();
    nstop=tstop.size();
    
    if(nstart!=nstop){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Size of tstart and tstop vectors do not match.";
        throw errHandler;
    }
    
    GTI.clear();
    GTI.resize(nstart);
    for(i=0; i<nstart; i++){
        GTI[i].tstart = tstart[i];
        GTI[i].tstop = tstop[i];
    }
}


//Function to remove reduntant entries from the GTI vector
//Mithun(08/12/15)

int GTIhandler::cleanup_GTI()
{
    long i;

    
    vector <double> new_tstart;
    vector <double> new_tstop;
    for(i=0;i<GTI.size();i++)
    {
	    if(GTI[i].tstart>=GTI[i].tstop) // THIS SHOULD ACTUALLY BE ==. For bypassing errors it is set to >=
	    {
		    GTI.erase(GTI.begin()+i);
		    i-=1;
	    }
    }

    long nrows=GTI.size();


    new_tstart.push_back(GTI[0].tstart);

    for(i=0;i<nrows-1;i++)
    {
        if(GTI[i].tstop!=GTI[i+1].tstart) 
        {
            new_tstop.push_back(GTI[i].tstop);
            new_tstart.push_back(GTI[i+1].tstart);
        }

    }

    new_tstop.push_back(GTI[nrows-1].tstop);

    set_GTI(new_tstart,new_tstop);

    return(EXIT_SUCCESS);    
}

// Function to write to GTI file

int GTIhandler::write_gti_file(fitsfile *fgti,char *extname)
{
    int status=0;
    long nrows=GTI.size(),i;

    fits_movnam_hdu(fgti, BINARY_TBL, extname, 0, &status);
    
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in moving to extensions in GTI file";
        return (EXIT_FAILURE);
    }


    double exposure_time=0;
   
    for(i=0;i<nrows;i++)
    {
        fits_write_col(fgti,TDOUBLE,1,i+1,1,1,&GTI[i].tstart,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in writing start column";
            return (EXIT_FAILURE);
        }

        fits_write_col(fgti,TDOUBLE,2,i+1,1,1,&GTI[i].tstop,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in writing stop column***";
            return (EXIT_FAILURE);
        }

    	exposure_time+=GTI[i].tstop-GTI[i].tstart;
	}

    fits_update_key(fgti, TLONG, "NAXIS2",&nrows, NULL, &status);
    if(status){
        fits_report_error(stderr,status);
        LOG(ERROR)<<"***Error in updating naxis2 for GTI***";
        return (EXIT_FAILURE);
    }

	//Update exposure time

    fits_update_key(fgti,TDOUBLE,"EXPOSURE",&exposure_time, NULL, &status);
    if(status) {fits_report_error(stderr,status); return(EXIT_FAILURE);}

    cztHeaderParam headKey;
    headKey.writeTimekey(GTI[0].tstart,GTI[GTI.size()-1].tstop,fgti);

    return(EXIT_SUCCESS); 
}*/


//Independent functions

/* Function re-written by Rakesh (13/12/15)
   based on the earlier version of the function
 * */
 
 //commented by mayuri
/*int find_intersecting_range(vector<GTIrecord> gti1, vector<GTIrecord> gti2, vector<GTIrecord>* gtiout)
{
    long nsize1=0,nsize2=0;
    
    nsize1=gti1.size();
    nsize2=gti2.size();

    bool gtiValidityFlag=true;
    vector<GTIrecord> mergedVec;

    GTIrecord intersec;
    int i,j,k,l,i1,i2;
    double index[4],subgti[4];
    

    if(nsize1>0){
        for(i1=0; i1<nsize1; i1++){
            if(gti1[i1].tstart >gti1[i1].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            LOG(ERROR)<<"GTI1 tstart should be less than tstop";
            return(EXIT_FAILURE);
        }
    }
    if(nsize2>0){
        for(i2=0; i2<nsize2; i2++){
            if(gti2[i2].tstart >gti2[i2].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            LOG(ERROR)<<"GTI2 tstart should be less than tstop";
            return(EXIT_FAILURE);
        }
    }


    if(nsize1==0 && nsize2==0){
        LOG(ERROR)<<"Empty GTIs.... exiting..";
        return(EXIT_FAILURE);
    }
    else if(nsize1==0)
    {
        *gtiout = gti2;
        return(EXIT_SUCCESS);
    } else if(nsize2==0)
    {
        *gtiout=gti1;
        return(EXIT_SUCCESS);
    } 

    for(i=0;i<gti1.size();i++)
       {
            for(j=0;j<gti2.size();j++)
            {
                subgti[0]=gti1[i].tstart;
                subgti[1]=gti1[i].tstop;
                subgti[2]=gti2[j].tstart;
                subgti[3]=gti2[j].tstop;
                
                index[0]=0;
                index[1]=1;
                index[2]=0;
                index[3]=1;

                for(k=0;k<4;k++)
                {
                    double gtiTemp=0;
                    int indTemp=0;

                    for(l=k+1;l<4;l++)
                    {
                        if(subgti[k]>subgti[l])
                        {
                        gtiTemp=subgti[k];
                        indTemp=index[k];
                        subgti[k]=subgti[l];
                        index[k]=index[l];
                        subgti[l]=gtiTemp;
                        index[l]=indTemp;
                        }
                    }
                }

                    if(index[1]==0 && index[2]==1)
                    {
                        intersec.tstart=subgti[1];
                        intersec.tstop=subgti[2];
                        mergedVec.push_back(intersec);
                    }                        

            }
        
       }

        
    //FIX to Remove the reduntant entries with tstart==tstop
        for(i=0;i<mergedVec.size();i++)
        {
            if(mergedVec[i].tstart==mergedVec[i].tstop)
            {
                mergedVec.erase(mergedVec.begin()+i);
                i-=1;
            }
        }


    *gtiout=mergedVec;

    return(EXIT_SUCCESS);
}    
*/
//already commented
/*

void find_intersecting_range(vector<GTIrecord> gti1, vector<GTIrecord> gti2, vector<GTIrecord>* gtiout){
    long nsize1=0; //size of gti record vector 1
    long nsize2=0; //size of gti record vector 2
    long i1=0, i2=0, i=0;
    ErrorHandler errHandler;
    GTIrecord intersection;
    GTIrecord merged;
    vector<GTIrecord> intersectionVec;
    vector<GTIrecord> mergedVec;
    bool stopFlag=true;
    bool gtiValidityFlag=true;
    nsize1=gti1.size();
    nsize2=gti2.size();

    //Checking whether tstart and tstop of gti1 are in proper format i.e. tstart <=tstop
    if(nsize1>0){
        for(i1=0; i1<nsize1; i1++){
            if(gti1[i1].tstart >gti1[i1].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            errHandler.severity = errERROR;
            errHandler.errorStatus =IMPROPER_GTI;
            errHandler.errorMsg = "GTI1tstart value is greater than tstop.";
            throw errHandler;
        }
    }
    if(nsize2>0){
        for(i2=0; i2<nsize2; i2++){
            if(gti2[i2].tstart >gti2[i2].tstop){
                gtiValidityFlag=false;
                break;
            }
        }
        if(gtiValidityFlag==false){
            errHandler.severity = errERROR;
            errHandler.errorStatus =IMPROPER_GTI;
            errHandler.errorMsg = "GTI2 tstart value is greater than tstop.";
            throw errHandler;
        }
    }
    if(nsize1==0 && nsize2==0){
        errHandler.severity = errERROR;
        errHandler.errorStatus = NO_OUTPUT_POSSIBLE;
        errHandler.errorMsg = "Error in generating intersecting GTI range as both GTI inputs are empty";
        throw errHandler;
    } else if(nsize1==0){
        *gtiout = gti2;
    } else if(nsize2==0){
        *gtiout=gti1;
    } else{
        //Finding GTI intersection range
    
        intersectionVec.clear();
        for(i1=0; i1<nsize1; i1++){
            for(i2=0; i2<nsize2; i2++){
                intersection.tstart = max(gti1[i1].tstart, gti2[i2].tstart);
                intersection.tstop = min(gti1[i1].tstop, gti2[i2].tstop);
                if(intersection.tstop<intersection.tstart){
                    continue;
                }
                intersectionVec.push_back(intersection);
            }
        }
        sort(intersectionVec.begin(), intersectionVec.end(), min_tstart_value);

        for(i=0; i<(intersectionVec.size()-1); i++){
            if(stopFlag==true){
                merged.tstart=intersectionVec[i].tstart;
                stopFlag=false;
            }
            if(intersectionVec[i].tstop<intersectionVec[i+1].tstart){
                merged.tstop=intersectionVec[i].tstop;
                mergedVec.push_back(merged);
                stopFlag=true;
            }
        }
        
        merged.tstop = intersectionVec[intersectionVec.size()-1].tstop;
        mergedVec.push_back(merged);

       //FIX to Remove the reduntant entries with tstart==tstop
        for(i=0;i<mergedVec.size();i++)
        {
            if(mergedVec[i].tstart==mergedVec[i].tstop)
            {
                mergedVec.erase(mergedVec.begin()+i);
                i-=1;
            }
        }
        // End of fix

        *gtiout = mergedVec;
    }
    
}

//compare function to sort gti
bool min_tstart_value(GTIrecord r1, GTIrecord r2){
    return r1.tstart < r2.tstart;
}
*/

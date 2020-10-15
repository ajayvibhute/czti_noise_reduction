
#include"cztgtigenv2.h"

using namespace std;

cztgtigen::cztgtigen(){
    
    strcpy(modulename,"cztgtigen_v");
    strcat(modulename,VERSION);
    strcpy(thresholdfile,"");
    strcpy(mkffile,"");
}

int cztgtigen::read(int argc, char** argv){
    int status=0;
    
    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error Initializing PIL***";
        return status;
    }
/*    
    LOG(INFO)<<"--------------------------------------";
    LOG(INFO)<<"Description for run mode parameter for CZTGTIGEN";
    LOG(INFO)<<"--------------------------------------";
    LOG(INFO)<<"It can have 4 values- copy/mkf/custom/custom-refined";
    LOG(INFO)<<"1) 'mkf'          - In mkf mode, it creates a gti based on the valid ranges given in ";
    LOG(INFO)<<"                    threshold file for each mkf parameter";
    LOG(INFO)<<"2) 'custom'       - In custom mode, it creates gti based on intervals given by user";
    LOG(INFO)<<"3) 'mkf-refined'  - In this mode it refines the gti produced in mkf mode with level2 gti";
    LOG(INFO)<<"4) 'custom-refined'- In this mode it refines the gti produced in custom mode with level2 gti";
    LOG(INFO)<<"5) 'copy' - copies input file to output file";
    LOG(INFO);

*/
    if(PIL_OK!=(status=PILGetFname("eventfile",eventfile))){
        LOG(ERROR)<<"***Error Reading eventfile name***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("mkffile",mkffile))){
        LOG(ERROR)<<"***Error reading mkffile file:"<<mkffile<<"***";
        return status;
    }
    if(PIL_OK!=(status=PILGetFname("thresholdfile",thresholdfile))){
        LOG(ERROR)<<"***Error reading threshold file:"<<thresholdfile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("outfile",outfile))){
        LOG(ERROR)<<"***Error reading outfile:"<<outfile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetFname("usergtifile",usergtifile))){
        LOG(ERROR)<<"***Error reading threshold file:"<<usergtifile<<"***";
        return status;
    } 
    
    if(PIL_OK!=(status=PILGetBool("history",&history))){
        LOG(ERROR)<<"***Error reading history parameter***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetBool("clobber",&clobber))){
        LOG(ERROR)<<"***Error reading clobber parameter***";
        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

int cztgtigen::read(char *eventfile,char *mkffile,char *thresholdfile,char *outfile, char *usergtifile,int clobber,int history)
{    
    strcpy(this->eventfile,eventfile);
    strcpy(this->mkffile,mkffile);
    strcpy(this->thresholdfile,thresholdfile);
    strcpy(this->outfile,outfile);
	strcpy(this->usergtifile,usergtifile);
    this->clobber=clobber;
    this->history=history;
    return (EXIT_SUCCESS);
}

void cztgtigen::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                           CZTGTIGEN PARAMETERS                               ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
	LOG(INFO)<<"Modulename         : "<<modulename;
    LOG(INFO)<<"Eventfile          : "<<eventfile;
    LOG(INFO)<<"MKF file           : "<<mkffile;     
    LOG(INFO)<<"MKF Threshold file : "<<thresholdfile;
    LOG(INFO)<<"Output gti file    : "<<outfile;      
	LOG(INFO)<<"User gti file      : "<<usergtifile; 
    LOG(INFO)<<"Clobber            : "<<clobber;
    LOG(INFO)<<"History            : "<<history;
    LOG(INFO)<<"----------------------------------------------------------------------------";
}

int cztgtigen::cztgtigenProcess(){

    if(clobber==YES){
        unlink(outfile);
    }
    else{
        LOG(ERROR)<<""<<outfile<<"  already exists";
        LOG(ERROR)<<" Use the option 'clobber=yes' to overwrite the file";
        return (EXIT_FAILURE);
    }
    
    vector<GTIrecord> mkfgti; 


    fitsfile *fgti,*fout,*fevt,*fuser;
    int status=0;
    
        GTIrecord objgti;
        mkfgti.clear();
        Mkf mkf;
        vector <double> tstart, tstop;
        tstart.clear(); tstop.clear();
        if(mkf.read_mkf_file((string) mkffile)){
            LOG(ERROR) << "Error in reading mkf file: " << mkffile;
            return EXIT_FAILURE;
        }

	LOG(INFO)<<"Using threadshold file:  "<<thresholdfile;
        if(mkf.get_filtered_gti(thresholdfile, tstart, tstop)){
            LOG(ERROR) << "Error in generating mkf filtered gti";
            return EXIT_FAILURE;
        }

        for(int itime=0; itime<tstart.size(); itime++){
            objgti.tstart= tstart[itime];
            objgti.tstop= tstop[itime];
            mkfgti.push_back(objgti);
        }
    
        LOG(INFO) << "Completed GTI creation based on mkf";

    	if(tstart.size() == 0) 
	{
		LOG(ERROR)<<"No Good time intervals found with MKF thresholds... EXITING....";
		return(EXIT_FAILURE);
	}	


        fits_open_file(&fevt, eventfile, READONLY, &status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in opening event file";
            return (EXIT_FAILURE);
        }

        fits_create_file(&fgti,outfile,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in creating GTI file***";
            return (EXIT_FAILURE);
        }

        
        fits_copy_file(fevt,fgti,1,0,0,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in copying file***";
            return (EXIT_FAILURE);
        }
        

        char extname[20];
        int i,j;
		GTIhandler usergti;       
		int isusergti=0; 
		
		if(strcmp(usergtifile,"-")==0)
			isusergti=0;
		else
		{
			if (FileExists(usergtifile)) 	
				isusergti=1;
			else
			{
				LOG(ERROR)<<"*****User defined GTI file doesn't exist. Enter - if not required*****";
				return EXIT_FAILURE;
			}
		}


        if(isusergti)
        {
            LOG(INFO)<<"User GTI is also being read";
			sprintf(extname,"GTI");
			usergti.read_gti_file((string)usergtifile,(string)extname);
        }
        
        for(i=-1;i<4;i++)
        {
            vector <GTIrecord> gtiout;
            gtiout.clear();
            GTIhandler l2gti;
            GTIhandler GTIfinal;

            if(i==-1)
                sprintf(extname,"GTI");
            else
                sprintf(extname,"Q%d_GTI",i);

            fits_movnam_hdu(fevt, BINARY_TBL, extname, 0, &status);
                if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in moving to EXTNAME GTI in event file***";
                return (EXIT_FAILURE);
            }
  
            fits_copy_header(fevt, fgti,  &status);
            if(status){
                fits_report_error(stderr,status);
                LOG(ERROR)<<"***Error in moving to EXTNAME GTI in event file***";
                return (EXIT_FAILURE);
            }
                
            l2gti.read_gti_file((string)eventfile,(string)extname);
           
			if(isusergti)
			{
	            vector <GTIrecord> gtiout1;
    	        gtiout1.clear();		
				find_intersecting_range(mkfgti,l2gti.get_GTI(),&gtiout1);
				find_intersecting_range(gtiout1,usergti.get_GTI(),&gtiout);
			}
			else
				find_intersecting_range(mkfgti,l2gti.get_GTI(),&gtiout);

            GTIfinal.set_GTI(gtiout);

            GTIfinal.write_gti_file(fgti,extname);
        }

        fits_close_file(fgti,&status);
        if(status){
            fits_report_error(stderr,status);
            LOG(ERROR)<<"***Error in closing GTI ***";
            return (EXIT_FAILURE);
        }


        if(history==YES){
            vector<string> vhistory;
            getHistory(vhistory);
            if(writeHistory(outfile,vhistory)){
              LOG(ERROR)<<"Error writing history";
              return (EXIT_FAILURE);
            }
        }

     updateKeywords(outfile,modulename);

    return (EXIT_SUCCESS);
}

/*
int refineGTI(vector<GTI> &gti2refine,vector<GTI> &ingti){
    int s1=gti2refine.size();
    int s2=ingti.size();
    if(s1<=0){
        LOG(ERROR)<<"Number of intervals are zero/less than zero in GTI";
        return (EXIT_FAILURE);
    }
    if(s2<=0){
        LOG(INFO)<<"Number of intervals in level2 GTI is zero";
        return (EXIT_SUCCESS);
    }
    int i,j;
   
    vector<GTI> outgti;
    
    GTI objgti;
    
    for(i=0;i<s1;i++){
        for(j=0;j<s2;j++){
            if(gti2refine[i].tstart>=ingti[j].tstart && gti2refine[i].tstop<=ingti[j].tstop)
                objgti.tstart=gti2refine[i].tstart;  
            //cout<<gti2refine[i].tstart<<" "<<ingti[i].tstart;
            if(gti2refine[i].tstart>ingti[j].tstart)
                objgti.tstart=gti2refine[i].tstart;
            else
                objgti.tstart=ingti[j].tstart;
            //cout<<gti2refine[i].tstop<<" "<<ingti[i].tstop;
            if(gti2refine[i].tstop<ingti[j].tstop)
                objgti.tstop=gti2refine[i].tstop;
            else
                objgti.tstop=ingti[j].tstop;
            outgti.push_back(objgti);
            //cout<<objgti.tstart<<" - "<<objgti.tstop;
            
        }
    }
    
    gti2refine.clear();
    vector<GTI>::iterator iter=outgti.begin();
    while(iter!=outgti.end()){
        if((*iter).tstart<=(*iter).tstop){
            gti2refine.push_back((*iter));
            LOG(INFO)<<(*iter).tstart<<"-"<<(*iter).tstop;
        }
        iter++;    
    }
   return (EXIT_SUCCESS);
}  

*/
void cztgtigen::getHistory(vector<string> &vhistory){
    //char *user=getlogin();
    strcpy(modulename,"cztgtigen_v");
    strcat(modulename,VERSION);

    char *user = getenv("USER");
	string str="Module run by "+(string)user;
    vhistory.push_back(str);
    vhistory.push_back("Parameter List START for "+(string)modulename);
    vhistory.push_back("P1 eventfile="+(string)eventfile);
    vhistory.push_back("P2 mkffile="+(string)mkffile);
    vhistory.push_back("P3 thresholdfile="+(string)thresholdfile);
    vhistory.push_back("P4 outfile="+(string)outfile);
    if(clobber==YES) 
        vhistory.push_back("P5 clobber=yes");
    else
        vhistory.push_back("P5 clobber=no");
    if(history==YES)
        vhistory.push_back("P6 history=yes");
    else
        vhistory.push_back("P6 history=no");
    vhistory.push_back("Parameter List END");
}





#include "cztflagbadpix.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztrspgen.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);
    		
    
    int status=0;
    time_t st, et;
    checkPFILESenv();
    checkParFile("cztflagbadpix");
    
    cztflagbadpix par;
    if(par.read(argc,argv)){
       LOG(ERROR)<<"***Error reading input parameters***";
       return (EXIT_FAILURE);        
    }
    
    par.display();
    LOG(INFO)<<"CZTFLAGBADPIX STARTED....................";

    st = time(NULL);
    status=par.cztflagbadpix_process();
    if(status){
        LOG(ERROR)<<"***Error in cztflagbadpix process***";
        exit(EXIT_FAILURE);
    }
    LOG(INFO)<<"CZTFLAGBADPIX COMPLETED SUCCESSFULLY";
    et = time(NULL);
    LOG(INFO) << "Time Elapsed:" << et - st << " sec";
    
    return (EXIT_SUCCESS);
 }


#include <cstdlib>
#include <ctime>
#include "cztgtigenv2.h"
#include "glog/logging.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztgtigen.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);
    		

    time_t st, et;
    checkPFILESenv();
    checkParFile("cztgtigen");
    
    cztgtigen par;
    if(par.read(argc,argv)){
        LOG(ERROR)<<"***Error reading parameters***";
        return (EXIT_FAILURE);
    }
    
    par.display();
    st = time(NULL);

    LOG(INFO)<<"CZTGTIGEN STARTED....................";
    int status=0;
    status=par.cztgtigenProcess();
    if(status){
        LOG(ERROR)<<"***Error in cztgtigen_process***";
        exit(EXIT_FAILURE);
    }
   LOG(INFO)<<"CZTGTIGEN COMPLETED SUCCESSFULLY";
   et=time(NULL);
   LOG(INFO)<<"Time Elapsed:"<< et-st << " seconds";
   return (EXIT_SUCCESS);
 }


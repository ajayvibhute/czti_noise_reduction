
#include <cstdlib>
#include "cztdpigen_v2.h"

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztdpigen.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);
    		

    time_t st, et;
    
    checkPFILESenv();
    checkParFile("cztdpigen");
    
    Cztdpigen par;
    par.read(argc,argv);
    par.display();
    st = time(NULL); 
    LOG(INFO)<<"CZTDPIGEN STARTED....................";
    int status=0;
    status=par.cztdpigenProcess();
    if(status){
        LOG(ERROR)<<"***Error in cztdpigen process***";
        exit(EXIT_FAILURE);
    }
    LOG(INFO)<<"CZTDPIGEN COMPLETED SUCCESSFULLY";
    et = time(NULL);
    LOG(INFO) << "Time Elapsed---" << et - st << " seconds" << endl;
    return (EXIT_SUCCESS);
 }


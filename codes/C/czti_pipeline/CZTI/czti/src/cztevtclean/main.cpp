
#include <cstdlib>
#include "cztevtclean.h"

using namespace std;

/*Wrapper function for CZTEVTCLEAN module
 * 
 */
int main(int argc, char** argv) {
        char const* log_file=std::getenv("logfile_env");
	time_t st,et;

    if ( log_file == NULL){
		log_file="cztevtclean.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);
    		

    
    checkPFILESenv();
    checkParFile("cztevtclean");
    
    cztevtclean param_var;
    param_var.read(argc,argv);
    param_var.display();
	 st=time(NULL);
    LOG(ERROR)<<"CZTEVTCLEAN STARTED....................";
    int status=0;
    status=param_var.cztevtcleanProcess();
    if(status){
        LOG(ERROR)<<"***Error in cztevtclean_process()***";
        exit(EXIT_FAILURE);
    }
    LOG(INFO)<<"CZTEVTCLEAN COMPLETED SUCCESSFULLY";
    et=time(NULL);
    LOG(INFO)<<"Time Elapsed:"<<et-st<<" seconds";

    return (EXIT_SUCCESS);
 }


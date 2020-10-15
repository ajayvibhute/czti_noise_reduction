/* main.cpp cztpixclean


Mithun NPS
(30/11/15)
*/

#include <ctime>
#include "cztpixclean.h"

using namespace std;

int main(int argc, char** argv) 
{
        char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztpixclean.log";
    }    
    google::InitGoogleLogging(log_file);
google::SetStderrLogging(google::INFO);
    		
			
    time_t st,et;

    checkPFILESenv();
    checkParFile("cztpixclean");

    cztpixclean par;

    if(par.read(argc,argv)){
        LOG(ERROR)<<"***Error reading parameters***";
        return (EXIT_FAILURE);
    }
    par.display();
	
    LOG(INFO)<<"CZTPIXCLEAN STARTED....................";
	st=time(NULL);
				
    int status=0;
    status=par.cztpixcleanProcess();

    if(status){
        LOG(ERROR)<<"***Error in cztpixclean process***";
        exit(EXIT_FAILURE);
    }

    LOG(INFO)<<"CZTPIXCLEAN COMPLETED SUCCESSFULLY";
    et=time(NULL);
    LOG(INFO)<<"Time Elapsed:"<<et-st<<" seconds";
    return (EXIT_SUCCESS);

}

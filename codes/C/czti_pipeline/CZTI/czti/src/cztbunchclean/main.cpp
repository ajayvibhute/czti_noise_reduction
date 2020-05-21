/* main.cpp cztbunchclean


Mithun NPS
(30/11/15)
*/

#include <ctime>
#include "cztbunchcleanv2.h"

using namespace std;

int main(int argc, char** argv) 
{
  
    //log_file=std::getenv("GLOG_log_dir")+std::getenv("logfile_env");
    //cout<<"log_file: "<<log_file;
    //google::SetLogDestination(google::INFO,log_file.c_str() );
    char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztbunchclean.log";
    }    
    google::InitGoogleLogging(log_file);

    google::SetStderrLogging(google::INFO);
		
    time_t st,et;

    checkPFILESenv();
    checkParFile("cztbunchclean");

    cztbunchclean par;

    if(par.read(argc,argv)){
        LOG(ERROR)<<"***Error reading parameters***";
        return (EXIT_FAILURE);
    }
    par.display();


    LOG(INFO)<<"CZTBUNCHCLEAN STARTED....................";
    st=time(NULL);
			
    int status=0;
    status=par.cztbunchcleanProcess();

    if(status){
        LOG(ERROR)<<"***Error in cztbunchclean process***";
        exit(EXIT_FAILURE);
    }

    LOG(INFO)<<"CZTBUNCHCLEAN COMPLETED SUCCESSFULLY\n";
    et=time(NULL);
    LOG(INFO)<<"Time Elapsed:"<<et-st<<" seconds\n";
    return (EXIT_SUCCESS);

}

/* main.cpp cztrspgen


Mithun NPS
(10/12/15)
*/

#include <ctime>
#include "cztrspgen.h"

using namespace std;

int main(int argc, char** argv) 
{
//    string log_file;
    //log_file=std::getenv("GLOG_log_dir")+std::getenv("logfile_env");
    //cout<<"log_file: "<<log_file;
    //google::SetLogDestination(google::INFO,log_file.c_str() );
    char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztrspgen.log";
    }    
    google::InitGoogleLogging(log_file);
    google::SetStderrLogging(google::INFO);
    		
    time_t st,et;

    checkPFILESenv();
    checkParFile("cztrspgen");

    cztrspgen par;

    if(par.read(argc,argv)){
        LOG(ERROR)<<"***Error reading parameters***";
        return (EXIT_FAILURE);
    }

	par.display();

    LOG(INFO)<<"CZTRSPGEN STARTED....................";
    st=time(NULL);
			
    int status=0;
    status=par.cztrspgenProcess();

    if(status){
        LOG(ERROR)<<"***Error in cztrspgen process***";
        exit(EXIT_FAILURE);
    }

    LOG(INFO)<<"CZTRSPGEN COMPLETED SUCCESSFULLY";
    et=time(NULL);
    LOG(INFO)<<"Time Elapsed:"<<et-st<<" seconds";
    return (EXIT_SUCCESS);

}

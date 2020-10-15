/* 
 * File:   main.cpp
 * Author: preeti
 *
 * Created on October 15, 2012, 6:11 PM
 */

#include <cstdlib>
#include <ctime>
#include "cztbindata.h"


using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    //log_file=std::getenv("GLOG_log_dir")+std::getenv("logfile_env");
    //cout<<"log_file: "<<log_file;
    //google::SetLogDestination(google::INFO,log_file.c_str() );
    char const* log_file=std::getenv("logfile_env");
    if ( log_file == NULL){
		log_file="cztbindata.log";
    }    
    google::InitGoogleLogging(log_file);

    google::SetStderrLogging(google::INFO);
    time_t st,et;
        
    checkPFILESenv();
    checkParFile("cztbindata");
    
    cztbindata par;
    if(par.read(argc,argv)){
        LOG(ERROR)<<"***Error reading parameters***";
        return (EXIT_FAILURE);
    }
    par.display();
    st=time(NULL);
	
    LOG(INFO)<<"CZTBINDATA STARTED....................";
    int status=0;
    status=par.cztbindataProcess();
    if(status){
        LOG(ERROR)<<"***Error in cztbindata_process***";
        exit(EXIT_FAILURE);
    }
    LOG(INFO)<<"CZTBINDATA COMPLETED SUCCESSFULLY";
    et=time(NULL);
    LOG(INFO)<<"Time Elapsed:"<<et-st<<" seconds";
    return (EXIT_SUCCESS);
 }


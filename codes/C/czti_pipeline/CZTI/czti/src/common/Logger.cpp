/* 
 * File:   Logger.cpp
 * Author: ajay
 * 
 * Created on August 2, 2012, 4:14 PM
 */

#include "Logger.h"
#include<iostream>
#include<string>
#include <stdlib.h>
#include "log4cpp/Category.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/PropertyConfigurator.hh"

using namespace std;

Logger::Logger() {}
log4cpp::Category& Logger::getLogger(string log4jcpppropfile)
{
        try 
        {
                log4cpp::PropertyConfigurator::configure(log4jcpppropfile);
        }    
        catch(log4cpp::ConfigureFailure& f)
        {
                cout << "Configure File " << f.what() <<endl;
                exit(0);
        }
        log4cpp::Category& root  = log4cpp::Category::getRoot();
        log4cpp::Category& log = log4cpp::Category::getInstance(std::string("log"));
        return log;
}
void Logger::shutdown()
{
    log4cpp::Category::shutdown();
}
Logger::Logger(const Logger& orig) {
}

Logger::~Logger() {
}


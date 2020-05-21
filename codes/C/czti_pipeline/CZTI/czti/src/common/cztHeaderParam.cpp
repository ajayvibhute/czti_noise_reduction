
#include "cztHeaderParam.h"

//cztHeaderParam (Header parameters of science,event,spec,lc,img etc)

//Set values of header parameters
int cztHeaderParam::set_default_values()
{
    this->mjdrefi=55197;
    this->mjdreff=0;
    this->timesys="UTC";
    this->equinox=2000;
    this->radecsys="ICRS";
    this->timeunit="s";
    this->mission="ASTROSAT";
    this->telescop="ASTROSAT";
    this->instrume="CZTI";
 
    return(EXIT_SUCCESS);
}


int cztHeaderParam::getTelapse(double &telapse)
{
	telapse=tstop-tstart;
	return(EXIT_SUCCESS);
}


int cztHeaderParam::writeTimekey(double tstart, double tstop,fitsfile *fptr)
{
	
    int status=0;
    
    this->tstart=tstart;
    this->tstop=tstop;

    this->tstarti=(long)(tstart);
    this->tstopi=(long)(tstop);

    this->tstartf=tstart-tstarti;
    this->tstopf=tstop-tstopi;


   fits_update_key(fptr, TLONG, "TSTARTI",&tstarti, "Start time of observation Integer part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
 //           return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TDOUBLE, "TSTARTF",&tstartf, "Start time of observation Fractional part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
   //         return (EXIT_FAILURE);
        }


    fits_update_key(fptr, TLONG, "TSTOPI",&tstopi, "Stop time of observation Integer part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
     //       return (EXIT_FAILURE);
        }


    fits_update_key(fptr, TDOUBLE, "TSTOPF",&tstopf, "Stop time of observation Fractional part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
       //     return (EXIT_FAILURE);
        }

        fits_update_key(fptr, TDOUBLE, "TSTART",&tstart, "Start time of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TDOUBLE, "TSTOP",&tstop,"Stop time of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
	
	double telapse=tstop-tstart;

        fits_update_key(fptr, TDOUBLE, "TELAPSE",&telapse, "Elapsed time", &status);
	if (status) {
	           LOG(ERROR) <<"Error in updating key.";
	              fits_report_error(stderr, status);
	}

    return(EXIT_SUCCESS);
}


int cztHeaderParam::readFromHeader(fitsfile *fptr)
{
    int status=0;
    this->mjdrefi=55197;
    this->mjdreff=0;
    this->timesys="UTC";
    this->equinox=2000;
    this->radecsys="ICRS";
    this->timeunit="s";
    this->mission="ASTROSAT";
    this->telescop="ASTROSAT";
    this->instrume="CZTI";

/*
    readFitsKeyStr(fptr,"MISSION",&mission,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

    readFitsKeyStr(fptr,"TELESCOP",&telescop,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

    readFitsKeyStr(fptr,"INSTRUME",&instrume,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

    readFitsKeyStr(fptr,"ORIGIN",&origin,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }
    */
	status=0;
    fits_read_key(fptr,TLONG,"TSTARTI",&tstarti,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	tstarti=0;
	}

	status=0;    
    fits_read_key(fptr,TLONG,"TSTOPI",&tstopi,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	tstopi=0;
	}

	status=0;    
    fits_read_key(fptr,TDOUBLE,"TSTARTF",&tstartf,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	tstartf=0.;
	}

    status=0; 
    fits_read_key(fptr,TDOUBLE,"TSTOPF",&tstopf,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	tstopf=0.;
	}

    status=0; 
	fits_read_key(fptr,TDOUBLE,"TSTART",&tstart,NULL, &status);
    if (status) {LOG(ERROR) <<"Error in reading key.";
    
    //return (EXIT_FAILURE);
    	tstart=0.0;
	}
    
	status=0; 
    fits_read_key(fptr,TDOUBLE,"TSTOP",&tstop,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	tstop=0.0;
	}

/*
    readFitsKeyStr(fptr,"TIMESYS",&timesys,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

    readFitsKeyStr(fptr,"TIMEUNIT",&timeunit,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }
*/
	status=0; 
    readFitsKeyStr(fptr,"OBJECT",&object,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }
   
	status=0; 
	fits_read_key(fptr,TFLOAT,"RA_PNT",&ra_pnt,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	ra_pnt=-999.0;
	}

    status=0;
    fits_read_key(fptr,TFLOAT,"DEC_PNT",&dec_pnt,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	dec_pnt=-999;
	}

	status=0; 
	readFitsKeyStr(fptr,"OBS_ID",&obs_id,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
	}

 	status=0; 
 	readFitsKeyStr(fptr,"OBS_MODE",&obs_mode,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

	status=0; 
    readFitsKeyStr(fptr,"DATE-OBS",&date_obs,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

	status=0; 
	readFitsKeyStr(fptr,"TIME-OBS",&time_obs,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

	status=0; 
    readFitsKeyStr(fptr,"DATE-END",&date_end,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }

	status=0;     
	readFitsKeyStr(fptr,"TIME-END",&time_end,NULL);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    }
    
	status=0; 
    fits_read_key(fptr,TDOUBLE,"TIMEDEL",&timedel,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	timedel=2.0e-5;
	}
	
	status=0; 
	fits_read_key(fptr,TDOUBLE,"EXPOSURE",&exp_time,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key.";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	exp_time=0.0;
	}

	status=0;    
   	fits_read_key(fptr,TDOUBLE,"TELAPSE",&telapse,NULL, &status);
    if (status) {
        LOG(ERROR) <<"Error in reading key telapse";
        fits_report_error(stderr, status);
     //   return (EXIT_FAILURE);
    	telapse=0.0;
	}

	//LOG(INFO)<<"EXP TIME"<<exp_time<<"\t"<<"TIME DEL"<<timedel<<"QUADID\t"<<quadid;
    LOG(INFO)<<"Read all the header keywords"<<endl; 
    return(EXIT_SUCCESS);
	
}

int cztHeaderParam::writeToHeader(fitsfile *fptr)
{
    int status=0;

    fits_update_key(fptr, TSTRING, "MISSION",(char *)mission.c_str(), "Name of the mission/satellite", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "TELESCOP",(char *)telescop.c_str(), "Name of the mission/satellite", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }
    fits_update_key(fptr, TSTRING, "INSTRUME",(char *)instrume.c_str(), "Name of the instrument/detector", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "ORIGIN",(char *)origin.c_str(), "Source of FITS FILE", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TLONG, "TSTARTI",&tstarti, "Start time of observation Integer part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
 //           return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TDOUBLE, "TSTARTF",&tstartf, "Start time of observation Fractional part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
   //         return (EXIT_FAILURE);
        }


    fits_update_key(fptr, TLONG, "TSTOPI",&tstopi, "Stop time of observation Integer part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
     //       return (EXIT_FAILURE);
        }


    fits_update_key(fptr, TDOUBLE, "TSTOPF",&tstopf, "Stop time of observation Fractional part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
       //     return (EXIT_FAILURE);
        }


        fits_update_key(fptr, TSTRING, "TIMESYS",(char *)timesys.c_str(),"Time is UTC", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

        fits_update_key(fptr, TSTRING, "TIMEUNIT",(char *)timeunit.c_str(), "Time is in seconds", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

       fits_update_key(fptr, TLONG, "MJDREFI",&mjdrefi, "MJDREF Integer part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }


        fits_update_key(fptr, TLONG, "MJDREFF",&mjdreff,"MJDREF Fractional part", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }
     fits_update_key(fptr, TDOUBLE, "TSTART",&tstart, "Start time of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TDOUBLE, "TSTOP",&tstop,"Stop time of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "OBJECT",(char *)object.c_str(), "Target name", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }


        fits_update_key(fptr, TSTRING, "RADECSYS",(char *)radecsys.c_str(),"Reference frame", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }


        fits_update_key(fptr, TINT, "EQUINOX",&equinox, "J2000", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TFLOAT, "RA_PNT",&ra_pnt, "Nominal Pointing RA", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TFLOAT, "DEC_PNT",&dec_pnt, "Nominal Pointing DEC", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }
    fits_update_key(fptr, TSTRING, "OBS_ID",(char *)obs_id.c_str(), "Observation ID", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }


    fits_update_key(fptr, TSTRING, "OBS_MODE",(char *)obs_mode.c_str(), NULL, &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "DATE-OBS",(char *)date_obs.c_str(), "Start date of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "TIME-OBS",(char *)time_obs.c_str(), "Start time of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "DATE-END",(char *)date_end.c_str(), "End date of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

    fits_update_key(fptr, TSTRING, "TIME-END",(char *)time_end.c_str(), "End time of observation", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }


        fits_update_key(fptr, TDOUBLE, "TIMEDEL",&timedel, "Time resolution", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }


        fits_update_key(fptr, TDOUBLE, "EXPOSURE",&exp_time, "Exposure time", &status);
        if (status) {
            LOG(ERROR) <<"Error in updating key.";
            fits_report_error(stderr, status);
         //   return (EXIT_FAILURE);
        }

        fits_update_key(fptr, TDOUBLE, "TELAPSE",&telapse, "Elapsed time", &status);
    	if (status) {
               LOG(ERROR) <<"Error in updating key Telapse";
                  fits_report_error(stderr, status);
    	}

  LOG(INFO)<<"Copied all the header keywords sucessfully"<<endl;      
    return(EXIT_SUCCESS);
}

// cztHeaderParam ends



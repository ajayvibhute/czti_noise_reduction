
#include <libgen.h>


#include "utils.h"

#include"utils.h"

#define PI 3.1415926535


/*Function to update livetime based on start and stop time of an interval.
* This sssumes that live time is computed for 1s bins.
* It can add or subtract from the the livetime bins based on add flag
* If tstart >= tstop, this routine does nothing
*
* Mithun N P S(14/12/15)
* 
* Edits:

	* Generalized for any binsize for live time calculation (30/08/16)
	 
**/
int updatelivetime(double gtstart, double gtstop,double *livetime,long tstarti, float fraction,bool add,long ntbins,double livetime_binsize)
{

    long gtstartbin=floor((gtstart-tstarti)/livetime_binsize);
    long gtstopbin=floor((gtstop-tstarti)/livetime_binsize);

	//LOG(INFO)<<gtstartbin<<"  "<<gtstopbin;

	long i;
   
    if(gtstart==0||gtstop==0||gtstartbin<0||gtstopbin<0||gtstartbin>=ntbins||gtstopbin>=ntbins)
    {
        //LOG(WARNING)<<"Negative bins. time less than tstarti";
        return(EXIT_SUCCESS);
    }
    else if(gtstop<=gtstart)
    {
        //LOG(INFO)<<"tstop less than or equal to tstart";
        return EXIT_SUCCESS;
    }
    else if(gtstartbin==gtstopbin)
    {
        if(add)
            livetime[gtstartbin]+=(gtstop-gtstart)*fraction;//*(livetime[gtstartbin]/livetime_binsize);
        else
            livetime[gtstartbin]-=(gtstop-gtstart)*fraction;//*(livetime[gtstartbin]/livetime_binsize);
    }
    else
    {
		//LOG(INFO)<<"DTIME "<<gtstop-gtstart;

		//LOG(INFO)<<gtstartbin<<"  BINS "<<gtstopbin; 
        if(add)
        {
            livetime[gtstartbin]+=((ceil((gtstart-tstarti)/livetime_binsize)*livetime_binsize)-(gtstart-tstarti))*fraction;//*(livetime[gtstartbin]/livetime_binsize);
            livetime[gtstopbin]+=(gtstop-tstarti-(floor((gtstop-tstarti)/livetime_binsize)*livetime_binsize))*fraction;//*(livetime[gtstopbin]/livetime_binsize);

            for(i=gtstartbin+1;i<gtstopbin;i++)
                livetime[i]+=livetime_binsize*fraction;
        }
        else 
        {
            livetime[gtstartbin]-=((ceil((gtstart-tstarti)/livetime_binsize)*livetime_binsize)-(gtstart-tstarti))*fraction;//*(livetime[gtstartbin]/livetime_binsize);
            livetime[gtstopbin]-=(gtstop-tstarti-(floor((gtstop-tstarti)/livetime_binsize)*livetime_binsize))*fraction;//*(livetime[gtstopbin]/livetime_binsize);
        
            for(i=gtstartbin+1;i<gtstopbin;i++)
                livetime[i]-=livetime_binsize*fraction;
        }
    
    }

    return(EXIT_SUCCESS);
}

/*
Function to update livetime based on bad time interval start and stop times 
and prior value of livetime. It accepts livetime values of all bins and 
correct the output livetime taking into account the prior livetime as well.

This would be used in cztpixclean

Mithun (31/08/16) 
*/

int re_updatelivetime(double gtstart, double gtstop,double *inlivetime,double *outlivetime,long tstarti, float fraction,bool add,long ntbins,double livetime_binsize)
{

    long gtstartbin=floor((gtstart-tstarti)/livetime_binsize);
    long gtstopbin=floor((gtstop-tstarti)/livetime_binsize);

    //LOG(INFO)<<gtstartbin<<"  "<<gtstopbin;

    long i;

    if(gtstart==0||gtstop==0||gtstartbin<0||gtstopbin<0||gtstartbin>=ntbins||gtstopbin>=ntbins)
    {
        //LOG(WARNING)<<"Negative bins. time less than tstarti";
        return(EXIT_SUCCESS);
    }
    else if(gtstop<=gtstart)
    {
        //LOG(INFO)<<"tstop less than or equal to tstart";
        return EXIT_SUCCESS;
    }
    else if(gtstartbin==gtstopbin)
    {
        if(add)
            outlivetime[gtstartbin]+=(gtstop-gtstart)*fraction*(inlivetime[gtstartbin]/livetime_binsize);
        else
            outlivetime[gtstartbin]-=(gtstop-gtstart)*fraction*(inlivetime[gtstartbin]/livetime_binsize);
    }
    else
    {
        //LOG(INFO)<<"DTIME "<<gtstop-gtstart;

        //LOG(INFO)<<gtstartbin<<"  BINS "<<gtstopbin; 
        if(add)
        {
            outlivetime[gtstartbin]+=((ceil((gtstart-tstarti)/livetime_binsize)*livetime_binsize)-(gtstart-tstarti))*fraction*(inlivetime[gtstartbin]/livetime_binsize);
            outlivetime[gtstopbin]+=(gtstop-tstarti-(floor((gtstop-tstarti)/livetime_binsize)*livetime_binsize))*fraction*(inlivetime[gtstopbin]/livetime_binsize);

            for(i=gtstartbin+1;i<gtstopbin;i++)
			{
                outlivetime[i]+=livetime_binsize*fraction*(inlivetime[i]/livetime_binsize);
				if(outlivetime[i]<0) outlivetime[i]=0;
			}
        }
        else
        {
            outlivetime[gtstartbin]-=((ceil((gtstart-tstarti)/livetime_binsize)*livetime_binsize)-(gtstart-tstarti))*fraction*(inlivetime[gtstartbin]/livetime_binsize);
            outlivetime[gtstopbin]-=(gtstop-tstarti-(floor((gtstop-tstarti)/livetime_binsize)*livetime_binsize))*fraction*(inlivetime[gtstopbin]/livetime_binsize);

            for(i=gtstartbin+1;i<gtstopbin;i++)
			{
                outlivetime[i]-=livetime_binsize*fraction*(inlivetime[i]/livetime_binsize);
				if(outlivetime[i]<0) outlivetime[i]=0;
			}
        }

    }


	if(outlivetime[gtstartbin]<0) outlivetime[gtstartbin]=0;
	if(outlivetime[gtstopbin]<0) outlivetime[gtstopbin]=0;
    return(EXIT_SUCCESS);
}



/* Updates the livetime with a given GTI 
 * Assumes that the given live time fraction is uniform over one 
 * second time bin and applies GTI to it.
 *
 * Mithun(15/12/15)
 * */
int updatelivetime_gti(fitsfile *fgti,double *timearray, double *livetime,long ntbins,double livetime_binsize)
{

    int status=0;
    long i,nrows;
    long startbin,stopbin;
    int k;
    long btstartbin,btstopbin;

    fits_get_num_rows(fgti, &nrows, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    double gti_tstart[nrows],gti_tstop[nrows];

    fits_read_col(fgti, TDOUBLE, 1, 1, 1, nrows, NULL, gti_tstart,NULL, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    fits_read_col(fgti, TDOUBLE, 2, 1, 1, nrows, NULL, gti_tstop,NULL, &status);
    if(status) { fits_report_error(stderr,status);  return (EXIT_FAILURE); }

    long tstart=(long)timearray[0];
    long tstop=(long)timearray[ntbins-1]+1;
    
    double btstart=tstart;
    double btstop;
    LOG(INFO)<<"Reached gti livtime";

    for(i=0;i<=nrows;i++)
    {
        if(i<nrows) btstop=gti_tstart[i];
		else btstop=tstop;

        btstartbin=floor((btstart-tstart)/livetime_binsize);
        btstopbin=floor((btstop-tstart)/livetime_binsize);

        if(btstartbin<0||btstopbin<0|| btstartbin>=ntbins||btstopbin>=ntbins)
            continue;

        if(btstartbin==btstopbin)
		    livetime[btstartbin]=livetime[btstartbin]-(btstop-btstart);
        else
        {
            for(k=btstartbin+1;k<btstopbin;k++)
                livetime[k]=0;
            livetime[btstartbin]-=((ceil((btstart-tstart)/livetime_binsize)*livetime_binsize)-(btstart-tstart));
            livetime[btstopbin]-=((btstop-tstart)-(floor((btstop-tstart)/livetime_binsize)*livetime_binsize));  
        }
        btstart=gti_tstop[i];
    
		if(livetime[btstartbin]<0) livetime[btstartbin]=0;
		if(livetime[btstopbin]<0) livetime[btstopbin]=0;

	}

    return(EXIT_SUCCESS);
}


//Maths vector class end
int copyKeywords(fitsfile *fptr1, fitsfile *fptr2, int n, ...) {
    //cout<<"\n----------copyKeywords----------\n";
    va_list keywords;
    va_start(keywords, n);

    int status = 0;
    char record[FLEN_CARD], *key;

    for (int i = 0; i < n; i++) {
        key = va_arg(keywords, char *);

        fits_read_card(fptr1, key, record, &status);
        if (status == 202) {
            LOG(ERROR) << "***" << key << " not found in header***";
        }
        if (status) {
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

        fits_write_record(fptr2, record, &status);
        if (status) {
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }

    }
    va_end(keywords);
    return (EXIT_SUCCESS);
}
//------------------------------------------------------------------------------

int copyUserKeyWords(fitsfile *fin, fitsfile *fout) {
    int status = 0, keyexist;
    char record[FLEN_CARD];
    char str[3];
    char filename[FLEN_FILENAME];
    fits_file_name(fin, filename, &status);
    if (status) return (EXIT_FAILURE);
    fits_get_hdrspace(fin, &keyexist, NULL, &status);
    if (status) {
        LOG(ERROR) << "***Could not find number of keywords in file " << filename << "***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    int keyclass;
    DLOG(INFO) << "Number of keywords found: "<<keyexist;
    for (int i = 1; i <= keyexist; i++) {
        fits_read_record(fin, i, record, &status);
        if (status) {
            LOG(ERROR) << "***Error in reading record number " << i << " in file " << filename << "***";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        keyclass = fits_get_keyclass(record);
        if (keyclass == TYP_COMM_KEY)
            continue;
        else if (keyclass == TYP_USER_KEY || keyclass == TYP_REFSYS_KEY || keyclass == TYP_WCS_KEY) {
            fits_write_record(fout, record, &status);
            if (status) {
                LOG(ERROR) << "***Error in writing keyword " << record << "***";
                fits_report_error(stderr, status);
                return (EXIT_FAILURE);
            }
        }
    }
    return (EXIT_SUCCESS);
}

int copyUserKeyWords(fitsfile *fin, fitsfile *fout, string inHDU, 
        string outHDU, vector <string> vecKeywords){
    int status=0;
    int i=0; //counter variable
    int nkeys=0;
    char record[FLEN_CARD];
    char filenameIn[FLEN_FILENAME];
    char filenameOut[FLEN_FILENAME];
    
    //Reading filenameIn
    fits_file_name(fin, filenameIn, &status);
    if (status) return (EXIT_FAILURE);
    //Reading filenameOut
    fits_file_name(fout, filenameOut, &status);
    if (status) return (EXIT_FAILURE);
    
    nkeys = vecKeywords.size();
    
    //Moving to required HDU in input file
    if(inHDU=="Primary"){
        fits_movabs_hdu(fin, 1, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << inHDU << " in file " << filenameIn;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    else {
        fits_movnam_hdu(fin, ANY_HDU, (char*) inHDU.c_str(), NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << inHDU << " in file " << filenameIn;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    //Moving to required HDU in output file
    if(inHDU=="Primary"){
        fits_movabs_hdu(fout, 1, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << outHDU << " in file " << filenameOut;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }
    else {
        fits_movnam_hdu(fout, ANY_HDU, (char*) outHDU.c_str(), NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << outHDU << " in file " << filenameOut;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    }

    //Reading records from input file and updating in output file
    for(i=0; i<nkeys; i++){
        fits_read_card(fin, (char*) vecKeywords[i].c_str(), record, &status);
        if(status){
            LOG(WARNING)<<"Error in reading record with key " << vecKeywords[i] << " from " << filenameIn;
        }
        fits_update_card(fout, (char*) vecKeywords[i].c_str(), record, &status);
        if(status){
            LOG(WARNING)<<"Error in updating record with key " << vecKeywords[i] << " in " << filenameOut;
        }
    }
    return (EXIT_SUCCESS);
}

//-----------------------------------------------------------------------------

void joinStrings(char *outstring, int n, ...) {
    va_list strings;
    va_start(strings, n);
    char *str;
    int len = 0;
    for (int i = 0; i < n; i++) {
        str = va_arg(strings, char *);
        len = len + strlen(str);
        if (i == 0) strcpy(outstring, str);
        else strcat(outstring, str);
    }
    outstring[len] = '\0';
    va_end(strings);
}
//-------------------------------------------------------------------------

int deleteDir(char *dir) {

    vector<string> directories;
    //cout<<"\nDir:"<<dir;
    struct dirent **namelist;
    int n, n_copy;
    char temp[2048];
    n = scandir(dir, &namelist, defaultfilter, alphasort);
    n_copy = n;
    //cout<<endl<<"n:"<<n;
    if (n < 0)
        return (EXIT_FAILURE);
    else if (n == 0) {
        LOG(INFO) << "1.Removing " << dir;
        if (remove(dir)) {
            LOG(ERROR) << "***1. Error in removing " << dir << "***";
            return (EXIT_FAILURE);
        }
    } else {
        while (n--) {
            strcpy(temp, dir);
            strcat(temp, "/");
            strcat(temp, namelist[n]->d_name);
            //cout<<endl<<temp;
            if (namelist[n]->d_type == DT_REG) {
                //cout<<endl<<"Deleting "<<temp;
                if (unlink(temp)) {
                    LOG(ERROR) << "***Error in removing " << temp << "***";
                    return (EXIT_FAILURE);
                }
            }/*else if((namelist[n]->d_type==DT_DIR) && 
                    (strcmp(namelist[n]->d_name,".")!=0) && 
                    (strcmp(namelist[n]->d_name,"..")!=0)){*/
            else if (namelist[n]->d_type == DT_DIR)
                deleteDir(temp);
            //cout<<endl<<"Namelist N: "<<namelist[n]->d_name;

            free(namelist[n]);
        }
        free(namelist);
    }

    if (DirExists(dir)) {
        LOG(INFO) << "2.Removing " << dir;
        if (remove(dir)) {
            LOG(ERROR) << "***2. Error in removing " << dir << "***";
            return (EXIT_FAILURE);
        }
    }
    return (EXIT_SUCCESS);
}
//------------------------------------------------------------------------------

int defaultfilter(const struct dirent *dptr) {
    int retval = 1;
    if (strcmp(dptr->d_name, ".") == 0) retval = 0;
    if (strcmp(dptr->d_name, "..") == 0) retval = 0;
    return retval;
}
//------------------------------------------------------------------------------

int parseString(string str, char delim, vector<string> &substr) {
    int pos = 0;
    string temp;
    int len = str.size();
    if (len <= 0) {
        LOG(ERROR) << "***Invalid input string" << str << "***";
        return (-1);
    }
    //cout<<"\nString Passed is "<<str;
    //cout<<"\nLength:"<<len;
    while (len > 0) {
        //cout<<endl<<"Length:"<<len;
        pos = str.find(delim, 0);
        //cout<<endl<<"pos:"<<pos;
        if (pos < 0) {
            substr.push_back(str);
            break;
        } else {
            temp = str.substr(0, pos);
            str = str.substr(pos + 1);
            substr.push_back(temp);
            len = len - pos;
            //cout<<"\nLength:"<<len;
        }
    }
    return (substr.size());
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

int getbins(char *bins, double *start, double *end, int nebins) {
    //cout<<"\nInside getbins\n";
    int i = 0;
    int temp;
    if (strcmp(bins, "-") == 0) {
        return 2;
    } else {
        string ebinranges = (string) bins;
        //string sets[MAX_BINS];
        string temp, temp2;
        //cout<<endl<<ebinranges<<endl;
        i = 0;
        int pos = 0;
        for (i = 0; i < nebins; i++) {
            //cout<<endl<<ebinranges;
            pos = ebinranges.find('-', 0);
            temp = ebinranges.substr(0, pos);
            temp2 = ebinranges.substr(pos + 1);
            ebinranges.clear();
            ebinranges.assign(temp2);
            start[i] = atof(temp.c_str());
            temp.clear();
            temp2.clear();

            pos = ebinranges.find(',', pos);
            temp = ebinranges.substr(0, pos);
            temp2 = ebinranges.substr(pos + 1);
            ebinranges.clear();
            ebinranges.assign(temp2);
            end[i] = atof(temp.c_str());
            temp.clear();
            temp2.clear();
        }
    }
    return (EXIT_SUCCESS);
}


int parseCommaSeparatedData(char *data, int *dataArray, int dataElements) {
    string Sdata = (string) data;
    return 0;
}
//------------------------------------------------------------------------------

int get_nbins(char *binstr) {
    int nbins = 0;
    if (strcmp(binstr, "-") == 0) {
        nbins = 1;
    } else {
        char bins[2000];
        strcpy(bins, binstr);
        char *token = strtok(bins, ",");
        while (token != NULL) {
            token = strtok(NULL, ",");
            nbins++;
        }
    }
    return nbins;
}

//--------------Quaternion related-----------

/*
 * 
 */

//Member Functions of struct Q

Q::Q() {
    q1 = 0;
    q2 = 0;
    q3 = 0;
    q4 = 0;
    norm = false;
    AxisAngle = false;
    qflag = false;
    theta = 0;
    x = 0;
    y = 0;
    z = 0;
}

Q::Q(double a, double b, double c, double d) {
    q1 = a;
    q2 = b;
    q3 = c;
    q4 = d;
    norm = false;
    AxisAngle = false;
}

Q::Q(const Q &q) {
    //cout<<endl<<"\nInside COPY CONSTRUCTOR\n";
    q1 = q.q1;
    q2 = q.q2;
    q3 = q.q3;
    q4 = q.q4;
    theta = q.theta;
    x = q.x;
    y = q.y;
    z = q.z;
    norm = false;
    AxisAngle = false;
}

void Q::readQ() {
    LOG(INFO) << "Enter q1,q2,q3,q4 [q1 is scalar]:";
    cin >> q1 >> q2 >> q3 >> q4;
}

void Q::readAxisAngle() {
    LOG(INFO) << "Enter theta, x, y and z:";
    cin >> theta >> x >> y >> z;
}

int Q::getAxisAngle() {
    //cout<<"\nInside getAxisAngle\n";
    if (AxisAngle == false) {
        double denominator = sqrt(1 - q1 * q1);
        theta = 2 * acos(q1);
        //theta=ceil(theta-0.5);  //rounding the theta angle
        if (denominator != 0) {
            x = q2 / denominator;
            y = q3 / denominator;
            z = q4 / denominator;
        } else {
            x = 1;
            y = 0;
            z = 0;
        }
    }
    AxisAngle = true;
    return (EXIT_SUCCESS);
}

int Q::getQuat() {
    if (qflag == false) {
        q1 = cos(theta / 2);
        q2 = x * sin(theta / 2);
        q3 = y * sin(theta / 2);
        q4 = z * sin(theta / 2);
        qflag = true;
    }

}

void Q::display() {
    LOG(INFO) << "(" << q1 << ")+(" << q2 << ")i+(" << q3 << ")j+(" << q4 << ")k";
    if (AxisAngle == true) {
        LOG(INFO) << "Angle=" << theta * (180 / M_PI) << " deg\t";
        LOG(INFO) << "Axis: x=" << x << "   y=" << y << "  z=" << z;
    }

}

void Q::normalize() {

    if (norm == false) {
        mod = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);
        //cout<<"Denominator:"<<denominator;
        if (mod != 0) {
            q1 = q1 / mod;
            q2 = q2 / mod;
            q3 = q3 / mod;
            q4 = q4 / mod;
            norm = true;
        }
    }
}
//-------------------------------------------------------------------
//Member functions of axis

void Axis::normalize() {
    mod = sqrt(x * x + y * y + z * z);
    if (mod != 0) {
        x = x / mod;
        y = y / mod;
        z = z / mod;
    }
    norm = true;
}

void Axis::display() {
    LOG(INFO) << "(" << x << ")i + (" << y << ")j + (" << z << ")k";
}

double Axis::getMod() {
    return (sqrt(x * x + y * y + z * z));
}
//------------------------------------------------------------------------------

void matrix_product(float **A, float **B, float **C, int a, int b, int c) {
    int i, j, k;
    float sum;
    for (i = 0; i < a; i++) {
        for (j = 0; j < c; j++) {
            sum = 0;
            for (k = 0; k < b; k++) {
                sum = sum + (A[i][k] * B[k][j]);
                C[i][j] = sum;
            }
        }
    }
}
//------------------------------------------------------------------------------

void quaternion_product(Q &Q1, Q &Q2, Q &Q3) {
    // LOG(INFO) << "Q1:";
    //Q1.display();
    // LOG(INFO) << "Q2:";
    //Q2.display();
    Q3.q1 = Q1.q1 * Q2.q1 - Q1.q2 * Q2.q2 - Q1.q3 * Q2.q3 - Q1.q4 * Q2.q4;
    Q3.q2 = Q1.q1 * Q2.q2 + Q1.q2 * Q2.q1 + Q1.q3 * Q2.q4 - Q1.q4 * Q2.q3;
    Q3.q3 = Q1.q1 * Q2.q3 - Q1.q2 * Q2.q4 + Q1.q3 * Q2.q1 + Q1.q4 * Q2.q2;
    Q3.q4 = Q1.q1 * Q2.q4 + Q1.q2 * Q2.q3 - Q1.q3 * Q2.q2 + Q1.q4 * Q2.q1;
    //LOG(INFO) << "Q3 : ";
    //Q3.display();
}
//------------------------------------------------------------------------------

void Inverse(Q &q, Q &inverseq) {
    //cout<<"\nInside inverse\n";
    inverseq.q1 = q.q1;
    inverseq.q2 = (-1) * q.q2;
    inverseq.q3 = (-1) * q.q3;
    inverseq.q4 = (-1) * q.q4;
}

//---------------------------------------------------------------------------------

int rotate(Axis &axis, Q &q, Axis &axisrot) {
    
    // axisrot = q*axis*qinv
    LOG(INFO) << "Transforming CZTI axis to inertial frame using the attitude quaternion obtained from attitude file.";
    Q qinv;
    Inverse(q, qinv);
    Q temp;
    Q v1(0, axis.x, axis.y, axis.z), v2;

    quaternion_product(v1, qinv, temp);

    quaternion_product(q, temp, v2);
    axisrot.x = v2.q2;
    axisrot.y = v2.q3;
    axisrot.z = v2.q4;
    return (EXIT_SUCCESS);
}
//-----------------------------------------------------------------------------------------

int getCC(char *aspectfile, double RA, double DEC, double *thetax, double *thetay) {
    // Converting RA & DEC of the source into thetaX and thetaY [Camera Co-ordinates]
    fitsfile *fptr;
    int status = 0;
    Q qc;
    //cout<<"\nInertial Coordinates-----RA:"<<RA<<"  DEC:"<<DEC;
    fits_open_file(&fptr, aspectfile, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC1", &qc.q1, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC2", &qc.q2, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC3", &qc.q3, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC4", &qc.q4, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //cout<<"\nQC from average aspect file : ";
    //qc.display();
    //cout<<endl;
    Q qc_inv;
    Inverse(qc, qc_inv);
    Axis inertial(sin(RA) * cos(DEC), sin(DEC), cos(RA) * cos(DEC)); // to be changed to inertial(cos(RA)*cos(DEC),sin(RA)*cos(DEC),sin(DEC))
    // if CZTI pointing is along y axis
    //cout<<"\nInertial Vector:";  inertial.display();
    Axis Nczti;
    rotate(inertial, qc_inv, Nczti);

    //cout<<"\nRotated vector:";  Nczti.display();

    *thetax = atan(Nczti.x / Nczti.z);
    *thetay = atan(Nczti.y / Nczti.z);
    //cout<<"\nCamera Coodinates-------Theta_x:"<<*thetax<<"  Theta_y:"<<*thetay;
    //cout<<"\nNczti:";   Nczti.display();
    //cout<<"\ntan theta_x:"<<(Nczti.x/Nczti.z)<<"\ntan theta_y:"<<(Nczti.y/Nczti.z)<<"\n";
    //cout<<"\n("<<RA<<","<<DEC<<")------>("<<*thetax<<","<<*thetay<<")";
    return (EXIT_SUCCESS);
}
//------------------------------------------------------------------------------------------------------------------------------

int getRaDec(char *aspectfile, double thetax, double thetay, double *RA, double *DEC) {
    //cout<<"\nInside getCC\n";
    fitsfile *fptr;
    int status = 0;
    Q qc;
    //cout<<"\nInertial Coordinates-----RA:"<<RA<<"  DEC:"<<DEC;
    fits_open_file(&fptr, aspectfile, READONLY, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC1", &qc.q1, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC2", &qc.q2, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC3", &qc.q3, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_read_key(fptr, TDOUBLE, "QC4", &qc.q4, NULL, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //cout<<"\nQC from average aspect file : ";
    //qc.display();
    //cout<<endl;
    double x = tan(thetax) / (sqrt(tan(thetax) * tan(thetax) + tan(thetay) * tan(thetay) + 1));
    double y = 1 / (sqrt(tan(thetax) * tan(thetax) + tan(thetay) * tan(thetay) + 1));
    double z = tan(thetay) / (sqrt(tan(thetax) * tan(thetax) + tan(thetay) * tan(thetay) + 1));

    Axis Nczti(x, y, z);
    Axis inertial;
    rotate(Nczti, qc, inertial);
    LOG(INFO) << "z:" << z;
    *RA = (sqrt(inertial.x * inertial.x + inertial.y * inertial.y) == 0) ? (0) : (acos(inertial.y / (sqrt(inertial.x * inertial.x + inertial.y * inertial.y)))); //in radians
    *DEC = asin(inertial.z);
    LOG(INFO) << "Inertial.z : " << inertial.z << "   " << asin(inertial.z);
    //in radians
    //cout<<endl<<"RA:"<<*RA<<"   DEC:"<<*DEC;
    return (EXIT_SUCCESS);
}


//------------------------------------------------------------------------------

int writeHistory(char *filename, vector<string> &vhistory) {
    int status = 0, numhdu;
    fitsfile *fptr;
    fits_open_file(&fptr, filename, READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file-writeHistory()";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_num_hdus(fptr, &numhdu, &status);
    if (status) {
        LOG(ERROR) << "Error in getting HDU number-writeHistory()";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    for (int i = 1; i <= numhdu; i++) {
        fits_movabs_hdu(fptr, i, NULL, &status);
        if (status) {
            LOG(ERROR) << "***Error moving in file " << filename << "-writeHistory()***";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        for (int j = 0; j < vhistory.size(); j++) {
            fits_write_history(fptr, (char *) vhistory[j].c_str(), &status);
            if (status) {
                LOG(ERROR) << "Error writing History-writeHistory()";
                fits_report_error(stderr, status);
                return (EXIT_FAILURE);
            }
        }
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error closing file " << filename << " - writeHistory()";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}
//------------------------------------------------------------------------------

int updateKeywords(char *filename, char *creator) {
    int status = 0, numhdu;
    fitsfile *fptr;
   
    fits_open_file(&fptr, filename, READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in updateKeywords()";
        fits_report_error(stderr, status);
        return status;
    }
    fits_get_num_hdus(fptr, &numhdu, &status);
    if (status) {
        LOG(ERROR) << "Error in updateKeywords()";
        fits_report_error(stderr, status);
        return status;
    }
    for (int i = 1; i <= numhdu; i++) {
        fits_movabs_hdu(fptr, i, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in updateKeywords()";
            fits_report_error(stderr, status);
            return status;
        }
        fits_write_date(fptr, &status);
        fits_update_key(fptr, TSTRING, "ORIGIN", (char *) ORIGIN, NULL, &status);
        fits_update_key(fptr, TSTRING, "CREATOR", creator, NULL, &status);
        fits_update_key(fptr, TSTRING, "FILENAME", filename, NULL, &status);
        fits_write_chksum(fptr, &status);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in updateKeywords()";
        fits_report_error(stderr, status);
        return status;
    }
    return (EXIT_SUCCESS);
}


int updateHdrTime(string filename, string extname, string colname){
    int status=0;
    fitsfile *fptr;
    vector <double> vecTime;
    double minTime=0.0;
    double maxTime=0.0;
    long i=0;
    double *timeArray;
    long tstarti=0, tstopi=0;
    double tstartf=0.0, tstopf=0.0;
    long nrows;
    string errorMsg="";
    int colnum=0;

    fits_open_file(&fptr, (char*) filename.c_str(), READWRITE, &status);
    if (status) {
        LOG(ERROR) << "Error in opening file " << filename;
        fits_report_error(stderr, status);
        return status;
    }

    fits_movnam_hdu(fptr, BINARY_TBL, (char*) (extname).c_str(), NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in moving to HDU " << extname << " of L2 event file " <<
                filename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    fits_get_num_rows(fptr, &nrows, &status);
    errorMsg = "Error in getting number of rows in " + extname + "of Level-2 Event File.";
    if (report_error(status, errorMsg)) {
        return EXIT_FAILURE;
    }

    if (nrows > 0) {
        fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
        errorMsg = "Error in getting column number for TIME column in " + extname + " of Level-2 Event File";
        if (report_error(status, errorMsg)) {
            return EXIT_FAILURE;
        }
        timeArray = new double[nrows];
        for (i = 0; i < nrows; i++) {
            timeArray[i] = 0.0;
        }
        fits_read_col(fptr, TDOUBLE, colnum, 1, 1, nrows, NULL, timeArray, NULL, &status);
        errorMsg = "Error in reading TIME extension  of " + extname + " of Level-2 Event File";
        if (report_error(status, errorMsg)) {
            return EXIT_FAILURE;
        }

        vecTime.clear();
        for (i = 0; i < nrows; i++) {
            vecTime.push_back(timeArray[i]);
        }

        compute_min_and_max(vecTime.begin(), vecTime.end(), minTime, maxTime);
        tstarti = (long) minTime;
        tstopi = (long) maxTime;
        tstartf = minTime - tstarti;
        tstopf = maxTime - tstopi;
        fits_update_key(fptr, TLONG, "TSTARTI", &tstarti, NULL, &status);
        fits_update_key(fptr, TLONG, "TSTOPI", &tstopi, NULL, &status);
        fits_update_key(fptr, TDOUBLE, "TSTARTF", &tstartf, NULL, &status);
        fits_update_key(fptr, TDOUBLE, "TSTOPF", &tstopf, NULL, &status);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file " << filename;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    return EXIT_SUCCESS;
}
//------------------------------------------------------------------------------

//size of q0/q1/q2/q3 is [XPIX_QUAD*YPIX_QUAD]
//size of fullImage is [XSIZE*YSIZE]

int makeFullImage(long *q0, long *q1, long *q2, long *q3, long *fullImage) {
    /**
     * Starting co-ordinates for quadrants are-
     * Q0   (0,0)     (0,0)
     * Q1   (64,0)    (XPIX_QUAD,0)
     * Q2   (64,64)   (XPIX_QUAD,YPIX_QUAD) 
     * Q3   (0,64)    (0,YPIX_QUAD)
     */
    int i,j,k,l=0; //counter variables.
    //for 1st quadrant
    //loop variables k,l for fullImage
    //loop variables i,j for quadrants
    for (i = 0; i < XSIZE * YSIZE; i++)
        fullImage[i] = 0;

    for (i = 0, k = 0; i < YPIX_QUAD; i++, k++) {
        for (j = 0, l = 0; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q0[i * XPIX_QUAD + j];
        }
    }
    //for 2nd quadrant
    for (i = 0, k = 0; i < YPIX_QUAD; i++, k++) {
        for (j = 0, l = XPIX_QUAD; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q1[i * XPIX_QUAD + j];
        }
    }
    //for 3rd quadrant
    for (i = 0, k = YPIX_QUAD; i < YPIX_QUAD; i++, k++) {
        for (j = 0, l = XPIX_QUAD; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q2[i * XPIX_QUAD + j];
        }
    }
    //for 4th quadrant
    for (i = 0, k = YPIX_QUAD; i < YPIX_QUAD; i++, k++) {
        for (j = 0, l = 0; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q3[i * XPIX_QUAD + j];
        }
    }
    return (EXIT_SUCCESS);
}



//------------------------------------------------------------------------------

int makeFullDPI(float *q0, float *q1, float *q2, float *q3, float *fullImage) {
    /**
     * Starting co-ordinates for quadrants are-
     * Q0   (0,0)     (0,0)
     * Q1   (64,0)    (XPIX_QUAD,0)
     * Q2   (64,64)   (XPIX_QUAD,YPIX_QUAD) 
     * Q3   (0,64)    (0,YPIX_QUAD)
     */
    //for 1st quadrant
    //loop variables k,l for fullImage
    //loop variables i,j for qs
    for (int i = 0; i < XSIZE * YSIZE; i++)
        fullImage[i] = 0;

    for (int i = 0, k = 0; i < YPIX_QUAD; i++, k++) {
        for (int j = 0, l = 0; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q0[i * XPIX_QUAD + j];
        }
    }
    //for 2nd quadrant
    for (int i = 0, k = 0; i < YPIX_QUAD; i++, k++) {
        for (int j = 0, l = XPIX_QUAD; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q1[i * XPIX_QUAD + j];
        }
    }
    //for 3rd quadrant
    for (int i = 0, k = YPIX_QUAD; i < YPIX_QUAD; i++, k++) {
        for (int j = 0, l = XPIX_QUAD; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q2[i * XPIX_QUAD + j];
        }
    }
    //for 4th quadrant
    for (int i = 0, k = YPIX_QUAD; i < YPIX_QUAD; i++, k++) {
        for (int j = 0, l = 0; j < XPIX_QUAD; j++, l++) {
            fullImage[k * XSIZE + l] = q3[i * XPIX_QUAD + j];
        }
    }
    return (EXIT_SUCCESS);
}

int mergeImages(long *outimage, int size, int n, ...) {
    va_list images;
    va_start(images, n);
    long *im;
    for (int i = 0; i < n; i++) {
        im = va_arg(images, long *);
        for (int j = 0; j < size; j++)
            outimage[j] = outimage[j] + im[j];
    }

    va_end(images);

    return 0;
}


#define MINQUADID 0
#define MAXQUADID 3
#define MINDETID 0
#define MAXDETID 15
#define MINPIXID 0
#define MAXPIXID 255

int getxy(int quadID, int detID, int pixID, unsigned char *x, unsigned char *y) { //x=column; y=row;
    int quad_x, quad_y, det_x, det_y, pix_x, pix_y, x_inquad, y_inquad;

    if (quadID > MAXQUADID || quadID < MINQUADID || detID > MAXDETID || detID < MINDETID || pixID > MAXPIXID || pixID < MINPIXID) {
        LOG(ERROR) << "Invalid quadid/detid/pixid";
        return (EXIT_FAILURE);
    }

    /*making quadrant id in the form
     *--------                    ---------
     *|0 | 1 |                    | 0 | 1 |
     *--------  where as given is ---------
     *|2 | 3 |                    | 3 | 2 |
     * -------                    ---------
     */

    if (quadID == 2) quadID = 3;
    else if (quadID == 3) quadID = 2;

    //cout<<"\n"<<quadID<<"  "<<detID<<"  "<<pixID;
    //changing quadrant_id to 2D x,y wrt to complete frame
    quad_y = (quadID) / NO_QUAD_X;
    quad_x = (quadID) % NO_QUAD_X;

    //changing detector_id to 2D x,y within the quadrant
    det_y = (detID) / NO_DET_X_PER_QUAD;
    det_x = (detID) % NO_DET_X_PER_QUAD;

    //changing pixel id to 2D x,y within the detector
    pix_y = (pixID) / NO_PIX_X_PER_DET;
    pix_x = (pixID) % NO_PIX_X_PER_DET;

    *y = pix_y + det_y * NO_PIX_Y_PER_DET + quad_y * NO_DET_Y_PER_QUAD*NO_PIX_Y_PER_DET;
    *x = pix_x + det_x * NO_PIX_X_PER_DET + quad_x * NO_DET_X_PER_QUAD*NO_PIX_X_PER_DET;

    return (EXIT_SUCCESS);
}

int getRaDec(char *file, int hdunum, double *RA, double *DEC) {
    int status = 0;
    fitsfile *fptr;
    fits_open_file(&fptr, file, READONLY, &status);
    if (status) return (EXIT_FAILURE);
    fits_movabs_hdu(fptr, hdunum, NULL, &status);
    if (status) return (EXIT_FAILURE);
    fits_read_key(fptr, TDOUBLE, "RA", &RA, NULL, &status);
    if (status) return (EXIT_FAILURE);
    fits_read_key(fptr, TDOUBLE, "DEC", &DEC, NULL, &status);
    if (status) return (EXIT_FAILURE);
    fits_close_file(fptr, &status);
    if (status) return (EXIT_FAILURE);
    return (EXIT_SUCCESS);
}

int writeImg(char *file, long *img, int m, int n) //just to check;
{
    int status = 0;
    fitsfile *fptr;
    fits_create_file(&fptr, file, &status);
    if (status) {
        LOG(ERROR) << "Error in creating file:" << file;
        return (EXIT_FAILURE);
    }
    int bitpix = LONG_IMG;
    int naxis = 2;
    long naxes[2];
    naxes[0] = m;
    naxes[1] = n;
    fits_create_img(fptr, bitpix, naxis, naxes, &status);
    if (status) {
        LOG(ERROR) << "Error in creating image:" << file;
        return (EXIT_FAILURE);
    }
    long fpixel[2];
    fpixel[0] = fpixel[1] = 1;
    fits_write_pix(fptr, TLONG, fpixel, m*n, img, &status);
    if (status) {
        LOG(ERROR) << "Error in writing pixels:" << file;
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file:" << file;
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int writeImg(char *file, float *img, int m, int n) //just to check;
{
    int status = 0;
    fitsfile *fptr;
    fits_create_file(&fptr, file, &status);
    if (status) {
        LOG(ERROR) << "Error in creating file:" << file;
        return (EXIT_FAILURE);
    }
    int bitpix = FLOAT_IMG;
    int naxis = 2;
    long naxes[2];
    naxes[0] = m;
    naxes[1] = n;
    fits_create_img(fptr, bitpix, naxis, naxes, &status);
    if (status) {
        LOG(ERROR) << "Error in creating image:" << file;
        return (EXIT_FAILURE);
    }
    long fpixel[2];
    fpixel[0] = fpixel[1] = 1;
    fits_write_pix(fptr, TFLOAT, fpixel, m*n, img, &status);
    if (status) {
        LOG(ERROR) << "Error in writing pixels:" << file;
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file:" << file;
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int writeImgD(char *file, double *img, int m, int n) //just to check;
{
    int status = 0;
    fitsfile *fptr;
    fits_create_file(&fptr, file, &status);
    if (status) {
        LOG(ERROR) << "Error in creating file:" << file;
        return (EXIT_FAILURE);
    }
    int bitpix = FLOAT_IMG;
    int naxis = 2;
    long naxes[2];
    naxes[0] = m;
    naxes[1] = n;
    fits_create_img(fptr, bitpix, naxis, naxes, &status);
    if (status) {
        LOG(ERROR) << "Error in creating image:" << file;
        return (EXIT_FAILURE);
    }
    long fpixel[2];
    fpixel[0] = fpixel[1] = 1;
    fits_write_pix(fptr, TDOUBLE, fpixel, m*n, img, &status);
    if (status) {
        LOG(ERROR) << "Error in writing pixels:" << file;
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "Error in closing file:" << file;
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int getqxqy(int x, int y, int *qx, int *qy, int *quad) {
    //cout<<"\n"<<x<<" , "<<y;
    if (x < XPIX_QUAD && y < YPIX_QUAD) {
        *qx = x;
        *qy = y;
        *quad = 0;
    } else if (x >= XPIX_QUAD && y < YPIX_QUAD) {
        *qx = x - XPIX_QUAD;
        *qy = y;
        *quad = 1;
    } else if (x >= XPIX_QUAD && y >= YPIX_QUAD) {
        *qx = x - XPIX_QUAD;
        *qy = y - YPIX_QUAD;
        *quad = 2;
    } else if (x < XPIX_QUAD && y >= YPIX_QUAD) {
        *qx = x;
        *qy = y - YPIX_QUAD;
        *quad = 3;
    } else {
        LOG(ERROR) << "***Invalid (x,y)***";
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int writeArray(char *file, float **array, int m, int n) {
    fitsfile *fptr;
    int status = 0;
    float data[m * n];
    int index = 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            data[index++] = array[i][j];
    fits_create_file(&fptr, file, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    int bitpix = FLOAT_IMG;
    int naxis = 2;
    long naxes[2];
    naxes[0] = n;
    naxes[1] = m;
    fits_create_img(fptr, bitpix, naxis, naxes, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    long fpixel[2];
    fpixel[0] = 1;
    fpixel[1] = 1;
    fits_write_pix(fptr, TFLOAT, fpixel, m*n, data, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int getxy4mqxqy(int qid, int qx, int qy, int *x, int *y) {

    if (qx > XPIX_QUAD || qy > YPIX_QUAD) {
        LOG(ERROR) << "***Invalid (qx,qy)***";
        return (EXIT_FAILURE);
    }
    switch (qid) {
        case 0: *x = qx;
            *y = qy;
            break;
        case 1: *x = qx + XPIX_QUAD;
            *y = qy;
            break;
        case 2: *x = qx + XPIX_QUAD;
            *y = qy + YPIX_QUAD;
            break;
        case 3: *x = qx;
            *y = qy + YPIX_QUAD;
            break;
        default:
            LOG(ERROR) << "Invalid Quadrant ID";
            return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

void checkPFILESenv() {

/*	
    if (getenv(ENV_PFILES) == NULL) {
        LOG(ERROR) << "Environment variable PFILES not set";
        exit(EXIT_FAILURE);
    }
*/
}

void checkParFile(char *modulename) {
 /*
    char *pardir, parfilename[NAMESIZE];
    pardir = getenv(ENV_PFILES);
    if (pardir == NULL) {
        LOG(ERROR) << "***Environment variable PFILES not set***";
        exit(EXIT_FAILURE);
    }
    strcpy(parfilename, pardir);
    strcat(parfilename, "/");
    strcat(parfilename, modulename);
    strcat(parfilename, ".par");
    //cout<<"\n"<<parfilename;
    if (!FileExists(parfilename)) {
        LOG(ERROR) << "***Parameter file not found for the module '" << modulename << "'***";
        exit(EXIT_FAILURE);
    }
*/
}

int readQE(string qefile, float *qearray) {
    fitsfile *fptr;
    int status = 0;
    fits_open_file(&fptr, (char*) qefile.c_str(), READONLY, &status);
    if (status) {
        LOG(ERROR) << "***Error opening file " << qefile;
        return (EXIT_FAILURE);
    }
    fits_movnam_hdu(fptr, IMAGE_HDU, "QE", 0, &status);
    if (status) {
        LOG(ERROR) << "***Error in moving to HDU 'QE' in file " << qefile;
        return (EXIT_FAILURE);
    }
    long fpixel[2];
    fpixel[0] = fpixel[1] = 1;
    fits_read_pix(fptr, TFLOAT, fpixel, XSIZE*YSIZE, NULL, qearray, NULL, &status);
    if (status) {
        LOG(ERROR) << "***Error reading pixel values in file " << qefile;
        return (EXIT_FAILURE);
    }
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "***Error closing file " << qefile;
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}



void quadToHDU(int quad_no, char* HDUname) {
    switch (quad_no) {
        case 0:
            strcpy(HDUname, "Q0");
            break;
        case 1:
            strcpy(HDUname, "Q1");
            break;
        case 2:
            strcpy(HDUname, "Q2");
            break;
        case 3:
            strcpy(HDUname, "Q3");
            break;
    }
    return;
}

void HDUtoQuad(string HDUname, unsigned char &quadNo) {
    if(HDUname=="Q0"){quadNo=0;}
    else if(HDUname=="Q1"){quadNo=1;}
    else if(HDUname=="Q2"){quadNo=2;}
    else if(HDUname=="Q3"){quadNo=3;}
    else {quadNo =-1;}
    
    return;
}

void DPIquadToHDU(int quad_no, char* HDUname) {
    switch (quad_no) {
        case 0:
            strcpy(HDUname, "Q0_DPI");
            break;
        case 1:
            strcpy(HDUname, "Q1_DPI");
            break;
        case 2:
            strcpy(HDUname, "Q2_DPI");
            break;
        case 3:
            strcpy(HDUname, "Q3_DPI");
            break;
    }
    return;
}

/******************************************************************************************************
 Function to break DPI array into 4 arrays based on quadrant
 ******************************************************************************************************/
int breakDPIarray(float *dpi, float **quad_dpi, int full_XSIZE, int full_YSIZE, int quad_XSIZE, int quad_YSIZE) {
    int index = 0, i = 0, j = 0;

    //creating quad0_dpi
    for (i = 0; i < quad_YSIZE; i++) {
        for (j = 0; j < quad_XSIZE; j++) {
            quad_dpi[0][index] = dpi[full_XSIZE * i + j];
            index++;
        }
    }
    index = 0;
    //creating quad1_dpi
    for (i = 0; i < quad_YSIZE; i++) {
        for (j = quad_XSIZE; j < full_XSIZE; j++) {
            quad_dpi[1][index] = dpi[full_XSIZE * i + j];
            index++;
        }
    }

    //creating quad2_dpi
    index = 0;
    for (i = quad_YSIZE; i < full_YSIZE; i++) {
        for (j = 0; j < quad_XSIZE; j++) {
            quad_dpi[2][index] = dpi[full_XSIZE * i + j];
            index++;
        }
    }

    //creating quad3_dpi
    index = 0;
    for (i = quad_YSIZE; i < full_YSIZE; i++) {
        for (j = quad_XSIZE; j < full_XSIZE; j++) {
            quad_dpi[3][index] = dpi[full_XSIZE * i + j];
            // LOG(INFO)<< quad_dpi[3][index];
            // LOG(INFO)<< dpi[full_XSIZE * i + j];
            index++;
        }
    }

    return 0;

}

int report_error(int status, string error_string){
    if (status){
        LOG(ERROR) << error_string;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    return status;
}

string quadExtName(int quadNo) {
    string tempExtName = "";
    switch (quadNo) {
        case 0:
            tempExtName = "Q0";
            break;
        case 1:
            tempExtName = "Q1";
            break;
        case 2:
            tempExtName = "Q2";
            break;
        case 3:
            tempExtName = "Q3";
            break;
    }
    
    return tempExtName;
}

//int find_interpolated_value(float x, float x1, float x2, float y1, float y2, float &y){
//    float slope;
//    slope = (y2-y1)/(x2-x1);
//    y = slope*(x-x1) + y1;
//    return (EXIT_SUCCESS);
//}

bool is_pix_present_in_quads(unsigned char locx, unsigned char locy, vector<int> quadsToProcess, int &quadID, int &pixx, int &pixy){
    int iquad=0;
    int i=0;
    vector<bool> pixPresent; //tells about pixel is present in which quad
    bool pixelPresent=false; //pixel is present in quads
    pixPresent.resize(4, false);
    quadID=-1;
    
    for(iquad=0; iquad<NUMQUAD; iquad++){
        if(iquad==0){
            if(locx>=0 && locx<=63 && locy>=64 && locy<=127){
                pixx = locx-0;
                pixy = locy-64;
                pixPresent[0]=true;
                quadID=0;
                break;
            }
        }
        else if(iquad==1){
            if(locx>=64 && locx<=127 && locy>=64 && locy<=127) {
                pixx = locx - 64;
                pixy = locy - 64;
                pixPresent[1]=true;
                quadID=1;
                break;
            }
        }
        else if(iquad==2){
            if(locx>=64 && locx<=127 && locy>=0 && locy<=63) {
                pixx = locx - 64;
                pixy = locy - 0;
                pixPresent[2]=true;
                quadID=2;
                break;
            }
        }
        else if(iquad==3){
            if(locx>=0 && locx<=63 && locy>=0 && locy<=63) {
                pixx = locx - 0;
                pixy = locy - 0;
                pixPresent[3]=true;
                quadID=3;
                break;
            }
        }
    }
    
    for(i=0; i<quadsToProcess.size(); i++){
        iquad=quadsToProcess[i];
        if(pixPresent[iquad]==true){
            pixelPresent=true;
            break;
        }
    }
    
    return pixelPresent;
}
int create_empty_fitsfile(string outFilename, string outTemplate){
    int status=0; //status variable
    fitsfile *fptr;
    string templateFileName = ""; //full path to template file.
    
    templateFileName = template_full_path_generator(outTemplate);
    if(templateFileName==""){
        LOG(ERROR)<< "Not able to generate Event file template path.";
        LOG(ERROR)<< "Probably Environment Variables are not declared properly.";
        return(EXIT_FAILURE);
    }
    LOG(INFO) << "Template file used : " << templateFileName;
    
    fits_create_template(&fptr, (char*) outFilename.c_str(), (char*) templateFileName.c_str(), &status);
    if(status){
        LOG(ERROR) << "Error in creating file : " << outFilename;
        fits_report_error(stderr, status);
        return(EXIT_FAILURE);
    }
    
    fits_close_file(fptr, &status);
    if (status) {
        LOG(ERROR) << "***Error in closing file : " << outFilename << "***";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    
    DLOG(INFO) << "Empty " << outFilename << " created successfully.";
    return (EXIT_SUCCESS);
}

string template_full_path_generator(string filename){
    string baseTemplateDirPath="";
    string TemplateFileFullPath="";
    try{
        string baseTemplateDirPath = getenv("CZTI_templates");

        string TemplateFileFullPath = baseTemplateDirPath + "/" + filename;
        return TemplateFileFullPath;
    }
    catch (std::exception& e) {
        LOG(ERROR) << "Error: " << e.what();
        return TemplateFileFullPath="";
    }

}

int read_fits_string_column(fitsfile *fptr, string colname,
        int datatype, long firstrow, long firstelem, long nelements,
        vector <string> &vec_data) {
    int status = 0;
    long i = 0;
    long nrows = 0;
    int colnum = 0;
    int typecode = 0;
    long repeat = 0;
    long width = 0;
    char **temp_data;
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        LOG(ERROR) << "Error in getting column number for " << colname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }
    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) {
        LOG(ERROR) << "Error in getting type code for column " << colname;
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            LOG(ERROR) << "Error in getting number of rows for " << colname << "COLUMN.";
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
        nelements = nrows;//* repeat*width;
    } 
    temp_data = new char* [nelements];
    for(i=0; i<nelements; i++){
        temp_data[i]=new char[width];
    }
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            temp_data, NULL, &status);
    if (status) {
        LOG(ERROR) << "Error in reading " << colname << " COLUMN.";
        fits_report_error(stderr, status);
        return (EXIT_FAILURE);
    }

    //Assigning values to output vector
        for(i=0; i<nelements; i++){
            vec_data.push_back((string) (temp_data[i]));
        }

    delete[] temp_data;
    return EXIT_SUCCESS;
}

int read_fits_string_columnN(fitsfile *fptr, string colname,
        int datatype, long firstrow, long firstelem, long nelements,
        vector <string> &vec_data) {
    int status = 0;
    long i = 0;
    long nrows = 0;
    int colnum = 0;
    int typecode = 0;
    long repeat = 0;
    long width = 0;
    char **temp_data;
    ErrorHandler errHandler;
    
    fits_get_colnum(fptr, CASEINSEN, (char*) colname.c_str(), &colnum, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in getting column number for " + colname;
        throw errHandler;
    }

    fits_get_coltype(fptr, colnum, &typecode, &repeat, &width, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in getting type code for column " + colname;
        throw errHandler;
    }

    if (nelements == -1) {
        fits_get_num_rows(fptr, &nrows, &status);
        if (status) {
            fits_read_errmsg(errHandler.fitsErrMsg);
            fits_get_errstatus(status, errHandler.fitsErrTxt);
            errHandler.fitsflag = true;
            errHandler.fitsErrorStatus = status;
            errHandler.severity = errERROR;
            errHandler.errorMsg = "Error in getting number of rows for " + colname + " COLUMN.";
            throw errHandler;
        }
        nelements = nrows;//* repeat*width;
    } 
    temp_data = new char* [nelements];
    for(i=0; i<nelements; i++){
        temp_data[i]=new char[width];
    }
    fits_read_col(fptr, datatype, colnum, firstrow, firstelem, nelements, NULL,
            temp_data, NULL, &status);
    if (status) {
        fits_read_errmsg(errHandler.fitsErrMsg);
        fits_get_errstatus(status, errHandler.fitsErrTxt);
        errHandler.fitsflag = true;
        errHandler.fitsErrorStatus = status;
        errHandler.severity = errERROR;
        errHandler.errorMsg = "Error in reading " + colname + " COLUMN.";
        throw errHandler;
    }

    //Assigning values to output vector
        for(i=0; i<nelements; i++){
            vec_data.push_back((string) (temp_data[i]));
        }

    delete[] temp_data;
    return EXIT_SUCCESS;
}


//ERROR HANDLING
void logError(ErrorHandler errHandler){
    if(errHandler.severity==errINFO){
        LOG(INFO) << errHandler.errorMsg;
        if(errHandler.fitsflag==true){
            DLOG(INFO) << "**** FITS STATUS CODE:" << errHandler.fitsErrorStatus << " MESSAGE: "<< errHandler.fitsErrMsg << " | " << errHandler.fitsErrTxt << "****"; 
        }    
    }
    else if(errHandler.severity==errWARNING){
        LOG(WARNING) << errHandler.errorMsg;
        if(errHandler.fitsflag==true){
            LOG(INFO) << "**** FITS STATUS CODE: " << errHandler.fitsErrorStatus << " MESSAGE: "<< errHandler.fitsErrMsg << " | " << errHandler.fitsErrTxt << "****"; 
        }    
    }
    if(errHandler.severity==errERROR){
        LOG(ERROR) << errHandler.errorMsg;
        if(errHandler.fitsflag==true){
            LOG(ERROR) << "**** FITS STATUS CODE:" << errHandler.fitsErrorStatus << " MESSAGE: "<< errHandler.fitsErrMsg << " | " << errHandler.fitsErrTxt << "****"; 
        }    
    }
}


//STATISTICS

bool max_pair_value(pair<float, long> pair1, pair<float, long> pair2) {
    return pair1.second < pair2.second;
}

//Getting quadrant vector

int get_quadsToProcessVec(string quadsToProcess, vector<int> *quadVec){
    int tempNoOccurance=0;
    int no_quads=0;
    int *quadArray;
    stringFinder((char*)quadsToProcess.c_str(), ",", 0, &tempNoOccurance);
    no_quads = tempNoOccurance+1;
    quadArray = new int[no_quads];

    if (csvStrToInt((char*)quadsToProcess.c_str(), ",", quadArray, &no_quads)) {
        LOG(ERROR) << "***Error in converting quadrant string array into integer array***";
        return EXIT_FAILURE;
    } 
    
    (*quadVec).insert((*quadVec).begin(),quadArray, quadArray+no_quads);
    return EXIT_SUCCESS;
}

//File/Direcory functions


bool FileExists(char *filename) {
    //cout<<"\n----------FileExists----------\n";
    struct stat filestat;
    int filest;
    filestat.st_mode = -1;
    bool exist = false;
    filest = stat(filename, &filestat);
    if (S_ISREG(filestat.st_mode)) {
        //cout<<endl<<filename<<":File exists:"<<endl;
        exist = true;
    } else exist = false;

    //cout<<"\n----------FileExists End----------\n";
    return exist;
}

void files_exist(vector<string> filenames, vector <string> *nonExistentFiles){
    int nfiles=filenames.size();
    int ifile=0;
    
    for(ifile=0; ifile<nfiles; ifile++){
        if(FileExists((char*)(filenames[ifile]).c_str())){
           continue; 
        } else {
            (*nonExistentFiles).push_back(filenames[ifile]);
        }
    }
}
//---------------------------------------------------------------------

bool DirExists(char *dirname) {
    struct stat dirstat;
    int t;
    dirstat.st_mode = -1;
    t = stat(dirname, &dirstat);
    if (S_ISDIR(dirstat.st_mode)) {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------

//Delete file if it exists and clobber is yes

int deleteFile(vector<string> filenames, int clobber) {
    int i = 0;
    int status = 0;
    for (i = 0; i < filenames.size(); i++) {
        if (FileExists((char*) (filenames[i]).c_str())) {
            if (clobber == YES) {
                unlink((char*) (filenames[i]).c_str());
            } else {
                LOG(ERROR) << "" << filenames[i] << " already exists";
                LOG(ERROR) << "Use clobber=yes for overwriting the file";
                return (EXIT_FAILURE);
            }
        }
    }

    return status;
}

vector <string> find_files_in_filelist(vector<string> filelist, string extname, vector <string> *extraChecks) {
    int ifile = 0;
    int icheck =0;
    size_t pos;
    string substring = "";
    string filename = "";
    bool badFileFlag = false;
    vector <string> filenameVec;
    for (ifile = 0; ifile < filelist.size(); ifile++) {
        substring = filelist[ifile].substr(filelist[ifile].find_last_of(".") + 1);
        if (substring == extname) {
            filename = filelist[ifile];
            if (extraChecks != NULL) {
                badFileFlag = false;
                substring = filelist[ifile].substr(filelist[ifile].find_last_of("/")+1);
                for(icheck=0; icheck<(*extraChecks).size(); icheck++){
                    pos=substring.rfind((*extraChecks)[icheck]);
                    if(pos==std::string::npos){
                        badFileFlag=true;
                        break;
                    }
                }
            }
            if (badFileFlag == false) {
                filenameVec.push_back(filename);
            }
        }
    }
    return filenameVec;
}

int getFiles(string dirName, vector<string> *filelist) {
    bool darkFlag = FALSE;
    

    struct dirent **namelist;
    int num = scandir(dirName.c_str(), &namelist, defaultfilter, alphasort);
    
    if (num <= 0) {
        LOG(WARNING) << "No file present in directory " << dirName << ".";
    }

    bool fileflag = true;

    while (num--) {
        string tempname = dirName + (string) "/" + (string) namelist[num]->d_name;
        if (namelist[num]->d_type == DT_REG) {
            (*filelist).push_back(tempname);
        } else if (namelist[num]->d_type == DT_DIR) {


            fileflag = false;
            getFiles(tempname, filelist);
        }
    }

    if (fileflag == true) {
        return (EXIT_SUCCESS);
    }

    return (EXIT_SUCCESS);
}

int getDirs(string dirName, vector<string> *subDirs) {

    struct dirent **namelist;
    int num = scandir(dirName.c_str(), &namelist, defaultfilter, alphasort);
    if (num <= 0) {
        //LOG(ERROR)<<endl<<"No files found in directory "<<dirName<<endl;
        return (EXIT_FAILURE);
    }

    bool fileflag = true;
    while (num--) {
        string tempname = dirName + (string) "/" + (string) namelist[num]->d_name;
        if (namelist[num]->d_type == DT_REG) {

        } else if (namelist[num]->d_type == DT_DIR) {
            fileflag = false;
            (*subDirs).push_back(tempname);
            getDirs(tempname, subDirs);
        }
    }
    if (fileflag == true) {
        return (EXIT_SUCCESS);
    }

    return (EXIT_SUCCESS);
}
int searchDir(string dirname, string pattern, vector<string> *fileList, vector<string> *dirList){
    int status=0;
    DIR *currentDir;
    vector <string> content;
    ErrorHandler errHandler;
    string tempcontent;
    int size_temp=0;
    int p = 0;
    int count = 0;
    struct dirent *dirp;
    if ((currentDir = opendir((const char *) dirname.c_str())) == NULL) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = ERROR_OPENING_DIRECTORY;
        errHandler.errorMsg = "Error in opening" +(string)dirname +" directory.";
        throw errHandler;
    }
    while (dirp = readdir(currentDir)) {
        if (strcmp(dirp->d_name, ".") == 0 || strcmp(dirp->d_name, "..") == 0) {
            continue;
        }
        tempcontent = dirp->d_name;
        size_temp = tempcontent.size();
        int npos = tempcontent.find(pattern, 0);
        if (npos >= 0 && npos < size_temp) {
            if (dirp->d_type == DT_DIR) {
                if (dirList != NULL) {
                    (*dirList).push_back(tempcontent);
                }
            } else if (dirp->d_type == DT_REG) {
                if (fileList != NULL) {
                    (*fileList).push_back(tempcontent);
                }
            }
        }

    }
    closedir(currentDir);

    return EXIT_SUCCESS;
}
int getHeaderKey(const char* infile, int hdunum,const char*key,char* value)
{
	fitsfile *fptr;
	char chr_value[100];
    	int status = 0;
    	fits_open_file(&fptr, infile, READONLY, &status);
	if (status) 
	{
        	fits_report_error(stderr, status);
        	LOG(ERROR) << "***Error in opening file : " << infile << "***";
       		return (EXIT_FAILURE);
    	}
        fits_movabs_hdu(fptr, hdunum, NULL, &status);
        if (status) {
            LOG(ERROR) << "Error in moving to HDU " << hdunum << " in file " << infile;
            fits_report_error(stderr, status);
            return (EXIT_FAILURE);
        }
    	fits_read_key(fptr, TSTRING, key, value, NULL, &status);
    	if (status) 
	{

        	LOG(ERROR) << "***Error in reading key : " << key << "***";
        	fits_report_error(stderr, status);
        	return (EXIT_FAILURE);
    	}

	fits_close_file(fptr, &status);
    	if (status) 
	{
        	LOG(ERROR) << "***Error in closing file : " << infile << "***";
        	fits_report_error(stderr, status);
        	return (EXIT_FAILURE);
    	}
}



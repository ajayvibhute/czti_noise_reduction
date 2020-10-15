

#include "cztevtclean.h"

using namespace std;

cztevtclean::cztevtclean(){
     strcpy(modulename,"cztevtclean_v");
     strcat(modulename,VERSION);
    
}

int cztevtclean::read(int argc,char **argv){
    int status=0;
    
    if(PIL_OK!=(status=PILInit(argc,argv))){
        LOG(ERROR)<<"***Error Initializing PIL***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("infile",infile))){
        LOG(ERROR)<<"***Error Reading Event data file:"<<infile<<"***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetFname("outfile",outfile))){
        LOG(ERROR)<<"***Error Reading outfile:"<<outfile<<"***";
        return status;
    }

    if(PIL_OK!=(status=PILGetString("alphaval",alphaval))){
            LOG(ERROR)<<"***Error Reading alphaval parameter***";
            return status;
        }
    if(PIL_OK!=(status=PILGetString("vetorange",vetorange))){
        LOG(ERROR)<<"***Error Reading vetorange parameter***";
        return status;
    }

    if(PIL_OK!=(status=PILGetBool("isdoubleEvent",&isdoubleEvent))){
        LOG(ERROR)<<"***Error Reading isdoubleevent parameter"<<isdoubleEvent<<"***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetBool("history",&history))){
        LOG(ERROR)<<"***Error Reading history parameter"<<history<<"***";
        return status;
    }
    
    if(PIL_OK!=(status=PILGetBool("clobber",&clobber))){
        LOG(ERROR)<<"***Error Reading clobber:"<<clobber<<"***";
        return status;
    }

    PILClose(status);
    return (EXIT_SUCCESS);
}

 int cztevtclean::read(char *infile,char *outfile,char *alphaval,
         char *vetorange,int isdoubleEvent,int clobber,int history){
     strcpy(this->infile,infile);
     strcpy(this->outfile,outfile);
     strcpy(this->alphaval,alphaval);
     strcpy(this->vetorange,vetorange);
     this->isdoubleEvent=isdoubleEvent;
     this->clobber=clobber;
     this->history=history;
     return (EXIT_SUCCESS);
 }
 
 void cztevtclean::display(){
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"                           CZTEVTCLEAN PARAMETERS                             ";
    LOG(INFO)<<"----------------------------------------------------------------------------";
    LOG(INFO)<<"Taskname                   : "<<modulename;
    LOG(INFO)<<"Input file                 : "<<infile;         
    LOG(INFO)<<"Output file                : "<<outfile;
    LOG(INFO)<<"Alpha Value                : "<<alphaval;
    LOG(INFO)<<"Veto Range                 : "<<vetorange;
    LOG(INFO)<<"IsdoubleEvent              : "<<isdoubleEvent;
    if(isdoubleEvent==YES)
	LOG(INFO)<<"Is this double event file  : YES";
    else
	LOG(INFO)<<"Is this double event file  : NO";
    if(clobber==YES)
	LOG(INFO)<<"Clobber                    : YES";
    else
	LOG(INFO)<<"Clobber                    : NO";
    if(history==YES)
	LOG(INFO)<<"History                    : YES";
    else
	LOG(INFO)<<"History                    : NO";	
	LOG(INFO)<<"----------------------------------------------------------------------------";
 }
 
 int cztevtclean::cztevtcleanProcess(){
     int status=0;
     char errstr[1024];
     if(strcmp(outfile,infile)==0){  }
     else if(strcmp(outfile,"")==0 || strcmp(outfile,"-")==0){
         strcpy(outfile,infile);
     }
     else {                           //infile and outfile are different
        if(FileExists(outfile)){
            if(clobber==YES){
                    if(unlink(outfile)!=0){
                       // cerr<<"Error in deleting "<<outfile; 
                        //return (EXIT_FAILURE);
                    }
            }
        else{
            LOG(ERROR)<<""<<outfile<<" already exists";
            LOG(ERROR)<<"Use clobber=yes for overwriting the file";
            return (EXIT_FAILURE);
           }
        }
     
         //copying input to output file
         fitsfile *fin,*fout;
         //fits_open_file(&fin,infile,READWRITE,&status);  //commented by mayuri
         fits_open_file(&fin,infile,READONLY,&status);
         
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<infile<<"***";
             return (EXIT_FAILURE);
         }
         fits_create_file(&fout,outfile,&status);
          if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in creating file "<<outfile<<"***";
             return (EXIT_FAILURE);
         }
         LOG(INFO) << "Output file created" << endl; 
         fits_copy_file(fin,fout,1,1,1,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in copying input to output file***";
             return (EXIT_FAILURE);
         }
         LOG(INFO)<<"File copied";
         fits_close_file(fout,&status); if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }
         fits_close_file(fin,&status);  if(status) { fits_report_error(stderr,status); return (EXIT_FAILURE); }
     }
     LOG(INFO)<<"Input file closed";
     //checking for valid input for alphaval
    if(!(strcmp(alphaval,"0")==0 || strcmp(alphaval,"1")==0 || 
            strcmp(alphaval,"0,1")==0 || strcmp(alphaval,"1,0")==0)){
        LOG(ERROR)<<"***Invalid input for alphaval***";
        LOG(ERROR)<<"***Use '0'/'1'/'0,1'/'1,0'***";
        return (EXIT_FAILURE);
    }
    //checking for valid input for vetorange
     
     
    fitsfile *fptr;
    fits_open_file(&fptr,outfile,READWRITE,&status);
    if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in opening file "<<outfile<<"***";
             return (EXIT_FAILURE);
         }
       
    int hdutype,numhdu=0;
  
    //for all four quadrants  - quadrant data is from HDU 2 to 5
    long i=0,j=0; 
    long nrows;
    int optionval,alphacolnum,vetocolnum;
    unsigned char *alpha;
    unsigned int *veto;
    //loop for each hdu for a quadrant  
    for(i=2;i<=5;i++){              
        LOG(INFO) << "Move into HDU" << i << endl;
        fits_movabs_hdu(fptr,i,&hdutype,&status);
        sprintf(errstr,"***Error moving to HDU %d***",i);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<errstr;
             return (EXIT_FAILURE);
         }
        //LOG(ERROR)<<"Moved to hdu number:"<<hdunum;
        if(hdutype!=BINARY_TBL){
            LOG(ERROR)<<"*** HDU "<<i<<" is not binary table that contains quadrant data***";
            return (EXIT_FAILURE);
        }
        
        fits_get_num_rows(fptr,&nrows,&status);
        sprintf(errstr,"***Error getting number of rows in quadrant %d - cztevtclean_process()***",i);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<errstr;
             return (EXIT_FAILURE);
         }
        LOG(INFO) << "Number of rows : " << nrows << endl;
        if(nrows==0) continue;
        alpha=new unsigned char[nrows];  
        veto=new unsigned int[nrows];
        if(alpha==NULL || veto==NULL){
            LOG(ERROR)<<"***Out of memory error***";
            return (EXIT_FAILURE);
        }
        
        //getting column number for alpha column
        fits_get_colnum(fptr,CASEINSEN, (char *)"alpha",&alphacolnum,&status);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Alpha column not found in file***";
             return (EXIT_FAILURE);
         }
        
        //column number for veto column
        fits_get_colnum(fptr,CASEINSEN,(char *)"veto",&vetocolnum,&status);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Veto column not found in file***";
             return (EXIT_FAILURE);
         }
                
        fits_read_col(fptr,TBYTE,alphacolnum,1,1,nrows,NULL,alpha,NULL,&status);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in reading alpha column***";
             return (EXIT_FAILURE);
         }
       // cout<<"alpha column number :"<<alphacolnum;
       //for(int i=0;i<nrows;i++) {   cout<<"value of alpha : " << (int)alpha[i]; } 
        
        fits_read_col(fptr,TUINT,vetocolnum,1,1,nrows,NULL,veto,NULL,&status);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in reading veto column***";
             return (EXIT_FAILURE);
         }
        vector<long> rows;   //stores row numbers which are to be deleted
       
        
        if(isdoubleEvent==0)
        {
        for(j=1;j<=nrows;j++){
            //cout<<(bool)alpha[j-1];
           // cout<<checkVeto(veto[j-1],vetorange)<<"\t"<<checkAlpha((bool)alpha[j-1],alphaval)<<"\t"<<j;
            if(checkVeto(veto[j-1],vetorange)==1 || checkAlpha(alpha[j-1],alphaval)==1){
                rows.push_back(j);
                //cout<<checkVeto(veto[j-1],vetorange)<<"\t"<<checkAlpha((bool)alpha[j-1],alphaval)<<"\t"<<j;
            }

            if(checkAlpha(alpha[j-1],alphaval)==-1){
               LOG(ERROR)<<"***Error in filtering for veto and alpha***";
               return (EXIT_FAILURE);
            }
          }
        }
        else //If it is a double event
        {
            for(j=1;j<nrows;j+=2){
                if(checkVeto(veto[j-1],vetorange)==1 || checkAlpha(alpha[j-1],alphaval)==1||checkVeto(veto[j],vetorange)==1||checkVeto(veto[j],vetorange)==1 ){
                rows.push_back(j);
                rows.push_back(j+1);
                }

            }
        }

	LOG(INFO)<<" Number of rows to be deleted "<<rows.size();
        fits_delete_rowlist(fptr,rows.data(),rows.size(),&status);
        if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in deleting rows from file***";
             return (EXIT_FAILURE);
         }
        LOG(INFO) << rows.size() << " rows deleted"<< endl ;
        rows.clear();
       delete[] alpha;  LOG(INFO)<<"Alpha deleted";
       delete[] veto;    LOG(INFO)<<"Veto deleted";
     }
      
    //loop for keyword writing/updating and writing history
    fits_get_num_hdus(fptr,&numhdu,&status);
    if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<"***Error in getting number of hdus in file"<<outfile<<"***";
             return (EXIT_FAILURE);
         }
    vector<string> vhistory;
    getHistory(vhistory);
    for(i=1;i<=numhdu;i++){
        fits_movabs_hdu(fptr,i,NULL,&status);
        sprintf(errstr,"***Error moving to HDU %d***",i);
         if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<errstr;
             return (EXIT_FAILURE);
         }
        fits_update_key(fptr,TSTRING,"CREATOR",modulename,NULL,&status);
        fits_write_date(fptr,&status);
        fits_write_chksum(fptr,&status);
        if(history==YES){
            for(int j=0;j<vhistory.size();j++){
                fits_write_history(fptr,vhistory[j].c_str(),&status);
                if(status){
                    fits_report_error(stderr,status);
                    LOG(ERROR)<<"***Error writing history***";
                    return (EXIT_FAILURE);
               }
             }
         }
     }
    fits_close_file(fptr,&status);
    sprintf(errstr,"***Error closing file %s ***",outfile);
    if(status){
             fits_report_error(stderr,status);
             LOG(ERROR)<<errstr;
             return (EXIT_FAILURE);
     }
    return (EXIT_SUCCESS);
 }

 //returns 0 if alpha satisfies the given value
 //returns 1 if alpha does not satisfies the range
 //returns -1 if input is incorrect
 int checkAlpha(unsigned char alpha,char *alphaval){
     if(strcmp(alphaval,"")==0 || strcmp(alphaval,"-")==0) return 0;
     else if(strcmp(alphaval,"0,1")==0 || strcmp(alphaval,"1,0")==0) return 0;
     else if(strcmp(alphaval,"1")==0){
         if(alpha==1) return 0;                 
         else return 1;
     }
     else if(strcmp(alphaval,"0")==0){
         if(alpha==0)  return 0;
         else return 1;
     }
     else{
         LOG(ERROR)<<"***alphaval not correct***";
         LOG(ERROR)<<"Valid inputs are '<space>','0','1','0,1'";
         return -1;
     } 
     return -1;
 }
 
 //returns 0 if veto value is within some range
 //returns 1 if veto value is not within any range
 int checkVeto(unsigned char veto,char *vetorange){
     if(strcmp(vetorange,"")==0 || strcmp(vetorange,"-")==0 || strcmp(vetorange,"0")==0   ) return 0;
     string vetorangestr=(string)vetorange;
     vector<string> range;
     parseString(vetorange,',',range);
     int min,max;
     for(int i=0;i<range.size();i++){
         vector<string> minmax;
         parseString(range[i],'-',minmax);
         min=atoi(minmax[0].c_str());
         max=atoi(minmax[1].c_str());
         //cout<<range[i]<<"\t"<<min<<"\t"<<max;
         if((int)veto>=min && (int)veto<=max){
             return 0;
         }
         else
             continue;
      }
     return 1;
 }
 
 int cztevtclean::getHistory(vector<string> &vhistory){
    //char *user=getlogin();
    strcpy(modulename,"cztevtclean_v");
    strcat(modulename,VERSION);

    char *user = getenv("USER");
	vhistory.push_back("Module run by "+(string)user);
    vhistory.push_back("Parameter List START for "+(string)modulename);
    vhistory.push_back("P1 infile="+(string)infile);
    vhistory.push_back("P2 outfile="+(string)outfile);
    vhistory.push_back("P3 alphaval="+(string)alphaval);
    vhistory.push_back("P4 vetorange="+(string)vetorange);
       
    if(clobber==YES) 
        vhistory.push_back("P5 clobber=yes");
    else
        vhistory.push_back("P5 clobber=no");
    
    if(history==YES)
        vhistory.push_back("P6 history=yes");
    else
        vhistory.push_back("P6 history=no");
    
    vhistory.push_back("Parameter List END");
    return (EXIT_SUCCESS);
}

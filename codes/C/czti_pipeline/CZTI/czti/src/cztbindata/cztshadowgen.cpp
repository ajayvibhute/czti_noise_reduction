/*
 * cztshadowgen
 *
 * Computes mask open fraction for each pixel for given quadrant
 *
 * Ajay Vibhute
 *
 * Edits:
 *
 * Maskelements recieved as input
 * Setting values in exposure table structure as used in bindata
 *
 * Mithun N P S
 * */
//Coordinate type

#include <cztbindata.h>

int calculate_renormalized_weights(ExposureTable &expTable,Badpix &badpix,int badpixThreshold,int quadID,string effareaFilename,string evtfilename)
{

	/* IMPORTANT NOTE: 
	 * In Badpix class badpixMap is indexed with [locy][locx].
	 * quad_badpixMap defined in this funciton is indexed with 
	 * [detx][dety] for that given quadrant
	 *
	 * Mithun N P S
	 * */

    int detx,dety,i,j;
    int index;
    EffArea effarea;
    char extname[FLEN_VALUE];

	vector < vector <unsigned char> > full_badpixMap;
	vector < vector <unsigned char> > quad_badpixMap(64);
	for(i=0;i<64;i++) quad_badpixMap[i].resize(64);

	full_badpixMap=badpix.get_badpix_map();


    //Reading effective area file
    quadToHDU(quadID, extname);
    if (effarea.read_effarea_file(effareaFilename, (string) extname)) {
        LOG(ERROR) << "Error in reading effective area file: " << effareaFilename <<
                " for quadrant " << quadID;
        return EXIT_FAILURE;
    }

	vector <float> effectiveArea;
	vector <float> temparea;
	vector <float> tempenergy;
    tempenergy.resize(1, 30.0);
	effectiveArea.resize(4096,0);


	//Assign the badpixMap to quadrant badpixMap for easier indexing and 
	//set weight for badpixels to zero
	for(detx=0;detx<64;detx++)
	{
		for(dety=0;dety<64;dety++)
		{
			if(quadID==0) quad_badpixMap[detx][dety]=full_badpixMap[64+dety][detx];
            if(quadID==1) quad_badpixMap[detx][dety]=full_badpixMap[64+dety][64+detx];
            if(quadID==2) quad_badpixMap[detx][dety]=full_badpixMap[dety][64+detx];
            if(quadID==3) quad_badpixMap[detx][dety]=full_badpixMap[dety][detx];			

			index=detx+dety*64;
			if(quad_badpixMap[detx][dety]>badpixThreshold) expTable.weights[index]=0;

			temparea = effarea.get_effective_area((unsigned char) detx, (unsigned char) dety,tempenergy);
			effectiveArea[index]=temparea[0];
		}
	}


	//Renormalization of maskweights. Maskweights for each detector is normalized
	//such that open and closed area in each detector is equal. Indexes i and j 
	//run over the detector number in x and y directions respectively
    
	double detRenormD,detRenormA;
	double denominator;
	double f2asum,fasum,pixarea;

	//for each detector in the quadrant
	for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
			detRenormD=0;
			detRenormA=0;
			f2asum=0;
			fasum=0;
			pixarea=0;

			//Compute renormalization factor and other factors for this detector
			for(detx=i*16;detx<(i+1)*16;detx++)
			{
				for(dety=j*16;dety<(j+1)*16;dety++)
				{
					index=detx+dety*64;
					if(quad_badpixMap[detx][dety]<=badpixThreshold)// for good pixels
					{
						pixarea=effectiveArea[index];
						detRenormD+=expTable.weights[index]*pixarea;
						detRenormA+=pixarea;
						f2asum+=expTable.openfrac[index]*expTable.openfrac[index]*pixarea;
						fasum+=expTable.openfrac[index]*pixarea;	
						//LOG(INFO)<<"DETNORM "<<expTable.weights[index];
					}
				}
			}

			detRenormD/=detRenormA;
			denominator=((2.0*f2asum/fasum)-(1+detRenormD))*(fasum);

			//if(fasum==0) LOG(INFO)<<i<<" "<<j<<" detRenormD "<<detRenormD<<" denominator "<<denominator<<" detRenormA "<<detRenormA;

			//Apply renormalization to mask-weights of this detector
            for(detx=i*16;detx<(i+1)*16;detx++)
            {
                for(dety=j*16;dety<(j+1)*16;dety++)
                {
                    index=detx+dety*64;
                    if(quad_badpixMap[detx][dety]<=badpixThreshold)
                        expTable.weights[index]=(expTable.weights[index]-detRenormD)/(denominator*16.0);
                }
            }

        
        }
    }

	return EXIT_SUCCESS;
}

int getShadow(float tx,float ty,int qid,float *shadow_pixels,int *maskElements, ExposureTable &exptable)
{
	int totalElements=TOTALROWS,status=0;
	int i=0,j=0,ii=0,jj=0,height=478;	
	float *temp_exp;
	float pix_height,pix_width,ref_pix=2.46;
	float frac,weight;
//	int *maskElements;

	temp_exp=(float*)malloc(sizeof(float)*DET_NUMPIXELS*DET_NUMPIXELS);
/*
	maskElements=(int*)malloc(sizeof(int)*TOTALROWS*COMP_COLS);

	getMaskPattern(MASKFILENAME,maskElements,qid+2);
	if(maskElements==NULL||temp_exp==NULL)
	{
		printf("Error(%s:%d): Unable to allocate memory\n",__FILE__,__LINE__);
		exit(0);
	}
*/	


	calculateTrans(tx,ty,height,maskElements,temp_exp);
	for(ii=0;ii<DET_NUMPIXELS;ii++)
	{
		if(ii%16==0||(ii+1)%16==0)
		{
			pix_width=2.28;
		}
		else
		{
			pix_width=2.46;
		}
		pix_width/=ref_pix;
		for(jj=0;jj<DET_NUMPIXELS;jj++)
		{
			if(jj%16==0||(jj+1)%16==0)
			{
				pix_height=2.28;
			}
			else
			{
				pix_height=2.46;
			}
			pix_height/=ref_pix;
			//temp_exp[ii*DET_NUMPIXELS+jj]*=(pix_height*pix_width);
		}
	}

	for(j=0;j<DET_NUMPIXELS;j++)
	for(i=0;i<DET_NUMPIXELS;i++)
	{
		shadow_pixels[(j*DET_NUMPIXELS)+i]=temp_exp[(j*DET_NUMPIXELS)+i];
		frac=shadow_pixels[(j*DET_NUMPIXELS)+i];
		weight=2*frac-1;
		status=exptable.set_exposure((unsigned char) i, (unsigned char)j, (unsigned char)i, (unsigned char)j, frac, weight);
		if(status){LOG(ERROR)<<"Unable to set exposure table";return EXIT_FAILURE;}
		//printf("%f\n",frac);
	}

	return EXIT_SUCCESS;
//	free(maskElements);
}



void calculateTrans(float tx,float ty,float height,int *maskPattern,float * transValues)
{
	double xmin=0,ymin=0,x_open=0,y_open=0,totalOpen=0,trans=0,thetaX=0.0,thetaY=0.0,thickness=0.5;
	double mask_lower_left_x=0,mask_lower_left_y=0,x=0,y=0,x_mask=0,y_mask=0;;
	int i=0,j=0;
	char temp[100];
	thetaX=ty*(MATHPI/180.0);
	thetaY=tx*(MATHPI/180);
	int totalClose=0;
	int ii=0,jj=0;	
	long index=0;
	int yy=0,xx=0;
	float yapitch=0,xpitch=0;
	float starty=0,startx=0,endy=0,endx=0;
	double lower_left_x=0,lower_left_y=0;
	int rowCounter=0,colCounter=0;	

	//All y's are X and all X's are y
	for(yy=0;yy<64;yy+=16,yapitch+=2.5)//this is for x
	{
	for(j=yy;j<(yy+16);j++)
	{
		
		starty=(COLS*(yy/16))+(yapitch/0.02);
		endy=(COLS*((yy/16)+1))+(yapitch/0.02);

		if(j==yy || j==yy+15)
		{
			rowCounter=114;
		}
		else
		{
			rowCounter=123;
		}
		lower_left_y=0;	
		
		for(xx=0,lower_left_y=0,xpitch=0;xx<64;xx+=16,xpitch+=2.5) //this is for y
		{
			for(i=xx;i<xx+16;i++)
			{

				
				startx=(COLS*(xx/16))+(xpitch/0.02);
				endx=(COLS*((xx/16)+1))+(xpitch/0.02);
				
				if(i==xx || i==xx+15)
				{
					colCounter=114;
				}
				else
				{
					colCounter=123;
				}

				x=lower_left_y;
				y=lower_left_x;
				x_mask=x+(height*(tan(thetaX)));
				y_mask=y+(height*(tan(thetaY)));
				x_mask/=0.02;
				y_mask/=0.02;	
				totalOpen=0;
				totalClose=0;
				for(jj=(int)y_mask;jj<(int)y_mask+rowCounter;jj++,index++)//this is for x
				{
				for(ii=(int)x_mask;ii<(int)x_mask+colCounter;ii++)// this is for y
				{
						if(jj>=starty && jj<endy && ii>=startx && ii<endx)
						{
							if(checkbit(maskPattern,ii,jj))
							{
								totalOpen++;
							}
							else
							{
								totalClose++;
							}
						}
					}
				}	
				transValues[i*DET_NUMPIXELS+j]=totalOpen/(colCounter*rowCounter);
				float thickBlockage=(thickness*tan(thetaY));//(thickness*tan(thetaY)));
				
				if(thickBlockage<0)
				{
					thickBlockage*=-1;	
				}	
				if(transValues[i*DET_NUMPIXELS+j]-thickBlockage<=0)
				{
					transValues[i*DET_NUMPIXELS+j]=0;
				}				
				else
				{
					transValues[i*DET_NUMPIXELS+j]-=thickBlockage;
				}
				thickBlockage=(thickness*tan(thetaX));
				if(thickBlockage<0)
				{
					thickBlockage*=-1;	
				}	
				if(transValues[i*DET_NUMPIXELS+j]-thickBlockage<=0)
				{
					transValues[i*DET_NUMPIXELS+j]=0;
				}				
				else
				{
					transValues[i*DET_NUMPIXELS+j]-=thickBlockage;
				}
				if(transValues[i*DET_NUMPIXELS+j]>1 || transValues[i*DET_NUMPIXELS+j]<0)
				{
					printf("Error:(%s:%d):Abnormal Trans value=%f\n",__FILE__,__LINE__,transValues[i*DET_NUMPIXELS+j]);
					exit(0);
				}
				if(i==xx || i==xx+15)
				{
					lower_left_y+=2.28;
				}
				else
				{
					lower_left_y+=2.46;
				}
			}
			lower_left_y+=2.5;
		}//end loop for y
		if(j==yy || j==yy+15)
		{
			lower_left_x+=2.28;
		}
		else
		{
			lower_left_x+=2.46;
		}
	}
		lower_left_x+=2.5;
	}//end loop for x
}



int checkbit(int *maskpattern,int i,int j)
{
	
	
	int temp_val=0,jj=0,bit_no=0;
	jj=j/25;
	bit_no=j%25;
	temp_val=maskpattern[i*COMP_COLS+jj];	
	if(temp_val&=(1<<bit_no))
	{
		return 1;
	}
	else
	{
		return 0;
	}	
}



void getMaskPattern(char*filename,int *buffer,int hduNo)
{
	fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    	int status,  nfound, anynull,hdutype,i,j;
    	long naxes[2], fpixel, nbuffer, npixels, ii;
    	float datamin, datamax, nullval;
    	status = 0;

    	if ( fits_open_file(&fptr, filename, READONLY, &status) )
    		printerror( status );

   	if ( fits_movabs_hdu(fptr, hduNo, &hdutype, &status) ) 
   		printerror( status );
	
    	if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
        	printerror( status );
	 
    	npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    	fpixel   = 1;
    	nullval  = 0;                /* don't check for null values in the image */
	while (npixels > 0)
    	{
      		nbuffer = npixels;
      		if ( fits_read_img(fptr, TINT, fpixel, nbuffer, &nullval,buffer, &anynull, &status) )
           		printerror( status ); 
      		npixels -= nbuffer;    
      		fpixel  += nbuffer;    
    	}
    	if ( fits_close_file(fptr, &status) )
        	printerror( status );
    	return;
}


void printerror( int status)
{
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */
       exit( status );    /* terminate the program, returning error status */
    }
    return;
}


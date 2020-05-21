/* 
 * @file fft.h
 * @author Tanul Gupta
 * @date Created on Oct 21, 2015, 6:08 PM
 * @brief FFT (radix 2)
 * @details 
 */

#include "errorHandler.h"
#include <iostream>
#include "utils.h"
#include <cmath>

using namespace std;

template <class T> vector <vector <T> > swap_correlation_matrix(vector <T> *data, int oversamplingFactor) ;
template<class T>
inline void SWAP(T &a, T &b) {
    T dum = a;
    a = b;
    b = dum;
}
/**
 * Function to evaluate fourier and inverse fourier transform using radix-2
 * algorithm.
 * @param data: data for which fourier transform needs to be taken
 * @param n
 * @param isign
 */
template <class T>
void four1(vector <T> *data, const int isign) {
    int n; //number of values in data
    int nn, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
    ErrorHandler errHandler;
    
    n = (*data).size()/2;
    if (n < 2 || n & (n - 1)) {
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Number of data points should be power of 2.";
        throw errHandler;
    }
    
    nn = n << 1; //actual size of data (real as well as complex parts)
    //Bit reversal started
    j = 1;
    for (i = 1; i < nn; i += 2) {
        if (j > i) {
            SWAP((*data)[j - 1], (*data)[i - 1]);
            SWAP((*data)[j], (*data)[i]);
        }
        m = n;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    //END Bit reversal ends here.
    
    mmax = 2;
    while (nn > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp*wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= nn; i += istep) {
                j = i + mmax;
                tempr = wr * ((*data)[j - 1]) - wi * ((*data)[j]);
                tempi = wr * (*data)[j] + wi * (*data)[j - 1];
                (*data)[j - 1] = (*data)[i - 1] - tempr;
                (*data)[j] = (*data)[i] - tempi;
                (*data)[i - 1] += tempr;
                (*data)[i] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

/**
 * Generates correlation matrix for oversampled data1 and data2;
 * @param data1
 * @param data2
 * @param correlationMatrix
 * @param nrows
 * @param ncols
 * @param oversamplingFactor
 */
template <class T>
vector <vector <T> > generate_correlation_matrix2D(vector <T> data1, vector <T> data2,
        vector <T> *correlationMatrix, int oversamplingFactor,float xshift,float yshift){
    ErrorHandler errHandler;
    vector <T> data1complex;
    vector <T> data2complex;
    vector <T> tmp_complex;
    vector <T> phase_complex;
    vector <T> tempCorrelationMatrix;
    vector < vector <T> > image;
    long i=0, j=0,k=0,ii=0,jj=0;
    int n_=0;
    float n02=0;
    float shiftfactor=2*M_PI*oversamplingFactor;
  //  float xshift,yshift;
    float fy=0,fx=0;
    long data1size=(data1).size();
    long data2size=(data2).size();
    long complexMatrixSize = 0;
/*
    if(xshift<0)
	yshift+=1;
*/
    xshift*=shiftfactor;
    yshift*=shiftfactor;
    if(data1size!=data2size){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Size of two datasets (to be correlated) is not equal.";
        throw errHandler;
    }
    
    complexMatrixSize = 2*data1size;
	LOG(INFO)<<"Data1size   "<<data1size;
    (*correlationMatrix).resize(complexMatrixSize+1, 0);
    data1complex.resize(complexMatrixSize+1, 0);
    data2complex.resize(complexMatrixSize+1, 0);
    phase_complex.resize(complexMatrixSize+1, 0);
    tmp_complex.resize(complexMatrixSize+1, 0);

    for (i = 0, j = 0; i < complexMatrixSize; i++) {
        if (i % 2 == 0) {
            //assigning real values at indices 1,3,5...
            data1complex[i] = data1[j];
            data2complex[i] = data2[j];
            j++;
        } else {
            //assigning complex value 0 at indices 0,2,4...
            data1complex[i] = 0;
            data2complex[i] = 0;
        }
		phase_complex[i]=0;
    }


    try{
        //taking fourier transform of data1 and data2.
        LOG(INFO) << "taking fourier transform of both data 1 and data2";
        four1(&data1complex, 1);
        four1(&data2complex, 1);
    } catch (ErrorHandler errHandler){
        throw errHandler;
    }
    
    n02 = (float) ((*correlationMatrix).size()/2);
	LOG(INFO)<<"N02  "<<n02;

	n_=(int)sqrt(data1size);
	LOG(INFO)<<"N:"<<n_<<"  Oversampling Factor: "<<oversamplingFactor;
	k=0;	
	for(i=0;i<n_;i++)
	{

	if(i<n_/2)
		ii=i;
	if(i>n_/2)
		ii=(i-n_);		

	for(j=0;j<n_;j++)
	{
		if(j<n_/2)
			jj=j;
		if(j>n_/2)
			jj=(j-n_);		

		fx=(ii/(n_*1.0))*xshift;
		fy=(jj/(n_*1.0))*yshift;

		phase_complex[k]=(cos(fx)*cos(fy)) - (sin(fx)*sin(fy)) ;
		phase_complex[k+1]=(sin(fx)*cos(fy)) + (cos(fx)*sin(fy));
		k+=2;


	}
	}


    for(i=0; i < complexMatrixSize; i+=2){
        (tmp_complex)[i] = data1complex[i]*data2complex[i] + data1complex[i+1]*data2complex[i+1];  
        (tmp_complex)[i+1] = data1complex[i+1]*data2complex[i] - data1complex[i]*data2complex[i+1];

  		(*correlationMatrix)[i]=((tmp_complex)[i]*phase_complex[i])-((tmp_complex)[i+1]*phase_complex[i+1]);
		(*correlationMatrix)[i+1] =((tmp_complex)[i+1]*phase_complex[i])+((tmp_complex)[i]*phase_complex[i+1]);


        (*correlationMatrix)[i] /= n02;
        (*correlationMatrix)[i+1] /= n02;
    }
    
    
    
    //taking inverse fourier transform of correlation matrix
    four1(correlationMatrix, -1);
    //getting absolute values
    tempCorrelationMatrix.resize(complexMatrixSize/2, 0);
    for(i=0, j=0; i<complexMatrixSize; i+=2, j++){
        tempCorrelationMatrix[j] = sqrt((*correlationMatrix)[i]*(*correlationMatrix)[i] +
                (*correlationMatrix)[i+1]*(*correlationMatrix)[i+1]);
    }
    (*correlationMatrix).clear();
    *correlationMatrix = tempCorrelationMatrix;
    image = swap_correlation_matrix(correlationMatrix, oversamplingFactor);
    return image;
    //keepCentral(correlationMatrix);
}

template <class T>
void oversample(vector <T> *data, int ysize, int xsize, int oversamplingFactor){
    int size = xsize*oversamplingFactor;
    vector <T> tempData; //to temporarily store data while changing the size of vector
    ErrorHandler errHandler;
    long i=0, j=0, k=0, l=0;
    tempData = *data;
    (*data).clear();
    if((tempData).size()!=ysize*xsize){
        errHandler.severity = errERROR;
        errHandler.errorStatus = IMPROPER_INPUT;
        errHandler.errorMsg = "Size of data: " + itoa((*data).size()) + " and ysize x xsize [" + itoa(ysize) + "x" + itoa(xsize) + "] do not match.";
        throw errHandler;
    }
    //Making an oversampled size vector
    (*data).resize(size*size, 0);
    for(i=0; i<ysize; i++){
        for(j=0; j<xsize; j++){
            for(k=i*oversamplingFactor; k<(i+1)*oversamplingFactor; k++){
                for(l=j*oversamplingFactor; l<(j+1)*oversamplingFactor; l++){
                    (*data)[k*size+l] = tempData[i*xsize+j];
                }
            }
        }
    }
}

template <class T>
void swap_correlation_matrix(vector <T> *data){
    long xsize = sqrt((*data).size());
    long ysize = sqrt((*data).size());
    int i=0, ii=0, j=0, jj=0;
    long n=ysize/2;
    T tempValue=0;
    
    for(i=0, ii=n; i<n; i++, ii++){
        for(j=0, jj=n; j<n; j++, jj++){
            tempValue = (*data)[i*xsize + j];
            (*data)[i*xsize + j] = (*data)[ii*xsize + jj];
            (*data)[ii * xsize + jj] = tempValue;
        }
    }
    
    for(i=n, ii=0; i<n+n; i++, ii++) {
        for (j = 0, jj = n; j < n; j++, jj++) {
            tempValue = (*data)[i*xsize+j];
            (*data)[i*xsize + j] = (*data)[ii*xsize+jj];
            (*data)[ii*xsize+jj]= tempValue;
        }
    }
}

template <class T>
vector <vector <T> > swap_correlation_matrix(vector <T> *data, int oversamplingFactor) {
    long xsize = sqrt((*data).size());
    long ysize = sqrt((*data).size());
    long imsize = 2*(15*oversamplingFactor) + 1; //both x and y axis size
    int i=0, ii=0, j=0, jj=0;
    int imid=15*oversamplingFactor;
    int jmid = imid;
    int yoff=0, xoff=0;
    vector <T> tempimg;
    vector < vector <T> > image;
    tempimg.resize(imsize, 0.0);
    image.resize(imsize, tempimg); // to store final image which is 2*(15*oversamp + 1)
    for(j=jmid; j>=0; j--){
        yoff = (jmid - j)*xsize;
        for(i=imid; i>=0; i--){
            image[j][i] = (*data)[(imid-i)+yoff];
        }
        for(i=imid+1; i<imsize; i++){
            image[j][i] = (*data)[(xsize+imid-i) + yoff];
        }
    }
    
    for(j=jmid+1; j<imsize; j++){
        yoff = (ysize + jmid-j)*xsize;
        for(i=imid; i>=0; i--){
            image[j][i] = (*data)[(imid-i)+yoff];
        }
        for (i = imid + 1; i < imsize; i++) {
            image[j][i] = (*data)[(xsize + imid - i) + yoff];
        }
    }
    
    return image;
    
}
template <class T>
int keepCentral(vector <T> *data){
    vector <T> temp;
    int i=0, j=0, ii=0, jj=0;
    int totalElements = sqrt((*data).size());
    LOG(INFO) << "totalElements : " << totalElements;
    int centralCount = totalElements/4;
    int startIndex = (totalElements/2) - (centralCount/2);
    int endIndex = startIndex + centralCount;
    
    temp.resize((*data).size(),0);
    temp = (*data); //assigning all (*data) values to temp;
    (*data).clear();
    (*data).resize((endIndex-startIndex)*(endIndex-startIndex), 0.0);
    for(i=startIndex, ii=0; i<endIndex; i++, ii++){
        for(j=startIndex, jj=0; j<endIndex; j++, jj++)
            (*data)[ii*centralCount + jj] = temp[i * totalElements + j];
    }
    LOG(INFO) << "dATA SIZE: " << (*data).size();
}





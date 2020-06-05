#Author : Ajay Ratheesh 

import os
import numpy as np
from astropy.io import fits,ascii
from matplotlib import pyplot as plt
import sys
import argparse

#How to run the code
#python bunchfile_fits2txt.py --inputdir NAME_OF_INPUT_DIR --bunchfile NAME_OF_BUNCHFILE --outdir NAME_OF_THE_OUTPUT_DIR


formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=formatter)
parser.add_argument("--inputdir",  help='Path to the input directory', 
                    type=str)
parser.add_argument("--bunchfile",  help='Name of the bunchfile (eg : *level2_bunch.fits)', 
                    type=str)
parser.add_argument("--outdir",  help='Path to the output directory', 
                    type=str)                  



def bunchfits2txt(**kwargs): 
	
	input_dir 			= kwargs['inputdir']
	bunch_file_name 	= kwargs['bunchfile']
	output_directory 	= kwargs['outdir']
	bunch_file = fits.open(input_dir+"/"+bunch_file_name)
	header_bunch = 'TIME'    +"\t"	+ 'NumEvent'+"\t" + 'Time_dfs'+"\t" + 'Time_dsl'+"\t" + 'DetID1'  +"\t" + 'DetID2'  +"\t" + 'DetID3'  +"\t" + 'DetID4'  +"\t" + 'DetID_fevt'+"\t"+'PixID_fevt'+"\t" +'DetID_sevt'+"\t"+'PixID_sevt'+"\t"+'DetID_tevt'+"\t"+'PixID_tevt'+"\t"+'RowNum'  +"\n"
	fmt_bunch  = ['%1.6f','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d','%d']
	

	for Q_n in range(1,5):
		outfilename = bunch_file_name.split('.')[0]+ "Q"+str(Q_n)+".txt"
		#~ outfile = open(output_directory+"/"+outfilename,"w")
		
		
		bunchlist = bunch_file[Q_n].data
		time_bunch   	= np.array(bunchlist['TIME'])
		bunch_length 	= np.array(bunchlist['NumEvent']).astype(int)
		bunch_dfs    	= np.array(bunchlist['Time_dfs']).astype(int)
		bunch_dsl    	= np.array(bunchlist['Time_dsl']).astype(int)
		bunch_detid1    = np.array(bunchlist['DetID1']).astype(int)
		bunch_detid2    = np.array(bunchlist['DetID2']).astype(int)
		bunch_detid3    = np.array(bunchlist['DetID3']).astype(int)
		bunch_detid4    = np.array(bunchlist['DetID4']).astype(int)
		bunch_detidfevt    = np.array(bunchlist['DetID_fevt']).astype(int)
		bunch_pixidfevt    = np.array(bunchlist['PixID_fevt']).astype(int)
		bunch_detidsevt    = np.array(bunchlist['DetID_sevt']).astype(int)
		bunch_pixidsevt    = np.array(bunchlist['PixID_sevt']).astype(int)
		bunch_detidtevt    = np.array(bunchlist['DetID_tevt']).astype(int)
		bunch_pixidtevt    = np.array(bunchlist['PixID_tevt']).astype(int)
		rownum    		   = np.array(bunchlist['RowNum']).astype(int)
		

		np.savetxt(output_directory+"/"+outfilename, np.transpose([time_bunch, bunch_dfs, bunch_dsl, bunch_length, bunch_detid1, bunch_detid2, bunch_detid3, bunch_detid4, rownum,  bunch_detidfevt, bunch_pixidfevt, bunch_detidsevt, bunch_pixidsevt, bunch_detidtevt, bunch_pixidtevt ]), delimiter='\t', header= header_bunch, fmt=fmt_bunch)


if __name__ == '__main__':
	

	args = parser.parse_args()
	bunchfits2txt(**args.__dict__)

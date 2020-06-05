'''
Author: Ajay Ratheesh
		University of Rome "Tor Vergata"
		27/01/2020
'''
from __future__ import division
import argparse
import os
import numpy as np
from astropy.io import fits,ascii
import math
import sys
#~ from matplotlib import pyplot as plt
#~ import scipy.stats.distributions as ps
#~ from scipy.optimize import curve_fit


formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=formatter)
parser.add_argument("--indir",  help='Path to the input directory', 
                    type=str)  
parser.add_argument("--Eventfile",  help='Path to the cleaned single event file (eg : *quad_clean.evt)', 
                    type=str)
parser.add_argument("--outdir",  help='Path to the output directory', 
                    type=str)  


def create_modified_eventfile(**kwargs): 
	#reading files 	 
	input_dir = kwargs['indir']
	evt_file_name 	= kwargs['Eventfile']
	output_directory 	= kwargs['outdir']
	
	evt_file = fits.open(input_dir+evt_file_name)
	header =[None]*13
	data =[None]*13
	name = [None]*13
	for i in range(0,13):
		header[i] = evt_file[i].header
		data[i]	  = evt_file[i].data
		name[i]	  = evt_file[i].name

	filename_new = evt_file_name.split(".")
	new_evtfilename = filename_new[0]+"_cztntickCorrected."+filename_new[1]
	table = [None]*4
	new_table = [None]*13


	for Q_n in range(1,5):
	#############################Reading Event file################################################################
		time_evt 		= evt_file[Q_n].data['TIME']
		pha_evt 		= evt_file[Q_n].data['PHA']
		cztseccnt_evt 	= evt_file[Q_n].data['CZTSECCNT']
		cztntick_evt 	= evt_file[Q_n].data['CZTNTICK']
		PI_evt 			= evt_file[Q_n].data['PI']
		energy_evt 		= evt_file[Q_n].data['ENERGY']
		det_id 			= evt_file[Q_n].data['DetID']
		pix_id 			= evt_file[Q_n].data['pixID']
		alpha			= evt_file[Q_n].data['alpha']
		veto 			= evt_file[Q_n].data['veto']
		detx			= evt_file[Q_n].data['DETX']
		dety 			= evt_file[Q_n].data['DETY']
		
		#~ print (type(cztntick_evt[0]))
		cztseccnt_evt_corr = cztseccnt_evt.astype(float) - cztseccnt_evt.astype(int) 
		

		
		cztntick_evt_temp = cztseccnt_evt_corr *50000 #10E6 / 20
		cztntick_evt_temp = np.round(cztntick_evt_temp,0)
		
		
		cztntick_evt_corr = cztntick_evt_temp.astype('I')
		#~ cztntick_evt = cztntick_evt.astype('I')

		
		#~ print (str(min(cztntick_evt_corr)) + '\t' + str(max(cztntick_evt_corr)) + '\n')
		#~ print (str(min(cztntick_evt)) + '\t' + str(max(cztntick_evt)) + '\n')
		

		
		#~ for i in range(0,len(cztntick_evt)):
			#~ if cztntick_evt_new[i] != cztntick_evt[i]:
			#~ if cztseccnt_evt_new[i] ==0.0:
				#~ print (str(i)+"\t"+str(cztseccnt_evt[i-1])+"\t"+str(cztseccnt_evt[i])+"\t"+str(cztntick_evt_temp[i])+"\t"+str(cztntick_evt_new[i]) + "\t" + str(cztntick_evt[i]) + "\n")

		#~ print (time_evt[0])
		#~ time_evt = np.round(time_evt, 5)
		#~ print (time_evt[0])
		
		
		time_evt_new			= fits.Column(name='TIME', format='D',array = time_evt)		
		pha_evt_new				= fits.Column(name='PHA', format='I',bzero= 32768, bscale =1.0, array = pha_evt)
		cztseccnt_evt_new		= fits.Column(name='CZTSECCNT', format='D',array = cztseccnt_evt)
		cztntick_evt_new		= fits.Column(name='CZTNTICK', format='I', bzero= 32768, bscale =1.0, array = cztntick_evt_corr)
		PI_evt_new				= fits.Column(name='PI', format='I',array = PI_evt)
		energy_evt_new			= fits.Column(name='ENERGY', format='E',array = energy_evt)
		det_id_new				= fits.Column(name='DetID', format='B',array = det_id)
		pix_id_new 				= fits.Column(name='pixID', format='B',array = pix_id)
		alpha_new 				= fits.Column(name='alpha', format='B',array = alpha)
		veto_new	 			= fits.Column(name='veto', format='I',bzero= 32768, bscale =1.0, array = veto)
		detx_new	 			= fits.Column(name='DETX', format='B',array = detx)
		dety_new	 			= fits.Column(name='DETY', format='B',array = dety)
	                   		
		table[Q_n-1] = fits.ColDefs([time_evt_new,cztseccnt_evt_new,cztntick_evt_new,pha_evt_new,det_id_new,pix_id_new,detx_new,dety_new,veto_new,alpha_new,PI_evt_new,energy_evt_new])
		new_table[Q_n-1] = fits.BinTableHDU.from_columns(table[Q_n-1], name=name[Q_n],header=header[Q_n])

	for i in range(4,12):
		new_table[i] = fits.BinTableHDU.from_columns(data[i+1],name=name[i+1],header=header[i+1])
	

	hdu_primary = fits.PrimaryHDU(header=header[0])


	hdulist = fits.HDUList([hdu_primary,new_table[0],new_table[1],new_table[2],new_table[3],new_table[4],new_table[5],new_table[6],new_table[7],new_table[8],new_table[9],new_table[10],new_table[11]])


	hdulist.writeto(output_directory + new_evtfilename,clobber='Yes')



if __name__ == '__main__':
	
	args = parser.parse_args()
	create_modified_eventfile(**args.__dict__)

'''
table = fits.ColDefs([new_time,new_rate,new_error,new_fracexp])
new_table = fits.BinTableHDU.from_columns(table, name='RATE',header=header_1)
hdu_primary = fits.PrimaryHDU(header=header_0)

hdulist = fits.HDUList([hdu_primary, new_table])
hdulist.writeto(new_lcfilename,clobber='Yes')
'''

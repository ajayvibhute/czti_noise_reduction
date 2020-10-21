'''
Author: Ajay Ratheesh
		University of Rome "Tor Vergata"
		25/01/2019
'''
from __future__ import division
import os
import numpy as np
import numpy.ma as ma
from astropy.io import fits,ascii
import math
import argparse
import sys
from matplotlib import pyplot as plt
#import scipy.stats.distributions as ps
#from scipy.optimize import curve_fit



formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=formatter)
parser.add_argument("--indir",  help='Path to the input directory', 
                    type=str)  
#~ parser.add_argument("--Eventfile",  help='cleaned single event file name(eg : *quad_clean.evt)', 
                    #~ type=str)
parser.add_argument("--doubleEventfile",  help='cleaned double event file name (eg : *quad_clean.dblevt)', 
                    type=str)
#~ parser.add_argument("--livetimefile",  help='livetime file name (eg : *quad_livetime.fits)', 
                    #~ type=str)                   
parser.add_argument("--badpixfile",  help='badpix file name (eg : *quad_badpix.fits)', 
                    type=str)  
parser.add_argument("--caldb_path",  help='Path to the CALDB bcf folder', 
                    type=str, default="/home/user/software/caldb/data/as1/czti/bcf/")      
                                 
parser.add_argument("--caldb_badpixfilename",  help='CALDB Badpix file name(eg : AS1cztbadpix20160908v01.fits', 
                    type=str, default="AS1cztbadpix20160908v01.fits")

parser.add_argument("--outdir",  help='Path to the output directory', 
                    type=str)       
parser.add_argument("--sigma_good",  help='Sigma for DPH flagging for good pixels ', 
                    type=int, default=4) 
parser.add_argument("--sigma_banana",  help='Sigma for DPH flagging for banana pixels ', 
                    type=int, default=3)                
parser.add_argument("--badpix_outputfiletype",  help='Either txt or fits',
                    type=str, default="txt")                  
                    
                              
                    
#Function read the start and stop point from the image
def onclick(event):
	global ix, iy
	ix, iy = event.xdata, event.ydata
	global interval
	interval.append((ix, iy))
	if len(interval) == 2:
		fig.canvas.mpl_disconnect(cid)
		plt.close(1)
		return
		
def select_events(time_dblevt_,detx_dbl_, dety_dbl_, energy_dblevt_, flag, tstart, tstop):
	
	#Selecting events from flag
	time_dblevt_	=	time_dblevt_[flag]
	detx_dbl_		=   detx_dbl_[flag]
	dety_dbl_		=   dety_dbl_[flag]
	energy_dblevt_	=	energy_dblevt_[flag]

	#Selecting time
	time_index_cut		= [(time_dblevt_ > tstart)&(time_dblevt_ < tstop)][0]
	time_dblevt_curr  	= time_dblevt_[time_index_cut]
	detx_dbl_curr		= detx_dbl_[time_index_cut]
	dety_dbl_curr		= dety_dbl_[time_index_cut]
	energy_dblevt_curr 	= energy_dblevt_[time_index_cut]
	
	return time_dblevt_curr, detx_dbl_curr, dety_dbl_curr, energy_dblevt_curr

def calculate_distance(x1,y1,x2,y2):
    d = np.sqrt(  np.power( (x1-x2) ,2) + np.power( (y1-y2) ,2)) 
  
    return abs(d)

def calculate_doubleEventrates(detx_,dety_,pixelid_,energy_,pixflag_,exposure_,comment_):
	
	detx_1	= detx_[::2]
	dety_1	= dety_[::2]
	pixelid_1 = pixelid_[::2]
	energy_1	= energy_[::2]
	
	detx_2	    = detx_[1::2]
	dety_2	    = dety_[1::2]
	pixelid_2   = pixelid_[1::2]
	energy_2	= energy_[1::2]
	
	index_1_noise = np.array([False]*len(pixelid_1))
	index_2_noise = np.array([False]*len(pixelid_2))
	compton_index = np.array([False]*len(pixelid_2))
	double_index = np.array([False]*len(pixelid_2))
	
	#~ print (len(detx_))
	#~ print (len(detx_1))
	#~ print (len(detx_2))

	
	for i in range(0,len(pixelid_1)):
		
		distance = calculate_distance(int(detx_1[i]), int(dety_1[i]), int(detx_2[i]), int(dety_2[i]))
		
		if distance<1.5:
			double_index[i] = True
		if (100.0 < energy_1[i]+energy_2[i] < 360.0) and distance<1.5:
			compton_index[i] = True
		if pixflag_[pixelid_1[i]]>1:
			index_1_noise[i] = True
		if pixflag_[pixelid_2[i]]>1:
			index_2_noise[i] = True
			
	index_noise = np.add(index_1_noise,index_2_noise)
	index_genuine = ~index_noise
	
	#~ print ("Rate of noise events = "+ str(np.sum(index_noise)/exposure_) )
	index_double  = np.multiply(double_index, index_genuine)
	index_compton = np.multiply(compton_index, index_genuine)
	print ("---------"+comment_+"--------")
	print ("Rate of total double events = "+ str(np.sum(index_genuine)/exposure_) )
	print ("Rate of neighbouring double events = "+ str(np.sum(index_double)/exposure_))
	print ("Rate of compton events = "+ str(np.sum(index_compton)/exposure_))
	
	return

def DPH_cut(detx_dbl_, dety_dbl_, sigma):
	#########################DPH analysis of good pixels##########################
	DPH = np.histogram2d(detx_dbl_, dety_dbl_,bins=64)[0]
	num_pix_dead = 4096-len(np.nonzero(DPH)[0])
	DPH_nonzero = DPH[np.nonzero(DPH)]
	#~ print (np.mean(DPH_nonzero))
	#~ print (np.std(DPH_nonzero))
	DPH_ind 			= np.where(	DPH > (np.mean(DPH_nonzero)	+	sigma *np.std(DPH_nonzero))	)
	#~ DPH_outlier 	= DPH_ind
	#~ print ("Number of DPH outliers in good pixels = "+str(len(DPH_ind[0])))
	#~ print ("good DPH mean = "+str(np.mean(DPH)))
	#~ print ("good DPH mean non zero = "+str(np.mean(DPH_nonzero)))
	#~ print ("good DPH std  = "+str(np.std(DPH)))
	#~ print ("good DPH std non zero = "+str(np.std(DPH_nonzero)))
	#print (DPH[DPH_ind])
	#~ plt.imshow(DPH)
	#~ plt.legend()
	#~ plt.colorbar()
	#~ plt.show()
	#~ plt.savefig(str(sr_n)+"_"+str(Q_n)+"_good_dph.png")
	#~ plt.clf()
	#~ neighbouring_spectra(time_dblevt_curr, DPH_outlier_good, DPH[DPH_outlier_good], energy_dblevt_curr, detx_dbl_curr, dety_dbl_curr)
		
	outlier_index= np.array([False]*len(detx_dbl_))
	for j in range(0,len(DPH_ind[0])):
		outlier_index = np.add(	outlier_index, np.multiply(	[ detx_dbl_==DPH_ind[0][j] ][0], [ dety_dbl_==DPH_ind[1][j] ][0]	)	)
	
	return DPH_ind, outlier_index
	
	#selecting the noise events
	#time_dph_outlier_good 	= time_dblevt_curr[outlier_index_good]
	#energy_dph_outlier_good 	= energy_dblevt_curr[outlier_index_good]
	#time_dph_normal_good	= time_dblevt_curr[~outlier_index_good]
	#energy_dph_normal_good 	= energy_dblevt_curr[~outlier_index_good]		
		
		
		
def runDPHanalysis(**kwargs): 
	#reading files 	 
	
	ind_dir = kwargs['indir']
	#~ evt_file_name = kwargs['Eventfile']
	dblevt_file_name = kwargs['doubleEventfile']
	#~ livetime_file_name = kwargs['livetimefile']
	badpix_file_name = kwargs['badpixfile']
	caldb_badpix_file_name = kwargs['caldb_badpixfilename']
	caldb_path = kwargs['caldb_path']
	output_directory = kwargs['outdir']
	sigma_good   = kwargs['sigma_good']
	sigma_banana = kwargs['sigma_banana']
	output_file_type = kwargs['badpix_outputfiletype']
	print("\n\n")
	print ("OUTPUT FILES WILL BE SAVED IN DIRECTORY = "+output_directory)
	#~ print (evt_file_name)
	#~ print (dblevt_file_name)
	#~ print (livetime_file_name)
	#~ print (badpix_file_name)
	#~ print (caldb_path)
	#~ print (caldb_badpix_file_name)
	#~ print (output_directory)
	#~ print (sigma_good)
	#~ print (sigma_banana)
	
	
	table = [None]*4
	new_table = [None]*14
	
	#~ path = "/home/user/Ajay/tifr/GRB_pol/GRBs/161218B/"
	#~ evt_file  = fits.open(ind_dir+'/'+evt_file_name)
	#~ livetime_file = fits.open(ind_dir+'/'+livetime_file_name)
	dblevt_file = fits.open(ind_dir+'/'+dblevt_file_name)
	badpix_file  = fits.open(ind_dir+'/'+badpix_file_name)
	#~ caldb_badpix = fits.open(caldb_path+caldb_badpix_file_name)
	
	
	header =[None]*5
	data_badpix =[None]*5
	name = [None]*5
	for j in range(0,5):
		header[j] = badpix_file[j].header
		data_badpix[j]	  = badpix_file[j].data
		name[j]	  = badpix_file[j].name
	
	#~ sr_n = 1
	#~ time_live = livetime_file[3].data['TIME']
	#~ lc_nrows = len(time_live)
	#~ lc_total = [0]*lc_nrows
	initial_time = dblevt_file[1].data['TIME']
	tstart = initial_time[0]
	tstop = initial_time[len(initial_time)-1]
	#~ print (tstart)
	#~ print (tstop)

	#CALDB file 
	caldb_badpix_file = caldb_path + caldb_badpix_file_name 
	#~ "/home/user/software/caldb/data/as1/czti/bcf/AS1cztbadpix20160908v01.fits"
	caldb_badpix = fits.open(caldb_badpix_file)
	
	for Q_n in range(1,5):
		print("\n\n")
		print ("##################################  Q"+str(Q_n-1)+" DPH analysis ######################################")
		
	
		#Reading the data from event file
		time_dblevt = dblevt_file[Q_n].data['TIME']
		pha_dblevt = dblevt_file[Q_n].data['PHA']
		cztseccnt_dblevt = dblevt_file[Q_n].data['CZTSECCNT']
		cztntick_dblevt = dblevt_file[Q_n].data['CZTNTICK']
		PI_dblevt = dblevt_file[Q_n].data['PI']
		energy_dblevt = dblevt_file[Q_n].data['ENERGY']
		det_id_dbl = dblevt_file[Q_n].data['DetID']
		pix_id_dbl = dblevt_file[Q_n].data['pixID']
		alpha_dbl = dblevt_file[Q_n].data['alpha']
		veto_dbl = dblevt_file[Q_n].data['veto']
		detx_dbl = dblevt_file[Q_n].data['DETX']
		dety_dbl = dblevt_file[Q_n].data['DETY']
		pixel_id_dbl = det_id_dbl*256+pix_id_dbl
		pixel_flag = caldb_badpix[Q_n].data['PIX_FLAG']
		detx_dbl = detx_dbl.astype(int)
		dety_dbl = dety_dbl.astype(int)
		pixel_flag_currdat = badpix_file[Q_n].data['PIX_FLAG']
		detid_flag_currdat = badpix_file[Q_n].data['DETID']
		pixid_flag_currdat = badpix_file[Q_n].data['PIXID']
		exposure = dblevt_file[Q_n].header['EXPOSURE']
	
		banana_flag = np.array([False]*len(time_dblevt))
		good_flag   = np.array([False]*len(time_dblevt))
		neigh_flag = np.array([False]*len(time_dblevt))
		side_flag   = np.array([False]*len(time_dblevt))
		corner_flag   = np.array([False]*len(time_dblevt))
		
		#Printing the double event and Compton event rates
		calculate_doubleEventrates(detx_dbl,dety_dbl,pixel_id_dbl,energy_dblevt,pixel_flag_currdat,exposure,"Before DPH cut")
	
		#Getting the flag for banana and good pixels and correcting the energy for banana pixels
		for i in range(0,len(time_dblevt)):
			if i >0:
				del_time = time_dblevt[i]-time_dblevt[i-1]
				if abs(del_time) < 30.0E-6:
					distance = np.sqrt( float( (detx_dbl[i]-detx_dbl[i-1])**2    +  (dety_dbl[i]-dety_dbl[i-1])**2   ) )
					if distance < 1.415:
						neigh_flag[i]   = True
						neigh_flag[i-1] = True
						if distance == 1.0:
							
							side_flag[i]   = True
							side_flag[i-1] = True
						else: 
							corner_flag[i]   = True
							corner_flag[i-1] = True
			
			
			if pixel_flag[pixel_id_dbl[i]]==1:
				banana_flag[i] = True
				energy_dblevt[i] = energy_dblevt[i]*3.7
				
			if pixel_flag[pixel_id_dbl[i]]==0:
				good_flag[i] = True
		
		#good_flag   = np.multiply(good_flag, neigh_flag)
		#banana_flag = np.multiply(banana_flag, neigh_flag)
		#~ print (np.sum(side_flag))
		#~ print (np.sum(corner_flag))
		
		good_side_flag 		= np.multiply(good_flag, side_flag)
		banana_side_flag 	= np.multiply(banana_flag, side_flag)
		good_corner_flag 	= np.multiply(good_flag, corner_flag)
		banana_corner_flag 	= np.multiply(banana_flag, corner_flag)
		#~ print (np.sum(good_side_flag))
		#~ print (np.sum(good_corner_flag))
		
		
		#~ time_dblevt_good,detx_dbl_good,dety_dbl_good,energy_dblevt_good =  select_events(time_dblevt,detx_dbl, dety_dbl, energy_dblevt, good_flag, tstart, tstop)
		time_dblevt_banana, detx_dbl_banana, dety_dbl_banana, energy_dblevt_banana =  select_events(time_dblevt,detx_dbl, dety_dbl, energy_dblevt, banana_flag, tstart, tstop)
		
		time_dblevt_good_side,detx_dbl_good_side,dety_dbl_good_side,energy_dblevt_good_side =  select_events(time_dblevt,detx_dbl, dety_dbl, energy_dblevt, good_side_flag, tstart, tstop)
		time_dblevt_good_corn,detx_dbl_good_corn,dety_dbl_good_corn,energy_dblevt_good_corn =  select_events(time_dblevt,detx_dbl, dety_dbl, energy_dblevt, good_corner_flag, tstart, tstop)
		DPH_ind_gs,outlier_index_good_side = DPH_cut(detx_dbl_good_side, dety_dbl_good_side, sigma_good)
		DPH_ind_gc,outlier_index_good_corner = DPH_cut(detx_dbl_good_corn, dety_dbl_good_corn, sigma_good)
		DPH_ind =  np.append(DPH_ind_gs, DPH_ind_gc, axis =1 )

		print ("No of outliers in side pixels (good) = " + str(len(DPH_ind_gs[0])))
		print ("No of outliers in corner pixels (good) = " + str(len(DPH_ind_gc[0])))
		outlier_index_good = np.append(outlier_index_good_side, outlier_index_good_corner )
		#~ DPH_ind,outlier_index_good = DPH_cut(detx_dbl_good, dety_dbl_good, sigma_good)
		DPH_outlier_good = DPH_ind
		print ("Number of DPH outliers in good pixels = "+str(len(DPH_ind[0])))

		
		DPH_ind , outlier_index_banana = DPH_cut(detx_dbl_banana, dety_dbl_banana, sigma_banana)
		DPH_outlier_banana = DPH_ind
		print ("Number of DPH outliers in banana pixels = "+str(len(DPH_ind[0])))

	 
	
		new_pixel_flag = pixel_flag_currdat
	
		detx_new = np.array([None]*64)
		dety_new = np.array([None]*64)
	
		outfile = open(output_directory+"/badpixlist_Q"+str(Q_n-1)+".txt", "w")
	
	
		for i in range(0,4096):
			
			detx_temp =((detid_flag_currdat[i]%4)*16)+(pixid_flag_currdat[i]%16)
			dety_temp =(int(detid_flag_currdat[i]/4)*16)+int(pixid_flag_currdat[i]/16)
				
			if Q_n==1 or Q_n==4:
				dety_temp=63-dety_temp
			else:
				detx_temp=63-detx_temp
				
	
			if new_pixel_flag[i]==0:
				for j in range(0,len(DPH_outlier_good[0])): 
					#~ print (str(i)+ " "+ str(j)+"\n")
					#~ print (str(DPH_outlier_good[0][j])+ "\t"+ str(detx_temp)+"\t"+str(DPH_outlier_good[0][j])+ "\t"+ str(dety_temp)+"\n")
					
					if DPH_outlier_good[0][j]==detx_temp and DPH_outlier_good[1][j]==dety_temp:
						#print ("----------------------------------------------_REACHED HERE-------------------------------------")
						new_pixel_flag[i]=5
						break
			if new_pixel_flag[i]==1:
				for j in range(0,len(DPH_outlier_banana[0])): 
					if DPH_outlier_banana[0][j]==detx_temp and DPH_outlier_banana[1][j]==dety_temp:
						#print ("----------------------------------------------_REACHED HERE-------------------------------------")
						new_pixel_flag[i]=5
						break
			
			

		
			
			if output_file_type == "txt":
				outfile.write("	   "+str(detid_flag_currdat[i])+"	   "+str(pixid_flag_currdat[i])+"	   "+str(new_pixel_flag[i])+"\n")
				
				
		if output_file_type == "fits":
			detid_new		= fits.Column(name='DETID', format='B',array = np.array(detid_flag_currdat))
			pixid_new		= fits.Column(name='PIXID', format='B',array = np.array(pixid_flag_currdat))
			flag_new		= fits.Column(name='PIX_FLAG', format='B',array = np.array(new_pixel_flag))
			table[Q_n-1] = fits.ColDefs([detid_new,pixid_new,flag_new])
			new_table[Q_n-1] = fits.BinTableHDU.from_columns(table[Q_n-1], name=name[Q_n],header=header[Q_n])
	
		calculate_doubleEventrates(detx_dbl,dety_dbl,pixel_id_dbl,energy_dblevt,pixel_flag_currdat,exposure,"After DPH cut")
	
	if output_file_type == "fits":
		hdu_primary = fits.PrimaryHDU(header=header[0])
		hdulist = fits.HDUList([hdu_primary,new_table[0],new_table[1],new_table[2],new_table[3]])
		hdulist.writeto(output_directory+"/badpixlist.fits",overwrite='Yes')

if __name__ == '__main__':
    args = parser.parse_args()
    runDPHanalysis(**args.__dict__)


'''
Author:Ajay Ratheesh, TIFR
Date : 22-06-2020
'''
from __future__ import division
import os
import numpy as np
from astropy.io import fits,ascii
import math
import sys
import collections
from matplotlib import pyplot as plt
import scipy.stats.distributions as ps
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=formatter)
parser.add_argument("--indir",  help='Path to the input directory', 
                    type=str)  
parser.add_argument("--eventfile",  help='double event file name(eg : *.dblevt)', 
                    type=str)
parser.add_argument("--timebin", help='To be entered in milli seconds', 
                    type=float)                                      
parser.add_argument("--Quad",  help='Quadrant number to be processed (0,1,2,3)', 
                    type=int)                    
parser.add_argument("--outdir",  help='event file name(eg : *.evt)', 
                    type=str) 
parser.add_argument("--ifplotornot",  help='y to plot, n for not', 
                    type=str)                                      
                   
def get_total_exposure(gtitstart_, gtitstop_, tstart, tstop ): 
	expo= 0.0
	for i in range(len(gtitstart_)):
		if ( tstart <= gtitstart_[i] <= tstop ) and ( tstart <= gtitstop_[i] <= tstop ): 
			expo = expo + (gtitstop_[i]-gtitstart_[i])
			
		elif (gtitstart_[i] < tstart) and ( tstart <= gtitstop_[i] <= tstop ):
			expo = expo + (gtitstop_[i]-tstart)
		
		elif ( tstart <= gtitstart_[i] <= tstop ) and (gtitstop_[i] > tstop):
			expo = expo + (tstop-gtitstart_[i])
		elif (gtitstart_[i] < tstart) and (gtitstop_[i] > tstop): 
			expo = expo + tstop-tstart
	
	return expo
			 
			
			

def run_poisson_data_test(**kwargs): 
	
	#~ plt.rc('axes',linewidth=2)
	#~ plt.rc('text',usetex=True)
	#~ plt.rcParams['text.latex.preamble']=[r'\boldmath']
	padding=8
	
	
	ind_dir 	  		= kwargs['indir']
	evt_file_name 		= kwargs['eventfile']
	binsize 			= kwargs['timebin']/1000.0
	if_plot 			= kwargs['ifplotornot']
	output_directory 	= kwargs['outdir']
	Q_n					= kwargs['Quad']+1
	
	outfile_name = evt_file_name.split('.')[0]+"_poissontest_Q"+str(Q_n)+".txt"
	outfile = open(output_directory+"/"+outfile_name, "w")

	evt_file 		= fits.open(ind_dir+"/"+evt_file_name)
	time_evt		= evt_file[Q_n].data['TIME']
	energy_evt 		= evt_file[Q_n].data['ENERGY']
	det_id 			= evt_file[Q_n].data['DetID']
	alpha 			= evt_file[Q_n].data['alpha']
	veto 			= evt_file[Q_n].data['veto']
	exposure 		= evt_file[Q_n].header['exposure']
	gti_tstart 		= evt_file[Q_n+8].data['START']
	gti_tstop 		= evt_file[Q_n+8].data['STOP']
	
	time_evt_1		= time_evt[::2]
	energy_evt_1 	= energy_evt[::2]
	det_id_1 		= det_id[::2]	
	alpha_1			= alpha[::2]	
	veto_1			= veto[::2]		
	
	time_evt_2		= time_evt[1::2]
	energy_evt_2 	= energy_evt[1::2]
	det_id_2 		= det_id[1::2]	
	alpha_2			= alpha[1::2]	
	veto_2			= veto[1::2]		

	index_veto_1  = [alpha_1==0][0]
	index_alpha_1 = [veto_1 ==0][0]
	index_veto_2  = [alpha_2==0][0]
	index_alpha_2 = [veto_2 ==0][0]
	index_veto = np.add(index_veto_1,index_veto_2)
	index_alpha = np.add(index_alpha_1,index_alpha_2)
	
	
	index_veto_alpha = np.add(index_veto, index_alpha)
	time_evt_1		= time_evt_1		[index_veto_alpha]
	energy_evt_1 	= energy_evt_1 	[index_veto_alpha]	
	det_id_1 		= det_id_1 		[index_veto_alpha]
	#~ alpha 			= alpha 		[index_veto_alpha]
	#~ veto 			= veto 			[index_veto_alpha]
    
	
	UT_0 = int(time_evt_1[0])
	UT_last = int(time_evt[len(time_evt_1)-1])

	hist_total = np.array([0]*39)

	T1=float(UT_0)
	T2=float(UT_0)+100.0
	hist_exp =  np.array([0]*39)
	hist_exp_total =  np.array([0]*39)
	
	x = 1.0/binsize # number of bins per second.
	#print x
	while True:
		
		if T2>=UT_last:
			T2 = float(UT_last)
			
		total_expo = get_total_exposure(gti_tstart, gti_tstop, T1, T2 )
		
		total_time = T2-T1
		
		#~ print (str(total_time)+ "\t" + str(total_expo) + "\n")
		
		if  (total_time < total_expo) or (total_expo<=0.0) :
			T1=T2
			T2= T2+100.0 
			if T1>= UT_last: 
				break
			continue
			
		
		num_of_expected_zeros = (total_time - total_expo)/binsize
		index_curr = [(time_evt_1>=T1)&(time_evt_1<=T2)][0]
		time_evt_curr = time_evt_1[index_curr]
		det_id_curr = det_id_1[index_curr]
		sum_0 = 0
		for j in range(0,16):
			index_det = [det_id_curr==j][0]
			time_evt_curr_det = time_evt_curr[index_det]
			time_evt_lc, bins = np.histogram(time_evt_curr_det, bins= np.arange(T1,T2,binsize))
		
			binhist = np.arange(0,40,1)
			hist_curr ,bin_e = np.histogram(time_evt_lc, bins=binhist)
			hist_curr[0] = hist_curr[0] - num_of_expected_zeros
			hist_total = np.add(hist_total,hist_curr)
			bin_e = np.delete(bin_e,len(bin_e)-1)
			sum_0 = sum_0+np.sum(time_evt_lc)
		

		no_of_points = total_time*16.0*x
		average = float(sum_0/no_of_points)

		for i in bin_e:
			hist_exp[i] = (np.exp(-average)*(average**i)/float(math.factorial(i)))*(no_of_points)
		#print hist_exp[i]
		hist_exp_total = hist_exp_total+hist_exp
	#print hist_exp
		T1=T2
		T2= T2+100.0
	
	
	#print str(T1)+"\t"+str(T2)+"\t"+str(UT_last)
		if T1>= UT_last: 
			break 
	

	
	#~ plt.title("1ms Histogram of counts per detector")

	plt.xlabel("Counts",fontsize='15')
	plt.ylabel("N", fontsize='15')
	plt.yscale("log")
	plt.xlim(0,10)
	plt.ylim(0,1e7)
	#~ print bin_e
	#~ bin_e = bin_e+0.5
	#~ plt.plot(bin_e,hist_total, 'r', lw=1.0, label='Observed')
	print (bin_e)
	plt.plot(bin_e,hist_exp_total,color='blue', label='Expected')

	plt.errorbar(bin_e,hist_total, yerr = np.sqrt(hist_total),color='red', lw=1.0,label='Observed')
	#~ plt.errorbar(bin_e-0.5,hist_exp, yerr = np.sqrt(hist_total), fmt='.r', label='gvjg', lw=0.5)
	plt.legend()
	#~ plt.tight_layout()
	
	if if_plot == 'y':
		#~ plt.show()
		plt.savefig(output_directory+"/"+"poisson_deviation_doubleevent_after_noise_clean.png")
		plt.clf()

	diff = 0
	#~ outfile.write("READ SERR 2\n")
	for i in range(0,len(bin_e)):
		outfile.write(str(bin_e[i])+"\t"+str(hist_total[i])+"\t"+str(np.sqrt(hist_total[i]))+"\t"+str(hist_exp_total[i])+"\n")
		if i==0:
			continue
		if float(hist_total[i])-float(hist_exp_total[i]) > 0:
			diff = diff + float(hist_total[i])-float(hist_exp_total[i])
	
		


	print ("Noise counts = "+str(diff))
	print ("exposure = "+str(exposure))
	print ("Noise countrate = "+ str(np.round(diff/exposure,2)))




if __name__ == '__main__':
    args = parser.parse_args()
    run_poisson_data_test(**args.__dict__)

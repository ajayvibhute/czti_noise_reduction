#Author : Ajay Ratheesh 

import os
import numpy as np
from astropy.io import fits,ascii
from matplotlib import pyplot as plt
import sys
import argparse

#How to run the code
#python --inputdir NAME_OF_INPUT_DIR --Eventfile NAME_OF_EVTFILE --bunchfile NAME_OF_BUNCHFILE --outdir NAME_OF_THE_OUTPUT_DIR



formatter = argparse.ArgumentDefaultsHelpFormatter
parser = argparse.ArgumentParser(formatter_class=formatter)
parser.add_argument("--inputdir",  help='Path to the input directory', 
                    type=str)
parser.add_argument("--Eventfile",  help='Name of the event file (eg : *level2_bc_nc.evt)', 
                    type=str)
parser.add_argument("--bunchfile",  help='Name of the bunchfile (eg : *level2_bunch.fits)', 
                    type=str)
parser.add_argument("--outdir",  help='Path to the output directory', 
                    type=str)                  


def Dphstructures(detx, dety): 
	
	thresh = 0.707
	allowable_thresh = 3.0
	DPH = np.histogram2d(detx,dety,bins = 64)[0]
	#~ plt.imshow(DPH)
	#~ plt.colorbar()
	#~ plt.legend()
	#~ plt.show()
	
	detx_nonzero, dety_nonzero = np.where(DPH>0)
	#~ print (DPH)
	#~ print (np.shape(DPH))
	#~ print (detx_nonzero)
	#~ print (dety_nonzero)
	DPH_nonzero = DPH[np.where(DPH>0)]
	#~ print (DPH_nonzero)
	flag_hot_pix = np.array([False]*len(detx_nonzero))
	hot_sum 	 = 0.0
	OUTCOME = False
	for i in range(0,len(detx_nonzero)-1):
		for j in range(i+1,len(dety_nonzero)):
			distance = np.sqrt(  ( float(detx_nonzero[i])- float(detx_nonzero[j]) )**2   +  ( float(dety_nonzero[i])- float(dety_nonzero[j]) )**2 )
			cij = DPH_nonzero[i]*DPH_nonzero[j]/distance
			
			neigh_sum = 0
			if distance<1.5 and DPH_nonzero[i]==1 and DPH_nonzero[i]==1:
				
				low_x = detx_nonzero[i]-1 if detx_nonzero[i]>0 else 0
				high_x = detx_nonzero[i]+1 if detx_nonzero[i]<63 else 63
				low_y = dety_nonzero[i]-1 if dety_nonzero[i]>0 else 0
				high_y = dety_nonzero[i]+1 if dety_nonzero[i]<63 else  63

				for k in [low_x, detx_nonzero[i], high_x]:
					for m in [low_y, dety_nonzero[i], high_y]:
						neigh_sum = neigh_sum + DPH[k][m]
				
				low_x = detx_nonzero[j]-1  if detx_nonzero[j]>0 else  0
				high_x = detx_nonzero[j]+1 if detx_nonzero[j]<63 else 63
				low_y = dety_nonzero[j]-1  if dety_nonzero[j]>0 else  0
				high_y = dety_nonzero[j]+1 if dety_nonzero[j]<63 else 63	
					
				for k in [low_x, detx_nonzero[j] , high_x]:
					for m in [low_y, dety_nonzero[j], high_y]:
						neigh_sum = neigh_sum + DPH[k][m]
						
				if neigh_sum >2:
					continue
			
			if cij>thresh:
				hot_sum = hot_sum + cij
				flag_hot_pix[i]=True
				flag_hot_pix[j]=True
			
	if float(np.sum(flag_hot_pix)) >0:
		allowable = hot_sum/float(np.sum(flag_hot_pix))
	else: 
		allowable = 0
	if allowable > allowable_thresh: 
		OUTCOME = True
		#~ DPH_NOISE = "FLAGGED"
		#~ plt.imshow(DPH)
		#~ plt.title(DPH_NOISE)
		#~ plt.colorbar()
		#~ plt.savefig("Noise_dph"+str(dph_counter)+".png")
	else:
		OUTCOME = False
		#~ DPH_NOISE = "NOT FALGGED"
	

	#~ plt.legend()
	#~ plt.show()
	
	
	
	return OUTCOME, allowable, float(np.sum(flag_hot_pix))

def read_data(path, file_name, Q_n):
	event_file = fits.open(path+file_name)
	event_list = event_file[Q_n].data
	return event_list
	
def runprobeDPHstructure(**kwargs): 
	
	input_dir 			= kwargs['inputdir']
	evt_file_name 		= kwargs['Eventfile']
	bunch_file_name 	= kwargs['bunchfile']
	output_directory 	= kwargs['outdir']
	
	outfilename = evt_file_name.split('.')[0]+ ".txt"
	outfile = open(output_directory+"/"+outfilename,"w")
	outfile.write("Quad_no"+"\t"+"Time_DPHstructure"+"\t"+"Allowable"+"\t"+"NoOfhotPixels"+"\n")
	for Q_n in range(1,5):
		print ("Quadrant ----"+str(Q_n))
		eventlist = read_data(input_dir, evt_file_name, Q_n)
		bunchlist = read_data(input_dir, bunch_file_name, Q_n)
		
		time_evt = eventlist['TIME']
		detx_evt = eventlist['DETX']
		dety_evt = eventlist['DETY']
		
		bin_edge = np.arange(time_evt[0], time_evt[len(time_evt)-1], 0.1 )
		dph_counter = 0
		for i in bin_edge: 
			index_curr = np.where(   (i < time_evt)     &  (i+0.1 > time_evt)   )
			detx_curr = detx_evt[index_curr]
			dety_curr = dety_evt[index_curr] 
			flag_noise, allowable, Npixels = Dphstructures(detx_curr, dety_curr)
			
			if flag_noise == True:
				outfile.write(str(Q_n-1)+"\t"+str(i+0.05)+"\t"+str(allowable)+"\t"+str(Npixels)+"\n")
				DPH = np.histogram2d(detx_curr,dety_curr,bins = 64)[0]
				dph_counter = dph_counter + 1 
				DPH_NOISE = "Q"+str(Q_n-1)+",DPHstructureNO=" + str(dph_counter) +",time ="+str(i+0.05)
				plt.imshow(DPH)
				plt.title(DPH_NOISE)
				plt.colorbar()
				plt.savefig("Noise_dph_Q"+str(Q_n-1)+"_num"+str(dph_counter)+".png")
				plt.clf()
		print ("The number of noise dph in Q"+str(Q_n-1)+" is "+str(dph_counter))
		
		

if __name__ == '__main__':
	

	args = parser.parse_args()
	runprobeDPHstructure(**args.__dict__)
	
	
	

	
	
	

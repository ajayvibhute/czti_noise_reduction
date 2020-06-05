from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys

evtorg=sys.argv[1]

org=fits.open(evtorg)

obsid=evtorg.split('_')[2]
orbitno=evtorg.split('_')[3].split('c')[0]

noisy_detx=int(sys.argv[2])
noisy_dety=int(sys.argv[3])
fig=plt.figure()
for qid in range(1, 5):
	pixdata=org[qid].data['Time'][np.where((org[qid].data['DETX']==noisy_detx) & (org[qid].data['DETY']==noisy_dety))]
	print pixdata
	oh,b1=np.histogram(pixdata,bins=(int)(pixdata[-1]-pixdata[0]))

	plt.plot(b1[:-1],oh)
	plt.show()
	#fig.savefig(obsid+"_"+orbitno+"_Q"+str(qid)+".png")
org.close()

#! /usr/bin/env python
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys

file1=sys.argv[1]
hdu=fits.open(file1)
tstart=float(sys.argv[2])
tstop=float(sys.argv[3])
print(file1+" "+str(tstart)+" "+str(tstop) )
fig=plt.figure()
for qid in range(1, 5):
	time=hdu[qid].data['Time']
	print time[0]
	time=time[(time > tstart) & (time < tstop)]
	print time[0]
	oh,b1=np.histogram(time,bins=(int)(time[-1]-time[0]))

	plt.plot(b1[:-1],oh)
	plt.show()
hdu.close()

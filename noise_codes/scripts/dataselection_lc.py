from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys

evtorg=sys.argv[1]
evtds1=sys.argv[2]
#evtds2=sys.argv[3]
#evtds3=sys.argv[4]
#evtds4=sys.argv[5]

org=fits.open(evtorg)
ds1=fits.open(evtds1)
#ds2=fits.open(evtds2)
#ds3=fits.open(evtds3)
#ds4=fits.open(evtds4)

qid=1#int(sys.argv[6])

obsid=evtorg.split('_')[2]
orbitno=evtorg.split('_')[3].split('c')[0]

#for energy
orgt = org[qid].data["ENERGY"]
ds1t = ds1[qid].data["ENERGY"]
#for time
#orgt = org[qid].data["Time"]
#ds1t = ds1[qid].data["Time"]
#ds2t = ds2[qid].data.field(0)
#ds3t = ds3[qid].data.field(0)
#ds4t = ds4[qid].data.field(0)

fig=plt.figure()
#histogram of energy
oh,b1=np.histogram(orgt,bins=90,range=(10,100))#bins=(int)(orgt[-1]-orgt[0]))
ds1h,b2=np.histogram(ds1t,bins=90,range=(10,100))#(int)(ds1t[-1]-ds1t[0]))
#histogram of time
#oh,b1=np.histogram(orgt,bins=(int)(orgt[-1]-orgt[0]))
#ds1h,b2=np.histogram(ds1t,bins=(int)(ds1t[-1]-ds1t[0]))
#ds2h,b3=np.histogram(ds2t,bins=(int)(ds2t[-1]-ds2t[0]))
#ds3h,b4=np.histogram(ds3t,bins=(int)(ds3t[-1]-ds3t[0]))
#ds4h,b5=np.histogram(ds4t,bins=(int)(ds4t[-1]-ds4t[0]))

plt.plot(b1[:-1],oh,color="red")
plt.plot(b2[:-1],ds1h,color="blue")
#plt.plot(b3[:-1],ds2h)
#plt.plot(b4[:-1],ds3h)
#plt.plot(b5[:-1],ds4h)

plt.show()
#fig.savefig(obsid+"_"+orbitno+"_Q"+str(qid)+".png")
org.close()

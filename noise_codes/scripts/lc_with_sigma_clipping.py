#! /usr/bin/env python
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys

evtfile=sys.argv[1]

evt=fits.open(evtfile)

qid=(int)(sys.argv[2])

#for time
evttime = evt[qid].data["pixID"]
evttimearr = np.array(evttime)
mean = np.mean(evttimearr, axis=0)
sd = np.std(evttimearr, axis=0)
print "Mean : "+str(mean)+" SD : "+str(sd);
final_evttime = [time for time in evttime if (time > mean - 3 * sd)]
final_evttime = [time for time in final_evttime if (time < mean + 3 * sd)]
#print(final_evttime)

#histogram of time
oh,b1=np.histogram(final_evttime,bins=(int)(final_evttime[-1]-final_evttime[0]))
plt.plot(b1[:-1],oh)

plt.show()
evt.close()

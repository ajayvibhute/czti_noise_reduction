#! /usr/bin/env python2.7
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys
#newfile="/data2/czti/testarea/mayuri/20170927_A03_086T01_9000001568_level2/0.1sec/20170927_A03_086T01_9000001568_level2/czti/modeM0/AS1A03_086T01_9000001568cztM0_level2_quad_clean_1.0_Q0.lc"
#oldfile="/data2/czti/testarea/mayuri/20170927_A03_086T01_9000001568_level2/pipeline_1.0_result/lv_1.0_lc_1.0/20170927_A03_086T01_9000001568_level2/czti/modeM0/AS1A03_086T01_9000001568cztM0_level2_quad_clean_Q0.lc"
newfile=sys.argv[1]
oldfile=sys.argv[2]
newhdu=fits.open(newfile)
oldhdu=fits.open(oldfile)
newrate=newhdu[1].data["RATE"]
newfracexp=newhdu[1].data["FRACEXP"]
oldrate=oldhdu[1].data["RATE"]
oldfracexp=oldhdu[1].data["FRACEXP"]

oldfracexp=oldfracexp[:-1]
oldrate=oldrate[:-1]


plt.plot(newrate,color="red")
plt.plot(oldrate,color="blue")
plt.show()

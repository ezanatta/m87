# -*- coding: utf-8 -*-
"""
Created on Fri Sep 14 14:12:06 2018

@author: emilio
"""

import matplotlib.pyplot as plt
import numpy as np

hrv = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/strader_acs_old_hawki.cat'
acs = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/old_acs_hawki.cat'

RAacs = np.loadtxt(acs, usecols=(0,))
DECacs = np.loadtxt(acs, usecols=(1,))

RAhrv = np.loadtxt(hrv, usecols=(0,))
DEChrv = np.loadtxt(hrv, usecols=(1,))

v = np.loadtxt(hrv, usecols=(8,))
err_v = np.loadtxt(hrv, usecols=(9,))

plt.plot(RAacs, DECacs, marker='+', color='black', linestyle='none', ms=4, alpha=0.2, label='ACS+Oldham+HAWKI')
plt.scatter(RAhrv, DEChrv, c=v, vmin=1000, vmax=1800, cmap='rainbow', s=40, label='Strader+2011 matched GCs')
cb = plt.colorbar(fraction=0.05)
cb.set_label('Velocity $(km/s)$')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.legend(loc='lower left', numpoints=1)
plt.show()
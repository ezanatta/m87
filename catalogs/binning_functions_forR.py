def binning(RA, DEC, RAgal, DECgal, d, nbins, pa, incl):
  import numpy as np
  from astropy import units as u
  from astropy.coordinates import SkyCoord
  
  pa = pa*u.deg
  incl = incl*u.deg
  
  galcenter = SkyCoord(RAgal, DECgal)
  
  #print("%.4f" %RA[1])
  
  RA = RA*u.deg
  DEC = DEC*u.deg
  n = len(RA)
  
  gcs = list()
  dist = list()
  d = d*u.Mpc
  
  pa_rad = pa.to(u.rad)
  pa_rad = pa_rad.value

  i_rad = incl.to(u.rad)
  i_rad = i_rad.value

  n = len(RA)

  ra = np.linspace(0,0,n)
  dec = np.linspace(0,0,n)

  for i in range(0,n):
	  aux = SkyCoord(RA[i], DEC[i])             #reading RA and DEC from GC
	  ra[i] = aux.ra.arcsec
	  dec[i] = aux.dec.arcsec
       
   
  rad = np.zeros(n)  
  for i in range(0, n):
         gc = SkyCoord(RA[i], DEC[i])
         rad[i] = SkyCoord.separation(gc, galcenter).arcsec
	
  y_rad = galcenter.dec.radian   # galaxy center DEC coordinate in radians
    #print 'Y RAD', y_rad
        
  xmc = (ra - galcenter.ra.arcsec)
  xm = xmc*(np.cos(y_rad))                       # trasnform the coordinates from spherical to planar   
  ym = (dec - galcenter.dec.arcsec)             # setting the y coordinates of the GC from its DEC value relative to galactic center
	
  xs = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)            # 	rotating coordinates 
  ys = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)    
	
  cos_i = np.cos(i_rad)               # cosine of inclination angle, used to shrink coordinates

  xd = ((xs)*np.sqrt(cos_i))
  yd = ((ys)/np.sqrt(cos_i))
     
  distance=np.sqrt((xs**2)*cos_i+((ys**2)/cos_i))
  
  ###binning ####
  ###
  r = sorted(rad)
  n = len(r)   #lenght of r --- number of objects to iterate on the following loops
  rgal = list()      #rgal contains the limits of each bin
  
  nbins = int(nbins)
  
  for h in range(1, nbins):
    rgal.append(r[int(n*h/nbins)])     #binning
  rgal.insert(nbins, r[n-1])
  rgal.insert(0, r[0])
  
  with open('temp_radius_Rpy', 'w') as f:
    for i in range(0, len(r)):
      print >>f, rad[i], RA[i].value, DEC[i].value
      
  with open('temp_bins_Rpy', 'w') as f:
    for item in rgal:
      print >>f, item
      
#  with open('temp_radius_Rpy', 'w') as f:
#    for item in x1l:
#      print >>f, item
#  with open('temp_radius_Rpy', 'w') as f:
#    for item in y1l:
#      print >>f, item
#  with open('temp_radius_Rpy', 'w') as f:
#    for item in z1l:
#      print >>f, item  
      
#def density(nbin, r, rgal):
#    import numpy as np
#    
#    n = len(r)
##    NGC = np.zeros(nbin)    #number of GC in each bin
#    rmed = [[] for i in range(0,nbin)]
#    for h in range(0, nbin):                              
#         for i in range(0, n):
#             if (rgal[h+1] >= r[i] and r[i] > rgal[h]):
#                 NGC[h] = NGC[h] + 1
#                 rmed[h].append(r[i])
#
#    area = []
#    for h in range(0, nbin):
#        area.append((np.pi*rgal[h+1]**2)-(np.pi*rgal[h]**2)) #area per bin
#    
#    binsize = np.linspace(0,0,nbin)    
#    for h in range(0,nbin):
#        binsize[h] = rgal[h+1]-rgal[h]  
#    binsize = binsize/2
#    
#    return area, NGC, median, binsize, poi_err 
    
def binning2(RA, DEC, RAgal, DECgal, d, nbins, pa, incl):
  import numpy as np
  from astropy import units as u
  from astropy.coordinates import SkyCoord
  
  pa = pa*u.deg
  incl = incl*u.deg
  
  galcenter = SkyCoord(RAgal, DECgal)
  
  #print("%.4f" %RA[1])
  
  RA = RA*u.deg
  DEC = DEC*u.deg
  n = len(RA)
  
  gcs = list()
  dist = list()
  d = d*u.Mpc
  
  pa_rad = pa.to(u.rad)
  pa_rad = pa_rad.value

  i_rad = incl.to(u.rad)
  i_rad = i_rad.value

  n = len(RA)

  ra = np.linspace(0,0,n)
  dec = np.linspace(0,0,n)

  for i in range(0,n):
	  aux = SkyCoord(RA[i], DEC[i])             #reading RA and DEC from GC
	  ra[i] = aux.ra.arcsec
	  dec[i] = aux.dec.arcsec
       
   
  rad = np.zeros(n)  
  for i in range(0, n):
         gc = SkyCoord(RA[i], DEC[i])
         rad[i] = SkyCoord.separation(gc, galcenter).arcsec
	
  y_rad = galcenter.dec.radian   # galaxy center DEC coordinate in radians
    #print 'Y RAD', y_rad
        
  xmc = (ra - galcenter.ra.arcsec)
  xm = xmc*(np.cos(y_rad))                       # trasnform the coordinates from spherical to planar   
  ym = (dec - galcenter.dec.arcsec)             # setting the y coordinates of the GC from its DEC value relative to galactic center
	
  xs = xm*np.cos(pa_rad)-ym*np.sin(pa_rad)            # 	rotating coordinates 
  ys = xm*np.sin(pa_rad)+ym*np.cos(pa_rad)    
	
  cos_i = np.cos(i_rad)               # cosine of inclination angle, used to shrink coordinates

  xd = ((xs)*np.sqrt(cos_i))
  yd = ((ys)/np.sqrt(cos_i))
     
  distance=np.sqrt((xs**2)*cos_i+((ys**2)/cos_i))
  
  ###binning ####
  ###
  r = sorted(rad)
  n = len(r)   #lenght of r --- number of objects to iterate on the following loops
  rgal = list()      #rgal contains the limits of each bin
  
  nbins = int(nbins)
  
  for h in range(1, nbins):
    rgal.append(r[int(n*h/nbins)])     #binning
  rgal.insert(nbins, r[n-1])
  rgal.insert(0, r[0])
  
  result = [rad, distance, rgal]
  
  return result

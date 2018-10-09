require(mixtools)
require(ggplot2)
require(mclust)
require(fpc)
require(factoextra)
source('/home/emilio/HAWKI/products-inter-combine/functions_to_bimodal.r')

sspselect <- function(n){
  
  col1 = d[[n]]$V7 - d[[n]]$V10
  col2 = d[[n]]$V10 - (e[[n]]$V14+1.827)
  
  return(data.frame(col1, col2))
}

corr_cat = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/masked_old_acs_hawki.cat')

g = corr_cat$V3
gerr = corr_cat$V4
z = corr_cat$V5
zerr = corr_cat$V6
k = corr_cat$V7
kerr = corr_cat$V8

####### redenning
### must be added both for oldham and hawki data
# u, g, r, i, z, k
#av = c(0.097,0.076,0.053,0.039,0.029,0.007)

#removing outliers in the k band:

#bad_index = c()
#bad_index = which(kerr >= 0.2)

cols = data.frame(g, z, k, gerr, zerr, kerr)

#cols = remove_out(df_cols, bad_index)

av = c(0.097,0.076,0.053,0.039,0.029,0.007)

g = cols$g-av[2]
gerr = cols$gerr
z = cols$z-av[5]
zerr = cols$zerr
k = cols$k-av[6]
kerr = cols$kerr

#converting from VEGA to AB system:

k = k+1.827


#### colors

zk = z-k
gk = g-k
gz = g-z


#zkpure = zk[-mw]
#gkpure = gk[-mw]
#gzpure = gz[-mw]

#kpure = k[-mw]
#kerrpure = kerr[-mw]

########### error in colors

errzk = sqrt(zerr^2+kerr^2)
errgk = sqrt(gerr^2+kerr^2)
errgz = sqrt(gerr^2+zerr^2)

#errzkpure = errzk[-mw]
#errgkpure = errgk[-mw]
#errgzpure = errgz[-mw]

######### bimodality tests

estgz = normalmixEM(gz)
estzk = normalmixEM(zk)
estgk = normalmixEM(gk)

plot(estzk, which=2)
plot(estgk, which=2)
plot(estgz, which=2)

##### using ggplot to plot bimodality tests

gg.mixEM(estzk, '(z-k)')
gg.mixEM(estgk, '(g-k)')
gg.mixEM(estgz, '(g-z)')

#dat = data.frame(zk, gk, uk, gz)
#datpure = data.frame(zkpure, gkpure, ukpure, gzpure)
#ggplot(datpure, aes(x=zkpure)) + geom_histogram(aes(y=..density..), fill='darkgrey')+geom_density()+labs(title='M87 GC (z-k) color distribution', x='(z-k)', y='Density')

##### histograms according to k-band error 

Kerr_hist(gz, kerr, 0.03, 0.06, '(g-z)')
Kerr_hist(zk, kerr, 0.03, 0.06, '(z-k)')
Kerr_hist(gk, kerr, 0.03, 0.06, '(g-k)')

plot(gk, g, ylim=c(27,19),xlab='(g-k)', ylab='g ', pch=5)
arrows(gk, g-gerr, gk, g+gerr, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(gk-errgk, g, gk+errgk, g, length=0.05, angle=90, code=3, lwd= 0.3)

df = data.frame(g, z, k, gerr, zerr, kerr, gz, zk, gk, errgz, errzk, errgk)

ggplot(df, aes(x=gz, y=zk))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha=0.5)+
  scale_fill_distiller(palette='YlOrRd', direction=1, name='Density') +
  xlim(c(-1,5))+ylim(c(-4, 4))+
  #geom_point(size=0.5, alpha=0.2)+
  #geom_errorbar(alpha=0.5, size= 0.2,width=.1, aes(ymin=zk-errzk, ymax=zk+errzk)) +
  #geom_errorbarh(alpha=0.5, size= 0.2,height = .1, aes(xmin = gz-errgz,xmax = gz+errgz))+
  geom_density2d(color='black')+
  geom_abline(intercept = -0.2179487179, slope=1.2820512821, col='purple')+
  geom_hline(yintercept = 2.5, col='purple', linetype='dashed')+
  geom_abline(intercept = -0.7, slope=0.4, col='purple', linetype='dotted')+
  theme_bw()+labs(title=' ', x='(g-z)', y='(z-k)')

db <- fpc::dbscan(data.frame(gz, zk), eps = 0.07, MinPts = 4)
a <- fviz_cluster(db, data = data.frame(gz, zk), stand = FALSE,
             ellipse = FALSE, show.clust.cent = FALSE, outlier.color='red',xlim=c(-1, 5),ylim=c(-4.5, 4.5),
             geom = "point", palette = "jco", ggtheme = theme_classic())
plot(a)

#db_gmm(db, df, 1) # only outliers
#db_gmm(db, df, 0) # principal component of clustering

dfpure = iterdbscan(gz, zk)

ggplot(dfpure, aes(x=x, y=y))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(c(-1,5))+ylim(c(-4, 4))+
  scale_fill_distiller(palette='YlOrRd', direction=1, name='Density') +
  #geom_point(size=0.5, alpha=0.5)+
  geom_density2d(color='black')+
  geom_abline(intercept = -0.2179487179, slope=1.2820512821, col='purple')+
  geom_hline(yintercept = 2.5, col='purple', linetype='dashed')+
  geom_abline(intercept = -0.7, slope=0.4, col='purple', linetype='dotted')+
  theme_bw()+labs(title=' ', x='(g-z)', y='(z-k)')

powalka = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/2psf_sub_powalka.cat'
ngvs = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/NGVS.cat'

cols = readcat('uik', powalka)

up = cols$u-av[1]
uperr = cols$uerr
gp = cols$g-av[2]
gperr = cols$gerr
rp = cols$r-av[3]
rperr = cols$rerr
ip_band = cols$i_band-av[4]
ip_band_err = cols$i_band_err
zp = cols$z-av[5]
zperr = cols$zerr
kp = cols$k-av[6]
kperr = cols$kerr
kp_hawki = cols$k_hawki-av[6]
kperr_hawki = cols$kerr_hawki

#converting HAWKI data from VEGA to AB system:

kp_hawki = kp_hawki+1.827

#### colors

zkp = zp-kp
zkhp = zp-kp_hawki
gkp = gp-kp
gkhp = gp-kp_hawki
ukp = up-kp
ukhp = up-kp_hawki
ikp = ip_band - kp
ikhp = ip_band-kp_hawki
gzp = gp-zp
uip = up-ip_band

errzkp = sqrt(zperr^2+kperr^2)
errzkhp = sqrt(zperr^2+kperr_hawki^2)
errgkp = sqrt(gperr^2+kperr^2)
errgkhp = sqrt(gperr^2+kperr_hawki^2)
errukp = sqrt(uperr^2+kperr^2)
errukhp = sqrt(uperr^2+kperr_hawki^2)
errikp = sqrt(ip_band_err^2 + kperr^2)
errikhp = sqrt(ip_band_err^2 + kperr_hawki^2)
errgzp = sqrt(gperr^2+zperr^2)
erruip = sqrt(uperr^2+ip_band_err^2)


##### acs_old_hawki ra and dec:
RA = corr_cat$V1
DEC = corr_cat$V2
RAp = read.table(ngvs)$V2
DECp = read.table(ngvs)$V3

index = which(gz==dfpure$x)

radec = data.frame(RA, DEC)
radecp = data.frame(RAp, DECp)

ggplot(radec, aes(x=RA, y=DEC))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", bins=10)+
  xlim(c(187.8,187.6))+
  scale_fill_distiller(palette='YlOrRd', direction=1, name='Density') +
  geom_point(size=0.5, alpha=0.5)+
  #geom_density2d(color='black')+
  theme_bw()+labs(title=' ', x='RA', y='DEC')

ggplot(radecp, aes(x=RAp, y=DECp))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(c(187.85,187.56))+ylim(c(12.5, 12.25))+
  scale_fill_distiller(palette='YlOrRd', direction=1, name='Density') +
  geom_point(data=radec, aes(x=RA, y=DEC), col='grey')+
  geom_point(size=0.5)+
  geom_density2d(color='black')+
  theme_bw()+labs(title=' ', x='RA', y='DEC')

#### MILES + Padova+00 SSPs:
#### 

ssp = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/MILES/sdss_ch_iPp0.00.MAG')
ssp2 = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/MILES/UBBVRIJHK_ch_iPp0.00.MAG')
ssp = data.frame(ssp)
ssp2 = data.frame(ssp2)
d = split(ssp, rep(1:7, each=50))
e = split(ssp2, rep(1:7, each=50))

plot(gz, zk, xlim=c(0, 2), ylim=c(-2, 1), pch=19, col='grey')
points(dfpure$x, dfpure$y, pch=19)
points(gzp, zkp, col='orange', pch=2)
points(sspselect(1)$col1, sspselect(1)$col2, col='red', pch=19)
points(sspselect(2)$col1, sspselect(2)$col2, col='green', pch=19)
points(sspselect(3)$col1, sspselect(3)$col2, col='blue', pch=19)
points(sspselect(4)$col1, sspselect(4)$col2, col='yellow', pch=19)
points(sspselect(5)$col1, sspselect(5)$col2, col='purple', pch=19)
points(sspselect(6)$col1, sspselect(6)$col2, col='cyan', pch=19)
points(sspselect(7)$col1, sspselect(7)$col2, col='magenta', pch=19)

old_oldham = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/masked_sub_oldham.cat')

cols = readcat('uik', oldham)

av = c(0.097,0.076,0.053,0.039,0.029,0.007)

uo = cols$u-av[1]
uoerr = cols$uerr
go = cols$g-av[2]
goerr = cols$gerr
ro = cols$r-av[3]
roerr = cols$rerr
io_band = cols$i_band-av[4]
io_band_err = cols$i_band_err
zo = cols$z-av[5]
zoerr = cols$zerr
ko = cols$k-av[6]
koerr = cols$kerr

#converting from VEGA to AB system:

#ko = ko+1.827


#### colors

zko = zo-ko
gko = go-ko
uko = uo-ko       #with outliers
iko = io_band - k 
gzo = go-zo
uio = uo-io_band

########### error in colors

errzko = sqrt(zoerr^2+koerr^2)
errgko = sqrt(goerr^2+koerr^2)
erruko = sqrt(uoerr^2+koerr^2)
erriko = sqrt(io_band_err^2 + koerr^2)
errgzo = sqrt(goerr^2+zoerr^2)
erruio = sqrt(uoerr^2+io_band_err^2)

dfpure_oldham = iterdbscan(uio, iko)

ggplot(dfpure_oldham, aes(x=x, y=y))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(c(-1,5))+ylim(c(-4, 4))+
  scale_fill_distiller(palette='YlOrRd', direction=1, name='Density') +
  geom_point(size=0.5)+
  geom_point(data=data.frame(uip, ikp), aes(x=uip, y=ikp), col='green', alpha=0.2)+
  geom_density2d(color='black')+
  theme_bw()+labs(title=' ', x='(u-i)', y='(i-k)')

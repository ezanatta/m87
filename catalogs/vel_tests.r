source('/home/emilio/HAWKI/products-inter-combine/functions_to_bimodal.r')

library(mixtools)
library(stats4)
library(ggplot2)
library(deldir)
library(latticeExtra)
library(graphics)
library(ggvoronoi)
library(astro)

hrv = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/strader_acs_old_hawki.cat')
RA = hrv$V1
DEC = hrv$V2
g = hrv$V3
gerr = hrv$V4
z = hrv$V5
zerr = hrv$V6
k = hrv$V7
kerr = hrv$V8
v = hrv$V9
verr = hrv$V10
class=hrv$V11

for(i in 1:length(class)){
  if(class[i]==1){
    class[i]=21
  }
  if(class[i]==2){
    class[i]=24
  }
  if(class[i]==3){
    class[i]=22
  }
}

corr_cat = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/masked_sub_oldham.cat')

RAp = corr_cat$V1
DECp = corr_cat$V2
up = corr_cat$V3
uperr = corr_cat$V4
gp = corr_cat$V5
gperr = corr_cat$V6
rp = corr_cat$V7
rperr = corr_cat$V8
ip_band = corr_cat$V9
ip_band_err = corr_cat$V10
zp = corr_cat$V11
zperr = corr_cat$V12
kp = corr_cat$V13
kperr = corr_cat$V14
Vp = corr_cat$V23
Vperr = corr_cat$V24

photdb = data.frame(RAp, DECp, up, uperr, gp, gperr, rp, rperr, ip_band, ip_band_err, zp, zperr, kp, kperr, Vp, Vperr)

ggplot(data=data.frame(RA, DEC, v, levels=class))+
  #geom_voronoi(aes(x=RA,y=DEC, fill=v))+
  #geom_point(data=photdb, aes(x=photdb$RAp, y=photdb$DECp), pch=1)+
  geom_point(aes(x=RA, y=DEC, colour=v), pch=class, size=4)+
  scale_colour_gradient2(low="blue", midpoint=1350, high="red")+
  scale_fill_gradient2(low="blue", midpoint=1350, high="red")+
  #geom_density2d(aes(x=RA, y=DEC))+
  #geom_density2d(data=photdb, aes(x=photdb$RAp, y=photdb$DECp))+
  theme_bw()

plot1 = ggplot(data = data.frame(RA, DEC, v), aes(x=RA, y=DEC))+
  #geom_tile()+
  geom_point(data=photdb, aes(x=photdb$RAp, y=photdb$DECp), inherit.aes = FALSE)+
  geom_point(aes(x=RA, y=DEC, col=v), size=3)+
  scale_colour_gradient2(low='royalblue4',midpoint=1330, high='red4')+
  geom_density2d(aes(x=RA, y=DEC))+
  theme_bw()

print(plot1)

vsys = 1284.011278  #M87 systemic velocity

model = normalmixEM(v-vsys)
#gg.mclust(model, 'Velocities', 'V (km/s)')
gg.mixEM(model, 'Velocity', 'V (km/s)')

LL = function(mu, sigma){
  R = dnorm(v, mu, sigma)
  return(-sum(log(R)))
}

binormal = function(mu, sigma, sigma2){
  R = dnorm(v, mu, sigma)+dnorm(v, 0, sigma2)
  return(R)
}

model = mle(minuslogl = LL, start=list(mu=mean(v), sigma=sd(v)))

like = binormal(mu=as.numeric(coef(model)[1]), sigma=as.numeric(coef(model)[2]), sigma2=as.numeric(coef(model)[3]))

hist(like)

so = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/strader_oldham.cat')
so = data.frame(so)

newflag = c()

for(i in 1:length(so$V1)){
  if(so$V17[i]=='trans' | so$V17[i]=='ucd'){
    print(so$V17[i])
    newflag = append(newflag, 'UCD_cand')
  }
  else{
    newflag = append(newflag, 'GC')
  }
}

so$V18 = newflag

ucd = split(so, factor(so$V18))
gk_ucd = ucd$UCD_cand$V5-ucd$UCD_cand$V13

model_ucd = normalmixEM(gk_ucd)
gg.mixEM(model_ucd, 'ucd', '(g-k)')

rh_old = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/rh_old_acs.cat', header=TRUE)
rh_old=data.frame(rh_old)

obj_r = angsize(z=0.004283, r=rh_old$Rh, out = 'pc')
hist(obj_r, xlab='Rh (kpc)')
plot(rh_old$gmag-rh_old$kmag, obj_r, xlim=c(0, 5.0), ylim=c(0, 8), xlab='(g-k)', ylab='Rh (pc)')
abline(h=5)

ggplot(data = data.frame(rh_old$gmag-rh_old$kmag, obj_r), aes(x=rh_old$gmag-rh_old$kmag, y=obj_r))+
  geom_point()+geom_density2d()+theme_bw()+geom_abline(slope=0, intercept=5)+labs(title=' ', x='(g-k)', y='Rh (pc)')

ggplot(data = data.frame(rh_old$gmag-rh_old$zmag, obj_r), aes(x=rh_old$gmag-rh_old$zmag, y=obj_r))+
  geom_point()+geom_density2d()+theme_bw()+xlim(c(0, 2))+geom_abline(slope=0, intercept=5)+labs(title=' ', x='(g-z)', y='Rh (pc)')

ggplot(data = data.frame(rh_old$gmag-rh_old$imag, obj_r), aes(x=rh_old$gmag-rh_old$imag, y=obj_r))+
  geom_point()+geom_density2d()+theme_bw()+xlim(c(0, 2))+geom_abline(slope=0, intercept=5)+labs(title=' ', x='(g-i)', y='Rh (pc)')


library(reticulate)
source_python('/home/emilio/HAWKI/products-inter-combine/binning_functions_forR.py')
pa = 153.0
incl = 31.0
RAgal = '12h30m49.4s'
DECgal = '+12d23m28s'
d = 22.2
nbins=4

results = binning2(RA, DEC, RAgal, DECgal, d, nbins, pa, incl)

vbin = data.frame(radius = results[[1]], v=v)

velbins = replicate(nbins, c())
radbins = replicate(nbins, c())
rbins=results[[3]]

for(i in 1:nbins){
  for(j in 1:length(vbin$radius)){
    if((rbins[i+1]>=vbin$radius[j])&(vbin$radius[j]>rbins[i])){
      velbins[[i]] = append(velbins[[i]], j)
      radbins[[i]] = append(radbins[[i]], j)
    }
  }
}

rbins = angsize(z=0.004283, r=unlist(rbins), out = 'kpc')

mean_bin = c()

for(i in 1:nbins){
  mean_bin = append(mean_bin, mean(seq(rbins[[i]], rbins[[i+1]], length.out = 50)))
}

rad_err = c()

for(i in 1:nbins){
  rad_err = append(rad_err, rbins[[i+1]]-rbins[[i]])
}

vrms = c()
srms = c()
vrms_err = c()
srms_err = c()
vrms_err_min = c()
vrms_err_max = c()
srms_err_min = c()
srms_err_max = c()
vhe = v-vsys

for(i in 1:nbins){
  LL = function(mu, sigma){
    R = dnorm(vhe[velbins[[i]]], mu, sigma)
    return(-sum(log(R)))
  }
  
  model = mle(minuslogl = LL, start=list(mu=0, sigma=sd(vhe[velbins[[i]]])))
  #model = mle(minuslogl = LL)
  
  vrms = append(vrms, coef(model)[1])
  srms = append(srms, coef(model)[2])
  vrms_err = append(vrms, sqrt(vcov(model)[1]))
  srms_err = append(srms, sqrt(vcov(model)[4]))
  vrms_err_min = append(vrms_err_min, confint(model, level=0.9)[1])
  vrms_err_max = append(vrms_err_max, confint(model, level=0.9)[2])
  srms_err_min = append(srms_err_min, confint(model, level=0.9)[3])
  srms_err_max = append(srms_err_max, confint(model, level=0.9)[4])
}

plot(mean_bin, vrms, pch=19, xlab='R (kpc)', ylab='Vrot (km/s)', ylim=c(0, 350), xlim=c(0, max(mean_bin+rad_err/2)))
arrows(mean_bin, vrms - vrms_err[-length(vrms_err)], mean_bin, vrms + vrms_err[-length(vrms_err)], length=0.05, angle=90, code=3)
arrows(mean_bin - rad_err/2, vrms, mean_bin + rad_err/2, vrms, length=0.05, angle=90, code=3)

plot(mean_bin, srms, pch=19, xlab='R (kpc)', ylab=expression(paste(sigma, ' (km/s)')), ylim=c(0, 500), xlim=c(0, max(mean_bin+rad_err/2)))
arrows(mean_bin, srms - srms_err[-length(vrms_err)], mean_bin, srms_err+srms[-length(vrms_err)], length=0.05, angle=90, code=3)
arrows(mean_bin - rad_err/2, srms, mean_bin + rad_err/2, srms, length=0.05, angle=90, code=3)


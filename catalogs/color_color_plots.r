source('/home/emilio/HAWKI/products-inter-combine/functions_to_bimodal.r')

oldham = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/masked_sub_oldham.cat'
powalka = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/masked_sub_powalka.cat'
liris1 = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/2psf_sub_liris.cat'
liris2 = '/home/emilio/HAWKI/products-inter-combine/liris_oldham.dat'
acs = '/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/2psf_sub_acs.cat'

###### main catalog

print('Starting Oldham...')

cols = readcat('uik', oldham)

av = c(0.097,0.076,0.053,0.039,0.029,0.007)

u = cols$u-av[1]
uerr = cols$uerr
g = cols$g-av[2]
gerr = cols$gerr
r = cols$r-av[3]
rerr = cols$rerr
i_band = cols$i_band-av[4]
i_band_err = cols$i_band_err
z = cols$z-av[5]
zerr = cols$zerr
k = cols$k-av[6]
kerr = cols$kerr

#converting from VEGA to AB system:

k = k+1.827


#### colors

zk = z-k
gk = g-k
uk = u-k       #with outliers
ik = i_band - k 
gz = g-z
ui = u-i_band

########### error in colors

errzk = sqrt(zerr^2+kerr^2)
errgk = sqrt(gerr^2+kerr^2)
erruk = sqrt(uerr^2+kerr^2)
errik = sqrt(i_band_err^2 + kerr^2)
errgz = sqrt(gerr^2+zerr^2)
errui = sqrt(uerr^2+i_band_err^2)


############ plots

color_color(zk, gk, errzk, errgk, '(z-k)', '(g-k)', main='HAWKI_OLDHAM')
color_color(ui, ik, errui, errik, '(u-i)', '(i-k)', xlim=c(-1,7), ylim=c(-2,4), main='HAWKI_OLDHAM')

ggplot(data.frame(ui, ik), aes(x=ui, y=ik))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon")+
  xlim(c(-1,7))+ylim(c(-4, 4))+
  scale_fill_distiller(palette='YlOrRd', direction=1, name='Density') +
  geom_point(size=1, alpha=1)+
  geom_density2d(color='white')+
  theme_bw()+labs(title=' ', x='(u-i)', y='(i-k)')

color_color(gz, zk, errgz, errzk, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), main='HAWKI_OLDHAM')
abline(-0.2179487179,1.2820512821, col='purple')
abline(h=2.5, col='purple', lty=2)
abline(-0.7,0.4, col='purple', lty='dotted')

plot(gk, g, ylim=c(27,19),xlab='(g-k)', ylab='g ', pch=5, main='HAWKI_OLDHAM')
arrows(gk, g-gerr, gk, g+gerr, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(gk-errgk, g, gk+errgk, g, length=0.05, angle=90, code=3, lwd= 0.3)

op = readline('Deseja tentar clustering? y/n ')
cond = 0
while(cond==0){
  if(op == 'y'){
    op2 = readline('uik or gzk? 1/2 ')
    require(dbscan)
    
    df = data.frame(g, z, k,  gerr, zerr, kerr, gz, zk, gk, ui, ik,  errgz, errzk, errgk, errui, errik)
    
    if(op2==1){dfaux = data.frame(ui, ik)}else{dfaux = data.frame(gz, zk)}
    
    kNNdistplot(dfaux, k=4)
    abline(h=0.1, lty=2)
    abline(h=0.2, lty=2)
    abline(h=0.3, lty=2)
    abline(h=0.4, lty=2)
    
    eps = readline('Entre com o valor mais adequado de eps: ')
    
    db <- fpc::dbscan(dfaux, eps = eps, MinPts = 4)
    fviz_cluster(db, data = dfaux, stand = FALSE,
                 ellipse = FALSE, show.clust.cent = FALSE, outlier.color='red',xlim=c(-1, 5),ylim=c(-4.5, 4.5),
                 geom = "point", palette = "jco", ggtheme = theme_classic())
    
    
    cols2 = readcat('uik', powalka)
    
    u2 = cols2$u-av[1]
    u2err = cols2$uerr
    i2_band = cols2$i_band-av[4]
    i2_band_err = cols2$i_band_err
    k2 = cols2$k-av[6]
    k2err = cols2$kerr
    ui2 = u2-i2_band
    ik2 = i2_band - k2
    errui2 = sqrt(u2err^2+i2_band_err^2)
    errik2 = sqrt(i2_band_err^2 + k2err^2)
    color_color(ui2, ik2, 0, 0, '(u-i)', '(i-k)', xlim=c(-1,7), ylim=c(-2,4), m_cor='blue', p_cor='red', p_ch=1, alpha=0.2, axes=FALSE)
    par(new=TRUE)
    db_gmm(db, df, 0, op2)
    
    #db_gmm(db, df, seq(1:20), op2)
  }
  
  cond = readline('Enter 0 to redo clustering or 1 to keep going: ')
}

######## additional data
######## 
######## 

print('Starting Powalka...')
cols = readcat('uik', powalka)

u = cols$u-av[1]
uerr = cols$uerr
g = cols$g-av[2]
gerr = cols$gerr
r = cols$r-av[3]
rerr = cols$rerr
i_band = cols$i_band-av[4]
i_band_err = cols$i_band_err
z = cols$z-av[5]
zerr = cols$zerr
k = cols$k-av[6]
kerr = cols$kerr
k_hawki = cols$k_hawki-av[6]
kerr_hawki = cols$kerr_hawki

#converting HAWKI data from VEGA to AB system:

k_hawki = k_hawki+1.82

#### colors

zk = z-k
zkh = z-k_hawki
gk = g-k
gkh = g-k_hawki
uk = u-k
ukh = u-k_hawki
ik = i_band - k
ikh = i_band-k_hawki
gz = g-z
ui = u-i_band

errzk = sqrt(zerr^2+kerr^2)
errzkh = sqrt(zerr^2+kerr_hawki^2)
errgk = sqrt(gerr^2+kerr^2)
errgkh = sqrt(gerr^2+kerr_hawki^2)
erruk = sqrt(uerr^2+kerr^2)
errukh = sqrt(uerr^2+kerr_hawki^2)
errik = sqrt(i_band_err^2 + kerr^2)
errikh = sqrt(i_band_err^2 + kerr_hawki^2)
errgz = sqrt(gerr^2+zerr^2)
errui = sqrt(uerr^2+i_band_err^2)

############ plots again

color_color(zk, gk, errzk, errgk, '(z-k)', '(g-k)', main='HAWKI_POWALKA')
par(new=TRUE)
color_color(zkh, gkh, errzkh, errgkh, '(z-k)', '(g-k)', m_cor='blue', p_cor='red', p_ch=1, axes=FALSE)
leg.txt = c('Powalka','Hawki')
legend('bottomright',leg.txt, pch=c(19, 1), col=c('black', 'red'))

color_color(ui, ik, errui, errik, '(u-i)', '(i-k)', xlim=c(-1,7), ylim=c(-2,4), main='HAWKI_POWALKA')
par(new=TRUE)
color_color(ui, ikh, errui, errikh, '(u-i)', '(i-k)', xlim=c(-1,7), ylim=c(-2,4), m_cor='blue', p_cor='red',p_ch=1, axes=FALSE)
color_color(gz, zk, errgz, errzk, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), main='HAWKI_POWALKA')
par(new=TRUE)
color_color(gz, zkh, errgz, errzkh, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), m_cor='blue', p_cor='red', p_ch=1, axes=FALSE)
abline(-0.2179487179,1.2820512821, col='purple')
abline(h=2.5, col='purple', lty=2)
abline(-0.7,0.4, col='purple', lty='dotted')

ggplot(data.frame(gz, zk), aes(x=gz, y=zk))+
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

ploterr(k, k_hawki, kerr, kerr_hawki, xlab='K_POWALKA', ylab='K_HAWKI', pch=19)
abline(0, 1, col='blue', lwd=5)
mk = mean(kerr)
mk_hawki = mean(kerr_hawki)
leg.txt = c('x=y', sprintf('<err_powalka>: %f',mk), sprintf('<err_hawki>: %f ',mk_hawki))
legend('topleft',leg.txt, lty=c(1, NA, NA), lwd=c(5, NA, NA), col=c('blue', NA, NA))

plot(k_hawki, k-k_hawki, xlab='(K_HAWKI)', ylab='(K_POWALKA)-(K_HAWKI)')
abline(h=-0.8)
abline(h=0.5)

op = readline('Deseja tentar clustering? y/n ')


if(op == 'y'){
  op2 = readline('uik or gzk? 1/2 ')
  require(dbscan)
  
  k = k_hawki
  gk = gkh
  zk = zkh
  kerr = kerr_hawki
  errgk = errgkh
  errzk = errzkh
  errik = errikh
  
  
  df = data.frame(g, z, k,  gerr, zerr, kerr, gz, zk, gk, ui, ik,  errgz, errzk, errgk, errui, errik)
  
  if(op2==1){dfaux = data.frame(ui, ik)}else{dfaux = data.frame(gz, zk)}
  
  kNNdistplot(dfaux, k=4)
  abline(h=0.1, lty=2)
  abline(h=0.2, lty=2)
  abline(h=0.3, lty=2)
  abline(h=0.4, lty=2)
  
  eps = readline('Entre com o valor mais adequado de eps: ')
  
  db <- fpc::dbscan(dfaux, eps = eps, MinPts = 4)
  fviz_cluster(db, data = dfaux, stand = FALSE,
               ellipse = FALSE, show.clust.cent = FALSE, outlier.color='red',xlim=c(-1, 5),ylim=c(-4.5, 4.5),
               geom = "point", palette = "jco", ggtheme = theme_classic())
  
  db_gmm(db, df, 1, op2)
  db_gmm(db, df, 0, op2)
}

readline('Press ENTER to keep going ')

######### liris
######### 
######### 
print('Starting LIRIS...')
cols = readcat('gzk', liris1)

g = cols$g-av[2]
gerr = cols$gerr
z = cols$z-av[5]
zerr = cols$zerr
k = cols$k-av[6]
kerr = cols$kerr
k_hawki = cols$k_hawki-av[6]
kerr_hawki = cols$kerr_hawki

#converting from VEGA to AB system:

k_hawki = k_hawki+1.827
k = k+1.827

#### colors

zk = z-k
zkh = z-k_hawki
gk = g-k
gkh = g-k_hawki
gz = g-z

########### error in colors

errzk = sqrt(zerr^2+kerr^2)
errzkh = sqrt(zerr^2+kerr_hawki^2)
errgk = sqrt(gerr^2+kerr^2)
errgkh = sqrt(gerr^2+kerr_hawki^2)
errgz = sqrt(gerr^2+zerr^2)

color_color(zk, gk, errzk, errgk, '(z-k)', '(g-k)', main='HAWKI_LIRIS')
par(new=TRUE)
color_color(zkh, gkh, errzkh, errgkh, '(z-k)', '(g-k)', m_cor='red', p_cor='blue', axes=FALSE)

color_color(gz, zk, errgz, errzk, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), main='HAWKI_LIRIS')
par(new=TRUE)
color_color(gz, zkh, errgz, errzkh, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), m_cor='blue', p_cor='red', p_ch=1, axes=FALSE)
abline(-0.2179487179,1.2820512821, col='purple')
abline(h=2.5, col='purple', lty=2)
abline(-0.7,0.4, col='purple', lty='dotted')

ploterr(k, k_hawki, kerr, kerr_hawki, xlab='K_LIRIS', ylab='K_HAWKI', pch=19)
abline(0, 1, col='blue', lwd=5)
mk = mean(kerr)
mk_hawki = mean(kerr_hawki)
leg.txt = c('x=y', sprintf('<err_liris>: %f',mk), sprintf('<err_hawki>: %f ',mk_hawki))
legend('topleft',leg.txt, lty=c(1, NA, NA), lwd=c(5, NA, NA), col=c('blue', NA, NA))

op = readline('Deseja tentar clustering? y/n ')

if(op == 'y'){
  op2 = 2
  require(dbscan)
  
  k = k_hawki
  gk = gkh
  zk = zkh
  kerr = kerr_hawki
  errgk = errgkh
  errzk = errzkh
  errik = errikh
  
  
  df = data.frame(g, z, k,  gerr, zerr, kerr, gz, zk, gk, errgz, errzk, errgk)
  
  dfaux = data.frame(gz, zk)
  
  kNNdistplot(dfaux, k=4)
  abline(h=0.1, lty=2)
  abline(h=0.2, lty=2)
  abline(h=0.3, lty=2)
  abline(h=0.4, lty=2)
  
  eps = readline('Entre com o valor mais adequado de eps: ')
  
  db <- fpc::dbscan(dfaux, eps = eps, MinPts = 4)
  fviz_cluster(db, data = dfaux, stand = FALSE,
               ellipse = FALSE, show.clust.cent = FALSE, outlier.color='red',xlim=c(-1, 5),ylim=c(-4.5, 4.5),
               geom = "point", palette = "jco", ggtheme = theme_classic())
  
  db_gmm(db, df, 1, op2)
  db_gmm(db, df, 0, op2)
}

readline('Press ENTER to keep going ')

###### liris-hawki

print('Starting Liris w/Oldham...')

lir_old = read.table('/home/emilio/HAWKI/products-inter-combine/liris_oldham.dat')

g_old = lir_old$V11
g_olderr = lir_old$V12
z_old = lir_old$V17
z_olderr = lir_old$V18
k_old = lir_old$V7
k_old_err = lir_old$V8

k_old = k_old+1.827

#### colors

zko = z_old-k_old
gko = g_old-k_old
gzo = g_old-z_old

########### error in colors

errzko = sqrt(z_olderr^2+k_old_err^2)
errgko = sqrt(g_olderr^2+k_old_err^2)
errgzo = sqrt(g_olderr^2+z_olderr^2)

color_color(zko, gko, errzko, errgko, '(z-k)', '(g-k)', main='OLDHAM_LIRIS')
par(new=TRUE)
color_color(zkh, gkh, errzkh, errgkh, '(z-k)', '(g-k)', m_cor='red', p_cor='blue', axes=FALSE)

color_color(gzo, zko, errgzo, errzko, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), main='OLDHAM_LIRIS')
par(new=TRUE)
color_color(gz, zkh, errgz, errzkh, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), m_cor='blue', p_cor='red', p_ch=1, axes=FALSE)
abline(-0.2179487179,1.2820512821, col='purple')
abline(h=2.5, col='purple', lty=2)
abline(-0.7,0.4, col='purple', lty='dotted')


######### ACS
######### 
######### 
print('Starting ACS...')
cols = readcat('gzk', acs)

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

########### error in colors

errzk = sqrt(zerr^2+kerr^2)
errgk = sqrt(gerr^2+kerr^2)
errgz = sqrt(gerr^2+zerr^2)


############ plots
color_color(zk, gk, errzk, errgk, '(z-k)', '(g-k)', main='HAWKI_ACS', colour=pgc)

rb = colorRampPalette(c('red', 'blue'))
col_pgc = rb(10)[as.numeric(cut(pgc, breaks=10))]

color_color(gz, zk, errgz, errzk, '(g-z)', '(z-k)', xlim=c(-1, 5), ylim=c(-1.5,3.5), main='HAWKI_ACS', p_cor=col_pgc)
abline(-0.2179487179,1.2820512821, col='purple')
abline(h=2.5, col='purple', lty=2)
abline(-0.7,0.4, col='purple', lty='dotted')

hist(pgc, main='Probability of being a GC (ACS)', xlab='P_GC', breaks=50)

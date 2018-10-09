remove_out <- function(df, bad_index){

  u = df$u[-c(bad_index)]
  g = df$g[-c(bad_index)]
  r = df$r[-c(bad_index)]
  i_band = df$i_band[-c(bad_index)]
  z = df$z[-c(bad_index)]
  k = df$k[-c(bad_index)]
  uerr = df$uerr[-c(bad_index)]
  gerr = df$gerr[-c(bad_index)]
  rerr = df$rerr[-c(bad_index)]
  i_band_err = df$i_band_err[-c(bad_index)]
  zerr = df$zerr[-c(bad_index)]
  kerr = df$kerr[-c(bad_index)]
  
  df = data.frame(u, g, r, i_band, z, k, uerr, gerr, rerr, i_band_err, zerr, kerr)
  
  return(df)
}


get_err_vecs <- function(est, err, err_thres1, err_thres2){
  est_less_half = c()
  est_less_1 = c()
  
  for(i in 1:length(est)){
    if(err[i] < err_thres1){
      est_less_half = append(est_less_half, est[i])
    }
    if(err[i] < err_thres2){
      est_less_1 = append(est_less_1, est[i])
    }
  }
  
  outp = list(est_less_1, est_less_half)
  
  return(outp)
}

Kerr_hist <- function(est, err, err_thres1, err_thres2, title){
  
  #this function only works if Kerr and K_band magnitudes are already loaded in the main code
  
  est_less_half = c()
  est_less_1 = c()
  
  for(i in 1:length(est)){
    if(err[i] < err_thres1){
      est_less_half = append(est_less_half, est[i])
    }
    if(err[i] < err_thres2){
      est_less_1 = append(est_less_1, est[i])
    }
  }
  
  hist(est, breaks='Sturges', main=NA, xlab=title)
  hist(est_less_1, add=TRUE, color='grey', density=10, breaks='Sturges')
  hist(est_less_half, add=TRUE, color='black', density=100, breaks='Sturges')
  legend('topleft', c('full sample', paste0('K_err < ',err_thres2), paste0('K_err < ', err_thres1)), fill=c('white', 'black', 'black'), density=c(NA, 10, 100))

  #ests = data.frame(est_less_half, est_less_1)
  return(est_less_half)
}

radialbins <- function(RA, DEC, est, err, err_thres, title){
  
  #this function only works if Kerr and K_band magnitudes are already loaded in the main code
  
  est_less_half = c()
  est_less_1 = c()
  RA_less_half = c()
  DEC_less_half = c()
  
  for(i in 1:length(est)){
    if(err[i] < err_thres){
      est_less_half = append(est_less_half, est[i])
      RA_less_half = append(RA_less_half, RA[i])
      DEC_less_half = append(DEC_less_half, DEC[i])
    }
  }
  
  hist(est, breaks='Sturges', main=NA, xlab=title)
  hist(est_less_half, add=TRUE, color='black', density=100, breaks='Sturges')
  legend('topleft', c('full sample', paste0('K_err < ',err_thres)), fill=c('white', 'black'), density=c(NA, 10))
  
  halfgc = data.frame(RA_less_half, DEC_less_half, est_less_half)
  return(halfgc)
}


gg.mixEM <- function(est1, t, xl) {
  require(ggplot2)
  
  n = length(est1$x)
  bwd = (2 * IQR(est1$x) / length(est1$x)^(1/3))
  bszs = diff(range(est1$x)) / bwd
  scl = n*bwd
  
  
  xrand=c(rnorm(10000,8,2),rnorm(10000,17,4))
  sdnorm =
    function(xrand, med=0, std=1, lambda=1, n=length(est1$x), bw=scl){lambda*dnorm(xrand, mean=med, sd=std)*bw}
  mixnorm = 
    function(xrand, med1=0, med2=0, std1=1, std2=1, lambda1=1,lambda2=1, n=length(est1$x), bw=scl){
      (lambda1*(dnorm(xrand, mean=med1, sd=std1))+
         lambda2*(dnorm(xrand, mean=med2, sd=std2)))*scl}
  
  f1 = stat_function(fun=sdnorm,
                     args=list(med=est1$mu[1],
                               std=est1$sigma[1],
                               lambda=est1$lambda[1]),
                     colour="red",geom="area", fill='red', alpha=0.5)
  f2 = stat_function(fun=sdnorm,
                     args=list(med=est1$mu[2],
                               std=est1$sigma[2],
                               lambda=est1$lambda[2]),
                     colour="blue",geom="area", fill='blue', alpha=0.5)
  
  fmix = stat_function(fun=mixnorm,
                       args=list(med1=est1$mu[1],
                                 std1=est1$sigma[1],
                                 lambda1=est1$lambda[1],
                                 med2=est1$mu[2],
                                 std2=est1$sigma[2],
                                 lambda2=est1$lambda[2]),
                       colour="purple",geom="line", size=1)
  
  densis = ggplot(data=data.frame(est1$x))+geom_density(aes(x=est1$x))
  p = ggplot_build(densis)
  d = as.data.frame(p$data)
  
  d$y = d$y*scl
  
  a = ggplot(data=data.frame(est1$x), aes(x=est1$x)) + 
    geom_histogram(aes(x=est1$x,y=..count..),fill="grey",color="black", bins = bszs) +
    geom_line(data=d, aes(x=d$x, y=d$y), size=1) +
    fmix+f1+f2+
    #labs(title=t, x=xl, y='NGC')+theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  plot(a)
  
  #ggsave(fname)
}

gg.mclust <- function(est1, t, xl) {
  require(ggplot2)
  
  n = length(est1$data)
  bwd = (2 * IQR(est1$data) / length(est1$data)^(1/3))
  bszs = diff(range(est1$data)) / bwd
  scl = n*bwd
  
  
  xrand=c(rnorm(10000,8,2),rnorm(10000,17,4))
  sdnorm =
    function(xrand, med=0, std=1, lambda=1, n=length(est1$data), bw=bszs){lambda*dnorm(xrand, mean=med, sd=std)*n*bw}
  mixnorm = 
    function(xrand, med1=0, med2=0, std1=1, std2=1, lambda1=1,lambda2=1, n=length(est1$data), bw=bszs){
      (lambda1*(dnorm(xrand, mean=med1, sd=std1))+
         lambda2*(dnorm(xrand, mean=med2, sd=std2)))*n*bw}
  
  f1 = stat_function(fun=sdnorm,
                     args=list(med=est1$parameters$mean[1],
                               std=est1$parameters$variance$sigmasq[1],
                               lambda=est1$parameters$pro[1]),
                     colour="red",geom="area", fill='red', alpha=0.5)
  f2 = stat_function(fun=sdnorm,
                     args=list(med=est1$parameters$mean[2],
                               std=est1$parameters$variance$sigmasq[2],
                               lambda=est1$parameters$pro[2]),
                   colour="blue",geom="area", fill='blue', alpha=0.5)
  
  fmix = stat_function(fun=mixnorm,
                       args=list(med1=est1$parameters$mean[1],
                                 std1=est1$parameters$variance$sigmasq[1],
                                 lambda1=est1$parameters$pro[1],
                                 med2=est1$parameters$mean[2],
                                 std2=est1$parameters$variance$sigmasq[2],
                                 lambda2=est1$parameters$pro[2]),
                       colour="purple",geom="line", size=1)
  
  densis = ggplot(data=data.frame(est1$data))+geom_density(aes(x=est1$data))
  p = ggplot_build(densis)
  d = as.data.frame(p$data)
  
  d$y = d$y#*scl
  
  a = ggplot(data=data.frame(est1$data), aes(x=est1$data)) + 
    geom_histogram(aes(x=est1$data,y=..count..),fill="grey",color="black", bins = bszs) +
    geom_line(data=d, aes(x=d$x, y=d$y), size=1) +
    fmix+f1+f2+
    labs(title=t, x=xl, y='NGC')+theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  plot(a)
  
  #ggsave(fname)
}

ploterr <- function(x, y, xerr,yerr, ylim=NULL, col='black', ..., revy='N') {
  if (is.null(ylim))
    ylim <- c(min(y-yerr), max(y+yerr))
  if (is.null(yerr)){
    plot(x, y, col=col, ...)
    arrows(x-xerr, y, x+xerr, y, length=0.05, angle=90, code=3)
  }else{
    if (revy=='N'){
  plot(x, y, col=col, ...)
  arrows(x, y-yerr, x, y+yerr, length=0.05, angle=90, code=3, col=col)
  arrows(x-xerr, y, x+xerr, y, length=0.05, angle=90, code=3, col=col)
    }else{
      plot(x, y, col=col, ..., ylim=ylim)
      arrows(x, y-yerr, x, y+yerr, length=0.05, angle=90, code=3, col=col)
      arrows(x-xerr, y, x+xerr, y, length=0.05, angle=90, code=3, col=col)
    }
  }
}

pbcm <- function(est1){
  require(mclust)
  
  est1.gmm = Mclust(est1, G=2, modelNames = 'V')
  est1.gmm.1 = Mclust(est1, G=1, modelNames = 'V')
  
  g1 = 0
  g2 = 0
  
  for(item in est1.gmm$class){
    if(item == 1){
      g1 = g1+1
    }else{
      g2 = g2+1
    }
  }
  
  m1.1 = est1.gmm$parameters$mean[1]
  m1.2 = est1.gmm$parameters$mean[2]
  v1.1 = est1.gmm$parameters$variance$sigmasq[1]
  v1.2 = est1.gmm$parameters$variance$sigmasq[2]
  
  m2.1 = est1.gmm.1$parameters$mean[1]
  v2.1 = est1.gmm.1$parameters$variance$sigmasq
  
  print(summary(est1.gmm))
  print(summary(est1.gmm.1))
  
  
  set.seed(7809)
  B = 1000;    x2.d = vector(length=B);    x1.d = vector(length=B)
  for(i in 1:B){
    x2      = c(rnorm(g1, mean=m1.1, sd=sqrt(v1.1)), rnorm(g2, mean=m1.2, sd=sqrt(v1.2)))
    x1      = rnorm(length(est1), mean=m2.1, sd=sqrt(v2.1))
    x2.d[i] = Mclust(x2, G=2)$loglik - Mclust(x2, G=1)$loglik
    x1.d[i] = Mclust(x1, G=2)$loglik - Mclust(x1, G=1)$loglik
  }
  
  x2.d = sort(x2.d);  x1.d = sort(x1.d)

  print(summary(x1.d))
  print(summary(x2.d))

  plot(density(x1.d), xlab='difference in log-likelihood', xlim=c(min(x1.d), max(x2.d)))
  lines(density(x2.d), col='green')
  legend('topright', c('1 component', '2 component'), lty=c(1, 1), col=c('black', 'green'))
}

color_color <- function(col1, col2, errcol1, errcol2, xlabel, ylabel, m_cor='green', p_cor='black', p_ch=1, ...){
  
  plot(col1, col2, xlab=xlabel, ylab=ylabel, pch=p_ch, col=p_cor, ...)
  arrows(col1, col2-errcol2, col1, col2+errcol2, length=0.05, angle=90, code=3, lwd= 0.3)
  arrows(col1-errcol1, col2, col1+errcol1, col2, length=0.05, angle=90, code=3, lwd= 0.3)
  #### binning and median plot
  
  #number_of_obj_per_bin = floor(length(col1)/10)
  
  #col1_or = sort(col1, decreasing=FALSE)
  #col2_or = sort(col2, decreasing=FALSE)
  
  #col1_bins = split(col1_or, ceiling(seq_along(col1_or)/number_of_obj_per_bin))
  #col2_bins = split(col2_or, ceiling(seq_along(col2_or)/number_of_obj_per_bin))
  
  #col1avg = list()
  #col2avg = list()
  
  #for(i in 1:11){
  #  col1avg = append(col1avg, median(col1_bins[[i]]))
  #  col2avg = append(col2avg, median(col2_bins[[i]]))
  #}
  
  #lines(col1avg, col2avg, type='l', col=m_cor, lwd=2)
  #leg.txt = c('M87 GCs','Median')
  #legend('bottomright',leg.txt, pch=c(5, NA), lty=c(NA, 1), col=c('black', 'green') )
}

gzbins <- function(RA, DEC, nbins, gz, set_bins=FALSE, bins=list()){
  
  require(reticulate)
  require(mclust)
  source_python('/home/emilio/HAWKI/products-inter-combine/binning_functions_forR.py')
  pa = 153.0
  incl = 31.0
  RAgal = '12h30m49.4s'
  DECgal = '+12d23m28s'
  d = 22.2
  binning(RA, DEC, RAgal, DECgal, d, nbins, pa, incl)
  
  r = read.table('temp_radius_Rpy')$V1
  x = read.table('temp_radius_Rpy')$V2
  y = read.table('temp_radius_Rpy')$V3
  
  if(set_bins==TRUE){rbins = bins}else{rbins = read.table('temp_bins_Rpy')$V1}
  
  ibins = replicate(nbins, c())
  xbins = replicate(nbins, c())
  ybins = replicate(nbins, c())
  zbins = replicate(nbins, c())

  for(i in 1:nbins){
    for(j in 1:length(r)){
      if((rbins[i+1]>=r[j])&(r[j]>rbins[i])){
        ibins[[i]] = append(ibins[[i]], j)
        xbins[[i]] = append(xbins[[i]], j)
        ybins[[i]] = append(ybins[[i]], j)
        zbins[[i]] = append(zbins[[i]], j)
      }
    }
  }
   ########## verifying if (g-z) is bimodal in each bin
 
  colbins = replicate(nbins, c())
  estbins = replicate(nbins, list())
  
  require(propagate)
  for(i in 1:nbins){
    #pbcm(gzbins[[i]])
    pv = 10
    iter=0
    exit=0
    #set.seed(5678)
    print('calculating best chisq...')
    
    pvalues=c()
    
    while(exit!=1){
    colbins[[i]] = gz[ibins[[i]]]
    n = length(colbins[[i]])
    est1 = Mclust(colbins[[i]], G=2, modelNames='V', verbose=FALSE)
    estbins[[i]] = append(estbins[[i]], est1)
    titles = paste0(round(rbins[i], 2), ' < r (\'\') < ', round(rbins[i+1], 2))
    #plot(est1, which=2, xlab2='(z-k)', main2=title, col2=c('blue', 'red'), breaks=20)

    bwd = (2 * IQR(colbins[[i]]) / length(colbins[[i]])^(1/3))
    bszs = diff(range(colbins[[i]])) / bwd
    scl = n*bwd
    
    xrand=c(rnorm(10000,8,2),rnorm(10000,17,4))
    sdnorm =
      function(xrand, mean=0, sd=1, lambda=1, n=length(colbins[[i]]), bw=bwd){lambda*dnorm(xrand, mean=mean, sd=sd)*n*bw}
    mixnorm = 
      function(xrand, mean1=0, mean2=0, sd1=1, sd2=1, lambda1=1,lambda2=1, n=length(colbins[[i]]), bw=bwd){(lambda1*(dnorm(xrand, mean=mean1, sd=sd1))+lambda2*(dnorm(xrand, mean=mean2, sd=sd2)))*n*bw}
    mixnorm_alt = 
      function(xrand, mean1=0, mean2=0, sd1=1, sd2=1, lambda1=1,lambda2=1, n=length(colbins[[i]]), bw=bwd){
        f = (lambda1*(dnorm(xrand, mean=mean1, sd=sd1))+lambda2*(dnorm(xrand, mean=mean2, sd=sd2)))*n*bw
        return(f)}
    trimixnorm = 
      function(xrand, mean1=0, mean2=0, mean3=0, sd1=1, sd2=1,sd3=1, lambda1=1,lambda2=1,lambda3=1, n=length(colbins[[i]]), bw=bwd){(lambda1*(dnorm(xrand, mean=mean1, sd=sd1))+lambda2*(dnorm(xrand, mean=mean2, sd=sd2))+lambda3*(dnorm(xrand, mean=mean3, sd=sd3)))*n*bw}
    
    bw <- diff(range(colbins[[i]])) / (2 * IQR(colbins[[i]]) / length(colbins[[i]])^(1/3))
    
    densis = ggplot(data=data.frame(colbins[[i]]))+geom_density(aes(x=colbins[[i]]))
    p = ggplot_build(densis)
    d = as.data.frame(p$data)
    
    d$y = d$y*scl
    xchi=c(rnorm(length(d$y)/2,8,2),rnorm(length(d$y)/2,17,4))
  
    chi = chisq.test(d$y/scl, mixnorm_alt(xchi, mean1=est1$parameters$mean[1],
                             sd1=est1$parameters$variance$sigmasq[1],
                             lambda1=est1$parameters$pro[1],
                             mean2=est1$parameters$mean[2],
                             sd2=est1$parameters$variance$sigmasq[2],
                             lambda2=est1$parameters$pro[2])/scl)
    pv = chi$p.value
    pvalues = append(pvalues, pv)
    iter = iter+1
    
    if(iter==1000){
      
      pvalues_k1 = c()
      print(paste0('WARNING: not convergent for k=2, testing k=1... best p_value for k=2: ', 
                   round(min(pvalues), digits=4)))
      iter2 = 0
      exit2 = 0
      
      while(exit2!=1){
          colbins[[i]] = gz[ibins[[i]]]
          n = length(colbins[[i]])
          est1_k1 = Mclust(colbins[[i]], G=1, modelNames='V', verbose=FALSE)
          estbins[[i]] = append(estbins[[i]], est1_k1)
          titles = paste0(round(rbins[i], 2), ' < r (\'\') < ', round(rbins[i+1], 2))
          
          densis = ggplot(data=data.frame(colbins[[i]]))+geom_density(aes(x=colbins[[i]]))
          p = ggplot_build(densis)
          d = as.data.frame(p$data)
        
          d$y = d$y*scl
          xchi=c(rnorm(length(d$y)/2,8,2),rnorm(length(d$y)/2,17,4))
      
          chi = chisq.test(d$y/scl, sdnorm(xchi, mean=est1_k1$parameters$mean[1],
                             sd=est1_k1$parameters$variance$sigmasq[1],
                             lambda=est1_k1$parameters$pro[1])/scl)
          iter2 = iter2+1
          pv_k1 = chi$p.value
          pvalues_k1 = append(pvalues_k1, pv)
          
          if(pv_k1 <=0.05){
              
              exit2=1}
          if(iter2==1000){
            
            
            pv_k1 = min(pvalues_k1)
            print(paste0('WARNING: not convergent for k=1 too... best p_value for k=1: ', round(pv_k1, digits=4)))
            exit2=1}
      }
      
      if(round(pv_k1, digits=4) >= round(pv, digits=4)){
        est1 = est1_k1
        pv = min(pvalues_k1)
      }else{
        pv=min(pvalues)}
      
      exit = 1
    }
    if(pv <= 0.05){exit=1}
    }
    
    
    print(paste0('best p-value found: p = ', pv))
    
    f1 = stat_function(fun=sdnorm,
                    args=list(mean=est1$parameters$mean[1],
                             sd=est1$parameters$variance$sigmasq[1],
                             lambda=est1$parameters$pro[1]),
                    colour="red",geom="area", fill='red', alpha=0.5)
    f2 = stat_function(fun=sdnorm,
                    args=list(mean=est1$parameters$mean[2],
                             sd=est1$parameters$variance$sigmasq[2],
                             lambda=est1$parameters$pro[2]),
                    colour="blue",geom="area", fill='blue', alpha=0.5)
    f3 = stat_function(fun=sdnorm,
                    args=list(mean=est1$parameters$mean[3],
                             sd=est1$parameters$variance$sigmasq[3],
                             lambda=est1$parameters$pro[3]),
                    colour="green",geom="area", fill='green', alpha=0.5)
    
    fmix = stat_function(fun=mixnorm,
                    args=list(mean1=est1$parameters$mean[1],
                             sd1=est1$parameters$variance$sigmasq[1],
                             lambda1=est1$parameters$pro[1],
                             mean2=est1$parameters$mean[2],
                             sd2=est1$parameters$variance$sigmasq[2],
                             lambda2=est1$parameters$pro[2]),
                    colour="purple",geom="line", size=1)
    
    ftrimix = stat_function(fun=trimixnorm,
                    args=list(mean1=est1$parameters$mean[1],
                             sd1=est1$parameters$variance$sigmasq[1],
                             lambda1=est1$parameters$pro[1],
                             mean2=est1$parameters$mean[2],
                             sd2=est1$parameters$variance$sigmasq[2],
                             lambda2=est1$parameters$pro[2],
                             mean3=est1$parameters$mean[3],
                             sd3=est1$parameters$variance$sigmasq[3],
                             lambda3=est1$parameters$pro[3]),
                    colour="purple",geom="line", size=1)
    
       
    a = ggplot(data=data.frame(colbins[[i]]), aes(x=colbins[[i]])) + 
      geom_histogram(aes(x=colbins[[i]],y=..count..),fill="grey",color="black", bins = bszs) +
      geom_line(data=d, aes(x=d$x, y=d$y), size=1) +
      fmix+f1+f2+
      labs(title=titles, x='(z-k)', y='NGC')+theme(plot.title = element_text(hjust = 0.5))+
      theme_bw()
    
    plot(a)
    
    b = ggplot(data = data.frame(x, y))+
      geom_point(aes(x=x, y=y))+
      geom_point(data=data.frame(x[xbins[[i]]], y[ybins[[i]]]), aes(x=x[xbins[[i]]], y=y[ybins[[i]]], color='red'))+
      geom_density2d()+
      theme_bw()+guides(color=FALSE)+
      labs(x='RA (deg)', y='DEC (deg)')
    
    plot(b)
    
    print(paste0('nbin= ', i))
    #to aid gmm usage:
    gmmfile = paste0('/home/emilio/HAWKI/products-inter-combine/gmm-files-bins/gmm-gzbin', i)
    write.table(as.double(format(colbins[[i]], digits=4)), file=gmmfile, row.names = FALSE, col.names = FALSE)
    print(paste0('kurtosis= ', kurtosis(colbins[[i]])))
    print(paste0('Number of GC = ', length(colbins[[i]])))
    
  }
  return(list(estbins, xbins, ybins, rbins))
}

rungmm <- function(gmmfile, title){
  gmmfile = paste0('/home/emilio/HAWKI/products-inter-combine/gmm-files-bins/hawki/', gmmfile)
  f = read.table(gmmfile)$V1
  est = normalmixEM(f)
  m1 = est$mu[1]
  m2 = est$mu[2]
  print(m1)
  print(m2)
  system2('./gmm', args=c(gmmfile, '0', as.character(as.integer(m1)), as.character(as.integer(m2)), paste0('> /home/emilio/Downloads/gmm/gmm_output_', title)))
  system2('mv',  args=c('/home/emilio/Downloads/gmm/peakprob.out', paste0('/home/emilio/Downloads/gmm/peakprob_',title)))
  system2('./dip', args=c(as.character(length(f)), gmmfile, paste0('> /home/emilio/Downloads/gmm/dipout_',title)))
}

db_gmm <- function(db_object, dframe, n, op=2){
  cluster = which(db_object$cluster != n)
  
  if(op==2){
    col1 = dframe$zk[cluster]
    col2 = dframe$gz[cluster]
    col3 = dframe$gk[cluster]
    col1err = dframe$errzk[cluster]
    col2err = dframe$errgz[cluster]
    col3err = dframe$errgk[cluster]
  }else{
      
    col1 = dframe$ui[cluster]
    col2 = dframe$ik[cluster]
    col3 = dframe$gk[cluster]
    col1err = dframe$errui[cluster]
    col2err = dframe$errik[cluster]
    col3err = dframe$errgk[cluster]
    
    }
  
  if(op==2){
    color_color(col2, col1, col2err, col1err, '(g-z)', '(z-k)', xlim=c(-1,7), ylim=c(-2,4))
    abline(-0.2179487179, 1.2820512821, col='purple')
    abline(h=2.5, col='purple', lty=2)
    abline(-0.7, 0.4, col='purple', lty='dotted')
    
    estgz1 = normalmixEM(col2)
    estzk1 = normalmixEM(col1)
    estgk1 = normalmixEM(col3)
    
    plot(estgz1, which=2, xlab2='(g-z)', col2 = c('blue', 'red'))
    lines(density(col2))
    plot(estzk1, which=2, xlab2='(z-k)', col2 = c('blue', 'red'))
    lines(density(col1))
    plot(estgk1, which=2, xlab2='(g-k)', col2 = c('blue', 'red'))
    lines(density(col3))
  }else{
    color_color(col1, col2, col1err, col2err, '(u-i)', '(i-k)', xlim=c(-1,7), ylim=c(-2,4))
    abline(h=-0.4, col='purple', lty=2)
    abline(h=1.0, col='purple', lty=2)
  }
}

readcat <- function(op, cat){
  
  if(op=='uik'){
    corr_cat = read.table(cat)
    
    u = corr_cat$V3
    uerr = corr_cat$V4
    g = corr_cat$V5
    gerr = corr_cat$V6
    r = corr_cat$V7
    rerr = corr_cat$V8
    i_band = corr_cat$V9
    i_band_err = corr_cat$V10
    z = corr_cat$V11
    zerr = corr_cat$V12
    k = corr_cat$V13
    kerr = corr_cat$V14
    k_hawki = corr_cat$V15
    kerr_hawki = corr_cat$V16
    
    cols = data.frame(u, g, r, i_band, z, k, k_hawki, uerr, gerr, rerr, i_band_err, zerr, kerr, kerr_hawki)
    
    return(cols)
    
  }else if(op=='gzk'){
    corr_cat = read.table(cat)
  
    
    g = corr_cat$V3
    gerr = corr_cat$V4
    z = corr_cat$V5
    zerr = corr_cat$V6
    k = corr_cat$V7
    kerr = corr_cat$V8
    k_hawki = corr_cat$V9
    kerr_hawki = corr_cat$V10
    
    cols = data.frame(g, z, k, k_hawki, gerr, zerr, kerr, kerr_hawki)
    
    return(cols)
  }
}

iterdbscan <- function(x, y){
  control = 0
  
  while(control!=1){
    
    dfaux = data.frame(x, y)
    
    kNNdistplot(dfaux, k=4)
    abline(h=0.1, lty=2)
    abline(h=0.2, lty=2)
    abline(h=0.3, lty=2)
    abline(h=0.4, lty=2)
    
    eps = readline('Enter the most adequate eps value (approximatelly where the curve \'bends\' in the y-axis): ')
    
    db <- fpc::dbscan(data.frame(x, y), eps = eps, MinPts = 4)
    
    a = fviz_cluster(db, data = data.frame(x, y), stand = FALSE,
                     ellipse = FALSE, show.clust.cent = FALSE, outlier.color='red',xlim=c(-1, 5),ylim=c(-4.5, 4.5),
                     geom = "point", palette = "jco", ggtheme = theme_classic())
    plot(a)
    
    ms=which(db$cluster==1)
    xms = x[ms]
    yms = y[ms]
    
    message = sprintf('Number of points to consider for next run: %d', length(xms))
    print(message)
    
    
    control = readline('Input 1 to stop or 0 to keep going: ')
    
    x = xms
    y = yms
  }
  
  return(data.frame(x, y))
}

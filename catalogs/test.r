corr_cat = read.table('sub_oldham.cat')


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

zk = z-k
gk = g-k
uk = u-k      
gz = g-z

errzk = sqrt(zerr^2+kerr^2)
errgk = sqrt(gerr^2+kerr^2)
erruk = sqrt(uerr^2+kerr^2)
errgz = sqrt(gerr^2+zerr^2)

col1=zk
col2=gk
errcol1=errzk
errcol2=errgk
xlabel='zk'
ylabel='gk'
  
plot(col1, col2, xlab=xlabel, ylab=ylabel, pch=5)
arrows(col1, col2-errcol2, col1, col2+errcol2, length=0.05, angle=90, code=3, lwd= 0.3)
arrows(col1-errcol1, col2, col1+errcol1, col2, length=0.05, angle=90, code=3, lwd= 0.3)
#### binning and median plot
  
number_of_obj_per_bin = floor(length(col1)/10)
  
col1_or = sort(col1, decreasing=FALSE)
col2_or = sort(col2, decreasing=FALSE)
  
col1_bins = split(col1_or, ceiling(seq_along(col1_or)/number_of_obj_per_bin))
col2_bins = split(col2_or, ceiling(seq_along(col2_or)/number_of_obj_per_bin))
  
#col1_bins = list(col1_bins)
#col2_bins = list(col2_bins)
  
col1avg = list()
col2avg = list()
  
for(i in 1:11){
    col1avg = append(col1avg, median(col1_bins[[i]]))
    col2avg = append(col2avg, median(col2_bins[[i]]))
}
  
lines(col1avg, col2avg, type='l', col='green', lwd=2)
leg.txt = c('M87 GCs','Median')
legend('bottomright',leg.txt, pch=c(5, NA), lty=c(NA, 1), col=c('black', 'green') )






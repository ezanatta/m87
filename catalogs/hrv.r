library(ggplot2)

hrv = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/ko_roma_acs_hawki.cat')
acs = read.table('/home/emilio/Downloads/HAWKI-M87-NEW/processed_by_eso/combined/catalogs/old_acs_hawki.cat')

RAacs = acs$V1
DECacs = acs$V2

RAhrv = hrv$V1
DEChrv = hrv$V2

v = hrv$V9
err_v = hrv$V10

dfv = data.frame(RAhrv, DEChrv, v, err_v)
dfacs = data.frame(RAacs, DECacs)

ggplot(dfv, aes(x=RAhrv, y=DEChrv, color=v, size=2)) +
  geom_point(data=dfacs, aes(x=RAacs, y=DECacs), col='grey', size=1)+
  geom_point() + theme_bw()+
  scale_colour_gradientn(colours = c('blue', 'green', 'red'), space='Lab', 
                                                   na.value='red')




# Mengatur direktori kerja
setwd("C://Users//Hp//Downloads//Script Zee//downloadedScript")

# Memuat package
library(sp)
library(ggplot2)
library(sf)
library(ggspatial)
library(rgdal)
library(geoR)
library(Metrics)
library(raster)

########## EKSPLORASI DATA ##########

# Input data
dta0 <- read.csv("CHJabar.csv",sep=",")

# Membuat dataframe statistik deskriptif
eda_df <- data.frame(bulan="none",rata2=0,minimum=0,q1=0,median=0,q3=0,
				maksimum=0,jangkauan=0,varians=0,sim.baku=0)

list_bulan <- c("Januari","Februari","Maret","April","Mei","Juni","Juli",
			"Agustus","September","Oktober","November","Desember")

for (i in seq(6,17)) {
	nama <- list_bulan[i-5]
	rata2 <- mean(dta0[[i]])
	minm <- min(dta0[[i]])
	q1 <- quantile(dta0[[i]],0.25)
	q2 <- quantile(dta0[[i]],0.5)
	q3 <- quantile(dta0[[i]],0.75)
	maks <- max(dta0[[i]])
	jangkauan <- maks-minm
	varians <- var(dta0[[i]])
	sdv <- sd(dta0[[i]])

	eda_lst <- list(bulan=nama,rata2=rata2, minimum=minm, q1=q1, median=q2, 
			q3=q3,maksimum=maks, jangkauan=jangkauan, varians=varians, 
			sim.baku=sdv)
	eda_df = rbind(eda_df, eda_lst)
}

eda_df = eda_df[-1,]

########## VISUALISASI SEBARAN DATA SECARA SPASIAL ##########

dtaSP <- st_as_sf(dta0,coords=c("CURRENT.LONGITUDE","CURRENT.LATITUDE"),
                  crs=4326)

shp <- st_read("C://Users//Hp//Downloads//Script Zee//data/ProvJabar.shp")

# Fungsi visualisasi sebaran data secara spasial
sp_viz <- function(bulan, judul){
  	p <- ggplot()+
		geom_sf(data=shp)+
		geom_sf(data=dtaSP, aes(color=bulan), size=2)+
		labs(title=judul,
			color="Curah\nHujan (mm)")+
  		scale_color_gradientn(
    			colours = c('#FFFF00', '#FF9900', '#FF0000'),
    			rescaler = ~ scales::rescale_mid(.x, mid = 803),
    			limits = c(2, 1608))+
		xlab("Longitude")+
		ylab("Latitude")+
		ggspatial::annotation_north_arrow(location="tr")+
		ggspatial::annotation_scale(location="bl")
  	p <- p + theme_bw()
  	p <- p + theme(plot.title = element_text(hjust = 0.5),
          		legend.title = element_text(size = 10),
          		legend.text = element_text(size = 8))
  	return(p)
}

# Pemanggilan fungsi sp_viz
sp_viz(dtaSP$MAR,"Maret")
sp_viz(dtaSP$APR,"April")
sp_viz(dtaSP$MAY,"Mei")
sp_viz(dtaSP$JUN,"Juni")
sp_viz(dtaSP$JUL,"Juli")
sp_viz(dtaSP$AUG,"Agustus")
sp_viz(dtaSP$SEP,"September")
sp_viz(dtaSP$OKT,"Oktober")
sp_viz(dtaSP$NOV,"November")
sp_viz(dtaSP$DES,"Desember")

########## KONVERSI KOORDINAT KE UTM ##########

cord.dec <- SpatialPoints(cbind(dta0$CURRENT.LONGITUDE, -dta0$CURRENT.LATITUDE), 
                         proj4string = CRS("+proj=longlat"))

cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32748"))

dta <- cbind(x=cord.UTM@coords[,1], y=cord.UTM@coords[,2], dta0[,6:17])

########## PEMBAGIAN DATA ##########

# Memilih baris data untuk training set dan test set secara acak
test_ind <- ceiling(seq(1,214,length.out=42))
train_ind <- seq(1,214)[!seq(1,214) %in% test_ind]

train <- dta[train_ind,]
test <- dta[test_ind,]

dim(train)
dim(test)

# Visualisasi sebaran training dan test set
plot(cord.UTM@coords[train_ind,], col='green', pch=12, xlab="Koordinat X", 
	ylab="Koordinat Y")
points(cord.UTM@coords[test_ind,], col='red', pch=12)
legend(650000, 10860000, legend=c("Training", "Test"),  
       fill = c("green","red"), cex=0.7)

########## PERHITUNGAN SEMIVARIOGRAM EKSPERIMENTAL ##########

# Fungsi menghitung semivariogram eksperimental
smvg_eks <- function(cols, judul=NULL, show.plot=TRUE, set="train"){
	# Memilih training atau test set
	if (set=="test") {
	   dtaGEO <- as.geodata(test, coords.col=1:2,data.col=cols)
	} else {
	   dtaGEO <- as.geodata(train, coords.col=1:2,data.col=cols)
	}

	# Perhitungan semivariogram
	vario <- variog(dtaGEO, option="cloud", messages=FALSE)
	k <- 1+3.3*log10(length(vario$u))
	bin1 <- variog(dtaGEO, uvec = k, messages=FALSE)

	if (show.plot==FALSE){
	   # Mengembalikan jarak, banyak pasangan, nilai semivariogram tiap bin
	   return(list(cloud=vario,bin=bin1))
	} else {
	   # Menampilkan plot semivariogram eksperimental
	   plot(bin1, xlab="Jarak",ylab="Semivariogram")
	   title(judul)
	}
}

# Pemanggilan fungsi smvg_eks
print(smvg_eks(3,show.plot=FALSE)$bin$u) # jarak antar titik tiap bin
print(smvg_eks(3,show.plot=FALSE)$bin$n) # banyak pasangan tiap bin

for (i in 1:12) {
	# nilai semivariogram tiap bin
	print(list_bulan[i])
	print(smvg_eks(i+2,show.plot=FALSE)$bin$v)
}

# Plot semivariogram eksperimental 12 bulan
for (i in seq(1,12)){
  smvg_eks(i+2, judul=list_bulan[i], set="train")
}

########## PEMILIHAN MODEL SEMIVARIOGRAM TEORITIS ##########

# Fungsi menghitung skor AIC dan BIC suatu model semivariogram
fit_var <- function(kolom, ini.sig, ini.phi, ini.nug, teo.mod){

	dtaGEO <- as.geodata(train, coords.col=1:2,data.col=kolom)
	est <- likfit(dtaGEO, ini = c(ini.sig,ini.phi), nugget=ini.nug, 
			cov.model=teo.mod, messages=FALSE)
  
	return(list(aic=est$AIC,bic=est$BIC))
}

# Fungsi membandingkan model spherical, eksponensial, dan Gaussian
model_comp <- function(kolom, ini.sig, ini.phi, ini.nug){

	sph <- fit_var(kolom, ini.sig, ini.phi, ini.nug, "spherical")
	exp <- fit_var(kolom, ini.sig, ini.phi, ini.nug, "exponential")
	gaus <- fit_var(kolom, ini.sig, ini.phi, ini.nug, "gaussian")

	res_mat <- rbind(sph, exp, gaus)
	return(res_mat)
}

# Perbandingan skor AIC dan BIC untuk ketiga model
for (i in 1:12) {
	print(list_bulan[i])
	print(model_comp(i+2,20000,50000,10000))
}

########## ESTIMASI PARAMETER VARIOGRAM ##########

# Fungsi Bayesian Kriging
bayes_krig <- function(kolom, covm, sigsq, dfsig, 
                       phimin, phimax, taumin, taumax, loc) {
  dtaGEO = as.geodata(train, coords.col=1:2,data.col=kolom)
  
  bs <- krige.bayes(dtaGEO, locations=loc,
    	  model=model.control(cov.m=covm),
    	  prior=prior.control(beta.prior="flat",
    				sigmasq.prior="sc.inv.chisq",sigmasq=sigsq,df.sigmasq=dfsig,
    				phi.prior="uniform",phi.discrete=seq(phimin,phimax,l=100),
    				tausq.rel.prior="uniform",tausq.rel.discrete=c(taumin,taumax)),
    	  output=output.control(n.post=5000,messages=FALSE))
  
  return(bs)
}

# Fungsi estimasi parameter
smvg_par <- function(bs){
            return(list(
              meanBeta = mean(bs$posterior$sample$beta),
              sdBeta = sd(bs$posterior$sample$beta),
              meanSigma = mean(bs$posterior$sample$sigmasq),
              sdSigma = sd(bs$posterior$sample$sigmasq),
              meanPhi = mean(bs$posterior$sample$phi),
              sdPhi = sd(bs$posterior$sample$phi),
              meanTau = mean(bs$posterior$sample$tausq.rel),
              sdTau = sd(bs$posterior$sample$tausq.rel)
            ))
}

# Menghitung posterior parameter
jan.fix <- bayes_krig(3, "spherical", 17500, nrow(train), 
                      25000, 100000, 0.3, 0.999, loc=train_mat)
feb.fix <- bayes_krig(4, "spherical", 20000, nrow(train), 
                      10000, 100000, 0.5, 0.999, loc=train_mat)
mar.fix <- bayes_krig(5, "exponential", 22000, nrow(train), 
                      10000, 100000, 0.3, 0.5, loc=train_mat)
apr.fix <- bayes_krig(6, "exponential", 10000, nrow(train), 
                    10000, 100000, 0.3, 0.999, loc=train_mat)
mei.fix <- bayes_krig(7, "exponential", 20000, nrow(train), 
                      200000, 300000, 0.001, 0.1, loc=train_mat)
jun.fix <- bayes_krig(8, "exponential", var(train[[8]]), nrow(train), 
                      70000, 170000, 0.001, 0.999, loc=train_mat)
jul.fix <- bayes_krig(9, "spherical", 1050, nrow(train), 
                       50000, 100000, 0.1, 0.999, loc=train_mat)
aug.fix = bayes_krig(10, "spherical", 2800, nrow(train), 
                     270000, 290000, 0.1, 0.25, loc=train_mat)
sep.fix = bayes_krig(11, "exponential", 15000, nrow(train),
                     150000, 250000, 0.001, 0.999, loc=train_mat)
okt.fix = bayes_krig(12, "exponential", 35000, nrow(train), 
                     60000, 70000, 0.01, 0.5, loc=train_mat)
nov.fix = bayes_krig(13, "exponential",30000, nrow(train), 
                     100000, 150000, 0.01, 0.2, loc=train_mat)
des.fix = bayes_krig(14, "exponential", 75000, nrow(train),
                     25000, 130000, 0.1, 0.9, loc=train_mat)

# Hasil estimasi titik parameter dan ketidakpastiannya
smvg_par(jan.fix)
smvg_par(feb.fix)
smvg_par(mar.fix)
smvg_par(apr.fix)
smvg_par(mei.fix)
smvg_par(jun.fix)
smvg_par(jul.fix)
smvg_par(aug.fix)
smvg_par(sep.fix)
smvg_par(okt.fix)
smvg_par(nov.fix)
smvg_par(des.fix)

# Plot semivariogram teoritis dengan semivariogram eksperimental
## Januari
smvg_eks(3,"Januari")
lines(jan.fix, summ=mean)
## Februari
smvg_eks(4,"Februari")
lines(feb.fix, summ=mean)
## Maret
smvg_eks(5,"Maret")
lines(mar.fix, summ=mean)
## April
smvg_eks(6,"Apr")
lines(apr.fix, summ=mean)
## Mei
smvg_eks(7,"Mei")
lines(mei.fix, summ=mean)
## Juni
smvg_eks(8,"Juni")
lines(jun.fix, summ=mean)
## Juli
smvg_eks(9,"Juli")
lines(jul.fix, summ=mean)
## Agustus
smvg_eks(10,"Agustus")
lines(aug.fix, summ=mean)
## September
smvg_eks(11,"September")
lines(sep.fix, summ=mean)
## Oktober
smvg_eks(12,"Oktober")
lines(okt.fix, summ=mean)
## November
smvg_eks(13,"November")
lines(nov.fix, summ=mean)
## Desember
smvg_eks(14,"Desember")
lines(des.fix, summ=mean)

#Catatan: Jika semivariogram teoritis tidak mengikuti pola semivariogram
#eksperimental maka perlu dicoba kombinasi lain dari nilai hyperparameter
#pada fungsi bayes_krig

########## EVALUASI HASIL INTERPOLASI ##########

# Fungsi menghitung MAE, RMSE, MAPE
pred_acc <- function(bs, bulan, acc=TRUE, set) {
  
  pred_mean = bs$predictive$mean.simulations
  pred_var = bs$predictive$variance.simulations
  
  if (acc == FALSE) {
    return(list(pred_mean,pred_var))
  } else {
    if (set=="test"){
        return(list(mae=mae(c(test[[bulan]]), pred_mean),
                    rmse=rmse(c(test[[bulan]]), pred_mean),
                    mape=mape(c(test[[bulan]]), pred_mean)))
    } else if (set=="train") {
        return(list(mae=mae(c(train[[bulan]]), pred_mean),
                    rmse=rmse(c(train[[bulan]]), pred_mean),
                    mape=mape(c(train[[bulan]]), pred_mean)))
    } else {
      return("Error")
    }}
  }

# Training error
rbind(jan=pred_acc(jan.fix,3,set="train"),
      feb=pred_acc(feb.fix,4,set="train"),
      mar=pred_acc(mar.fix,5,set="train"),
      apr=pred_acc(apr.fix,6,set="train"),
      mei=pred_acc(mei.fix,7,set="train"),
      jun=pred_acc(jun.fix,8,set="train"),
      jul=pred_acc(jul.fix,9,set="train"),
      aug=pred_acc(aug.fix,10,set="train"),
      sep=pred_acc(sep.fix,11,set="train"),
      okt=pred_acc(okt.fix,12,set="train"),
      nov=pred_acc(nov.fix,13,set="train"),
      des=pred_acc(des.fix,14,set="train")
)

# Melakukan prediksi pada test set
jan.fix.test <- bayes_krig(3, "spherical", 17500,nrow(train),
                           25000, 100000, 0.3, 0.999,loc=test_mat)
feb.fix.test <- bayes_krig(4, "spherical", 20000, nrow(train),
                           10000, 100000, 0.5, 0.999, loc=test_mat)
mar.fix.test <- bayes_krig(5, "exponential", 22000, nrow(train),
                           10000, 100000, 0.3, 0.5, loc=test_mat)
apr.fix.test <- bayes_krig(6, "exponential", 10000, nrow(train), 
                           10000, 100000, 0.3, 0.999, loc=test_mat)
mei.fix.test <- bayes_krig(7, "exponential", 20000, nrow(train),
                           200000, 300000, 0.001, 0.1, loc=test_mat)
jun.fix.test <- bayes_krig(8, "exponential", var(train[[8]]), nrow(train), 
                           70000, 170000, 0.001, 0.999, loc=test_mat)
jul.fix.test <- bayes_krig(9, "spherical", 1050, nrow(train), 
                            50000, 100000, 0.3, 0.999, loc=test_mat)
aug.fix.test <- bayes_krig(10, "spherical", 2800, nrow(train), 
                           270000, 290000, 0.1, 0.25, loc=test_mat)
sep.fix.test <- bayes_krig(11, "exponential", 15000, nrow(train),
                           150000, 250000, 0.001, 0.999, loc=test_mat)
okt.fix.test <- bayes_krig(12, "exponential", 35000, nrow(train), 
                           60000, 70000, 0.01, 0.5, loc=test_mat)
nov.fix.test <- bayes_krig(13, "exponential",30000, nrow(train), 
                           100000, 150000, 0.01, 0.2, loc=test_mat)
des.fix.test <- bayes_krig(14, "exponential", 75000, nrow(train),
                           25000, 130000, 0.1, 0.9, loc=test_mat)

# Testing error
rbind(jan=pred_acc(jan.fix.test,3,set="test"),
      feb=pred_acc(feb.fix.test,4,set="test"),
      mar=pred_acc(mar.fix.test,5,set="test"),
      apr=pred_acc(apr.fix.test,6,set="test"),
      mei=pred_acc(mei.fix.test,7,set="test"),
      jun=pred_acc(jun.fix.test,8,set="test"),
      jul=pred_acc(jul.fix.test,9,set="test"),
      aug=pred_acc(aug.fix.test,10,set="test"),
      sep=pred_acc(sep.fix.test,11,set="test"),
      okt=pred_acc(okt.fix.test,12,set="test"),
      nov=pred_acc(nov.fix.test,13,set="test"),
      des=pred_acc(des.fix.test,14,set="test")
)

########## VISUALISASI PETA PREDIKSI ##########

# Menentukan grid prediksi 
## Batas derajat bujur dan derajat lintang pada shapefiles Jawa Barat
x_bnd <- c(106.3705,108.8338)
y_bnd <- c(-7.823398,-5.91377)
## Konversi ke koordinat UTM
bnd.dec <- SpatialPoints(cbind(x_bnd,-y_bnd),proj4string = CRS("+proj=longlat"))
bnd.UTM <- spTransform(bnd.dec, CRS("+init=epsg:32748"))
## Grid prediksi untuk wilayah Jawa Barat (resolusi 31x31)
pred.grid <- expand.grid(seq(bnd.UTM@coords[1,1],bnd.UTM@coords[2,1], l=31),
				seq(bnd.UTM@coords[2,2],bnd.UTM@coords[1,2], l=31))

# Melakukan interpolasi pada grid prediksi
jan.p <- bayes_krig(3, "spherical", 17500, nrow(train), 
                    25000, 100000, 0.3, 0.999, loc=pred.grid)
feb.p <- bayes_krig(4, "spherical", 20000, nrow(train), 
                    10000, 100000, 0.5, 0.999, loc=pred.grid)
mar.p <- bayes_krig(5, "exponential", 22000, nrow(train), 
                    10000, 100000, 0.3, 0.5, loc=pred.grid)
apr.p = bayes_krig(6, "exponential", 10000, nrow(train), 
                     10000, 100000, 0.3, 0.999, loc=pred.grid)
mei.p <- bayes_krig(7, "exponential", 20000, nrow(train), 
           200000, 300000, 0.001, 0.1, loc=pred.grid)
jun.p <- bayes_krig(8, "exponential", var(train[[8]]), nrow(train), 
                    70000, 170000, 0.001, 0.999, loc=pred.grid)
jul.p <- bayes_krig(9, "spherical", 1050, nrow(train), 
                       50000, 100000, 0.1, 0.999, loc=pred.grid)
aug.p <- bayes_krig(10, "spherical", 2800, nrow(train), 
                     270000, 290000, 0.1, 0.25, loc=pred.grid)
sep.p <- bayes_krig(11, "exponential", 15000, nrow(train), 
                     150000, 250000, 0.001, 0.999, loc=pred.grid)
okt.p <- bayes_krig(12, "exponential", 35000, nrow(train), 
                     60000, 70000, 0.01, 0.5, loc=pred.grid)
nov.p <- bayes_krig(13, "exponential",30000, nrow(train), 
                     100000, 150000, 0.01, 0.2, loc=pred.grid)
des.p <- bayes_krig(14, "exponential", 75000, nrow(train), 
                     25000, 130000, 0.1, 0.9, loc=pred.grid)

# Fungsi visualisasi peta prediksi
map_ch <- function(bk, judul) {
	   # Grid prediksi dalam derajat bujur dan derajat lintang
	   grid_dgr <- expand.grid(seq(106.3705,108.8338, l=31), 
			   seq(-7.823398,-5.91377, l=31))

	   # Mengubah hasil prediksi menjadi raster
	   raster.df <- data.frame(cbind(x=grid_dgr$Var1,y=grid_dgr$Var2,
			   id=c(bk$predictive$mean.simulations)))
	   rst <- rasterFromXYZ(raster.df)

	   # Potong raster sesuai bentuk shapefiles Jawa Barat
	   cropp <- crop(coba,extent(shp))
	   masked <- mask(cropp,shp)

	   # Menampilkan peta 
	   plot(masked,col=hcl.colors(15,palette="viridis",rev=TRUE), zlim=c(0,1300),  
		  main=sprintf("Peta Prediksi Curah Hujan Bulanan Jawa Barat\n%s 2020", judul))
	   plot(shp,add=TRUE,color="transparent")
}

# Memanggil fungsi map_ch
map_ch(jan.p, "Januari")
map_ch(feb.p, "Februari")
map_ch(mar.p, "Maret")
map_ch(apr.p, "April")
map_ch(mei.p, "Mei")
map_ch(jun.p, "Juni")
map_ch(jul.p, "Juli")
map_ch(aug.p, "Agustus")
map_ch(sep.p, "September")
map_ch(okt.p, "Oktober")
map_ch(nov.p, "November")
map_ch(des.p, "Desember")
 

##########stock assessment of S. putnamae######
##############set dictionary and load data########
setwd("D:/Vinoth")
mydata<-read.csv("Putnamae_data.csv") 
head(mydata1)

library(dplyr)
library(TropFishR)
library(ks)
library(ggplot2)

##########data arrangement#########
mydata1<- mydata[3:4]
colnames(mydata1)<-c("date","length")
head(mydata1)
ggplot(mydata1, aes(length)) +
  geom_histogram()

#jpeg("Histogram_length.jpg", res = 600,height = 6,width = 9,units = "in")
hist(mydata1$length, xlim = c(0,120),ylim = c(0,2000),xlab= "Length (cm)",main = "")
abline(v=38.26,col="blue")
text(x = mean(mydata1$length), y = 1900, label = paste("Mean:", round(mean(mydata1$length), 2)))
dev.off()

L_mean= mean(mydata1$length)
max(mydata1$length)
#create LF
mydata1$date<- as.Date(mydata1$date, format = "%d.%m.%Y")
putname5<- lfqCreate(data = mydata1, Lname = "length", Dname = "date", bin_size = 5,
                         species = "Sphyraena putnamae", stock = "SCI")
plot(putname5, Fname = "catch")

putname5<- lfqModify(putname5, bin_size = 5) #best one
plot(putname5, Fname = "catch")

putname10<- lfqModify(putname5, bin_size = 10) 
plot(putname10, Fname = "catch")

#LF restructuring with MA
putname5_ma3<-lfqRestructure(putname5, MA=3)
putname5_ma5<-lfqRestructure(putname5, MA=5) 
putname5_ma7<-lfqRestructure(putname5, MA=7)#best one
putname5_ma9<-lfqRestructure(putname5, MA=9)
plot(putname5_ma7,hist.sc = 0.75)
plot(putname5_ma7, Fname = "rcounts")


###############Linf and K#######
####### Powell Wetherall plot
par(mfcol = c(1,1))
PW <- powell_wetherall(param = putname5_ma7,
                       catch_columns = 1:ncol(putname5_ma7$catch),
                       reg_int = c(16,21))
PW$Linf_est
PW$confidenceInt_Linf

#########K scan
Lmax<- max(mydata1$length)
Linf_cal<-Lmax/0.95

KScan <- ELEFAN(putname5_ma7, method = "cross", Linf_fix = Linf_cal, K_range = seq(0.1,1.5,by=0.05),
                cross.date = putname5_ma7$dates[1],cross.midLength = putname5_ma7$midLengths[1])
#MA=7, addl.sqrt = TRUE, hide.progressbar = TRUE)
#K scan with optimise method
KScan1 <- ELEFAN(putname5_ma7, method = "optimise", Linf_range = seq(95,115,by=0.3), 
                 K_range = seq(0.1,1.5,by=0.05))
# show results
KScan1$par
KScan1$Rn_max

res_SA <- ELEFAN_SA(putname5, MA = 7, seasonalised = F, 
                     init_par = list(Linf = 105, K = 0.6, t_anchor = 0.5, C=0.5, ts = 0.5),
                     low_par = list(Linf = 100, K = 0.4, t_anchor = 0, C = 0, ts = 0),
                     up_par = list(Linf = 110, K = 0.9, t_anchor = 1, C = 1, ts = 1))

# show results
res_SA$par; res_SA$Rn_max

#ELEFAN_SA bootstrapped
lfq<-putname5
MA <- 7
init_par <- NULL
low_par <- list(Linf = 100, K = 0.4, t_anchor = 0)
up_par <- list(Linf = 110, K = 0.9, t_anchor = 1)
SA_time <- 15
SA_temp <- 1e5
nresamp <- 500

# parallel version
library(parallel)
t1 <- Sys.time()
res <- ELEFAN_SA_boot(lfq=lfq, MA = MA, seasonalised = F,
                      init_par = init_par, up_par = up_par, low_par = low_par,
                      SA_time = SA_time, SA_temp = SA_temp,
                      nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
                      seed = 1)
t2 <- Sys.time()
t2 - t1
res
# plot resulting distributions
univariate_density(res, use_hist = F)

# plot scatterhist of Linf and K
LinfK_scatterhist(res)

#jpeg("density_LinfK.jpg", res = 600,height = 6,width = 9,units = "in")
univariate_density(res, use_hist = F)
dev.off()

#plot the LF with SA1 result
#jpeg("LFcurve.jpg", res = 600,height = 6,width = 9,units = "in")
plot(res_SA,draw =F)
tmp <-lfqFitCurves(res_SA,par =res_SA$par,col=4,lty=2,draw =TRUE)
legend("top",ncol=1,legend =c("estimated"),col=4,lty=1)
dev.off()
############# assign estimates to the data list##############
putname5_1 <- c(putname5, res_SA$par)
#t0<-0 #not required in the TropfishR package as it gives any meaning
#coilia_wc1.0_1$t0=t0
###parameter updated from boostrap density
putname5_1$Linf<-102.52
putname5_1$K<-0.44
putname5_1$t_anchor<-0.98
putname5_1$phiL<-3.67
class(putname5_1) <- "lfq"

#################Estimate M,F,Z,E ################
# estimation of M
?M_empirical
tmax=t0+3/putname5_1$K
tmax1<- res_SA$agemax

Ms <- M_empirical(Linf = putname5_1$Linf, K_l = putname5_1$K, method = "Then_growth")
Ms.pauly <- M_empirical(Linf = putname5_1$Linf, K_l = putname5_1$K,temp = 29, 
                        method = "Pauly_Linf")
M.list<-M_empirical(Linf = res_SA1$par$Linf, K_l = res_SA1$par$K,temp = 28, tmax = tmax,
                    method = c("Then_growth","Pauly_Linf","AlversonCarney","Hoenig",
                               "Then_tmax"))
Ms
Ms.pauly
M.list
putname5_1$M <- as.numeric(Ms.pauly)
# show results M
paste("M =", as.numeric(Ms))

# run catch curve
# summarise catch matrix into vector and add plus group which is smaller than Linf
putname5_2 <- lfqModify(putname5_1, vectorise_catch = TRUE, plus_group = 102.5)
putname5_2$catch<-as.matrix(rowMeans(putname5_2$catch))

res_cc <- catchCurve(putname5_2,  catch_columns = 1:ncol(putname5_2$catch), 
                     reg_int = c(6,14), calc_ogive = T)
res_cc

#jpeg("Z_LCCC.jpg", res = 600,height = 6,width = 9,units = "in")
catchCurve(putname5_2,  catch_columns = 1:ncol(putname5_2$catch), 
                     reg_int = c(6,14), calc_ogive = F)
dev.off()
#conver age to length est. pop of catch at 50% Length

#jpeg("L50.jpg", res = 600,height = 6,width = 9,units = "in")
plot(res_cc$midLengths,res_cc$Sest, type = "l",lwd=2,col="navyblue",
     xlab = "Length (cm)",ylab = "Proportion")
abline(v=res_cc$L50,col="blue")
text(x = 30, y = 0.05, label = paste(30.43))
#abline(h=0.5,col="blue")
dev.off()

#calculate F and E from Z and M
f=res_cc$Z-Ms.pauly
f
res_cc$Z
expl.rate=f/res_cc$Z
expl.rate

# assign the estimated parameters to the data
putname5_2$M <- as.numeric(Ms.pauly)
putname5_2$Z <- res_cc$Z
putname5_2$currF <- as.numeric(putname5_2$Z - putname5_2$M)
putname5_2$E <- putname5_2$currF/putname5_2$Z
putname5_2$a<-0.018
putname5_2$b<-2.632
str(putname5_2)

##########VPA###################
set.seed(1000)

vpa <- VPA(param = putname5_2,catch_columns = 1:ncol(putname5_2$catch), 
           terminalF = putname5_2$currF, terminalE = 0.5,analysis_type = "VPA", plot=T)
vpa
#plot
#jpeg("VPA.jpg", res = 600,height = 6,width = 9,units = "in")
VPA(param = putname5_2,catch_columns = 1:ncol(putname5_2$catch), 
    terminalF = putname5_2$currF, terminalE = 0.5,analysis_type = "VPA", plot=T)
dev.off()
###########Thompsen and bell#############
vpa$FM_calc
putname5_2$FM<-vpa$FM_calc

#add market value per kg for each length group
#rlfdata1$meanValue<-c(30,30,50,80,80,130,130,160,160,160,180,180,180,180,200,200)
#code for adding repeating values as alternative to above one
putname5_2$meanValue<-c(rep(150,2),rep(200,4),rep(250,5),rep(300,8))
#coilia_wc1.0_2$meanValue<-NULL

# Thompson and Bell model with changes in F
TB1 <- predict_mod(putname5_2, type = "ThompBell",
                   FM_change = seq(0,4,0.05), stock_size_1 = 1,FM_relative = T, 
                   curr.E = putname5_2$E, plot = T, hide.progressbar = TRUE)

# Thompson and Bell model with changes in F and Lc
TB2 <- predict_mod(putname5_2, type = "ThompBell",
                   FM_change = seq(0,4,0.1), Lc_change = seq(10,80,1),
                   stock_size_1 = 1,
                   curr.E = putname5_2$E, curr.Lc = res_cc$L50,
                   s_list = list(selecType = "trawl_ogive",
                                 L50 = res_cc$L50, L75 = res_cc$L75),
                   plot = T, hide.progressbar = TRUE)

TB1_frm<-as.data.frame(cbind(TB1$meanB,TB1$totY,TB1$FM_change))
names(TB1_frm)<-c("Biomass","Yield","Fishing mortality")
TB1_frm
# plot results

#jpeg("TB_model1.jpg", res = 600,height = 7,width = 9,units = "in")
par(mfrow = c(2,1), mar = c(4,5,2,8)+0.1, oma = c(1,0,0,0))
plot(TB1, mark = TRUE)
mtext("(a)", side = 3, at = -1, line = 0.6)
plot(TB2, type = "Isopleth", xaxis1 = "FM", mark = TRUE, contour = 6)
mtext("(b)", side = 3, at = -0.1, line = 0.6)
dev.off()
# Biological reference levels
TB1$df_Es

# Current yield and biomass levels
TB1$currents


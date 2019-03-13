setwd("/Users/tareen/Desktop/second_paper_local/2A2D_Large_Run")
data_2K2P <- read.table(file="2A2D_2p7mill.txt",sep="",fill=TRUE)
head(data_2K2P)

# notation: kij -> jth protein acting on ith protein
k11 <- data_2K2P$V6
k12 <- data_2K2P$V7
k13 <- data_2K2P$V8
k14 <- data_2K2P$V9
k21 <- data_2K2P$V10
k22 <- data_2K2P$V11
k23 <- data_2K2P$V12
k24 <- data_2K2P$V13

k31 <- data_2K2P$V14
k32 <- data_2K2P$V15
k33 <- data_2K2P$V16
k34 <- data_2K2P$V17
k41 <- data_2K2P$V18
k42 <- data_2K2P$V19
k43 <- data_2K2P$V20
k44 <- data_2K2P$V21


library(ggplot2)
par(mfrow=c(1,2))
qplot((k31-k32)/(k31+k32),bins=30)
qplot((k13-k14)/((k13+k14)),bins=60)
hist(test1,breaks=40)

par(mfrow=c(2,1))
plot((k31[1:100000]-k32[1:100000])/(k31[1:100000]+k32[1:100000]),type='l',ylim=c(-1,1))
lines(phase_A1[1:100000],type='l',col='red')
lines(-phase_A2[1:100000],type='l',col='red')
plot((k13[1:100000]-k14[1:100000])/(k13[1:100000]+k14[1:100000]),type='l',ylim=c(-1,1))
lines(phase_D1[1:100000],type='l',col='red')
lines(-phase_D2[1:100000],type='l',col='red')

mean(k13[1:2000000])
mean(k11[1:2000000])
mean(k33[1:2000000])
mean(k31[1:2000000])
### symbiosis vs. parasitic

# rescale between 0 and 1
library("scales")
k21_norm <- k21[1:100000]
k21_norm <- rescale(k21_norm)

k12_norm <- k12[1:100000]
k12_norm <- rescale(k12_norm)


k22_norm <- k22[1:100000]
k22_norm <- rescale(k22_norm)

k_A1_D1 <- k13[1:100000]
k_A1_D1 <- rescale(k_A1_D1)



k_A1_A1 <- k11[1:2000000]
k_A1_A1 <- rescale(k_A1_A1)

k_A2_A2 <- k22[1:2000000]
k_A2_A2 <- rescale(k_A2_A2)

k_A2_A2 <- k22[1:2000000]
k_A2_A2 <- rescale(k_A2_A2)

k_D1_D1 <- k33[1:2000000]
k_D1_D1 <- rescale(k_D1_D1)

k_D2_D2 <- k44[1:2000000]
k_D2_D2 <- rescale(k_D2_D2)

k_A1_A2 <- k21[1:2000000]
k_A1_A2 <- rescale(k_A1_A2)

k_A2_A1 <- k12[1:2000000]
k_A2_A1 <- rescale(k_A2_A1)

plot(k_A1_A2[1:100000],type='l',lwd=3)
qplot((k_A1_A2-k_A2_A1)/(k_A1_A2+k_A2_A1),geom='histogram',bins=40)

k_D1_D2 <- k43[1:2000000]
k_D1_D2 <- rescale(k_D1_D2)

# more phases than we're acknowleding

phase_A1_Clean <- phase_A1[1:100000]
phase_D1_Clean <- phase_D1[1:100000]

# temporary creation of combined phase
combined_phase <- phase_A1_Clean
for(i in 1:100000)
{
  combined_phase[i] <- (phase_A1_Clean[i]+phase_D1_Clean[i])/2.0
}

plot(1.2*k_A1_D1[1:100000],type='l',col='black',lwd=1.5,ylim=c(0,1),ylab="k_A1_D1",main = "Phase A1 (Red), Phase D1 (Blue), k_A1_D1 (Black)",xlab='Accepted Mutations')
lines(combined_phase,lwd=2.5,col='red')
#lines(phase_A1_Clean,lwd=3.5,col='red')
#lines(phase_D1_Clean,lwd=3.5,col='blue')

#cor(k_A1_D1,combined_phase)

lines(phase_A1_Clean[1:50000],lwd=3.5,col='red')
lines(phase_D1_Clean[1:50000],lwd=3,col='blue')


# synergy vs. competition
par(mfrow=c(2,1),mar=c(5.1,5.1,4.1,2.1))

#plot(2.5*k_A1_A2[1:100000],type='l',ylim=c(0,1),col='darkgray',ylab='k_A1_A2',xlab='Accepted Mutations', main='Phase-A1 and k_A1_A2')
#lines(phase_A1[1:100000],col='red',lwd=2)

plot(2.5*k_A1_A2[1:100000],type='l',lwd=1.5,ylim=c(0,1),col='darkgray',ylab=expression(italic(k)[A1 %->% A2]),xlab='Accepted Mutations',cex.axis=1.8,cex.lab=2.0)
lines(phase_A1_Clean[1:100000],col='red',lwd=3)

#plot(2.5*k_D1_D2[1:100000],type='l',ylim=c(0,1),col='darkgray',ylab='k_D1_D2',xlab='Accepted Mutations', main='Phase-D1 and k_D1_D2')
#lines(phase_D1[1:100000],col='blue',lwd=2)

plot(2.5*k_D1_D2[1:100000],type='l',lwd=1.5,ylim=c(0,1),col='darkgray',ylab=expression(italic(k)[D1 %->% D2]),xlab='Accepted Mutations',cex.axis=1.8,cex.lab=2.0)
lines(phase_D1_Clean[1:100000],col='blue',lwd=3)

phase_A1_Clean <- phase_A1[1:100000]
phase_D1_Clean <- phase_D1[1:100000]
#plot(phase_D1_Clean[82000:100000],type='l')
plot(phase_A1_Clean[90000:92000],type='l')


phase_A2_Clean <- phase_A2[1:100000]
plot(phase_A2_Clean[90000:100000],type='l')

# clean up phase A2
#phase_A2_Clean[90000:92000] <- 1

# clean up for auto-activation
#phase_A1_Clean[44500:48000] <- 0
#phase_A1_Clean[90000:92000] <- 0

# clean up for auto-deactivation
#phase_D1_Clean[3000:4000] <- 0
#phase_D1_Clean[17000:18000] <- 0
#phase_D1_Clean[40000:41300] <- 0
#phase_D1_Clean[44000:45000] <- 0
#phase_D1_Clean[63800:65000] <- 0
#phase_D1_Clean[74000:78000] <- 0
#phase_D1_Clean[82000:100000] <- 0

#############################
#### Supplemental figure ####
#############################

# autoactivation vs. auto-deactivation
par(mfrow=c(2,1),mar=c(5.1,6.1,4.1,3.1))

plot(1.35*k_A1_A1[1:100000],ylim=c(0.,1),type='l',col='darkgray',ylab=expression(italic(k)[A1 %->% A1]),lwd=1.5,cex.axis=2,cex.lab=2.5,xlab="")
lines(phase_A1_Clean,col='red',lwd=2.5)

plot(2.25*k_D1_D1[1:100000],type='l',ylim=c(0.,1),col='darkgray',ylab=expression(italic(k)[D1 %->% D1]),xlab='Accepted Mutations',lwd=1.5,cex.axis=2,cex.lab=2.5)
lines(phase_D1_Clean,col='blue',lwd=2.5)

#plot(k_D2_D2,type='l',ylim=c(0,0.5),col='darkgray',ylab='k_D2_D2',xlab='Accepted Mutations', main='Phase-D1 and k_D2_D2')
#lines(0.5*phase_D1,col='red',lwd=2)

###############################
#### k_D1_A1 vs k_D1_vs_A2 ####
###############################

k_D1_A1 <- k13[1:100000]
k_D1_A1 <- rescale(k_D1_A1)

k_D1_A2 <- k23[1:100000]
k_D1_A2 <- rescale(k_D1_A2)

par(mfrow=c(2,1),mar=c(5.1,8.1,4.1,3.1))


k_D1_A1_vs_k_D1_A2 <- (k_D1_A1-k_D1_A2)/(k_D1_A1+k_D1_A2)
plot(k_D1_A1_vs_k_D1_A2,ylim=c(-1,1),type='l',col='darkgray',ylab=expression(frac(italic(k)[D1 %->% A1]-italic(k)[D1 %->% A2],italic(k)[D1 %->% A1]+italic(k)[D1 %->% A2])),lwd=1.5,cex.axis=2,cex.lab=2.5,xlab="Accepted Mutations")
#lines(2.25*k_D1_A2[1:100000],type='l',ylim=c(0.,1),col='blue',ylab=expression(italic(k)[D1 %->% A2]),xlab='Accepted Mutations',lwd=1.5)

#plot(2.25*k_D1_A2[1:100000],type='l',ylim=c(0.,1),col='darkgray',ylab=expression(italic(k)[D1 %->% D1]),xlab='Accepted Mutations',lwd=1.5,cex.axis=2,cex.lab=2.5)
#lines(-phase_D1_Clean,col='blue',lwd=2.5)
lines(phase_A1_Clean,col='red',lwd=2)
lines(-phase_A2_Clean,col='blue',lwd=2)

###############################
#### k_D2_A1 vs k_D2_vs_A2 ####
###############################

k_D2_A1 <- k14[1:100000]
k_D2_A1 <- rescale(k_D2_A1)

k_D2_A2 <- k24[1:100000]
k_D2_A2 <- rescale(k_D2_A2)

#par(mfrow=c(1,1),mar=c(5.1,6.1,4.1,3.1))


k_D2_A1_vs_k_D2_A2 <- (k_D2_A1-k_D2_A2)/(k_D2_A1+k_D2_A2)
plot(k_D2_A1_vs_k_D2_A2,ylim=c(-1,1),type='l',col='darkgray',ylab=expression(frac(italic(k)[D2 %->% A1]-italic(k)[D2 %->% A2],italic(k)[D2 %->% A1]+italic(k)[D2 %->% A2])),lwd=1.5,cex.axis=2,cex.lab=2.5,xlab="Accepted Mutations")
#lines(2.25*k_D1_A2[1:100000],type='l',ylim=c(0.,1),col='blue',ylab=expression(italic(k)[D1 %->% A2]),xlab='Accepted Mutations',lwd=1.5)

#plot(2.25*k_D1_A2[1:100000],type='l',ylim=c(0.,1),col='darkgray',ylab=expression(italic(k)[D1 %->% D1]),xlab='Accepted Mutations',lwd=1.5,cex.axis=2,cex.lab=2.5)
#lines(-phase_D1_Clean,col='blue',lwd=2.5)
lines(phase_A1_Clean,col='red',lwd=2)
lines(-phase_A2_Clean,col='blue',lwd=2)

###############################
#### k_A1_D1 vs k_A1_vs_D2 ####
###############################

k_A1_D1 <- k31[1:100000]
k_A1_D1 <- rescale(k_A1_D1)

k_A1_D2 <- k41[1:100000]
k_A1_D2 <- rescale(k_A1_D2)

k_A1_D1_vs_k_A1_D2 <- (k_A1_D1-k_A1_D2)/(k_A1_D1+k_A1_D2)
par(mfrow=c(1,1),mar=c(5.1,6.1,4.1,3.1))
plot(k_A1_D1_vs_k_A1_D2,ylim=c(-1,1),type='l',col='darkgray',ylab=expression(italic(k)[A1 %->% D1]-italic(k)[A1 %->% D2]),lwd=1.5,cex.axis=2,cex.lab=2.5,xlab="Accepted Mutations")
lines(2*phase_D1_Clean-1,col='red',lwd=2)

#min(k21_norm)

#########################
# Ranjan's hypothesis
#########################

par(mfrow=c(1,2))
hist(k13)
hist((k13+k23)/2)

library(ggplot2)

par(mfrow=c(1,2))
qplot(log(k14),bins=40)
qplot(((k13-k23)/(k13+k23)),bins=20)

require(gridExtra)
bin_variable <-12

#par(mfrow=c(2,1),mar=c(5.1,5.1,3.1,2.1))


plot1 <- qplot((k13-k23)/(k13+k23),bins=14,xlab=expression(frac(k[D1 %->% A1]-k[D1 %->% A2],k[D1 %->% A1]+k[D1 %->% A2])), 
               geom="histogram", 
               ylab="Count", 
               fill=I("grey"),
               col=I("black"))+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=22,color = 'black'),
        axis.text.y = element_text(size=24,color = 'black'),
        axis.title.x = element_text(size=24,color = 'black'),
        axis.title.y = element_text(size=24,color = 'black'),
        plot.margin=unit(c(1,1,2,2), "cm")
  )

plot2 <- qplot((k31-k41)/(k31+k41),bins=bin_variable,xlab=expression(frac(k[A1 %->% D1]-k[A1 %->% D2],k[A1 %->% D1]+k[A1 %->% D2])), 
               geom="histogram", 
               ylab="Count", 
               fill=I("grey"),
               col=I("black"))+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=22,color = 'black'),
        axis.text.y = element_text(size=24,color = 'black'),
        axis.title.x = element_text(size=24,color = 'black'),
        axis.title.y = element_text(size=24,color = 'black'),
        plot.margin=unit(c(1,1,2,2), "cm")
  )
grid.arrange(plot1, plot2, ncol=2)

###########################
## Deactivation asymmetry 2
###########################

deactivation_asym_2 <- ((k13+k23)-(k14+k24))/((k13+k23)+(k14+k24))

qplot(deactivation_asym_2,bins=14,xlab="Deactivation Asymmetry 2", 
      geom="histogram", 
      ylab="Count", 
      fill=I("grey"),
      col=I("black"))+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=22,color = 'black'),
        axis.text.y = element_text(size=24,color = 'black'),
        axis.title.x = element_text(size=24,color = 'black'),
        axis.title.y = element_text(size=24,color = 'black'),
        plot.margin=unit(c(1,1,2,2), "cm")
  )
        
#######################
## Self asymmetry 1 ###
#######################

k_A1_A2_m_A2_A1 <- (k21-k12)/(k21+k12)

plot1 <- qplot((k11-k22)/(k11+k22),bins=14,xlab=expression(frac(italic(k)[A1 %->% A1]-italic(k)[A2 %->% A2],italic(k)[A1 %->% A1]+italic(k)[A2 %->% A2])), 
      geom="histogram", 
      ylab="Count", 
      fill=I("grey"),
      col=I("black"))+theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=22,color = 'black'),
        axis.text.y = element_text(size=24,color = 'black'),
        axis.title.x = element_text(size=24,color = 'black'),
        axis.title.y = element_text(size=24,color = 'black'),
        plot.margin=unit(c(1,1,2,2), "cm")
  )

k_D1_D2_m_D2_D1 <- (k43-k34)/(k43+k34)

plot2 <- qplot((k33-k44)/(k33+k44),bins=14,xlab=expression(frac(italic(k)[D1 %->% D1]-italic(k)[D2 %->% D2],italic(k)[D1 %->% D1]+italic(k)[D2 %->% D2])), 
               geom="histogram", 
               ylab="Count", 
               fill=I("grey"),
               col=I("black"))+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=22,color = 'black'),
        axis.text.y = element_text(size=24,color = 'black'),
        axis.title.x = element_text(size=24,color = 'black'),
        axis.title.y = element_text(size=24,color = 'black'),
        plot.margin=unit(c(1,1,2,2), "cm")
  )
grid.arrange(plot1, plot2, ncol=2)




plot(phase_A1,type='l',lwd=2,col='red')
lines(k12_norm,type='l',lwd=2,col='blue')


# activator and deactivator asymmetry
#deactivation_asym <- ((k13+k23)-(k14+k24))/((k13+k23)+(k14+k24))

deactivation_asym <- (k13-k14)/(k13+k14)
activation_asym <- (k31-k32)/(k31+k32)

test_act    <- (k31-k41)/(k31+k41)
test_act2    <- (k32-k42)/(k32+k42)

test_deact <- (k13-k23)/(k13+k23)
test_deact2 <- (k14-k24)/(k14+k24)




#activation_asym <- ((k31+k32)-(k41+k42))/((k31+k32)+(k41+k42))
#activation_asym <- ((k31+k41)-(k32+k42))/((k31+k41)+(k32+k42))
#plot(deactivation_asym[1:100000],type='l')

#hist(k00_k22,breaks=40,xlab=expression(frac(k[A1 %->% A1]-k[D1 %->% D1],k[A1 %->% A1]+k[D1 %->% D1])),probability = TRUE,col='GRAY',main = '2A-2D')
#hist(k_A1A2_D1D2,breaks=100,xlab=expression(frac(k[A1 %->% A2]-k[D1 %->% D2],k[A1 %->% A2]+k[D1 %->% D2])),probability = TRUE,col='GRAY',main = '2A-2D')


ess_A1 <- data_2K2P$V2
#ess_A1 <- data$V2

ess_A1 <- data_2K2P$V2[1:2000000]
ess_D1 <- data_2K2P$V4[1:2000000]
ess_A2 <- data_2K2P$V3[1:2000000]
ess_D2 <- data_2K2P$V5[1:2000000]

n <- length(ess_A1)-500

phase_A1 <- make.phase(ess_A1,n)
phase_D1 <- make.phase(ess_D1,n)
phase_A2 <- make.phase(ess_A2,n)
phase_D2 <- make.phase(ess_D2,n)

plot(phase_D1[1:20000],type='l',col='red',lwd=2,ylim=c(-1,1))
lines(-phase_D2[1:20000],type='l',col='blue',lwd=2)

par(mfrow=c(1,2))
hist((k13-k14)/(k13+k14))
hist(k14[5000:10000])


#par(mfrow=c(3,1),mgp=c(2.5,1,0))
par(mfrow=c(3,1),mgp=c(2.5,0.5,0),mar = c(5.1, 4.1, 4.1, 2.1))

# deactivator phases
plot(ess_D1,type='l',ylim=c(-1,1),col='sienna1',cex.main=1.5,main="(a)",adj = 0, ylab="Deactivator phases",xlab='',xaxt='n',cex.axis=1.75,cex.lab=2.0)
lines(-ess_D2,type='l',col='sienna1')
lines(phase_D1[1:100000],col="red",type='l',lwd='2')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='2')


par(mfrow=c(2,1),mgp=c(2.5,0.5,0),mar = c(5.1, 6.1, 4.1, 2.1))
# deactivation asymmetry plot

plot(test_deact2[1:100000],type='l',cex.main=1.5,main="(b)",adj = 0,xaxt='n',ylab=expression(frac(k[D2 %->% A1]-k[D2 %->% A2],k[D2 %->% A1]+k[D2 %->% A2])),col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=1.5)
lines(phase_A1[1:100000],col="red",type='l',lwd='2')
lines(-phase_A2[1:100000],col="blue",type='l',lwd='2') 


###############
## new figure 3
###############
par(mfrow=c(2,1),mgp=c(3.5,1.5,0),mar = c(7.1, 8.1, 2.1, 2.1))



plot(test_deact[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab=expression(frac(italic(k)[D1 %->% A1]-italic(k)[D1 %->% A2],italic(k)[D1 %->% A1]+italic(k)[D1 %->% A2])),col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(phase_A1[1:100000],col="red",type='l',lwd='2')
lines(-phase_A2[1:100000],col="red",type='l',lwd='2')




# deactivation asymmetry plot
plot(test_act[1:100000],type='l',cex.main=1.75,main="(b)",adj = 0,xlab='Accepted Mutations',ylab=expression(frac(italic(k)[A1 %->% D1]-italic(k)[A1 %->% D2],italic(k)[A1 %->% D1]+italic(k)[A1 %->% D2])),col='darkgray',ylim=c(-1,1),cex.axis=2.25,cex.lab=2.25)
lines(phase_D1[1:100000],col="blue",type='l',lwd='2')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='2')

plot.new()
legend(x="bottom",inset =0, ncol=3,legend=c("Activator Phase","Deactivator Phase","Rate Constant Asymmetry"),
       fill=c("red","blue","darkgray"))


###############
## end new figure 3
###############

#plot(test_deact[1:100000],type='l',cex.main=1.5,main="(b)",adj = 0,xaxt='n',ylab=expression(frac(k[A2 %->% D1]-k[A2 %->% D2],k[A2 %->% D1]+k[A2 %->% D2])),col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=1.5)
#lines(phase_D1[1:100000],col="red",type='l',lwd='2')
#lines(-phase_D2[1:100000],col="blue",type='l',lwd='2') 

#par(mfrow=c(1,2))
#hist(deactivation_asym,bins=20)
#hist(activation_asym,bins=20)

#par(mfrow=c(1,2))
#qplot(test_deact,bins=12)
#qplot(test_act,bins=12)

# Activation asymmetry plot
plot(test_act[1:100000],type='l',cex.main=1.5,main="(c)",adj = 0,ylab=expression(frac(k[A1 %->% D1]-k[A1 %->% D2],k[A1 %->% D1]+k[A1 %->% D2])),col='darkgray',ylim=c(-1,1),xlab='Accepted Mutations',cex.axis=1.75,cex.lab=1.5)
lines(phase_D1[1:100000],col="red",type='l',lwd='2')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='2')

plot.new()
legend(x="bottom", ncol=4,legend=c("Phase A(D) 1","Phase A(D) 2","Asymmetry","Essentiality"),
       fill=c("red","blue","gray50","sienna1"))

legend("bottom", legend=c('Deactivation Asymmetry','Activation Asymmetry'),ncol=2,
col=c("red", "blue"),lty=1,lwd=2,cex=0.8)

#####

# https://stackoverflow.com/questions/12402319/centring-legend-below-two-plots-in-r

# run following section 2nd
layout(matrix(c(1,2,3,3), nrow=2, bycol=TRUE), heights=c(2, 1))

par(mai=rep(0.5, 4))
# plot 1
plot(deactivation_asym[1:100000],type='l',main='Deactivation Asymmetry',col='gray50',ylim=c(-1,1),xlab='Accepted Mutations',ylab='')
lines(phase_D1[1:100000],col="red",type='l',lwd='3')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='3')
# plot 2
plot(activation_asym[1:100000],type='l',main='Activation Asymmetry',col='gray50',ylim=c(-1,1),xlab='Accepted Mutations',ylab='')
lines(phase_A1[1:100000],col="red",type='l',lwd='3')
lines(-phase_A2[1:100000],col="blue",type='l',lwd='3')

# run following line 1st
par(mai=c(0,0,0,0),mgp=c(2,1,0))

# run whole thing 3rd
plot.new()
legend(x="center", ncol=3,legend=c("Phase A(D) 1","Phase A(D) 2","Asymmetry"),
       fill=c("red","blue","gray50"))

#####

head(phase)

# make a phase array
for(i in 1:n)
{
  a <- i+500
  if(sum(ess_A1[i:a])>475)
  {
    phase[i] <- 1      
  }
  else
  {
    phase[i] <- 0
  }
}

plot(ess_A1[1:100000],col="red",type='l')
lines(phase[1:100000],col="blue",lwd='3')
cor(phase[1:100000],ess_A1[1:100000],use = "complete.obs")

make.phase <- function(essentiality,n) {
  phase <- rep(0.0, n)
  for(i in 1:n)
  {
    a <- i+500
    if(sum(essentiality[i:a])>475)
    {
      phase[i] <- 1      
    }
    else
    {
      phase[i] <- 0
    }
  }
  phase
}

# function for distributions in various phases
make.distribution <- function(k_i_j, phase, n) {
  
  in_phase_list <- list()
  out_phase_list <- list()
  
  counter_in_phase <- 1
  counter_out_phase <- 1
  for(i in 1:n)
  {
    if(phase[i]==1)
    {
      in_phase_list[[counter_in_phase]] <- k_i_j[i]
      counter_in_phase <- counter_in_phase + 1
    }
    if(phase[i]==0)
    {
      out_phase_list[[counter_out_phase]] <- k_i_j[i]
      counter_out_phase <- counter_out_phase + 1
    }
    
  }
  return_phase_list <- list(in_phase_list,out_phase_list)
  return_phase_list
}

# function for distributions in various phases
make.distribution.inter_phase <- function(k_i_j, phase1,phase2, n) {
  
  in_phase1_list <- list()
  out_phase1_list <- list()
  
  in_phase2_list <- list()
  out_phase2_list <- list()
  
  counter_in_phase1  <- 1
  counter_out_phase1 <- 1
  
  counter_in_phase2  <- 1
  counter_out_phase2 <- 1
  
  in_1_in_2_phase_list <- list()
  counter_in1_in2 <- 1
  in_1_only_phase_list <- list()
  counter_in1_only <- 1
  in_2_only_phase_list <- list()
  counter_in2_only <- 1
  neihter_phase_list <- list()
  neither_counter <- 1
  
  for(i in 1:n)
  {
    if(phase1[i]==1)
    {
      in_phase1_list[[counter_in_phase1]] <- k_i_j[i]
      counter_in_phase1 <- counter_in_phase1 + 1
    }
    if(phase1[i]==0)
    {
      out_phase1_list[[counter_out_phase1]] <- k_i_j[i]
      counter_out_phase1 <- counter_out_phase1 + 1
    }
    
    if(phase2[i]==1)
    {
      in_phase2_list[[counter_in_phase2]] <- k_i_j[i]
      counter_in_phase2 <- counter_in_phase2 + 1
    }
    if(phase2[i]==0)
    {
      out_phase2_list[[counter_out_phase2]] <- k_i_j[i]
      counter_out_phase2 <- counter_out_phase2 + 1
    }
    
    # and
    if(phase1[i]==1 & phase2[i]==1)
    {
      in_1_in_2_phase_list[[counter_in1_in2]] <- k_i_j[i]
      counter_in1_in2 <- counter_in1_in2 + 1
    }
    # only 1
    if(phase1[i]==1 & phase2[i]==0)
    {
      in_1_only_phase_list[[counter_in1_only]] <- k_i_j[i]
      counter_in1_only <- counter_in1_only + 1
    }
    # only 2
    if(phase1[i]==0 & phase2[i]==1)
    {
      in_2_only_phase_list[[counter_in2_only]] <- k_i_j[i]
      counter_in2_only <- counter_in2_only + 1
    }
    # neither
    if(phase1[i]==0 & phase2[i]==0)
    {
      neihter_phase_list[[neither_counter]] <- k_i_j[i]
      neither_counter <- neither_counter + 1
    }
  }
  return_phase_list <- list(in_phase1_list,out_phase1_list,in_phase2_list,out_phase2_list,in_1_in_2_phase_list,in_1_only_phase_list,in_2_only_phase_list,neihter_phase_list)
  return_phase_list
}

k_A1_A2_phase_dists <- make.distribution(k_A1_A2,phase_A1,n)
k_D1_D2_phase_dists <- make.distribution(k_D1_D2,phase_D1,n)
k_A1_A1_phase_dists <- make.distribution(k_A1_A1,phase_A1,n)
k_A2_A2_phase_dists <- make.distribution(k_A2_A2,phase_A2,n)
k_D1_D1_phase_dists <- make.distribution(k_D1_D1,phase_D1,n)
k_D2_D2_phase_dists <- make.distribution(k_D2_D2,phase_D2,n)

k_A1_D1_phase_dists <- make.distribution.inter_phase(k_A1_D1,phase_A1,phase_D1,n)
k_D1_A1_phase_dists <- make.distribution.inter_phase(k_D1_A1,phase_A1,phase_D1,n)
k_A1_A2_phase_dists <- make.distribution.inter_phase(k_A1_A2,phase_A1,phase_A2,n)
k_D1_D2_phase_dists <- make.distribution.inter_phase(k_D1_D2,phase_D1,phase_D2,n)

# A1_A1, A2_A2, D1_D1, D2_D2

y <- qplot(unlist(k_A1_A1_phase_dists[[1]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A1',main='A1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot1 <- y

y <- qplot(unlist(k_A1_A1_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A1',main='A1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot2 <- y

y <- qplot(unlist(k_A2_A2_phase_dists[[1]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A2_A2',main='A2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot3 <- y

y <- qplot(unlist(k_A2_A2_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A2_A2',main='A2 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot4 <- y


y <- qplot(unlist(k_D1_D1_phase_dists[[1]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D1',main='D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot5 <- y

y <- qplot(unlist(k_D1_D1_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D1',main='D1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot6 <- y


y <- qplot(unlist(k_D2_D2_phase_dists[[1]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D2_D2',main='D2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot7 <- y

y <- qplot(unlist(k_D2_D2_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D2_D2',main='D2 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot8 <- y

require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol=2,nrow=4,top='Distributions of Auto acitvation/inhibition under various phases')



library(ggplot2)
# A1_A2
x <- qplot(unlist(k_A1_A2_phase_dists[[1]]),geom='histogram',binwidth=number_breaks, xlim = c(0,1),xlab='k_A1_A2',main='A1 dominant')
#plot1 <- qplot(unlist(k_A1_D1_phase_dists[[1]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot1 <- x
#x <- ggplot_build(x)

y <- qplot(unlist(k_A1_A2_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='A1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot2 <- y

y <- qplot(unlist(k_A1_A2_phase_dists[[3]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='A2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot3 <- y

y <- qplot(unlist(k_A1_A2_phase_dists[[4]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='A2 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot4 <- y
#y <- ggplot_build(y)

y <- qplot(unlist(k_A1_A2_phase_dists[[5]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='A1 and A2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot5 <- y

y <- qplot(unlist(k_A1_A2_phase_dists[[6]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='Only A1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot6 <- y

y <- qplot(unlist(k_A1_A2_phase_dists[[7]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='Only A2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot7 <- y

y <- qplot(unlist(k_A1_A2_phase_dists[[8]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_A2',main='Neither dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot8 <- y


require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol=2,nrow=4,top='Distributions of k_A1_A2 under various phases')



# D1_D2
x <- qplot(unlist(k_D1_D2_phase_dists[[1]]),geom='histogram',binwidth=number_breaks, xlim = c(0,1),xlab='k_D1_D2',main='D1 dominant')
#plot1 <- qplot(unlist(k_A1_D1_phase_dists[[1]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot1 <- x
#x <- ggplot_build(x)

y <- qplot(unlist(k_D1_D2_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='D1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot2 <- y

y <- qplot(unlist(k_D1_D2_phase_dists[[3]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='D2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot3 <- y

y <- qplot(unlist(k_D1_D2_phase_dists[[4]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='D2 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot4 <- y
#y <- ggplot_build(y)

y <- qplot(unlist(k_D1_D2_phase_dists[[5]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='D1 and D2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot5 <- y

y <- qplot(unlist(k_D1_D2_phase_dists[[6]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='Only D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot6 <- y

y <- qplot(unlist(k_D1_D2_phase_dists[[7]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='Only D2 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot7 <- y

y <- qplot(unlist(k_D1_D2_phase_dists[[8]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_D2',main='Neither dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot8 <- y


require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol=2,nrow=4,top='Distributions of k_D1_D2 under various phases')

# A1_D1
x <- qplot(unlist(k_A1_D1_phase_dists[[1]]),geom='histogram',binwidth=number_breaks, xlim = c(0,1),xlab='k_A1_D1',main='A1 dominant')
#plot1 <- qplot(unlist(k_A1_D1_phase_dists[[1]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot1 <- x
#x <- ggplot_build(x)

y <- qplot(unlist(k_A1_D1_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='A1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot2 <- y

y <- qplot(unlist(k_A1_D1_phase_dists[[3]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot3 <- y

y <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='D1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot4 <- y
#y <- ggplot_build(y)

y <- qplot(unlist(k_A1_D1_phase_dists[[5]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='A1 and D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot5 <- y

y <- qplot(unlist(k_A1_D1_phase_dists[[6]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='Only A1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot6 <- y

y <- qplot(unlist(k_A1_D1_phase_dists[[7]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='Only D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot7 <- y

y <- qplot(unlist(k_A1_D1_phase_dists[[8]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_A1_D1',main='Neither dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot8 <- y


require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol=2,nrow=4,top='Distributions of k_A1_D1 under various phases')

# D1_A1
x <- qplot(unlist(k_D1_A1_phase_dists[[1]]),geom='histogram',binwidth=number_breaks, xlim = c(0,1),xlab='k_D1_A1',main='A1 dominant')
#plot1 <- qplot(unlist(k_A1_D1_phase_dists[[1]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot1 <- x
#x <- ggplot_build(x)

y <- qplot(unlist(k_D1_A1_phase_dists[[2]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='A1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot2 <- y

y <- qplot(unlist(k_D1_A1_phase_dists[[3]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot3 <- y

y <- qplot(unlist(k_D1_A1_phase_dists[[4]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='D1 sub-dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot4 <- y
#y <- ggplot_build(y)

y <- qplot(unlist(k_D1_A1_phase_dists[[5]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='A1 and D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot5 <- y

y <- qplot(unlist(k_D1_A1_phase_dists[[6]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='Only A1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot6 <- y

y <- qplot(unlist(k_D1_A1_phase_dists[[7]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='Only D1 dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot7 <- y

y <- qplot(unlist(k_D1_A1_phase_dists[[8]]),geom='histogram',binwidth=number_breaks,xlim = c(0,1),xlab='k_D1_A1',main='Neither dominant')
#plot2 <- qplot(unlist(k_A1_D1_phase_dists[[4]]),geom='density',adjust=adjust_level, xlim = c(0,1))
plot8 <- y


require(gridExtra)
grid.arrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, ncol=2,nrow=4,top='Distributions of k_D1_A1 under various phases')



##### JS Divergence

x <- x$data[[1]][[1]]/sum(x$data[[1]][[1]])
y <- y$data[[1]][[1]]/sum(y$data[[1]][[1]])

# is this reasonable?
p <- x + 0.00000001
q <- y + 0.00000001


m <- 0.5 * (p + q)
#JS <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
JS <- 0.5 * (sum(p * log2(p / m)) + sum(q * log2(q / m)))
JS


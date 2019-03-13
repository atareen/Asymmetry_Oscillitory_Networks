########################################
#### Environment variables and data ####
########################################


# read data
setwd("/Users/tareen/Desktop/second_paper_local/2A2D_Large_Run")
data_2K2P <- read.table(file="2A2D_2p7mill.txt",sep="",fill=TRUE)
head(data_2K2P)

# rate constants
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
library("scales")
require(gridExtra)

# phase function
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

# essentiality variables
ess_A1 <- data_2K2P$V2[1:2000000]
ess_D1 <- data_2K2P$V4[1:2000000]
ess_A2 <- data_2K2P$V3[1:2000000]
ess_D2 <- data_2K2P$V5[1:2000000]

n <- length(ess_A1)-500

phase_A1 <- make.phase(ess_A1,n)
phase_D1 <- make.phase(ess_D1,n)
phase_A2 <- make.phase(ess_A2,n)
phase_D2 <- make.phase(ess_D2,n)

############################################
#### End Environment variables and data ####
############################################

#############################
#### Supplemental figure ####
#############################

k_A1_A1 <- k11[1:2000000]
k_A1_A1 <- rescale(k_A1_A1)

k_D1_D1 <- k33[1:2000000]
k_D1_D1 <- rescale(k_D1_D1)

# autoactivation vs. auto-deactivation
par(mfrow=c(2,2),mar=c(5.1,6.1,4.1,3.1))

plot(1.35*k_A1_A1[1:100000],ylim=c(0.,1),type='l',col='darkgray',ylab=expression(italic(k)[A1 %->% A1]),lwd=1.5,cex.axis=2,cex.lab=2.5,xlab="")
lines(phase_A1,col='red',lwd=2.5)

plot(2.25*k_D1_D1[1:100000],type='l',ylim=c(0.,1),col='darkgray',ylab=expression(italic(k)[D1 %->% D1]),xlab='Accepted Mutations',lwd=1.5,cex.axis=2,cex.lab=2.5)
lines(phase_D1,col='blue',lwd=2.5)


#################################
#### End Supplemental figure ####
#################################

####################################################
### Auto-activation and deactivation histograms  ###
####################################################

bins_variable <-20
plot1 <- qplot((k11-k22)/(k11+k22),bins=bins_variable,xlab=expression(frac(italic(k)[A1 %->% A1]-italic(k)[A2 %->% A2],italic(k)[A1 %->% A1]+italic(k)[A2 %->% A2])),
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

plot2 <- qplot((k33-k44)/(k33+k44),bins=bins_variable,xlab=expression(frac(italic(k)[D1 %->% D1]-italic(k)[D2 %->% D2],italic(k)[D1 %->% D1]+italic(k)[D2 %->% D2])),
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

########################################################
### End Auto-activation and deactivation histograms  ###
########################################################

#################################################################
### FIG 3: activation and deactivation asymmetry times series ###
#################################################################

activation_asym   <- (k31-k32)/(k31+k32)
deactivation_asym <- (k13-k14)/(k13+k14)

par(mfrow=c(3,1),mgp=c(3.5,1.5,0),mar = c(7.1, 8.1, 2.1, 2.1))

#plot.new()
#legend(x="bottom", ncol=4,legend=c("A Phase","D Phase","Rate Asymmetry", "Deactivator Essentiality"),
#       fill=c("red","blue","darkgray", "light green"))

#plot(activation_asym[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab=expression(frac(italic(k)[A1 %->% D1]-italic(k)[A2 %->% D1],italic(k)[A1 %->% D1]+italic(k)[A2 %->% D1])),col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
plot(ess_D1[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab="Deactivator \nEssentiality",col='lightgreen',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(-ess_D2[1:100000],type='l',col='lightgreen')
lines(phase_D1[1:100000],col="blue",type='l',lwd='2')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='2')

#plot(activation_asym[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab=expression(frac(italic(k)[A1 %->% D1]-italic(k)[A2 %->% D1],italic(k)[A1 %->% D1]+italic(k)[A2 %->% D1])),col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
plot(activation_asym[1:100000],type='l',cex.main=1.75,main="(b)",adj = 0,xaxt='n',ylab="  Activation \nAsymmetry",col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(phase_A1[1:100000],col="red",type='l',lwd='2')
lines(-phase_A2[1:100000],col="red",type='l',lwd='2')

# deactivation asymmetry plot
#plot(deactivation_asym[1:100000],type='l',cex.main=1.75,main="(b)",adj = 0,xlab='Accepted Mutations',ylab=expression(frac(italic(k)[D1 %->% A1]-italic(k)[D2 %->% A1],italic(k)[D1 %->% A1]+italic(k)[D2 %->% A1])),col='darkgray',ylim=c(-1,1),cex.axis=2.25,cex.lab=2.25)
plot(deactivation_asym[1:100000],type='l',cex.main=1.75,main="(c)",adj = 0,xlab='Accepted Mutations',ylab="Deactivation \n Asymmetry",col='darkgray',ylim=c(-1,1),cex.axis=2.25,cex.lab=2.25)
lines(phase_D1[1:100000],col="blue",type='l',lwd='2')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='2')


##### Fig 3 debugging #########

#plot(activation_asym[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab=expression(frac(italic(k)[A1 %->% D1]-italic(k)[A2 %->% D1],italic(k)[A1 %->% D1]+italic(k)[A2 %->% D1])),col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
plot(ess_D1[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab="Deactivator \nEssentiality",col='lightgreen',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(-ess_D2[1:100000],type='l',col='lightgreen')
lines(phase_D1[1:100000],col="blue",type='l',lwd='2')
lines(-phase_D2[1:100000],col="blue",type='l',lwd='2')

phase_D_clean <- phase_D1
phase_A_clean <- phase_A1


plot(ess_D1[1:100000],type='l',col='lightgreen',ylim=c(-1,1))
lines(-ess_D2[1:100000],type='l',col='lightgreen')
lines(2*phase_D_clean[1:100000]-1,type='l',lwd='2')

# build phases with consistent definitions
# +1, D1 essential, -1, D2 essential, 0, both essential
  n <-100000
  phase0 <- rep(0.0, n)
  for(i in 1:n)
  {
    if(ess_D1[i]==1 & ess_D2[i]==1)
    {
      phase0[i] <- 0      
    }
    else if(ess_D1[i]==1)
    {
      phase0[i] <- 1
    }
    
    else if(ess_D2[i]==1)
    {
      phase0[i] <- -1 
    }
    
  }
  

# slightly misaligned ACTIVATOR phase 
  
plot(phase1[1:10000],type='l',col='lightgreen',lwd=1.5)
lines(2*phase_A1[1:10000]-1,type='l',col='red',lwd=2)

phase_D_clean <- 2*phase_D1[1:100000]-1



# zoomed plot
plot(phase0[94000:96000],type='l',col='lightgreen',lwd=1.5)
lines(2*phase_D1[94000:96000]-1,type='l',col='red',lwd=2)
lines(phase_D_clean[94000:96000],type='l',col='black')

# cleaning of D phase

phase_D_clean[2000:2300]<-  1       # done
phase_D_clean[3000:4000]<- -1       # done
phase_D_clean[10500:11100]<- -1     # done
phase_D_clean[35000:35550] <- 1     # done
phase_D_clean[44500:45300] <-1      # done
phase_D_clean[64500:65300] <-1      # done
phase_D_clean[72000:73000] <- 1     # done
phase_D_clean[74000:74800] <- -1    # done 
phase_D_clean[74801:75480] <-  1    # done
phase_D_clean[77450:77600] <- -1    # done
phase_D_clean[80500:81400] <- 1     # done
phase_D_clean[86000:86150] <- -1
phase_D_clean[86500:87300] <- 1
phase_D_clean[93000:94000] <- -1
phase_D_clean[95000:95850] <- 1
# above is good from 1:100000


# same process of A phase
phase1 <- rep(0.0, n)
for(i in 1:n)
{
  if(ess_A1[i]==1 & ess_A2[i]==1)
  {
    phase1[i] <- 0      
  }
  else if(ess_A1[i]==1)
  {
    phase1[i] <- 1
  }
  
  else if(ess_A2[i]==1)
  {
    phase1[i] <- -1 
  }
  
}




### cleaning of A phase

phase_A_clean <- 2*phase_A1[1:100000]-1

# full plot, misaligned
plot(activation_asym[1:100000],type='l',cex.main=1.75,main="(b)",adj = 0,ylab="  Activation \nAsymmetry",col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(phase_A_clean[1:100000],type='l',col='red',lwd=2)

plot(activation_asym[80000:90000],type='l',cex.main=1.75,main="(b)",adj = 0,ylab="  Activation \nAsymmetry",col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(phase_A_clean[80000:90000],type='l',col='red',lwd=2)


# cleaning of A phase
phase_A_clean[44500:48000] <- -1
phase_A_clean[90000:92000] <- -1
phase_A_clean[5000:5800] <- 1
phase_A_clean[10001:11500] <- -1
phase_A_clean[82800:84000] <-1
phase_A_clean[84000:85700] <- 1
# done

##### End Fig 3 debugging #####


### FIG 3 PANEL (a) new 

activation_asym   <- (k31-k32)/(k31+k32)
deactivation_asym <- (k13-k14)/(k13+k14)

par(mfrow=c(3,1),mgp=c(3.5,1.5,0),mar = c(7.1, 8.1, 2.1, 2.1))

# full plot, clean
plot(phase0[1:100000],type='l',cex.main=1.75,main="(a)",adj = 0,xaxt='n',ylab="Deactivator \nEssentiality",col='lightgreen',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25,lwd=1.25)
#lines(2*phase_D1[1:100000]-1,type='l',col='red',lwd=2)
lines(phase_D_clean[1:100000],type='l',col='blue',lwd=2)

plot(activation_asym[1:100000],type='l',cex.main=1.75,main="(b)",adj = 0,xaxt='n',ylab="  Activation \nAsymmetry",col='darkgray',ylim=c(-1,1),xlab='',cex.axis=1.75,cex.lab=2.25)
lines(phase_A_clean[1:100000],type='l',col='red',lwd=2)

plot(deactivation_asym[1:100000],type='l',cex.main=1.75,main="(c)",adj = 0,xlab='Accepted Mutations',ylab="Deactivation \n Asymmetry",col='darkgray',ylim=c(-1,1),cex.axis=2.25,cex.lab=2.25)
lines(phase_D_clean[1:100000],type='l',col='blue',lwd=2)

####################################################################
### End FIG 3 activation and deactivation asymmetry times series ###
####################################################################

#############################################################
### FIGURE 2 activation and deactivation asymmetry hists  ###
#############################################################

bin_variable <-14

margin_var_1 <- 1 # squishes from top the larger 
margin_var_2 <- 1.1
margin_var_3 <- 3 # squishes from bottom the larger
margin_var_4 <- 2


plot1 <- qplot(activation_asym,bins=bin_variable,xlab="Activation asymmetry", 
               geom="histogram", 
               ylab="Counts", 
               fill=I("grey"),
               col=I("black"))+
  ggtitle("(a)")+
  ylim(0, 4e5)+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=23,color = 'black'),
        axis.text.y = element_text(size=25,color = 'black'),
        axis.title.x = element_text(size=27,color = 'black'),
        axis.title.y = element_text(size=25,color = 'black'),
        plot.title = element_text(size=25,color='black'),
        plot.margin=unit(c(margin_var_1,margin_var_2,margin_var_3,margin_var_4), "cm")
  )

plot2 <- qplot(deactivation_asym,bins=bin_variable,xlab="Deactivation Asymmetry", 
               geom="histogram", 
               ylab="",
               fill=I("grey"),
               col=I("black"))+
              ggtitle("(b)")+
  theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=23,color = 'black'),
        axis.text.y = element_text(size=25,color = 'black'),
        axis.title.x = element_text(size=27,color = 'black'),
        axis.title.y = element_text(size=25,color = 'black'),
        plot.title = element_text(size=25,color='black'),
        plot.margin=unit(c(margin_var_1,0.9,margin_var_3,margin_var_4), "cm")
  )
grid.arrange(plot1, plot2, ncol=2)

#######################################################
### End activation and deactivation asymmetry hists ###
#######################################################

###################################################
#### Make k_A_D dist in phase and out of phase ####
###################################################

plot(phase_A1[1:10000],type='l',lwd='3')
lines(activation_asym[1:10000])

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

# in  phase k A1->D1,D2
k_A1_D12_dom <- unlist(make.distribution((k31[1:1000000]+k41[1:1000000])/2,phase_A1,1000000)[[1]])

# out phase k A1->D1,D2
k_A1_D12_sub <-unlist(make.distribution((k31[1:1000000]+k41[1:1000000])/2,phase_A1,1000000)[[2]])

# in  phase k D1->A1,A2
k_D1_A12_dom <- unlist(make.distribution((k13[1:1000000]+k23[1:1000000])/2,phase_D1,1000000)[[1]])

# out phase k D1->A1,A2
k_D1_A12_sub <-unlist(make.distribution((k13[1:1000000]+k23[1:1000000])/2,phase_D1,1000000)[[2]])

#par(mfrow=c(2,2),mgp=c(3.5,1.5,0),mar = c(7.1, 8.1, 2.1, 2.1))

#breaks_var <- 28
#hist(k_A1_D12_dom,col='gray',breaks=breaks_var)
#hist(k_A1_D12_sub,col='gray',breaks=breaks_var)

#hist(k_D1_A12_dom,col='gray',breaks=breaks_var)
#hist(k_D1_A12_sub,col='gray',breaks=breaks_var)

k_A1_D12_dom <- rescale(k_A1_D12_dom)
k_A1_D12_sub <- rescale(k_A1_D12_sub)
k_D1_A12_dom <- rescale(k_D1_A12_dom)
k_D1_A12_sub <- rescale(k_D1_A12_sub)

bin_variable <-20

margin_var_1 <- 3 # squishes from top the larger 
margin_var_2 <- 1
margin_var_3 <- 1 # squishes from bottom the larger
margin_var_4 <- 1

# plot 3 d-e
dat1 = data.frame(x=k_A1_D12_sub, group="Sub-dominant")
dat2 = data.frame(x=k_A1_D12_dom, group="Dominant")

dat = rbind(dat1, dat2)
plotde <-ggplot(dat, aes(x, fill=group, colour=group)) +
  geom_histogram(aes(y=..density..), alpha=0.7,bins=16, 
                 position="identity", lwd=0.2) +
  xlab(expression(frac(italic(k)[A1 %->% D1]+italic(k)[A1 %->% D2],2))) +
  ylab("Density") +
  ggtitle("(d)")+
  xlim(0, 100)+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=23,color = 'black'),
        axis.text.y = element_text(size=25,color = 'black'),
        axis.title.x = element_text(size=25,color = 'black'),
        axis.title.y = element_text(size=25,color = 'black'),
        plot.title = element_text(size=25,color='black'),
        legend.position = c(.75, .6),
        plot.margin=unit(c(margin_var_1,margin_var_2,margin_var_3,margin_var_4), "cm")
  )+ theme(legend.text=element_text(size=20))


# plot 3 f-g
dat1 = data.frame(x=k_D1_A12_sub, group="Sub-dominant")
dat2 = data.frame(x=k_D1_A12_dom, group="Dominant")

dat = rbind(dat1, dat2)

#ggplot(dat, aes(x, fill=group, colour=group)) +
#  geom_histogram(breaks=seq(0,200,5), alpha=0.6, 
#                 position="identity", lwd=0.2) +
#  ggtitle("Unormalized")

plotfg <- ggplot(dat, aes(x, fill=group, colour=group)) +
  geom_histogram(aes(y=..density..), alpha=0.7,bins=22, 
                 position="identity", lwd=0.2) +
  xlab(expression(frac(italic(k)[D1 %->% A1]+italic(k)[D1 %->% A2],2))) +
  ylab("Density") +
  xlim(0, 300)+
  ggtitle("(e)")+
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=23,color = 'black'),
        axis.text.y = element_text(size=25,color = 'black'),
        axis.title.x = element_text(size=25,color = 'black'),
        axis.title.y = element_text(size=25,color = 'black'),
        plot.title = element_text(size=25,color='black'),
        legend.position = c(.75, .6),
        plot.margin=unit(c(margin_var_1,margin_var_2,margin_var_3,margin_var_4), "cm")
  )+ theme(legend.text=element_text(size=20))


# 3_de
# 3_fg
grid.arrange(plotde, plotfg,ncol=2, nrow=2)

###################################################
#### End Make k_A_D dist in phase and out of phase ####
###################################################

### JS entropy between pairs of dists


# pair 1
k_A1_D12_dom 
k_A1_D12_sub 

# pair 2
k_D1_A12_dom 
k_D1_A12_sub 


# pair 1 computation

k_A1_D12_dom_norm <- k_A1_D12_dom/sum(k_A1_D12_dom)
k_A1_D12_sub_norm <- k_A1_D12_sub/sum(k_A1_D12_sub)  

# is this reasonable?
p <- k_A1_D12_dom_norm + 0.00000001
q <- k_A1_D12_sub_norm + 0.00000001

m <- 0.5 * (p + q)
#JS <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
JS <- 0.5 * (sum(p * log2(p / m)) + sum(q * log2(q / m)))
JS # 0.1850717

# pair 2 computation
k_D1_A12_dom_norm <- k_D1_A12_dom/sum(k_D1_A12_dom)
k_D1_A12_sub_norm <- k_D1_A12_sub/sum(k_D1_A12_sub)  


# is this reasonable?
p <- k_D1_A12_dom_norm + 0.00000001
q <- k_D1_A12_sub_norm + 0.00000001

m <- 0.5 * (p + q)
#JS <- 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
JS <- 0.5 * (sum(p * log2(p / m)) + sum(q * log2(q / m)))
JS


#######################################################
#### End Make k_A_D dist in phase and out of phase ####
#######################################################

####################################################
### FIGURE 4 synergy vs competitive time series  ###
####################################################

k_A1_A2 <- k21[1:2000000]
#k_A1_A2 <- rescale(k_A1_A2)
k_A1_A2 <- 2*rescale(k_A1_A2)-1

k_D1_D2 <- k43[1:2000000]
#k_D1_D2 <- rescale(k_D1_D2)
k_D1_D2 <- 2*rescale(k_D1_D2)-1

# synergy vs. competition
par(mfrow=c(2,1),mar=c(5.1,5.1,4.1,2.1))

#plot(2.5*k_A1_A2[1:100000],type='l',ylim=c(0,1),col='darkgray',ylab='k_A1_A2',xlab='Accepted Mutations', main='Phase-A1 and k_A1_A2')
#lines(phase_A1[1:100000],col='red',lwd=2)

plot(2.5*k_A1_A2[1:100000],type='l',lwd=1.5,ylim=c(0,1),col='darkgray',ylab=expression(italic(k)[A1 %->% A2]),xlab='',cex.axis=1.8,cex.lab=2.0,xaxt='n')
lines((phase_A_clean[1:100000]+1)/2,col='red',lwd=3)

plot(2.5*k_D1_D2[1:100000],type='l',lwd=1.5,ylim=c(0,1),col='darkgray',ylab=expression(italic(k)[D1 %->% D2]),xlab='Accepted Mutations',cex.axis=1.8,cex.lab=2.0)
lines((phase_D_clean[1:100000]+1)/2,col='blue',lwd=3)

########################################################
### END FIGURE 4 synergy vs competitive time series  ###
########################################################





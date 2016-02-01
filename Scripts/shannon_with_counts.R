require(dplyr)
require(vegan)
require(reshape)

#load data
counts<-as.data.frame(read.csv("data/counts.csv"), header=FALSE)
quadrats<-as.data.frame(read.csv("data/occupancy.csv"))
sitename<-quadrats$X

#Have to get rid of unnecessary columns (the diversity function will not read) 
quadrats$sty<-NULL
quadrats$site<-NULL
quadrats$tx<-NULL
quadrats$yr<-NULL
quadrats$X<-NULL
sitename
#Now we create a table of exp shannon indices and corresponding site names
cbind(levels(sitename),apply(quadrats, MARGIN=1, FUN=function(x)
    {cbind(rownames(x),exp(diversity(x, "shannon")))}))
shannlist<-cbind(levels(sitename),apply(quadrats, MARGIN=1, FUN=function(x)
    {cbind(rownames(x),exp(diversity(x, "shannon")))}))

#Creating column names for more clarity
colnames(shannlist)[1]<- "Site Name"
colnames(shannlist)[2] <- "ENS"

str(shannlist)  
class(shannlist) 
shannlist

#Statistical comparison between controls and treatments
numers_cont<-shannlist[c(1,2,3,4,9,10,11,12,17,18,19,23,24,25,29,30,31,35,
                             36,37,41,42,43,44,49,50,51,55,56,57,58),2]
numers_treat<-shannlist[c(5,6,7,8,13,14,15,16,20,21,22,26,27,28,32,33,34,
                          38,39,40,45,46,47,48,52,53,54,59,60,61,62),2]
C<-as.numeric(as.character(numers_cont))
T<-as.numeric(as.character(numers_treat))
mean(C)
sd(C)
mean(T)
sd(T)

sum(T)/sum(C)*100 #treatment sites 28.8% more diverse than control sites
t.test(C,T) #p-value = 0.0381 --> statistically significant, reject null hypothesis

#Yearly comparisons
#Year 1
y1_cont<-shannlist[c(1,9,17,23,29,35,41,49,55),2]
y1_treat<-shannlist[c(5,13,20,26,32,38,45,52,59),2]

y1C<-as.numeric(as.character(y1_cont))
y1T<-as.numeric(as.character(y1_treat))
mean(y1C)
sd(y1C)
mean(y1T)
sd(y1T)

sum(y1T)/sum(y1C)*100 
t.test(y1C,y1T) #p-value --> not statistically significant
#Year 2
y2_cont<-shannlist[c(2,10,18,24,30,36,42,50,56),2]
y2_treat<-shannlist[c(6,14,21,27,33,39,46,53,60),2]

y2C<-as.numeric(as.character(y2_cont))
y2T<-as.numeric(as.character(y2_treat))
mean(y2C)
sd(y2C)
mean(y2T)
sd(y2T)

sum(y2T)/sum(y2C)*100 
t.test(y2C,y2T) #p-value --> not statistically significant
#Year 3
y3_cont<-shannlist[c(3,11,19,25,31,37,43,51,57),2]
y3_treat<-shannlist[c(7,15,22,28,34,40,47,54,61),2]

y3C<-as.numeric(as.character(y3_cont))
y3T<-as.numeric(as.character(y3_treat))
mean(y3C)
sd(y3C)
mean(y3T)
sd(y3T)

sum(y3T)/sum(y3C)*100 
t.test(y3C,y3T) #p-value = 0.1194 --> not statistically significant
#Year 4
y4_cont<-shannlist[c(4,12,44,58),2]
y4_treat<-shannlist[c(8,16,48,62),2]

y4C<-as.numeric(as.character(y4_cont))
y4T<-as.numeric(as.character(y4_treat))
mean(y4C)
sd(y4C)
mean(y4T)
sd(y4T)

sum(y4T)/sum(y4C)*100 
t.test(y4C,y4T) #p-value = 0.7949 --> not statistically significant



#Now to create graphs (definitely not the most efficient way to create them, but hey it works)

#Control Graph
x<-1:4
x2<-1:3
AB<-shannlist[,2][1:4] 
AL<-shannlist[,2][9:12] 
BW<-shannlist[,2][17:19] 
DR<-shannlist[,2][23:25]
FH<-shannlist[,2][29:31]
MA<-shannlist[,2][35:37]
RO<-shannlist[,2][41:44]
SH<-shannlist[,2][49:51]
URWA<-shannlist[,2][55:58]

plot(x,AB, col="red", xlim = c(1,4), ylim = c(0,22), main="Control Sites",
     xlab="Year", ylab="Exponential of Shannon-Wiener \nDiversity Index", 
     pch=".", cex=8, type="o", lwd=3, xaxt="n") 
points(x,AL, col="green", pch=".", cex=8, type="o", lwd=3)
points(x2,BW, col="blue", pch=".", cex=8, type="o", lwd=3)
points(x2,DR, col="purple", pch=".", cex=8, type="o", lwd=3)
points(x2,FH, col="midnightblue", pch=".", cex=8, type="o", lwd=3)
points(x2,MA, col="gold", pch=".", cex=8, type="o", lwd=3)
points(x,RO, col="deeppink",pch=".", cex=8, type="o", lwd=3)
points(x2,SH, col="yellowgreen", pch=".", cex=8, type="o", lwd=3)
points(x,URWA, col="darkorange1", pch=".", cex=8, type="o", lwd=3)
axis(side=1, at=c(1,2,3,4))

#Treatment Graph
x<-1:4
x2<-1:3
AB<-shannlist[,2][5:8] 
AL<-shannlist[,2][13:16] 
BW<-shannlist[,2][20:22] 
DR<-shannlist[,2][26:28]
FH<-shannlist[,2][32:34]
MA<-shannlist[,2][38:40]
RO<-shannlist[,2][45:48]
SH<-shannlist[,2][52:54]
URWA<-shannlist[,2][59:62]

plot(x,AB, col="red", xlim = c(1,4), ylim = c(0,22), main="Treatment Sites",
     xlab="Year", ylab="Exponential of Shannon-Wiener \nDiversity Index", 
     pch=".", cex=8, type="o", lwd=3, xaxt="n") 
points(x,AL, col="green", pch=".", cex=8, type="o", lwd=3)
points(x2,BW, col="blue", pch=".", cex=8, type="o", lwd=3)
points(x2,DR, col="purple", pch=".", cex=8, type="o", lwd=3)
points(x2,FH, col="midnightblue", pch=".", cex=8, type="o", lwd=3)
points(x2,MA, col="gold", pch=".", cex=8, type="o", lwd=3)
points(x,RO, col="deeppink",pch=".", cex=8, type="o", lwd=3)
points(x2,SH, col="yellowgreen", pch=".", cex=8, type="o", lwd=3)
points(x,URWA, col="darkorange1", pch=".", cex=8, type="o", lwd=3)
axis(side=1, at=c(1,2,3,4))

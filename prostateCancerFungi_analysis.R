setwd("C:/Microbime/Prostate cancer/fungi")
#install.packages("RSvgDevice")
library(vegan)
library(ecodist)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(scatterplot3d)
library(lattice)
library(reshape2)
library(RSvgDevice)
library(reshape)
library(indicspecies)
library(devEMF)
library(tidyr)
library(car)
library(Rmisc)
library(ggplot2)
#install.packages("car", dependencies=TRUE, INSTALL_opts = c('--no-lock'))


#remove the previous list
rm(list=ls())

#### first, lets read in data clinical and bacteria microbiome

otu.raw <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.OTU.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw)
otu.name <- otu.raw[,1]
#find the index of the column name '9715'
index.9715 <- grep("9715", colnames(otu.raw))
index.9992 <- grep("9992", colnames(otu.raw))
index.PP <- grep("PP", colnames(otu.raw))
index.PT <- grep("PT", colnames(otu.raw))
index.DA1 <- grep("DA1", colnames(otu.raw))
index.DA40 <- grep("DA40", colnames(otu.raw))

#get all otu counts
otu.raw <- otu.raw[,-c(1,2,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw)
colSums (otu.raw, na.rm = FALSE, dims = 1) 
dim(otu.raw)
#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw))
index.HC22 <- grep("HC22", colnames(otu.raw))
index.HC23 <- grep("HC23", colnames(otu.raw))
#get male otu counts
otu.raw <- otu.raw[,-c(index.HC14,index.HC22,index.HC23)]

#patient only
patients.info <- read.csv("patientsInfo.csv",sep=",", header=T,check.names=FALSE)
patients.info <- rename(patients.info,c("Patient_ID"="Sample_ID"))
patients.info <- patients.info[,-3] %>% data.frame()
patients.info$Sample_ID <- patients.info$Sample_ID %>% as.character()

#combine patients and health male
sample.all <- rbind(health.male.info,patients.info)
write.csv(sample.all, file="results/test.csv")

#get data frame
otu.raw <- as.data.frame(t(otu.raw))
dim(otu.raw)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw-negmean
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)



# calculate diversity
otu$simpson <- diversity(otu[,1:41],"simp")
otu$shannon <- diversity(otu[,1:41], index = "shannon")
otu$specnumb <- specnumber(otu[,1:41], MARGIN = 1)
otu$Sample_ID <- rownames(otu)
#typeof(otu$Sample_ID)

a.div <- otu[,c("Sample_ID","simpson","shannon","specnumb")]
write.csv(a.div, file="results/AlfaDiversity.csv")
#add group to a.div
a.div.all <- inner_join(a.div,sample.all[,c("Sample_ID","Group1")])

#a.div$Sample_ID <- a.div$Sample_ID %>% as.character()
#typeof(a.div$Sample_ID)
# calculate OTU percentage
otu.per <- prop.table(as.matrix(otu[,(1:41)]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/OtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#row.names(otu.per) <- a.div$Sample_ID




#### second, lets read in patient information data
#info <- read.csv("Samples information.csv",sep=",", header=T,check.names=FALSE)
#info <- info[,1:7]
#info$Sample_ID %>% as.character()
#info$Sample_ID <- as.character(info$Sample_ID)
#inner_john by Sample_ID
L <- inner_join(a.div[,c("Sample_ID","simpson","shannon","specnumb")],sample.all)
write.csv(L, file="results/test.csv")
#L[, 1] <- sapply(L[, 1], as.numeric)
## Beta diversity
## combine diversity data and percentage data

otu.per$Sample_ID <- row.names(otu.per)
otu.all <- inner_join(otu.per,L)
rownames(otu.all) <- otu.all$Sample_ID
write.csv(otu.all, file="results/test.csv")

### beta diversity
colnames(otu.all)
#discard the not numeric 
otu.all <- otu.all[,-c(42:49),drop=TRUE] 
otu.all <- subset(otu.all, rownames(otu.all)!=rownames(otu.all[apply(otu.all, MARGIN=1, sum)==0,]))
write.csv(otu.all, file="results/test.csv")
#delete NAs
otu.all <- otu.all[-60,]

# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA
# makes a gradient from red to blue
my.colors <- colorRampPalette(c('red','blue'))(10)

plot(d.chiq[,1], d.chiq[,2], col=my.colors, cex=3, pch=16)

#otu.all <- otu.all[,-c(42:45),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#L2.update <- rowSums(L2[1:61,] != 0) %>% data.frame()

#otu.all <- subset(otu.all, rownames(otu.all)!=rownames(otu.all[apply(otu.all, MARGIN=1, sum)==0,]))
#write.csv(otu.all, file="results/test.csv")
#delete NA in rows
#dim(otu.all)
#otu.all <- otu.all[-60,]






#Simpson box plots with dots
p <- ggplot(L, aes(x=Group1, y=simpson)) + labs(title="Simpson", y = "Simpson Diversity Index")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot()
p

# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))

#Shannon box plots with dots
p <- ggplot(L, aes(x=Group1, y=shannon)) + labs(title="Shannon", y = "Shannon Diversity Index")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot()
p

# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# Box plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))

#add p value
tg <- L
dlply(tg, .(Group1), summarise,
      p_value = (leveneTest(shannon ~ Group1 , center = mean))$`Pr(>F)` )
#??????dose???,Levene??????p????????????0.052???0.149???0.129,??????????????????????????????(var.equal = TRUE)???

p +???stat_compare_means(method = "t.test", method.args = list(var.equal = TRUE))

# u test
wilcox.test(L$specnumb~L$Group1)
wilcox.test(L$simpson~L$Group1)
wilcox.test(L$shannon~L$Group1)

# t test for simpson between groups
leveneTest(L$simpson,L$Group1)
t.test(L$simpson~L$Group1)
summary(lm(L$simpson~L$Group1))

### 1.2 beta diversity
otu.per$Sample_ID <- row.names(otu.per)
otu.all <- inner_join(otu.per, a.div)
dim(otu.all)

rownames(otu.all) <- otu.all$Sample_ID
write.csv(otu.all, file="results/test.csv")
colnames(otu.all)
map <- otu.all[otu.all$Sample_ID,]
#discard the not numeric 
otu.all <- otu.all[,-c(11:14),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#otu.rowsums.gt0 <- rowSums(otu.all[1:61,] != 0)%>% data.frame()
otu.all <- subset(otu.all, rownames(otu.all)!=rownames(otu.all[apply(otu.all, MARGIN=1, sum)==0,]))
#delete NA in rows
#rownames(otu.all)
#otu.all <- otu.all[-60,]

### beta diversity
# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA

phylum.color<-c("basidiomycota" = "#FFBBFF", 
                "ascomycota" = "#8A2BE2") 
plot(d.chiq[,1], d.chiq[,2], col=phylum.color, cex=3, pch=16)

## phylum otu and percentage
otu.raw.phylum <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.phylum.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw.phylum)
otu.name <- otu.raw.phylum[,1]


#find the index of patients
index.9715 <- grep("9715", colnames(otu.raw.phylum))
index.9992 <- grep("9992", colnames(otu.raw.phylum))
index.PP <- grep("PP", colnames(otu.raw.phylum))
index.PT <- grep("PT", colnames(otu.raw.phylum))
index.DA1 <- grep("DA1", colnames(otu.raw.phylum))
index.DA40 <- grep("DA40", colnames(otu.raw.phylum))

#get all otu counts
otu.raw.phylum <- otu.raw.phylum[,-c(1,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw.phylum)
colSums (otu.raw.phylum, na.rm = FALSE, dims = 1) 
dim(otu.raw.phylum)

#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw.phylum))
index.HC22 <- grep("HC22", colnames(otu.raw.phylum))
index.HC23 <- grep("HC23", colnames(otu.raw.phylum))
#get male otu counts
otu.raw.phylum <- otu.raw.phylum[,-c(index.HC14,index.HC22,index.HC23)]

#get data frame
otu.raw.phylum <- as.data.frame(t(otu.raw.phylum))
dim(otu.raw.phylum)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw.phylum[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw.phylum-negmean
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)
colnames(otu) <- otu.name

# calculate phylum OTU percentage
otu.per <- prop.table(as.matrix(otu[(1:61),]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/PhylumOtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#delete NA
otu.per[is.na(otu.per)] <- 0




#phylum percentage figure
### pt vs. healthy
otu.per$Sample_ID <- rownames(otu.per)
L1 <- inner_join(otu.per, a.div.all[,c("Sample_ID", "Group1")])

#loop through the column of the phylum table, testing each one
otu.per <- L1[,-c(3:4)]
pvals <- numeric(ncol(otu.per))
names(pvals) <- colnames(otu.per)
for (i in 1:ncol(otu.per)){
  fit <- lm(otu.per[,i] ~ L1$Group1)
  pvals[i] <- anova(fit)['L1$Group1','Pr(>F)']
}

# PDR: False discovery rate
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]

#L1$group <- revalue(L1$group, c("Patient"="HIV Patient"))
colMax <- function(data) sapply(data, max, na.rm=TRUE)
test <- cbind(sum=colSums(L1[,1:58]),mean=colMeans(L1[,1:58]),max=colMax(L1[,1:58]))%>%data.frame()

test <- data.frame(sapply(test, function(x) as.numeric(as.character(x))))
test$order <- colnames(L1[,1:58])
# classify any bacteria with a mean < 3% to "other bacteria"
colSums(test[,1:3])
test$order.other <-test$order
test$order.other[test$mean<=3] <- "bacteria_other"
test2 <- test[test$mean<=3,]
L1$bacteria_other <- rowSums(L4[,(colnames(L4)%in% test2$order)])
str(test2$order)
L1 <-L4[,-which(names(L4)%in%test2$order)]
str(L4)

phylum.level<- melt(L1, id.vars=c("Group1","Sample_ID"))


phylum.color<-c("basidiomycota" = "#FFBBFF", 
               "ascomycota" = "#8A2BE2")

phylum.level.2 <- aggregate(value~variable+Group1, data=phylum.level, FUN=mean) 

figureA <- ggplot(phylum.level.2, aes(x = factor(Group1), y = value)) +scale_fill_manual(values=phylum.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figureA


## class otu and percentage
otu.raw.class <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.class.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw.class)
otu.name <- otu.raw.class[,1]


#find the index of patients
index.9715 <- grep("9715", colnames(otu.raw.class))
index.9992 <- grep("9992", colnames(otu.raw.class))
index.PP <- grep("PP", colnames(otu.raw.class))
index.PT <- grep("PT", colnames(otu.raw.class))
index.DA1 <- grep("DA1", colnames(otu.raw.class))
index.DA40 <- grep("DA40", colnames(otu.raw.class))

#get all otu counts
otu.raw.class <- otu.raw.class[,-c(1,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw.class)
colSums (otu.raw.class, na.rm = FALSE, dims = 1) 
dim(otu.raw.class)

#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw.class))
index.HC22 <- grep("HC22", colnames(otu.raw.class))
index.HC23 <- grep("HC23", colnames(otu.raw.class))
#get male otu counts
otu.raw.class <- otu.raw.class[,-c(index.HC14,index.HC22,index.HC23)]

#get data frame
otu.raw.class <- as.data.frame(t(otu.raw.class))
dim(otu.raw.class)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw.class[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw.class-negmean
dim(otu2)
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)
colnames(otu) <- otu.name



# calculate class OTU percentage
otu.per <- prop.table(as.matrix(otu[(1:61),]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/ClassOtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#delete NA
otu.per[is.na(otu.per)] <- 0





#class percentage figure
### pt vs. healthy
otu.per$Sample_ID <- rownames(otu.per)
L2 <- inner_join(otu.per, a.div.all[,c("Sample_ID", "Group1")])
#L1$group <- revalue(L1$group, c("Patient"="HIV Patient"))

### 1.2 beta diversity


rownames(L2) <- L2$Sample_ID
write.csv(L2, file="results/test.csv")
colnames(L2)
#map <- otu.all[L2$Sample_ID,]
#discard the not numeric 
L2 <- L2[,-c(11,12),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#L2.update <- rowSums(L2[1:61,] != 0) %>% data.frame()

otu.all <- subset(L2, rownames(L2)!=rownames(L2[apply(L2, MARGIN=1, sum)==0,]))

#delete NA in rows



### beta diversity
# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA

class.color<-c("agaricomycetes" = "blue", 
               "eurotiomycetes" = "darkmagenta",
               "dothideomycetes" = "darkorange",
               "sordariomycetes" = "gold",
               "saccharomycetes" = "pink",
               "pezizomycetes" = "#458B00",
               "microbotryomycetes" = "#FF69B4",
               "tremellomycetes" = "#FF4040",
               "malasseziomycetes" = "darkolivegreen1",
               "leotiomycetes" = "#8A8A8A") 
plot(d.chiq[,1], d.chiq[,2], col=class.color, cex=3, pch=16)

#loop through the column of the class table, testing each one
colnames(L2)
otu.per <- L2[,-c(11:12)]
pvals <- numeric(ncol(otu.per))
names(pvals) <- colnames(otu.per)
for (i in 1:ncol(otu.per)){
  fit <- lm(otu.per[,i] ~ L2$Group1)
  pvals[i] <- anova(fit)['L2$Group1','Pr(>F)']
}

# PDR: False discovery rate
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]

class.level<- melt(L2, id.vars=c("Group1","Sample_ID"))


class.color<-c("agaricomycetes" = "blue", 
                "eurotiomycetes" = "darkmagenta",
                "dothideomycetes" = "darkorange",
                "sordariomycetes" = "gold",
                "saccharomycetes" = "pink",
                "pezizomycetes" = "#458B00",
                "microbotryomycetes" = "#FF69B4",
                "tremellomycetes" = "#FF4040",
                "malasseziomycetes" = "darkolivegreen1",
                "leotiomycetes" = "#8A8A8A")

class.level.2 <- aggregate(value~variable+Group1, data=class.level, FUN=mean) 

figureB <- ggplot(class.level.2, aes(x = factor(Group1), y = value)) +scale_fill_manual(values=class.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figureB

#order
## order otu and percentage
otu.raw.order <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.order.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw.order)
otu.name <- otu.raw.order[,1]


#find the index of patients
index.9715 <- grep("9715", colnames(otu.raw.order))
index.9992 <- grep("9992", colnames(otu.raw.order))
index.PP <- grep("PP", colnames(otu.raw.order))
index.PT <- grep("PT", colnames(otu.raw.order))
index.DA1 <- grep("DA1", colnames(otu.raw.order))
index.DA40 <- grep("DA40", colnames(otu.raw.order))

#get all otu counts
otu.raw.order <- otu.raw.order[,-c(1,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw.order)
colSums (otu.raw.order, na.rm = FALSE, dims = 1) 
dim(otu.raw.order)

#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw.order))
index.HC22 <- grep("HC22", colnames(otu.raw.order))
index.HC23 <- grep("HC23", colnames(otu.raw.order))
#get male otu counts
otu.raw.order <- otu.raw.order[,-c(index.HC14,index.HC22,index.HC23)]

#get data frame
otu.raw.order <- as.data.frame(t(otu.raw.order))
dim(otu.raw.order)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw.order[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw.order-negmean
dim(otu2)
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)
colnames(otu) <- otu.name

# calculate class OTU percentage
otu.per <- prop.table(as.matrix(otu[(1:61),]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/OrderOtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#delete NA
otu.per[is.na(otu.per)] <- 0




#class percentage figure
### pt vs. healthy
otu.per$Sample_ID <- rownames(otu.per)
L3 <- inner_join(otu.per, a.div.all[,c("Sample_ID", "Group1")])
#L1$group <- revalue(L1$group, c("Patient"="HIV Patient"))

#loop through the column of the order table, testing each one
colnames(L3)
otu.per <- L3[,-c(16:17)]
pvals <- numeric(ncol(otu.per))
names(pvals) <- colnames(otu.per)
for (i in 1:ncol(otu.per)){
  fit <- lm(otu.per[,i] ~ L3$Group1)
  pvals[i] <- anova(fit)['L3$Group1','Pr(>F)']
}

# PDR: False discovery rate
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]
#capnodiales 1.2x10-5

# order lever Beta diversity
L3 <- L3[,-c(16,17),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#L2.update <- rowSums(L2[1:61,] != 0) %>% data.frame()

otu.all <- subset(L3, rownames(L3)!=rownames(L3[apply(L3, MARGIN=1, sum)==0,]))

#delete NA in rows



### beta diversity
# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA

order.color<-c("agaricales" = "blue", 
               "hypocreales" = "darkmagenta",
               "saccharomycetales" = "darkorange",
               "filobasidiales" = "gold",
               "sporidiobolales" = "pink",
               "helotiales" = "#458B00",
               "capnodiales" = "#FF69B4",
               "onygenales" = "#FF4040",
               "xylariales" = "darkolivegreen1",
               "eurotiales" = "#8A8A8A",
               "pleosporales" = "#53868B",
               "pezizales" = "#000080",
               "tremellales" = "dodgerblue",
               "chaetothyriales" = "#B8860B",
               "malasseziales" = "#0000FF") 
plot(d.chiq[,1], d.chiq[,2], col=order.color, cex=3, pch=16)

order.level<- melt(L3, id.vars=c("Group1","Sample_ID"))


order.color<-c("agaricales" = "blue", 
               "hypocreales" = "darkmagenta",
               "saccharomycetales" = "darkorange",
               "filobasidiales" = "gold",
               "sporidiobolales" = "pink",
               "helotiales" = "#458B00",
               "capnodiales" = "#FF69B4",
               "onygenales" = "#FF4040",
               "xylariales" = "darkolivegreen1",
               "eurotiales" = "#8A8A8A",
               "pleosporales" = "#53868B",
               "pezizales" = "#000080",
               "tremellales" = "dodgerblue",
               "chaetothyriales" = "#B8860B",
               "malasseziales" = "#0000FF")

order.level.2 <- aggregate(value~variable+Group1, data=order.level, FUN=mean) 

figureC <- ggplot(order.level.2, aes(x = factor(Group1), y = value)) +scale_fill_manual(values=order.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figureC

#family
## family otu and percentage
otu.raw.family <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.famliy.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw.family)
otu.name <- otu.raw.family[,1]


#find the index of patients
index.9715 <- grep("9715", colnames(otu.raw.family))
index.9992 <- grep("9992", colnames(otu.raw.family))
index.PP <- grep("PP", colnames(otu.raw.family))
index.PT <- grep("PT", colnames(otu.raw.family))
index.DA1 <- grep("DA1", colnames(otu.raw.family))
index.DA40 <- grep("DA40", colnames(otu.raw.family))

#get all otu counts
otu.raw.family <- otu.raw.family[,-c(1,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw.family)
colSums (otu.raw.family, na.rm = FALSE, dims = 1) 
dim(otu.raw.family)

#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw.family))
index.HC22 <- grep("HC22", colnames(otu.raw.family))
index.HC23 <- grep("HC23", colnames(otu.raw.family))
#get male otu counts
otu.raw.family <- otu.raw.family[,-c(index.HC14,index.HC22,index.HC23)]

#get data frame
otu.raw.family <- as.data.frame(t(otu.raw.family))
dim(otu.raw.family)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw.family[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw.family-negmean
dim(otu2)
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)
colnames(otu) <- otu.name

# calculate class OTU percentage
otu.per <- prop.table(as.matrix(otu[(1:61),]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/FamilyOtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#delete NA
otu.per[is.na(otu.per)] <- 0




#famlily percentage figure
### pt vs. healthy
otu.per$Sample_ID <- rownames(otu.per)
L4 <- inner_join(otu.per, a.div.all[,c("Sample_ID", "Group1")])
#L1$group <- revalue(L1$group, c("Patient"="HIV Patient"))

#loop through the column of the family table, testing each one
colnames(L4)
otu.per <- L4[,-c(21:22)]
pvals <- numeric(ncol(otu.per))
names(pvals) <- colnames(otu.per)
for (i in 1:ncol(otu.per)){
  fit <- lm(otu.per[,i] ~ L4$Group1)
  pvals[i] <- anova(fit)['L4$Group1','Pr(>F)']
}

# PDR: False discovery rate
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]
# cladosporiaceae  pyronemataceae  filobasidiales   didymellaceae 
# P value = 0.0003, 0.05,0.05 and 0.05


L4 <- L4[,-c(21,22),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#L2.update <- rowSums(L2[1:61,] != 0) %>% data.frame()

otu.all <- subset(L4, rownames(L4)!=rownames(L4[apply(L4, MARGIN=1, sum)==0,]))

#delete NA in rows



### beta diversity
# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA


family.color<-c("arthrodermataceae" = "blue", 
                "pyronemataceae" = "darkmagenta",
                "filobasidiales" = "darkorange",
                "aspergillaceae" = "gold",
                "malasseziaceae" = "pink",
                "xylariaceae" = "#458B00",
                "nectriaceae" = "#FF69B4",
                "gymnoascaceae" = "#FF4040",
                "pleosporaceae" = "darkolivegreen1",
                "trichomeriaceae" = "#8A8A8A",
                "tricholomataceae" = "#53868B",
                "saccharomycetaceae" = "#000080",
                "helotiaceae" = "dodgerblue",
                "hydnangiaceae" = "#B8860B" ,
                "cladosporiaceae" = "#00FF00",
                "thermoascaceae" = "#C0FF3E",
                "clavicipitaceae" = "#FFB90F",
                "thermoascaceae" = "#A52A2A", 
                "sporidiobolales" = "#00FFFF",
                "didymellaceae" = "#C0FF3E")
plot(d.chiq[,1], d.chiq[,2], col=family.color, cex=3, pch=16)

family.level<- melt(L4, id.vars=c("Group1","Sample_ID"))


family.color<-c("arthrodermataceae" = "blue", 
               "pyronemataceae" = "darkmagenta",
               "filobasidiales" = "darkorange",
               "aspergillaceae" = "gold",
               "malasseziaceae" = "pink",
               "xylariaceae" = "#458B00",
               "nectriaceae" = "#FF69B4",
               "gymnoascaceae" = "#FF4040",
               "pleosporaceae" = "darkolivegreen1",
               "trichomeriaceae" = "#8A8A8A",
               "tricholomataceae" = "#53868B",
               "saccharomycetaceae" = "#000080",
               "helotiaceae" = "dodgerblue",
               "hydnangiaceae" = "#B8860B" ,
               "cladosporiaceae" = "#00FF00",
               "thermoascaceae" = "#C0FF3E",
               "clavicipitaceae" = "#FFB90F",
               "thermoascaceae" = "#A52A2A", 
               "sporidiobolales" = "#00FFFF",
               "didymellaceae" = "#C0FF3E")

family.level.2 <- aggregate(value~variable+Group1, data=family.level, FUN=mean) 

figureD <- ggplot(family.level.2, aes(x = factor(Group1), y = value)) +scale_fill_manual(values=family.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figureD

#genus
## genus otu and percentage
otu.raw.genus <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.genus.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw.genus)
otu.name <- otu.raw.genus[,1]


#find the index of patients
index.9715 <- grep("9715", colnames(otu.raw.genus))
index.9992 <- grep("9992", colnames(otu.raw.genus))
index.PP <- grep("PP", colnames(otu.raw.genus))
index.PT <- grep("PT", colnames(otu.raw.genus))
index.DA1 <- grep("DA1", colnames(otu.raw.genus))
index.DA40 <- grep("DA40", colnames(otu.raw.genus))

#get all otu counts
otu.raw.genus <- otu.raw.genus[,-c(1,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw.genus)
colSums (otu.raw.genus, na.rm = FALSE, dims = 1) 
dim(otu.raw.genus)

#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw.genus))
index.HC22 <- grep("HC22", colnames(otu.raw.genus))
index.HC23 <- grep("HC23", colnames(otu.raw.genus))
#get male otu counts
otu.raw.genus <- otu.raw.genus[,-c(index.HC14,index.HC22,index.HC23)]

#get data frame
otu.raw.genus <- as.data.frame(t(otu.raw.genus))
dim(otu.raw.genus)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw.genus[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw.genus-negmean
dim(otu2)
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)
colnames(otu) <- otu.name

# calculate class OTU percentage
otu.per <- prop.table(as.matrix(otu[(1:61),]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/GenusOtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#delete NA
otu.per[is.na(otu.per)] <- 0




#genus percentage figure
### pt vs. healthy
otu.per$Sample_ID <- rownames(otu.per)
L5 <- inner_join(otu.per, a.div.all[,c("Sample_ID", "Group1")])
#L1$group <- revalue(L1$group, c("Patient"="HIV Patient"))

#loop through the column of the genus table, testing each one
colnames(L5)
otu.per <- L5[,-c(22:23)]
pvals <- numeric(ncol(otu.per))
names(pvals) <- colnames(otu.per)
for (i in 1:ncol(otu.per)){
  fit <- lm(otu.per[,i] ~ L5$Group1)
  pvals[i] <- anova(fit)['L5$Group1','Pr(>F)']
}

# PDR: False discovery rate
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]
#cladosporium p value=0.0002
#genus beta diversity
L5 <- L5[,-c(22,23),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#L2.update <- rowSums(L2[1:61,] != 0) %>% data.frame()

otu.all <- subset(L5, rownames(L5)!=rownames(L5[apply(L5, MARGIN=1, sum)==0,]))

#delete NA in rows



### beta diversity
# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA


genus.color<-c("phoma" = "blue", 
               "trichophyton" = "darkmagenta",
               "cryptococcus" = "darkorange",
               "cladosporium" = "gold",
               "aspergillus" = "pink",
               "bipolaris" = "#458B00",
               "mycena" = "#FF69B4",
               "alternaria" = "#FF4040",
               "pyronema" = "darkolivegreen1",
               "metarhizium" = "#8A8A8A",
               "cudoniella" = "#53868B",
               "saccharomyces" = "#000080",
               "byssochlamys" = "dodgerblue",
               "gymnoascus" = "#B8860B" ,
               "rhodotorula" = "#00FF00",
               "knufia" = "#C0FF3E",
               "malassezia" = "#FFB90F",
               "fusarium" = "#A52A2A", 
               "xylaria" = "#00FFFF",
               "hydnangium" = "#C0FF3E",
               "penicillium" = "darkred")
plot(d.chiq[,1], d.chiq[,2], col=genus.color, cex=3, pch=16)


genus.level<- melt(L5, id.vars=c("Group1","Sample_ID"))


genus.color<-c("phoma" = "blue", 
                "trichophyton" = "darkmagenta",
                "cryptococcus" = "darkorange",
                "cladosporium" = "gold",
                "aspergillus" = "pink",
                "bipolaris" = "#458B00",
                "mycena" = "#FF69B4",
                "alternaria" = "#FF4040",
                "pyronema" = "darkolivegreen1",
                "metarhizium" = "#8A8A8A",
                "cudoniella" = "#53868B",
                "saccharomyces" = "#000080",
                "byssochlamys" = "dodgerblue",
                "gymnoascus" = "#B8860B" ,
                "rhodotorula" = "#00FF00",
                "knufia" = "#C0FF3E",
                "malassezia" = "#FFB90F",
                "fusarium" = "#A52A2A", 
                "xylaria" = "#00FFFF",
                "hydnangium" = "#C0FF3E",
                "penicillium" = "darkred")

genus.level.2 <- aggregate(value~variable+Group1, data=genus.level, FUN=mean) 

figureE <- ggplot(genus.level.2, aes(x = factor(Group1), y = value)) +scale_fill_manual(values=genus.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figureE

#species
## species otu and percentage
otu.raw.species <- read.csv(file="010318WJ515F-complete-pr.fasta.otus.fa.species.csv", sep=",", header=T,check.names=FALSE)
dim(otu.raw.species)
otu.name <- otu.raw.species[,1]


#find the index of patients
index.9715 <- grep("9715", colnames(otu.raw.species))
index.9992 <- grep("9992", colnames(otu.raw.species))
index.PP <- grep("PP", colnames(otu.raw.species))
index.PT <- grep("PT", colnames(otu.raw.species))
index.DA1 <- grep("DA1", colnames(otu.raw.species))
index.DA40 <- grep("DA40", colnames(otu.raw.species))

#get all otu counts
otu.raw.species <- otu.raw.species[,-c(1,index.DA1:index.DA40,index.PP,index.PT)]
dim(otu.raw.species)
colSums (otu.raw.species, na.rm = FALSE, dims = 1) 
dim(otu.raw.species)

#find male sample_ID in health
health.info <- read.csv("healthInfo.csv",sep=",", header=T,check.names=FALSE)
health.female.info <- health.info[which (health.info$gender=='F'),]
health.female.patient_ID <- health.female.info$Sample_ID
health.male.info <- health.info[which (health.info$gender=='M'),]
health.male.info <- health.male.info[,-c(3,4)]
health.male.info <- rename(health.male.info,c("Group_1"="Group1","age"="Age","gender"="Gender","race"="Race")) 
health.male.info$Sample_ID <- health.male.info$Sample_ID %>% as.character()
#captalize race
health.male.info$Race <- toupper(health.male.info$Race)
#find the index of HC14,HC22,HC23 (female health)
index.HC14 <- grep("HC14", colnames(otu.raw.species))
index.HC22 <- grep("HC22", colnames(otu.raw.species))
index.HC23 <- grep("HC23", colnames(otu.raw.species))
#get male otu counts
otu.raw.species <- otu.raw.species[,-c(index.HC14,index.HC22,index.HC23)]

#get data frame
otu.raw.species <- as.data.frame(t(otu.raw.species))
dim(otu.raw.species)
#sample value mins negative control means
#remove the absolute values in negative control
neg <- otu.raw.species[(62:64),]
negmean <- colMeans(neg)
negmean <- t(negmean)
otu2 <- otu.raw.species-negmean
dim(otu2)
#remove three rows water values
otu <- otu2[-(62:64),]
#replace negative value with zero
otu[otu < 0] <- 0 
dim(otu)
colnames(otu) <- otu.name

# calculate species OTU percentage
otu.per <- prop.table(as.matrix(otu[(1:61),]), 1) %>% data.frame()
colnames(otu.per) <- otu.name
write.csv(otu.per, file="results/SpeciesOtuPercentage.csv")
otu.per <- otu.per*100
rowSums(otu.per)
#delete NA
otu.per[is.na(otu.per)] <- 0




#genus percentage figure
### pt vs. healthy
otu.per$Sample_ID <- rownames(otu.per)
L6 <- inner_join(otu.per, a.div.all[,c("Sample_ID", "Group1")])
#L1$group <- revalue(L1$group, c("Patient"="HIV Patient"))

colnames(L6)
otu.per <- L6[,-c(28:29)]
pvals <- numeric(ncol(otu.per))
names(pvals) <- colnames(otu.per)
for (i in 1:ncol(otu.per)){
  fit <- lm(otu.per[,i] ~ L6$Group1)
  pvals[i] <- anova(fit)['L6$Group1','Pr(>F)']
}

# PDR: False discovery rate
qvals <- p.adjust(pvals,'fdr')
sort(qvals)[1:10]
# cladosporium cladosporioides,cryptococcus ater p value=0.00001 and 0.03
#species beta diversity
L6 <- L6[,-c(28,29),drop=TRUE] 

#write.csv(otu.all, file="results/test.csv")

#delete the rows that rowsums<0
#L2.update <- rowSums(L2[1:61,] != 0) %>% data.frame()

otu.all <- subset(L6, rownames(L6)!=rownames(L6[apply(L6, MARGIN=1, sum)==0,]))

#delete NA in rows



### beta diversity
# get Eucliden distance
d.euc <- dist(otu.all)
# get Bray-Curtis distance
d.bray <- vegdist(otu.all)
# # get Chi-square distance
my.ca <- cca(otu.all)
d.chiq <- as.matrix(dist(my.ca$CA$u[,1:2]))

#Plot Chi-square distances with gradient colors
# Plot Chi-square PCoA


species.color<-c("knufia capronia peltigerae" = "blue", 
                 "alternaria crivellia papaveracea" = "darkmagenta",
                 "aspergillus flavus" = "darkorange",
                 "saccharomyces cerevisiae" = "gold",
                 "cryptococcus ater" = "pink",
                 "cudoniella clavus" = "#458B00",
                 "fusarium cerealis" = "#FF69B4",
                 "byssochlamys nivea" = "#FF4040",
                 "penicillium commune" = "darkolivegreen1",
                 "mycena sp." = "#8A8A8A",
                 "alternaria sp." = "#53868B",
                 "rhodotorula mucilaginosa" = "#000080",
                 "malassezia restricta" = "dodgerblue",
                 "cryptococcus dejecticola" = "#B8860B" ,
                 "cladosporium sphaerospermum" = "#00FF00",
                 "phoma herbarum" = "#C0FF3E",
                 "aspergillus eurotium rubrum" = "#FFB90F",
                 "trichophyton rubrum" = "#A52A2A", 
                 "cladosporium cladosporioides" = "#00FFFF",
                 "hydnangium sp." = "#C0FF3E",
                 "bipolaris tetramera" = "darkred",
                 "xylaria hypoxylon" = "#00FF00",
                 "pyronema domesticum" = "#00CC00",
                 "metarhizium anisopliae" = "#009900",
                 "fusarium merismoides" = "#006600",
                 "gymnoascus reesii" = "#003300",
                 "fusarium gibberella pulicaris" = "#000000")
plot(d.chiq[,1], d.chiq[,2], col=species.color, cex=3, pch=16)
species.level<- melt(L6, id.vars=c("Group1","Sample_ID"))


species.color<-c("knufia capronia peltigerae" = "blue", 
               "alternaria crivellia papaveracea" = "darkmagenta",
               "aspergillus flavus" = "darkorange",
               "saccharomyces cerevisiae" = "gold",
               "cryptococcus ater" = "pink",
               "cudoniella clavus" = "#458B00",
               "fusarium cerealis" = "#FF69B4",
               "byssochlamys nivea" = "#FF4040",
               "penicillium commune" = "darkolivegreen1",
               "mycena sp." = "#8A8A8A",
               "alternaria sp." = "#53868B",
               "rhodotorula mucilaginosa" = "#000080",
               "malassezia restricta" = "dodgerblue",
               "cryptococcus dejecticola" = "#B8860B" ,
               "cladosporium sphaerospermum" = "#00FF00",
               "phoma herbarum" = "#C0FF3E",
               "aspergillus eurotium rubrum" = "#FFB90F",
               "trichophyton rubrum" = "#A52A2A", 
               "cladosporium cladosporioides" = "#00FFFF",
               "hydnangium sp." = "#C0FF3E",
               "bipolaris tetramera" = "darkred",
               "xylaria hypoxylon" = "#00FF00",
               "pyronema domesticum" = "#00CC00",
               "metarhizium anisopliae" = "#009900",
               "fusarium merismoides" = "#006600",
               "gymnoascus reesii" = "#003300",
               "fusarium gibberella pulicaris" = "#000000")

#my_palette <- colorRampPalette(c("green", "black", "red"))(n = 1000)

species.level.2 <- aggregate(value~variable+Group1, data=species.level, FUN=mean) 

figureG <- ggplot(species.level.2, aes(x = factor(Group1), y = value)) +scale_fill_manual(values=species.color, name="Key")+geom_bar(aes(fill=variable), position="fill", stat="identity") +labs( x="", y="Relative abundance (%)")+theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+theme(legend.text=element_text(face="italic"));figureG
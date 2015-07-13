#################################################
# Author: Robin Elahi
# Date: 150710
#################################################
# Question: Does body size influence the predictive capacity of 
# development time (T) vs lifetime fecundity (mx) and 
# survivorship (lx) on intrinsic rate of growth (r)?

library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size=12))

rm(list=ls(all=TRUE)) 

dat <- read.csv("./data/tdef_data.csv", header=TRUE, na.strings="NA")
bib <- read.csv("./data/tdef_biblio.csv", header=TRUE, na.strings="NA")

head(dat)
names(bib)

bib2 <- bib %>% select(study.no, genus, species, family, order, class, 
                       subphylum, phylum, lat, long, dry.mass.mg, 
                       author1, pub.yr)
names(bib2)

dat2 <- merge(dat, bib2, by = "study.no" )
dat2$study.no <- as.factor(dat2$study.no)

head(dat2)

# create study ID
studyF <- function(x) return(as.factor(paste(x$author1, x$pub.yr, sep = "_")))
dat2$studyID <- studyF(dat2)
unique(dat2$studyID)

# get the number of temperature measurements per critter
totalN <- dat2 %>% group_by(study.no) %>% 
  summarise(n = length(temp))

# subset data, and get number of temp measurements
datSub <- dat2 %>% filter(beyondTopt == "no") 

#
subsetN <- dat2 %>% filter(beyondTopt == "no") %>% group_by(study.no) %>%
  summarise(n = length(temp))
subsetN

sampleSizeDF <- inner_join(totalN, subsetN, by = "study.no") %>%
  rename(totalN = n.x, subsetN = n.y)
sampleSizeDF


fullDat <- inner_join(dat2, sampleSizeDF)
summary(fullDat)

# remove data points beyond thermal optimum
subDat <- fullDat %>% filter(beyondTopt == "no") %>% 
  filter(subsetN > 2)
summary(subDat)

###################################
# Plots

# Full dataset
qplot(temp, r, data = fullDat, geom = c("point", "smooth"), color = class, 
      ylab = "r - intrinsic rate of increase", 
      xlab = "temperature") + 
  facet_wrap(~ study.no, scales = "free")

# Remove points beyon thermal optimum, and remove studies with <3 points
qplot(temp, r, data = subDat, geom = c("point", "smooth"), color = class, 
      ylab = "r - intrinsic rate of increase", 
      xlab = "temperature") + 
  facet_wrap(~ study.no, scales = "free")

#######
# r vs time
qplot(time, r, data = subDat, geom = c("point", "smooth"), method = "lm", 
      se = FALSE, ylab = "r - intrinsic rate of increase", 
      xlab = "development time", color = class) + 
  facet_wrap(~ study.no, scales = "free")

# r vs fecundity
qplot(fecundity, r, data = subDat, geom = c("point", "smooth"), method = "lm", 
      se = FALSE, ylab = "r - intrinsic rate of increase", 
      xlab = "fecundity", color = class) + 
  facet_wrap(~ study.no, scales = "free")

############################################################
############################################################
############################################################
# Calculate correlations for full dataset
ggDat <- fullDat
summary(ggDat)

AllStudies <- unique(ggDat$study.no)
AllStudies
N <- length(AllStudies)
mat1 <- matrix(nrow = N, ncol = 6)
colnames(mat1) <- c("study.no", "temp.reps", "CCtime", "CCtimeP", 
                    "CCfec", "CCfecP")
mat1

study.i <- AllStudies[5]
study.i
ggDat.i <- ggDat[ggDat$study.no == study.i,]
ggDat.i
rep.i <- length(ggDat.i$time); rep.i
CCtime <- cor.test(ggDat.i$time, ggDat.i$r)$estimate; CCtime

for(i in 1:N) {
	study.i <- AllStudies[i]
	ggDat.i <- ggDat[ggDat$study.no == study.i,]
	rep.i <- length(ggDat.i$time)

	CCtime <- cor.test(ggDat.i$time, ggDat.i$r)$estimate
	CCtimeP <- cor.test(ggDat.i$time, ggDat.i$r)$p.value

	CCfec <- cor.test(ggDat.i$fecundity, ggDat.i$r)$estimate
	CCfecP <- cor.test(ggDat.i$fecundity, ggDat.i$r)$p.value

	mat1[i,] <- c(i, rep.i, CCtime, CCtimeP, CCfec, CCfecP)
}

mat1 <- as.data.frame(mat1)
summary(mat1)

mat1 <- inner_join(mat1, bib2)
mat1

###################################
# plot results
source("./R/multiplotF.R")

a <- qplot(dry.mass.mg, abs(CCtime), data = mat1, 
      geom = "blank",  ylim = c(0, 1), log = "x", 
      xlab = "Dry mass (log mg)", 
      ylab = "Abs[Pearson correlation (r ~ development time)]") + 
  geom_smooth(method = "lm", size = 1, color = "black") +
  geom_point(aes(color = class), size = 8, alpha = 0.5) +
  geom_text(aes(label = study.no, color = NULL), size = 3, face = "bold")


b <- qplot(dry.mass.mg, abs(CCfec), data = mat1, 
      geom = "blank",  ylim = c(0, 1), log = "x", 
      xlab = "Dry mass (log mg)", 
      ylab = "Abs[Pearson correlation (r ~ fecundity)]") + 
  geom_smooth(method = "lm", size = 1, color = "black") +
  geom_point(aes(color = class), size = 8, alpha = 0.5) +
  geom_text(aes(label = study.no, color = NULL), size = 3, face = "bold")

fullDataPlot <- multiplot(a,b)

############################################################
############################################################
############################################################
# Calculate correlations for reduced dataset
ggDat <- subDat
summary(ggDat)

AllStudies <- unique(ggDat$study.no)
AllStudies
N <- length(AllStudies)
mat1 <- matrix(nrow = N, ncol = 6)
colnames(mat1) <- c("study.no", "temp.reps", "CCtime", "CCtimeP", 
                    "CCfec", "CCfecP")
mat1

study.i <- AllStudies[5]
study.i
ggDat.i <- ggDat[ggDat$study.no == study.i,]
ggDat.i
rep.i <- length(ggDat.i$time); rep.i
CCtime <- cor.test(ggDat.i$time, ggDat.i$r)$estimate; CCtime

for(i in 1:N) {
  study.i <- AllStudies[i]
  ggDat.i <- ggDat[ggDat$study.no == study.i,]
  rep.i <- length(ggDat.i$time)
  
  CCtime <- cor.test(ggDat.i$time, ggDat.i$r)$estimate
  CCtimeP <- cor.test(ggDat.i$time, ggDat.i$r)$p.value
  
  CCfec <- cor.test(ggDat.i$fecundity, ggDat.i$r)$estimate
  CCfecP <- cor.test(ggDat.i$fecundity, ggDat.i$r)$p.value
  
  mat1[i,] <- c(i, rep.i, CCtime, CCtimeP, CCfec, CCfecP)
}

mat1 <- as.data.frame(mat1)
summary(mat1)

mat1 <- inner_join(mat1, bib2)
mat1

###################################
# plot results
source("./R/multiplotF.R")

aSub <- qplot(dry.mass.mg, abs(CCtime), data = mat1, 
           geom = "blank",  ylim = c(0, 1), log = "x", 
           xlab = "Dry mass (log mg)", 
           ylab = "Abs[Pearson correlation (r ~ development time)]") + 
  geom_smooth(method = "lm", size = 1, color = "black") +
  geom_point(aes(color = class), size = 8, alpha = 0.5) +
  geom_text(aes(label = study.no, color = NULL), size = 3, face = "bold")


bSub <- qplot(dry.mass.mg, abs(CCfec), data = mat1, 
           geom = "blank",  ylim = c(0, 1), log = "x", 
           xlab = "Dry mass (log mg)", 
           ylab = "Abs[Pearson correlation (r ~ fecundity)]") + 
  geom_smooth(method = "lm", size = 1, color = "black") +
  geom_point(aes(color = class), size = 8, alpha = 0.5) +
  geom_text(aes(label = study.no, color = NULL), size = 3, face = "bold")

subDataPlot <- multiplot(aSub, bSub)

  
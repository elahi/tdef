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
library(ggplot2)
library(grid)

a.base <- ggplot(data = mat1, aes(x = dry.mass.mg, y = CCtime))
a.base + geom_text(aes(label = study.no), size = 2.5, face = "bold") + 
	geom_smooth(stat = "smooth", method = "lm")

####################################################
# test hypotheses using lm (ignoring phylogenetic confounding)
fit1 <- lm(CCtime ~ dry.mass.mg, data = mat1)
summary(fit1)

fit2 <- lm(CCfec ~ dry.mass.mg, data = mat1)
summary(fit2) # not sig

################
# create 2 x 1 plot
# function to save typing
vplayout <- function (x, y)
	viewport(layout.pos.row = x , layout.pos.col = y)

fit1.text <- paste(
	"Adj R2 =", signif(summary(fit1)$adj.r.squared, 3),
	"; Intercept =", signif(fit1$coef[[1]], 3),
	"; Slope =", signif(fit1$coef[[2]], 3),
	"; P =", signif(summary(fit1)$coef[2,4], 3))
fit1.text 

fit2.text <- paste(
	"Adj R2 =", signif(summary(fit2)$adj.r.squared, 3),
	"; Intercept =", signif(fit2$coef[[1]], 3),
	"; Slope =", signif(fit2$coef[[2]], 3),
	"; P =", signif(summary(fit2)$coef[2,4], 3))
fit2.text 

a <- qplot(dry.mass.mg, abs(CCtime), data = mat1, 
	geom = c("smooth"), method = "lm", ylim = c(-1, 1), log = "x", 
	ylab = "Pearson correlation (r ~ development time)") +
	geom_text(aes(label = study.no), size = 2.5, face = "bold") + 
	theme(axis.title.y = element_text(size = 10)) + 
	geom_hline(yintercept = 0, color = "red", lty = 2) + 
	ggtitle(fit1.text) + theme(plot.title = element_text(size = 8, face = "bold"))
	
b <- qplot(dry.mass.mg, abs(CCfec), data = mat1, 
	geom = c("smooth"), method = "lm", ylim = c(-1, 1), log = "x", 
	ylab = "Pearson correlation (r ~ fecundity)") +
	geom_text(aes(label = study.no), size = 2.5, face = "bold") + 
	theme(axis.title.y = element_text(size = 10)) + 
	geom_hline(yintercept = 0, color = "red", lty = 2) +
	ggtitle(fit2.text) + theme(plot.title = element_text(size = 8, face = "bold"))

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(a, vp = vplayout (1, 1))
print(b, vp = vplayout (2, 1))
dev.off()


###################################
# calculate correlations between r and Tea, Tg, Fec
studyN
mat1 <- matrix(nrow = studyN, ncol = 5)
colnames(mat1) <- c("study.no", "temp.reps", "CC.Tea", "CC.Tg", "CC.fec")
mat1
ggDat

### write a for loop to automate this
for(i in 1:studyN) {
	dummy <- subset(dat, study.no == i)

rep <- length(dummy$Tea)
CC.Tea <- cor(dummy$Tea, dummy$r)
CC.Tg <- cor(dummy$Tg, dummy$r)
CC.fec <- cor(dummy$fecundity, dummy$r)

mat1[i,] <- c(i, rep, CC.Tea, CC.Tg, CC.fec )

}

mat1

###################################
# plot results
mat1
mat1 <- as.data.frame(mat1)
cat <- cat[1:16,]; dim(cat); cat <- cat[,1:7]
cat
mat1 <- merge(mat1, cat)
mat1

qplot(dry.mass.mg, CC.Tea,  data = mat1, geom = c("point", "smooth"), method = "lm", ylim = c(-1, 1), log = "x", main = "Fitness (r) ~ Development time (Tea)")
qplot(dry.mass.mg, CC.Tg,  data = mat1, geom = c("point", "smooth"), method = "lm", ylim = c(-1, 1), log = "x", main = "Fitness (r) ~ Generation time (Tg)")
qplot(dry.mass.mg, CC.fec,  data = mat1, geom = c("point", "smooth"), method = "lm", ylim = c(-1, 1), log = "x", main = "Fitness (r) ~ Lifetime fecundity")

# remove studies with 3 or less temps
mat1 <- mat1[mat1$temp.reps >= 4,]
mat1
qplot(dry.mass.mg, CC.Tea,  data = mat1, geom = c("point", "smooth"), method = "lm", ylim = c(-1, 1), log = "x", main = "Fitness (r) ~ Development time (Tea)")
qplot(dry.mass.mg, CC.Tg,  data = mat1, geom = c("point", "smooth"), method = "lm", ylim = c(-1, 1), log = "x", main = "Fitness (r) ~ Generation time (Tg)")
qplot(dry.mass.mg, CC.fec,  data = mat1, geom = c("point", "smooth"), method = "lm", ylim = c(-1, 1), log = "x", main = "Fitness (r) ~ Lifetime fecundity")





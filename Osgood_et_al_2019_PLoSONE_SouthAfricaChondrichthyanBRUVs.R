
#######################################################################################################
#This code is for running GLMs and analyzing MRT and boral models for Osgood et al. 2019 PLoS ONE BRUV
#######################################################################################################

library(boral) # for boral ordination
library(glmmADMB) # for running the GLMs
library(ade4) # for dummy coding for boral analysis
library(mvpart) # for multivariate regression tree
library(labdsv) # for DLI values
library(lmtest) # for likelihood ratio tests
library(spdep) #used to calculate spatial variables to account for spatial autocorrelation in the analyses
library(adespatial) #used to calculate spatial variables to account for spatial autocorrelation in the analyses

#####Read in and prepare data
#environmentals:
env<-read.csv("Osgood_et_al_2019_BRUV_sampling_metadata_environmentals.csv", header=TRUE, stringsAsFactors=FALSE) #environmental data

env$Vis_Code<-as.character(env$Vis_Code) #ensure visibility code is categorical character
env$Vis_Code[env$Vis_Code%in%c("2", "3", "4")]<-"2" #reduce to just low (1) and high (2) categories
env$Site<-factor(env$Site) #ensure site is a factor

#Set up a Julian date starting at 1 for first sampling day
env$StudyJulianDate<-julian(as.Date(env$Date, format = "%y-%m-%d"))-min(julian(as.Date(env$Date, format = "%y-%m-%d")))+1 

#Set up sine and cosine variables for modeling seasonality
env$SIN_TIME=sin(2*pi*env$StudyJulianDate/365.25)
env$COS_TIME=cos(2*pi*env$StudyJulianDate/365.25)

env$Month<-as.numeric(format(as.Date(env$Date), "%m"))
env$Year<-as.numeric(format(as.Date(env$Date), "%Y"))
env$Day<-as.numeric(format(as.Date(env$Date), "%d"))

#add categorical seasons:
env$Season<-env$Month

env$Season[env$Season%in%c("9", "10", "11", "12", "1", "2")]<-"Spring-Summer"
env$Season[env$Season%in%c("3", "4", "5","6", "7", "8")]<-"Fall-Winter"


#Ensure the factors have levels in the order I like (ie sand as baseline, etc.)
env$Habitat<-factor(env$Habitat, levels=c("Sand", "Rocky reef", "Kelp"))
env$Region<-factor(env$Region, levels=c("Walker Bay", "Betty's Bay"))
env$Protection<-factor(env$Protection, levels=c("Protected","Unprotected"))
env$Year<-factor(env$Year, levels=c("16","17", "18"))

#BRUV MaxN data:
BRUV.MaxN<-read.csv("Osgood_et_al_2019_chondrichthyans_MaxN_BRUV.csv", header=TRUE, stringsAsFactors=FALSE)

#Add Julian date to BRUV data for ordering below:
BRUV.MaxN$StudyJulianDate<-julian(as.Date(BRUV.MaxN$Date, format = "%d/%m/%Y"))-min(julian(as.Date(BRUV.MaxN$Date, format = "%d/%m/%Y")))+1 

#Ensure the two data frames are in the same order:
BRUV.MaxN<-BRUV.MaxN[order(BRUV.MaxN$Site, BRUV.MaxN$StudyJulianDate),]
env<-env[order(env$Site, env$StudyJulianDate),]

#Organize data frames by region
BRUV.MaxN<-BRUV.MaxN[order(env$Region), ]
env<-env[order(env$Region), ]

#Checking the two data frames line up (both should be zero)
sum(as.numeric(BRUV.MaxN$StudyJulianDate!=env$StudyJulianDate))
sum(as.numeric(BRUV.MaxN$Site!=env$Site))


########################################################
##Setting up spatial variables##########################
########################################################

LatLong<-data.frame(Longitude=BRUV.MaxN$Longitude, Latitude=BRUV.MaxN$Latitude) #make a data frame of lat long of each site
num_WB<-nrow(subset(env, Region=="Walker Bay")) #number of Walker Bay sites
num_BB<-nrow(env)-num_WB #number of Betty's Bay sites

#Define neighbours of each site using spdep package:
nb<-dnearneigh(as.matrix(LatLong), 0,5, longlat=TRUE) #Neighbours defined to be within 5 km

#Make spatial weighting matrix based on neighbours 
nb.dist<-nbdists(nb, as.matrix(LatLong)) # weights based on distance between neighbours
fdist <- lapply(nb.dist, function(x) 1-x/max(dist(as.matrix(LatLong)))) #Then, spatial weights are defined as a function of distance (1-x/max(dist(xyir))
#And the spatial weighting matrix is then created
listwBRUV <- nb2listw(nb, glist = fdist, style = "B")
groups<-c(num_WB, num_BB) # since I have two areas, I want a seperate set of variables for each matrix 
#I use groups so I don't have large values representing distances between sites in the different areas
#I just get spatial weighting for sites in the same area (that could be autocorrelated)

#spatial eigenfunctions using distance-based Moran's eigenvector maps using adespatial:
dbmem.out<-create.dbMEM.model(coord=as.matrix(LatLong), nsites=groups)
test <- moran.randtest(dbmem.out, listwBRUV, nrepet = 99) #run permutation test to find significant eigenfunctions
dbmem.out<-dbmem.out[,which(test$pvalue<0.05)] #just use significant eigenfunctions


#####################################
###Last bits of data preperation:
#####################################

#Ensure BRUV MaxN data contains only the abundance data (and not any metadata):
#by removing first four columns and the Julian Date Column I added:
BRUV.MaxN<-BRUV.MaxN[,-c(1:4,23)]
length(names(BRUV.MaxN)) #Should be 18 species

#total abundance at each bruv deployment:
sharks.abun<-rowSums(BRUV.MaxN)

#For richness analysis:
#Calculate species richness and presence-absence table
BRUV.Richness<-apply(BRUV.MaxN, 2, function(x) {as.numeric(x>0)})
sharks.rich<-rowSums(BRUV.Richness)


#######################################################
#####Setting up group specific data frames#############
#######################################################

catsharks.abun<-rowSums(BRUV.MaxN[,c("Dark.shyshark", "Puffadder.shyshark",
  "Leopard.catshark", "Pyjama.shark", "Tiger.catshark")])

bigsharks.abun<-rowSums(BRUV.MaxN[,c("Smooth.hound.shark", "Broadnose.sevengill.shark",
  "Soupfin.shark", "Spotted.gully.shark", "Bronze.whaler", "Smooth.hammerhead.shark")])

catsharks.Richness<-BRUV.Richness[,which(colnames(BRUV.Richness)%in%c("Dark.shyshark",
                                                                   "Puffadder.shyshark", "Leopard.catshark",
                                                                   "Pyjama.shark",
                                                                   "Tiger.catshark"))]
catsharks.rich<-rowSums(catsharks.Richness)

skates.Richness<-BRUV.Richness[,which(colnames(BRUV.Richness)%in%c("Spearnose.skate", "Biscuit.skate", "Eagle.ray", 
                                                                   "Lesser.guitarfish", "Short.tail.stingray"))]
skates.rich<-rowSums(skates.Richness)

########################################################
###Analyses#######################################
########################################################

#######################################################################
########GLMs on shark abundance########################################
#######################################################################
#I use the glmmADMB package
#Total shark abundance

#using poisson:
glm.shark.pois<-glmmadmb(sharks.abun~Region*Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+COS_TIME+SIN_TIME+dbmem.out+offset(log(Total_Time))+(1|Site), data=env, 
  family="poisson", extra.args="-ndi 30000")

summary(glm.shark.negbin) #p-values from Wald's z-tests

#Run other models without interaction, habitat, or region to check significance with likelihood ratio test

glm.shark.pois.NoInt<-update(glm.shark.pois, ~.-Region:Protection)
glm.shark.pois.NoHabitat<-update(glm.shark.pois, ~.-Habitat)
glm.shark.pois.NoRegion<-update(glm.shark.pois, ~.-Region)

#check diagnostic plots:
plot(glm.shark.pois$residuals~glm.shark.pois$fitted)

#Running likelihood ratio tests:
lrtest(glm.shark.pois.NoInt, glm.shark.pois) #Interaction significantly improves fit

lrtest(glm.shark.pois.NoRegion, glm.shark.pois) #Region significantly improves fit

lrtest(glm.shark.pois.NoHabitat, glm.shark.pois) #Habitat significantly improves fit

#######################
##Catshark abundance##
#######################
glm.catshark<-glmmadmb(catsharks.abun~Region+Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site) + offset(log(Total_Time)),
                       data=env, family="poisson")
summary(glm.catshark)

glm.catshark.NoProtection<-update(glm.catshark, ~.-Protection)
glm.catshark.NoHabitat<-update(glm.catshark, ~.-Habitat)
glm.catshark.NoRegion<-update(glm.catshark, ~.-Region)

lrtest(glm.catshark, glm.catshark.NoProtection) #Protection important
lrtest(glm.catshark, glm.catshark.NoHabitat) #Habitat important
lrtest(glm.catshark, glm.catshark.NoRegion) #Region important

##################################
##Large shark abundance###########
##################################

glm.bigsharks<-glmmadmb(bigsharks.abun~Region+Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site) + offset(log(Total_Time)), 
                         data=env, family="poisson")
summary(glm.bigsharks)

glm.bigsharks_NoProtection<-update(glm.bigsharks, ~.-Protection)
glm.bigsharks_NoRegion<-update(glm.bigsharks, ~.-Region)
glm.bigsharks_NoHabitat<-update(glm.bigsharks, ~.-Habitat)

lrtest(glm.bigsharks_NoProtection, glm.bigsharks)
lrtest(glm.bigsharks_NoRegion, glm.bigsharks)
lrtest(glm.bigsharks_NoHabitat, glm.bigsharks)

######################################################################
##GLMS on shark richness##############################################
######################################################################

glm.shark_rich<-glmmadmb(sharks.rich~Region*Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site), data=env, 
                         family="poisson")
summary(glm.shark_rich) #Wald z-test

plot(glm.shark_rich$residuals~glm.shark_rich$fitted) #diagnostic plots

glm.shark_rich.NoInt<-update(glm.shark_rich, ~.-Region:Protection)

glm.shark_rich.NoProtection<-update(glm.shark_rich, ~.-Protection)
glm.shark_rich.NoHabitat<-update(glm.shark_rich, ~.-Habitat)
glm.shark_rich.NoRegion<-update(glm.shark_rich, ~.-Region)

lrtest(glm.shark_rich, glm.shark_rich.NoInt) #Interaction significant 
lrtest(glm.shark_rich, glm.shark_rich.NoHabitat) #Habitat significant

####################################################
##############GLMs on FO############################
####################################################

glm.shark_FO<-glmmadmb(as.numeric(sharks.rich>0)~Region+Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site), data=env, extra.args="-ndi 30000",
                         family="binomial")
summary(glm.shark_FO) 

glm.shark_FO_NoRegion<-update(glm.shark_FO, ~.-Region)
glm.shark_FO_NoProtection<-update(glm.shark_FO, ~.-Protection)
glm.shark_FO_NoHabitat<-update(glm.shark_FO, ~.-Habitat)

#likelihood ratio tests
lrtest(glm.shark_FO_NoHabitat, glm.shark_FO)
lrtest(glm.shark_FO_NoRegion, glm.shark_FO)
lrtest(glm.shark_FO_NoProtection, glm.shark_FO)

#######################
##Catshark FO#########
#######################
glm.catshark_FO<-glmmadmb(as.numeric(catsharks.rich>0)~Region+Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site) + offset(log(Total_Time)), data=env, extra.args="-ndi 30000",
                         family="binomial")
summary(glm.catshark_FO) 

glm.shark_FO_NoHabitat<-update(glm.catshark_FO, ~.-Habitat)
glm.shark_FO_NoRegion<-update(glm.catshark_FO, ~.-Region)

lrtest(glm.shark_FO,glm.shark_FO_NoHabitat)
lrtest(glm.shark_FO,glm.shark_FO_NoRegion)

#######################
##Large sharks FO######
#######################

glm.bigsharksFO<-glmmadmb(as.numeric(bigsharks.abun>0)~Region+Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site) + offset(log(Total_Time)), 
                         data=env, family="binomial")
summary(glm.bigsharksFO)

glm.bigsharksFO_NoProtection<-update(glm.bigsharksFO, ~.-Protection)
glm.bigsharksFO_NoRegion<-update(glm.bigsharksFO, ~.-Region)
glm.bigsharksFO_NoHabitat<-update(glm.bigsharksFO, ~.-Habitat)

lrtest(glm.bigsharksFO_NoProtection, glm.bigsharksFO)
lrtest(glm.bigsharksFO_NoRegion, glm.bigsharksFO)
lrtest(glm.bigsharksFO_NoHabitat, glm.bigsharksFO)


#######################
##Skates FO############
#######################
glm.skates_FO<-glmmadmb(as.numeric(skates.rich>0)~Region+Protection+Habitat+Depth+Vis_Code+WaterTemp+
  Year+SIN_TIME+COS_TIME+dbmem.out+offset(log(Total_Time))+(1|Site) + offset(log(Total_Time)), data=env, family="binomial")
summary(glm.skates_FO)

glm.skates_FO_NoHabitat<-update(glm.skates_FO, ~.-Habitat)
lrtest(glm.skates_FO_NoHabitat, glm.skates_FO)


#################################################################
########Multivariate regression tree#############################
#################################################################

#calculate Bray-Curtis dissimilarity

g.out<-gdist(BRUV.MaxN, meth="bray", full=TRUE, sq=TRUE)

#Set NAs (ie comparisons between rows where row sum = 0) to 0 dissimilarity (ie completely similar)

g.out.sharks<-apply(g.out, 2, function(x) {x[is.na(x)==TRUE]<-0; return(x)})

#calculate the multivariate regression tree using Bray-Curtis similarity calculated above
mv.outSharkMaxN<-mvpart(g.out.sharks~Depth+Habitat+WaterTemp+Region+Protection+Season, 
                        data=env, xv="pick", xval=nrow(BRUV.MaxN), xvmult=100, plot=T)

#calculate indicator values
ind.out<-indval(BRUV.MaxN, mv.outSharkMaxN$where, numitr=1000) # total DLI of leaves
ind.out$pval
ind.out$maxcls[which(ind.out$pval<=0.05)] #which leaf they indicate
ind.out$indcls[which(ind.out$pval<=0.05)] # significant indicator values

###Determine individual values for first split - kelp and reef vs sand

first<-mv.outSharkMaxN$where
first[first%in%c(2)]<-1
first[first%in%c(4,5)]<-2

ind.out.first<-indval(BRUV.MaxN, first, numitr=1000) # total DLI of leaves
ind.out.first$pval
ind.out.first$maxcls[which(ind.out.first$pval<=0.05)] #which leaf they indicate
ind.out.first$indcls[which(ind.out.first$pval<=0.05)] # significant indicator values

###Determine individual values for second split - Betty's Bay vs Walker Bay in sand

second.sharks.MaxN<-subset(BRUV.MaxN, first==2)
second.sharks.MaxN<-second.sharks.MaxN[,colSums(second.sharks.MaxN)>0]
second<-subset(mv.outSharkMaxN$where, first==2)

ind.out.second<-indval(second.sharks.MaxN, second, numitr=1000) # total DLI of leaves
ind.out.second$pval
ind.out.second$maxcls[which(ind.out.second$pval<=0.05)]#which leaf they indicate
ind.out.second$indcls[which(ind.out.second$pval<=0.05)]# significant indicator values


###################################################
###Latent variable only models for paper###########
###################################################

#Boral on shark abundance - latent variable model (with row random effects)
boral.out_sharks_latent_roweff<-boral(sharks.MaxN, row.eff = "random",
  family="poisson", num.lv = 2,
   save.model=T, offset=env_sharks_TotalTime)
save(boral.out_sharks_latent_roweff, file="boral_sharks_latent_roweff.Rdata")

#Residual analysis
par(mfrow=c(2,2))
plot(boral.out_sharks_exp) 

#######################################################################################################
#This code is for the figures for Osgood et al. 2019 PLoS ONE BRUV
#######################################################################################################
library(reshape2) #for melt function

#####################
###Figure 2##########
#####################

#Plot total number observed by habitat

sharks.habitat<-aggregate(BRUV.MaxN, by=list(env$Habitat), sum)

SpeciesTotals<-data.frame(melt(sharks.habitat))  #Calculate Total number of individuals observed per species per habitat

names(SpeciesTotals)<-c("Habitat", "Species", "Total_MaxN") #rename column names

SpeciesTotals$Taxon<-as.character(SpeciesTotals$Species) #Create taxon column

#Put correct family into taxon column
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Soupfin.shark", "Spotted.gully.shark", "Smooth.hound.shark")]<-"Triakidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("St..Joseph.shark")]<-"Callorhinchidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Broadnose.sevengill.shark")]<-"Hexanchidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Short.tail.stingray")]<-"Dasyatidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Eagle.ray")]<-"Myliobatidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Spearnose.skate", "Biscuit.skate")]<-"Rajidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Shortnose.spurdog")]<-"Squalidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Lesser.guitarfish")]<-"Rhinobatidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Bronze.whaler")]<-"Carcharhinidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Dark.shyshark", "Leopard.catshark",
                           "Pyjama.shark", "Puffadder.shyshark", "Tiger.catshark")]<-"Scyliorhinidae"
SpeciesTotals$Taxon[SpeciesTotals$Taxon%in%c("Smooth.hammerhead.shark")]<-"Sphyrnidae"

# total MaxN by family by habitat
TaxonTotals<-aggregate(SpeciesTotals$Total_MaxN, by=list(SpeciesTotals$Habitat, SpeciesTotals$Taxon) , sum) 
names(TaxonTotals)<-c("Habitat", "Taxon", "Total_MaxN")

# Order families by total MaxN
TaxonTotals<-TaxonTotals[c(22:24, 31:33, 16:18, 7:12, 1:3, 28:30, 19:21, 13:15, 4:6, 25:27),] 

#Plot total MaxN versus family:

shark_order_names<-c("Scyliorhinidae", "Triakidae", "Rajidae", "Dasyatidae", "Hexanchidae",
  "Callorhinchidae", "Squalidae", "Rhinobatidae", "Myliobatidae", "Carcharhinidae",
  "Sphyrnidae") #Names for x axis labels of taxa in correct order

#Plot total number observed by taxonomic group and protection level

sharks.protection<-aggregate(BRUV.MaxN, by=list(env$Region, env$Protection), sum)

SpeciesTotals.prot<-data.frame(melt(sharks.protection))  #Calculate Total number of individuals observed per species per area of protection

names(SpeciesTotals.prot)<-c("Region", "Protection", "Species", "Total_MaxN")

SpeciesTotals.prot$Taxon<-as.character(SpeciesTotals.prot$Species) #Create taxon column

#Put correct family into taxon column
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Soupfin.shark", "Spotted.gully.shark", "Smooth.hound.shark")]<-"Triakidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("St..Joseph.shark")]<-"Callorhinchidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Broadnose.sevengill.shark")]<-"Hexanchidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Short.tail.stingray")]<-"Dasyatidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Eagle.ray")]<-"Myliobatidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Spearnose.skate", "Biscuit.skate")]<-"Rajidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Shortnose.spurdog")]<-"Squalidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Lesser.guitarfish")]<-"Rhinobatidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Bronze.whaler")]<-"Carcharhinidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Dark.shyshark", "Leopard.catshark",
                           "Pyjama.shark", "Puffadder.shyshark", "Tiger.catshark")]<-"Scyliorhinidae"
SpeciesTotals.prot$Taxon[SpeciesTotals.prot$Taxon%in%c("Smooth.hammerhead.shark")]<-"Sphyrnidae"

SpeciesTotals.prot$Protection<-factor(SpeciesTotals.prot$Protection, levels=c("Unprotected",
  "Protected"))

TaxonTotals.prot<-aggregate(SpeciesTotals.prot$Total_MaxN, by=list(SpeciesTotals.prot$Protection, SpeciesTotals.prot$Region, SpeciesTotals.prot$Taxon) , sum) # total MaxN by family by habitat
names(TaxonTotals.prot)<-c("Protection", "Region", "Taxon", "Total_MaxN")

TaxonTotals.prot<-TaxonTotals.prot[c(29:32, 41:44, 21:24, 9:16, 1:4, 37:40, 25:28, 17:20, 5:8, 33:36),] # Order families by total MaxN

####Code for figure below#####

#layout for the figure:
nf<-layout(matrix(c(rep(c(1,2,2,2,2,2,2),3),rep(c(3,4,4,4,4,4,4),3)), 6, 7, byrow=TRUE),respect=T)
par(mar=c(3,3,3,0.1), oma=c(4,2,2,2), xpd = NA)

#Making bar plot with catsharks seperated out due to higher abundance for them (a)

b<-barplot(TaxonTotals.prot$Total_MaxN[1:4],
        ylab="", cex.names = 0.8, xlab="", ylim=c(0, 300), xaxt='n', 
        col=c("#37AFA9", "#0C695D", "#FAAC77", "#F6893D"), 
        las=2)
text(2.5, par("usr")[3]-5,
     srt = 30, adj= 1, xpd = TRUE,
     labels = shark_order_names[1] , cex=1.3, xpd=NA)
mtext("Total Max N", side=2, line=2.5, adj=0.5)
mtext("(a)", side=3, line=0, adj=0, cex=1.3)

#Rest of species (b)
b<-barplot(TaxonTotals.prot$Total_MaxN[-c(1:4)], space=c(rep(c(1,0,0,0),10)),
        ylab="", cex.names = 1.3, xlab="", ylim=c(0, 30), 
        xaxt='n', col=c("#37AFA9", "#0C695D", "#FAAC77", "#F6893D"),
        las=2)
text(seq(3,48,by=5), par("usr")[3]-0.7,
     srt = 30, adj= 1, xpd = TRUE, 
     labels = shark_order_names[-1], cex=1.3, xpd=NA)
legend(19,30, c("Unprotected", "Protected"),
       fill=c("#37AFA9", "#0C695D"), title="Walker Bay", bty="n", cex=1.9)
legend(35,30, c("Unprotected", "Protected"), 
       fill=c("#FAAC77", "#F6893D"), title="Betty's Bay", bty="n", cex=1.9)
mtext("(b)", side=3, line=0, adj=0, cex=1.3)

#Making bar plot with catsharks seperated out due to their higher abundance (c)

b<-barplot(TaxonTotals$Total_MaxN[1:3],
        ylab="", cex.names = 0.8, xlab="", ylim=c(0, 600), xaxt='n', 
        col=c("#000000", "#e79f00", "#009E73"), 
        las=2)
text(1.9, par("usr")[3]-10,
     srt = 30, adj= 1, xpd = TRUE,
     labels = shark_order_names[1] , cex=1.3, xpd=NA)
mtext("Total Max N", side=2, line=2.5, adj=0.5)
mtext("(c)", side=3, line=0, adj=0, cex=1.3)

#Rest of species (d)
b<-barplot(TaxonTotals$Total_MaxN[-c(1:3)], space=c(rep(c(1,0,0),10)),
        ylab="", cex.names = 1.3, xlab="", ylim=c(0, 50), xaxt='n', 
        col=c("#000000", "#e79f00", "#009E73"),
        las=2)
text(seq(2.5,38.5,by=4), par("usr")[3]-0.7,
     srt = 30, adj= 1, xpd = TRUE, 
     labels = shark_order_names[-1], cex=1.3, xpd=NA)
legend(33,50, c("Sand", "Reef", "Kelp"),
       fill=c("#000000", "#e79f00", "#009E73"), title="Habitat", bty="n", cex=1.9)
mtext("(d)", side=3, line=0, adj=0, cex=1.3)
mtext("Family", side=1, line=5, adj=0.37, cex=1.3)

#######################################################
#############Figure 3##################################
#######################################################

#Setting up barplots

se<-function(x) {sd(x)/sqrt(length(x))} # Function for calculating standard error


env$Protection<-factor(env$Protection, levels=c("Unprotected","Protected"))

sharks.abun.mean<-aggregate(sharks.abun, by=list(env$Protection, env$Region), mean)
sharks.abun.se<-aggregate(sharks.abun, by=list(env$Protection, env$Region), se)

sharks.habitat.mean<-aggregate(sharks.abun, by=list(env$Habitat), mean)
sharks.habitat.se<-aggregate(sharks.abun, by=list(env$Habitat), se)

sharks.rich.mean<-aggregate(sharks.rich, by=list(env$Protection, env$Region), mean)
sharks.rich.se<-aggregate(sharks.rich, by=list(env$Protection, env$Region), se)

sharks.rich.habitat.mean<-aggregate(sharks.rich, by=list(env$Habitat), mean)
sharks.rich.habitat.se<-aggregate(sharks.rich, by=list(env$Habitat), se)

op <- par(mfrow=c(2,2), oma = c(5,3,1.5,0.2) + 0.1,
          mar = c(2.4,1,0.2,0.1) + 0.1)
sharks.abun.bar<-barplot(sharks.abun.mean$x, ylim=c(0,5),
                         ylab="", cex.axis=2, las=2, space=c(0.2,0,0.2,0), col="grey")
arrows(sharks.abun.bar, sharks.abun.mean$x-sharks.abun.se$x,
       sharks.abun.bar, sharks.abun.mean$x+sharks.abun.se$x,angle=90,code=3)
mtext("(a)", side=3, line=0, adj=0, cex=2)
mtext("Mean MaxN", side=2, line=2.2, adj=0.5, cex=1.5)
text(0.7,3, "a", cex=2)
text(1.7,3, "b", cex=2)
text(2.9,4.3, "c", cex=2)
text(3.9,4.3, "d", cex=2)


sharks.habitat.bar<-barplot(sharks.habitat.mean$x, ylim=c(0,5), yaxt="n", col="grey")
arrows(sharks.habitat.bar, sharks.habitat.mean$x-sharks.habitat.se$x,
       sharks.habitat.bar, sharks.habitat.mean$x+sharks.habitat.se$x,angle=90,code=3)
mtext("(b)", side=3, line=0, adj=0, cex=2)
text(1.9,4.3, "b", cex=2)
text(3.1,4.3, "b", cex=2)
text(0.7,2.1, "a", cex=2)


sharks.rich.bar<-barplot(sharks.rich.mean$x, ylim=c(0,3),
                         ylab="", yaxt="n", space=c(0.2,0,0.2,0), col="grey")
axis(2, labels=c(0,1,2,3), at=c(0,1,2,3), cex.axis=1.5, las=2)
text(c(0.7,1.7,2.9,3.9), par("usr")[3]-0.2,
     srt = 0, adj= 0.5, xpd = TRUE,
     labels = rep(c("Unprotected", "Protected"),2), cex=1.2)
arrows(sharks.rich.bar, sharks.rich.mean$x-sharks.rich.se$x,
       sharks.rich.bar, sharks.rich.mean$x+sharks.rich.se$x,angle=90,code=3)
mtext("(c)", side=3, line=0, adj=0, cex=2)
mtext("Mean species richness", side=2, line=2.2, adj=0.5, cex=1.5)
mtext(expression(paste("    Walker Bay\nWhale Sanctuary")), side=1, line=5, adj=0.2, cex=1.2)
mtext("Betty's Bay MPA", side=1, line=3, adj=0.85, cex=1.2)
text(0.7,1.8, "a", cex=2)
text(1.7,1.8, "b", cex=2)
text(2.9,2.9, "c", cex=2)
text(3.9,2.9, "d", cex=2)

sharks.rich.habitat.bar<-barplot(sharks.rich.habitat.mean$x, ylim=c(0,3), yaxt="n", col="grey")
text(c(0.7,1.9,3.1), par("usr")[3]-0.2, srt = 0, adj= 0.5, xpd = TRUE,
     labels = c("Sand", "Reef", "Kelp"), cex=1.5)
arrows(sharks.rich.habitat.bar, sharks.rich.habitat.mean$x-sharks.rich.habitat.se$x,
       sharks.rich.habitat.bar, sharks.rich.habitat.mean$x+sharks.rich.habitat.se$x,angle=90,code=3)
mtext("(d)", side=3, line=0, adj=0, cex=2)
text(1.9,2.7, "b", cex=2)
text(3.1,2.7, "a", cex=2)
text(0.7,1.5, "a", cex=2)

#######################################################
#############Figure 5b##################################
#######################################################

color_ellipse_shark<-c("#000000", "#e79f00", "#009E73") #set up color palette for colouring points based on habitat

#Get biplot scores (x and y values of points and species labels)
biplot.out_sharks<-lvsplot(boral.out_sharks_latent_roweff, biplot=TRUE, est = "mean", alpha=0.6, 
  return.val=TRUE)

##Figure 5b - boral ordination
par(mar=c(5,6,5,5))
plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
plot(biplot.out_sharks$scaled.lvs, col=color_ellipse_shark[env$Habitat], pch=16,  # have seperate colour for each habitat
     ylab="Latent variable 2", xlab="Latent variable 1", 
     cex.axis=1.5, cex.lab=2, xlim=c(-2,2), ylim=c(-2.5,2))
#this adds confidence ellipses for all points from the same habitat:
for (i in unique(as.numeric(env$Habitat))) ordiellipse (biplot.out_sharks$scaled.lvs[,], 
                                                        groups = as.numeric(env$Habitat), show.group = i, 
                                                        col = color_ellipse_shark[i], 
                                                        kind="sd", label = F, lwd=4)
#This adds the text labels for each species based on locations from the biplot saved above:
text(biplot.out_sharks$scaled.lv.coefs+0.01, labels=c("BS", "BG", "BW", "DS", "ER", "LC", "LG", "PS",
            "PJ","SR", "SD", "SH", "HH",  "SF", "SN", "SG", "SJ", "TC"), cex=1.6,
col=c("#0072B2"))
#the legend for the colours based on habitat:
legend(-2, 2, c("Sand", "Reef", "Kelp"), col=c("#000000", "#e79f00", "#009E73"), pch=16, cex=1.5)
#mtext("(b)", side=3, line=0.5, adj=0, cex=2) #put "b" on this panel

########################################################################
########Species accumulation curve######################################
########################################################################

sp.curve<-specaccum(BRUV.MaxN, method="random", permutations=100)
plot(sp.curve, ylab="Cumulative number of species", xlab="Number of BRUV samples")
#boxplot(sp.curve)
slopes <- with(sp.curve,diff(richness)/diff(sites)) #Reaches asymptote.
length(slopes)
which(slopes<0.05 & slopes>0.035)

sp.curveBettysBay<-specaccum(BRUV.MaxN[which(env$Region=="Betty's Bay"),], method="random", permutations=100)

with(sp.curveBettysBay,diff(richness)/diff(sites)) #Reaches asymptote.

sp.curveWalkerBay<-specaccum(BRUV.MaxN[which(env$Region=="Walker Bay"),], 
  method="random", permutations=100)

with(sp.curveWalkerBay,diff(richness)/diff(sites)) #Barely Reaches asymptote.

#Figure S1

par(mfrow=c(2,2), mar=c(5,5,3,3))
plot(sp.curve, ylab="Cumulative species richness", xlab="Number of samples", ci.col="grey")
mtext("(a)", side=3, line=0.3, adj=0)

plot(sp.curveBettysBay, ylab="", xlab="Number of samples", ci.col="grey")
mtext("(b)", side=3, line=0.3, adj=0)

plot(sp.curveWalkerBay, ylab="Cumulative species richness", xlab="Number of samples", ci.col="grey")
mtext("(c)", side=3, line=0.3, adj=0)

length(unique(env$Site))

##Table S1

data.frame(table(env$Region, env$Protection, env$Season, env$Year))$Freq

range(env$Date)
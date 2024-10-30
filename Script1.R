# SCRIPT 1: IMPORT AND FORMAT DATA
library(dplyr)
library(readr)
suppressWarnings(suppressMessages(library(tidyr)))

# Summarises data into tables:
# siteStrata lists stratum heights and cover for each site
# siteData lists environmental attributes and delegatensis composition in each site


# Fire regime data
siteData <- read.csv("siteData.csv") %>%
  mutate(solCool = sApr+sMay+sJun+sJul+sAug+sSep,
         solWarm = sOct+sNov+sDec+sJan+sFeb+sMar,
         sev91 = case_when(dNBR91 > 225 & lastFire03 == 1991 ~ 3,
                           dNBR91 > 25 & lastFire03 == 1991 ~ 2,
                           TRUE ~ 1),
         sev03 = case_when(dNBR03 > 225 & lastFire06 == 2003 ~ 3,
                           dNBR03 > 25 & lastFire06 == 2003 ~ 2,
                           TRUE ~ 1),
         sev06 = case_when(dNBR06 > 225 & lastFire07 == 2006 ~ 3,
                           dNBR06 > 25 & lastFire07 == 2006 ~ 2,
                           TRUE ~ 1),
         sev07 = case_when(dNBR07 > 225 & lastFire13 == 2007 ~ 3,
                           dNBR07 > 25 & lastFire13 == 2007 ~ 2,
                           TRUE ~ 1),
         short = as.numeric(lastFire91 >= 1986)+as.numeric(lastFire03 >= 1998)+
           as.numeric(lastFire06 >= 2001)+as.numeric(lastFire07 >= 2002),
         shortAsh = as.numeric(lastFire91 >= 1971)+as.numeric(lastFire03 >= 1983)+
           as.numeric(lastFire06 >= 1986)+as.numeric(lastFire07 >= 1987),
         TSF = surveyYear - lastFire13,
         obsReliability = (TSF < 10 & !is.na(sevObserved)),
         highSev = pmax((as.numeric(sev91 == 3)+as.numeric(sev03 == 3)+as.numeric(sev06 == 3)+as.numeric(sev07 == 3)), (obsReliability & sevObserved == 3)),
         TSFclass = ceiling(TSF/5)*5,
         lastSev = case_when(lastFire13 == 2007 ~ sev07,
                             lastFire13 == 2006 ~ sev06,
                             lastFire13 == 2003 ~ sev03,
                             TRUE ~ sev91))

# Create sev class for unknown sites
siteData$highSev[which(siteData$sev91 == 1 & siteData$sev03 == 1 & siteData$sev06 == 1 & siteData$sev07 == 1 & siteData$highSev == 0)] <- NA


# Floristics
smallPlants <- read.csv("under10cm.csv") %>%
  mutate(Numbers = New.number) %>%
  select(Site, Species, Numbers)

# Summarise delegatensis composition
smallDel <- smallPlants[which(smallPlants$Species=="Eucalyptus delegatensis"),] %>%
  mutate(del = Numbers) %>% select(Site, del) 
nSmall <- smallPlants %>%
  group_by(Site) %>%
  summarise_if(is.numeric, sum) %>%
  left_join(smallDel, by = "Site")
nSmall$del[which(is.na(nSmall$del))]<-0
nSmall$delCompY<-100*round(nSmall$del/nSmall$Numbers, 2)
nSmall <- nSmall %>%
  select(Site, delCompY)

largePlants <- read.csv("over10cm.csv") %>%
  select(Site, Species, DBH, Height) %>% # Ash ages From Fig. 2a in Mokany et. al. 2003 Tree Physiol.
  mutate(Age = (Species == "Eucalyptus delegatensis")*round(1.43*Height,0)) %>%
  left_join(siteData, by = "Site") %>%
  mutate(DBHdev = DBH - (2.39*Height))
largePlants$Age[which(largePlants$Age==0)]<-NA

# Summarise delegatensis composition
nLarge <- data.frame(matrix(ncol = 2, nrow = length(unique(largePlants$Site))))
colnames(nLarge) <- c('Site', 'delCompO')
n <- 1
for (site in unique(largePlants$Site)) {
  nLarge$Site[n] <- site
  nLarge$delCompO[n] <- round(length(which(largePlants$Species[largePlants$Site==site]=="Eucalyptus delegatensis"))/
                                length(largePlants$Species[largePlants$Site==site]),2)*100
  n<- n+1
}

delComp <- left_join(nSmall,nLarge, by = "Site")
delComp$delCompO[which(is.na(delComp$delCompO))]<-0

# Sort strata
siteStrataA <- read.csv("siteStrata.csv") %>%
  mutate(UpperHeight = pmax(LowerHeight, UpperHeight, na.rm = TRUE),
         Site = parse_number(SiteNumber),
         medianHeight = (LowerHeight + UpperHeight)/2)
siteStrata <- data.frame()
for (site in unique(siteStrataA$SiteNumber)) {
  s <- siteStrataA[siteStrataA$SiteNumber==site,]
  s <- s[order(s$UpperHeight, s$LowerHeight),]
  s$Stratum[1]<-"NS"
  s$Stratum[nrow(s)]<-"C"
  s$Stratum[s$Stratum!="NS" & s$Stratum!="C"][1]<-"E"
  s$Stratum[s$Stratum!="NS" & s$Stratum!="E" & s$Stratum!="C"][1]<-"M"
  s$Stratum[s$Stratum!="NS" & s$Stratum!="E" & s$Stratum!="M" & s$Stratum!="C"][1]<-"MU"
  siteStrata <- rbind(siteStrata, s)
}

siteStrata <- siteStrata %>%
  select(Site, Stratum, LowerHeight, medianHeight, UpperHeight, PercentCover)
siteData <- siteData %>%
  left_join(delComp, by = "Site")
ashData <- left_join(siteStrata, siteData, by = "Site")

# Divisions
siteMeans <- ashData %>%
  group_by(Site) %>%
  summarise_if(is.numeric, mean) %>%
  select(Elevation, Protection, Aspect, solCool, solWarm, BIO1, BIO5, BIO6, BIO12, BIO19)
altitude <- quantile(siteMeans$Elevation, probs = 0.5)
protection <- quantile(siteMeans$Protection, probs = 0.5)
solC <- quantile(siteMeans$solCool, probs = 0.5)
solW <- quantile(siteMeans$solWarm, probs = 0.5)
tDiv <- quantile(siteMeans$BIO1, probs = 0.5)
tmDiv <- quantile(siteMeans$BIO5, probs = 0.5)
tlDiv <- quantile(siteMeans$BIO6, probs = 0.5)
rDiv <- quantile(siteMeans$BIO12, probs = 0.5)
rcDiv <- quantile(siteMeans$BIO19, probs = 0.5)


for (row in 1:nrow(ashData)) {
  ashData$altClass[row] <- case_when(ashData$Elevation[row] <= altitude ~ round(mean(ashData$Elevation[ashData$Elevation<=altitude]),0),
                                     TRUE ~ round(mean(ashData$Elevation[ashData$Elevation>altitude]),0))
  ashData$protClass[row] <- case_when(ashData$Protection[row] <= protection ~ round(mean(ashData$Protection[ashData$Protection<=protection]),0),
                                      TRUE ~ round(mean(ashData$Protection[ashData$Protection>protection]),0))
  ashData$aspClass[row] <- case_when(ashData$Aspect[row] >= 90 & ashData$Aspect[row] <270 ~ 180,
                                      TRUE ~ 0)
  ashData$solCClass[row] <- case_when(ashData$solCool[row] <= solC ~ round(mean(ashData$solCool[ashData$solCool<=solC]),0),
                                      TRUE ~ round(mean(ashData$solCool[ashData$solCool>solC]),0))
  ashData$solWClass[row] <- case_when(ashData$solWarm[row] <= solW ~ round(mean(ashData$solWarm[ashData$solWarm<=solW]),0),
                                      TRUE ~ round(mean(ashData$solWarm[ashData$solWarm>solW]),0))
  ashData$tClass[row] <- case_when(ashData$BIO1[row] <= tDiv ~ round(mean(ashData$BIO1[ashData$BIO1<=tDiv]),0),
                                      TRUE ~ round(mean(ashData$BIO1[ashData$BIO1>tDiv]),0))
  ashData$tmClass[row] <- case_when(ashData$BIO5[row] <= tmDiv ~ round(mean(ashData$BIO5[ashData$BIO5<=tmDiv]),0),
                                   TRUE ~ round(mean(ashData$BIO5[ashData$BIO5>tmDiv]),0))
  ashData$tlClass[row] <- case_when(ashData$BIO6[row] <= tlDiv ~ round(mean(ashData$BIO6[ashData$BIO6<=tlDiv]),0),
                                   TRUE ~ round(mean(ashData$BIO6[ashData$BIO6>tlDiv]),0))
  ashData$rClass[row] <- case_when(ashData$BIO12[row] <= rDiv ~ round(mean(ashData$BIO12[ashData$BIO12<=rDiv]),0),
                                   TRUE ~ round(mean(ashData$BIO12[ashData$BIO12>rDiv]),0))
  ashData$rcClass[row] <- case_when(ashData$BIO19[row] <= rcDiv ~ round(mean(ashData$BIO19[ashData$BIO19<=rcDiv]),0),
                                   TRUE ~ round(mean(ashData$BIO19[ashData$BIO19>rcDiv]),0))
}

#####################
# Species composition
agesL <- ashData %>%
  select(Site, TSF, highSev, short, altClass, protClass, aspClass, solCClass, solWClass, Lithology, tClass, tmClass, tlClass, rClass, rcClass)

NS <- read.csv("NS.csv") %>%
left_join(agesL, by = "Site") %>%
  select(-1)
NS[is.na(NS)] <- 0
NS <- NS  %>%
  group_by(TSF, highSev, short, altClass, protClass, aspClass, solCClass, solWClass, Lithology, tClass, tmClass, tlClass, rClass, rcClass) %>%
  summarise_if(is.numeric, mean)
NSLong <- gather(NS, key = "Species", value = "Presence", -c(1,2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14))
NSLong$Stratum <- "NS"

E <- read.csv("E.csv") %>%
  left_join(agesL, by = "Site") %>%
  select(-1)
E[is.na(E)] <- 0
E <- E  %>%
  group_by(TSF, highSev, short, altClass, protClass, aspClass, solCClass, solWClass, Lithology, tClass, tmClass, tlClass, rClass, rcClass) %>%
  summarise_if(is.numeric, mean)
ELong <- gather(E, key = "Species", value = "Presence", -c(1,2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14))
ELong$Stratum <- "E"

M <- read.csv("M.csv") %>%
  left_join(agesL, by = "Site") %>%
  select(-1)
M[is.na(M)] <- 0
M <- M  %>%
  group_by(TSF, highSev, short, altClass, protClass, aspClass, solCClass, solWClass, Lithology, tClass, tmClass, tlClass, rClass, rcClass) %>%
  summarise_if(is.numeric, mean)
MLong <- gather(M, key = "Species", value = "Presence", -c(1,2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14))
MLong$Stratum <- "M"

C <- read.csv("C.csv") %>%
  left_join(agesL, by = "Site") %>%
  select(-1)
C[is.na(C)] <- 0
C <- C  %>%
  group_by(TSF, highSev, short, altClass, protClass, aspClass, solCClass, solWClass, Lithology, tClass, tmClass, tlClass, rClass, rcClass) %>%
  summarise_if(is.numeric, mean)
CLong <- gather(C, key = "Species", value = "Presence", -c(1,2,3,4,5,6,7,8, 9, 10, 11, 12, 13, 14))
CLong$Stratum <- "C"


specDynamics <- bind_rows(NSLong, ELong, MLong, CLong, .id = NULL) %>%
  mutate(Euc = grepl("Eucalyptus", Species))

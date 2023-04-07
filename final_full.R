#loading package ----
library(dplyr)        # A package for data manipulation
library(sf)           # Simple feature for R
library(spdep)        # Functions and tests for evaluating spatial patterns
library(tidyr)        # Tools to create tidy data
library(INLA)         # Integrated Nested Laplace Approximation package
library(ggplot2)      # A package for creating maps and graphs
library(viridis)      # A package providing color palettes
library(patchwork)    # A package to compose plots

# For tables in RMarkdown
library(knitr)
library(kableExtra)
#looking at the data
#first we will set the working directory
setwd("~/OneDrive/Documents/HDA/Advanced analytics/AA_project/DATA_SCIENTIFIC_PROJECTS/DS2_DengueBrazil")
#load in the csv
data <- read.csv("data_2015_2019.csv")
#dealing with the missing microregion ------
#The PDSI was missing for the island FERNANDO DE NORONHA with microregion code 26019,
##thus we excluded this microregion from the data set, so that you have now 557 microregions.
###Therefore, you should remove this microregion from the shapefile as well,
#before obtaining the graph for INLA
#(this is easy to do after importing the shapefile into R using the sf and dplyr packages).
#firs ned to install sf package
install.packages("sf")
#using the sf packages to read shapefile
library(sf)
shape <- read_sf(dsn = ".", layer = "shape_brazil")
#removing microregion code 26019
shape <- subset(shape, code!= 26019)
#now we hae removed that code as that microregion is unavaliable
#viewing map
library(ggplot2)
shape %>%
  ggplot() +
  geom_sf()
#creating neighbours and weights between the region
shape_nb = poly2nb(shape, snap=1000, queen=TRUE)
summary(shape_nb)
#adding the neighbour
nb2INLA("shape.graph",shape_nb)
shape.adj = paste(getwd(),"/shape.graph",sep="")



###pre processing ----
#looking for missingness
sapply(data,is.na)
colSums(is.na(data))
#there is no missingess (good)
## we need to include expected cases
data$E  <- data$population/10^5
#use log(E) within the inla function
#discovery analysis ----
#looking at the regions
ggplot() +
  geom_sf(data = shape, color = "blue", fill = "white") +
  coord_sf() +    #axis limits and CRS
  theme_bw() +    # dark-on-light theme
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
#we can confirm that the arrangment of the microegion fits the data presented 
##by inistry of Health Information Department (DATASUS)

#rebanding PDSI to categorical varaible for analysis
data$rebanded_pdsi <- cut(RESP_DATA_ST$pdsi, 
                                  breaks = c(-Inf, -3, -2, 1.9, 2.9, Inf),
                                  labels = c("severe_extreme_drought", "moderate_drought", 
                                             "normal", "unusually_moist", "very_extreme_moist"))
data$rebanded_pdsi <- relevel(RESP_DATA_ST$rebanded_pdsi, ref = "normal")




#looking at amount of average cases per month across the yeats
kable(data %>%
        group_by(month) %>%
        summarise(observed = sum(dengue_cases)/5, expected=sum(E)/5), booktabs = T,
      caption = "average monthly cases") %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center")
#looking at how each major regions differ in climatic stressors 
library(table1)
data1 <- data
#adding units to some of the varaibles 
units(data1$water_network) <- "percentage"
units(data1$dengue_cases) <- "per microregion result in region"
units(data1$E) <- "per microregion result in region"

# Subset the data to include only the relevant columns
sub_data <- data1 %>%
  select(year, region_name, dengue_cases, E, water_network, tmin, pdsi, rebanded_pdsi) %>%
  rename(observed_cases = dengue_cases, expected_cases = E)

# Create a table with counts and standard deviations
my_table <- table1(~ rebanded_pdsi + tmin + water_network + observed_cases + expected_cases | region_name, data=sub_data)

# Print the table
my_table
levels(RESP_DATA_ST$rebanded_pdsi)


#the vast fluctuation of the amount of cases per month suggest a poisson distribution

###Analysing SMR### 

#Aggregating data by microregion
data_agg = data%>% group_by(code) %>%
  summarize(observed = sum(dengue_cases),
            expected = sum(E)) %>%
  dplyr::rename(O = observed, E = expected)
#producing SMR using expected and outcome values
data_agg = data_agg %>% mutate(SMR = O/E)
#producing spatial map of aggregeated SMR only looking at SMR right now
data_agg$SMRcat = cut(data_agg$SMR,
                      breaks=c(min(data_agg$SMR),
                               0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6,
                               max(data_agg$SMR)), include.lowest = T)

map_SMR = left_join(shape, data_agg, by = c("code" = "code"))
#plotting map
SMR_plot<-ggplot() + geom_sf(data = map_SMR, col = NA) + aes(fill = SMRcat) +
  theme_bw() + scale_fill_viridis_d() +
  guides(fill=guide_legend(title="SMR"))
SMR_plot


###spatial analysis###

#adding covaraites to data agg 
data_agg <- data %>% 
  group_by(code) %>%
  summarize(observed = sum(dengue_cases),
            expected = sum(E),
            avg_tmin = mean(tmin),
            avg_pdsi = mean(pdsi),
            avg_water_network = mean(water_network)) %>%
  dplyr::rename(O = observed, E = expected)
#rebanding the avg_pdsi
data_agg$rebanded_pdsi <- cut(data_agg$avg_pdsi, 
                          breaks = c(-Inf, -3, -2, 1.9, 2.9, Inf),
                          labels = c("severe_extreme_drought", "moderate_drought", 
                                     "normal", "unusually_moist", "very_extreme_moist"))


ID = seq(1,557)
data_agg <- cbind(ID,data_agg)
library(INLA)
sapply(data_agg,class)
formula_BYM2 = O ~ f(ID,model="bym2", graph=shape.adj,
                     hyper=list(prec = list(
                       prior = "pc.prec",
                       param = c(0.5 / 0.31, 0.01)),
                       phi = list(
                         prior = "pc",
                         param = c(0.5, 2 / 3)))) + avg_pdsi + avg_tmin 

sBYM.model = inla(formula=formula_BYM2, family="poisson", data=data_agg, 
                  offset=log(E), control.compute=list(waic=TRUE),verbose=TRUE)
summary(sBYM.model)
summary(stIntI.BYM.model_1)
#getting posterior probabilities
#Relative risks
RR_sBYM = c()

for(i in 1:557){
  RR_sBYM[i] = inla.emarginal(function(x) exp(x),
                              sBYM.model$marginals.random$ID[[i]])
}

#Posterior probabilities
RR_sBYM_marg = sBYM.model$marginals.random$ID[1:557]
PP_sBYM = lapply(RR_sBYM_marg, function(x) {1-inla.pmarginal(0,x)})
#
resRR_PP = data.frame(resRR=RR_sBYM,
                      PP=unlist(PP_sBYM),
                      code=data_agg[,2])
#plotting posterior probabilirt of RR
resRR_PP$resRRcat = cut(resRR_PP$resRR, breaks=c(min(resRR_PP$resRR),
                                                 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6,
                                                 max(resRR_PP$resRR)),include.lowest = T)
# breakpoints
resRR_PP$PPcat = cut(resRR_PP$PP, c(0, 0.2, 0.8, 1.00), include.lowest = TRUE)#
#making map
map_RR_PP = left_join(shape, resRR_PP, by = c("code" = "code"))
#plotting it with gg plot
p1 = ggplot() + geom_sf(data = map_RR_PP) + aes(fill = resRRcat) +
  theme_bw() + scale_fill_brewer(palette = "PuOr") +
  guides(fill=guide_legend(title="RR")) + ggtitle("RR Spatial model") +
  theme(text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

p2 = ggplot() + geom_sf(data = map_RR_PP) + aes(fill = PPcat) +
  theme_bw() +
  scale_fill_viridis(
    option = "plasma", name="PP",
    discrete = T,
    direction = -1,
    guide = guide_legend(
      title.position = 'top',
      reverse = T
    )) +  ggtitle("PP Spatial model") + theme(text = element_text(size=15),
                                              axis.text.x = element_blank(),
                                              axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

(p1|p2) / (p9/p10)

#looking at spatial fraction
sBYM.model$summary.hyperpar
#time serier aabalysis -----
#0.8 of the spatial variability is explained here showing a strong spatial link
RESP_DATA_ST = left_join(shape, data, by="code")
##Rename the columns of Observed and Expected as we did before
RESP_DATA_ST = RESP_DATA_ST  %>% dplyr::rename(O = dengue_cases)
#Create the ID for year (time)
ID_time
#Create the ID for space
IS_space
#no covariates
formula_ST_noint_no_cov = O ~ f(ID_space, model="bym2", graph=shape.adj,
                         hyper=list(prec = list(
                           prior = "pc.prec",
                           param = c(0.5 / 0.31, 0.01)),
                           phi = list(
                             prior = "pc",
                             param = c(0.5, 2 / 3)))) + f(ID_time,model="rw1",
                                                          hyper=list(prec = list
                                                                     (
                                                                       prior = "pc.prec",
                                                                       param = c(0.5 / 0.31, 0.01)))) 
#mkaing inla
formula_ST_noint = O ~ f(ID_space, model="bym2", graph=shape.adj,
                         hyper=list(prec = list(
                           prior = "pc.prec",
                           param = c(0.5 / 0.31, 0.01)),
                           phi = list(
                             prior = "pc",
                             param = c(0.5, 2 / 3)))) + f(ID_time,model="rw1",
                                                          hyper=list(prec = list
                                                                     (
                               prior = "pc.prec",
                               param = c(0.5 / 0.31, 0.01)))) + tmin + rebanded_pdsi + water_network

stBYM.model = inla(formula=formula_ST_noint, family="poisson",
                   data=RESP_DATA_ST, offset = log(E),
                   control.compute=list(waic=TRUE),verbose = TRUE)
stBYM.model_nc = inla(formula=formula_ST_noint_no_cov, family="poisson",
                   data=RESP_DATA_ST, offset = log(E),
                   control.compute=list(waic=TRUE),verbose = TRUE)
#waic of both formula similar and covaraites selected 
#the interactions
#Spatial Relative risks
RR_stBYM = c()

for(i in 1:557){
  RR_stBYM[i] = inla.emarginal(function(x) exp(x),
                               stBYM.model$marginals.random$ID_space[[i]])
}

#Posterior probabilities (for spatial RR)
RR_stBYM_marg = stBYM.model$marginals.random$ID_space[1:557]
PP_stBYM = lapply(RR_stBYM_marg, function(x) {1-inla.pmarginal(0,x)})

#Temporal Relative risks and CI95
RR_stRW_RR = c()
RR_stRW_lo = c()
RR_stRW_hi = c()

for(i in 1:60){
  #Posterior mean
  RR_stRW_RR[i] = inla.emarginal(function(x) exp(x),
                                 stBYM.model$marginals.random$ID_time[[i]])
  #2.5% quantile
  RR_stRW_lo[i] = inla.qmarginal(0.025,inla.tmarginal(function(x) exp(x), stBYM.model$marginals.random$ID_time[[i]]))
  #97.5% quantile
  RR_stRW_hi[i] = inla.qmarginal(0.975, inla.tmarginal(function(x) exp(x), stBYM.model$marginals.random$ID_time[[i]]))
}

RR_stRW = data.frame(RR=RR_stRW_RR,low=RR_stRW_lo,high=RR_stRW_hi)
#Plot the temporal residual RRs (RR_stWR)
Temp1 = ggplot(RR_stRW, aes(seq(1,60), RR)) + geom_line() +
  ggtitle("ST model No Int") + geom_ribbon(aes(ymin=low,ymax=high), alpha=0.2) + labs(x="month")

Temp1
#plotting porper cycle
ggplot(RR_stRW, aes(seq(1,60), RR)) +
  geom_line() +
  ggtitle("ST model No Int") +
  geom_ribbon(aes(ymin=low,ymax=high), alpha=0.2) +
  labs(x="month/year") +
  scale_x_continuous(breaks = seq(1,60,12),
                     labels = c("Month 1\n2015", "Month 13\n2016", "Month 25\n2017", "Month 37\n2018", "Month 49\n2019"))

#using code to plot on the map
resRR_PP_st = data.frame(resRR=RR_stBYM,
                         PP=unlist(PP_stBYM),
                         code=data_agg[,2])
# breakpoints
resRR_PP_st$resRRcat = cut(resRR_PP_st$resRR, breaks=c(min(resRR_PP_st$resRR),
                                                       0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6,
                                                       max(resRR_PP_st$resRR)),include.lowest = T)

resRR_PP_st$PPcat = cut(resRR_PP_st$PP, c(0, 0.2, 0.8, 1.00), include.lowest = TRUE)

map_RR_ST = left_join(shape, resRR_PP_st, by = c("code" = "code"))
#plotting map
p3 = ggplot() + geom_sf(data = map_RR_ST) + aes(fill = resRRcat) +
  theme_bw() + scale_fill_brewer(palette = "PuOr") +
  guides(fill=guide_legend(title="RR")) +  ggtitle("RR Spatio-temporal model") +
  theme(text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

p4 = ggplot() + geom_sf(data = map_RR_ST) + aes(fill = PPcat) +
  theme_bw() +
  scale_fill_viridis(
    option = "plasma",
    name = "PP Spatio-temporal model",
    discrete = T,
    direction = -1,
    guide = guide_legend(
      title.position = 'top',
      reverse = T
    )) +  ggtitle("PP ST model") + theme(text = element_text(size=15),
                                         axis.text.x = element_blank(),
                                         axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

(p1|p2) / (p3|p4)
summary(stBYM.model)
summary(sBYM.model)
#the temporaral effect seems to explain very little of the varaince

### SPAtio-temporal model type 1 interaction ###
RESP_DATA_ST$ID.space.time = seq(1,dim(RESP_DATA_ST)[1])
formula_ST_intI_nc = O ~ f(ID_space, model="bym2", graph=shape.adj,
                        hyper=list(prec = list(
                          prior = "pc.prec",
                          param = c(0.5 / 0.31, 0.01)),
                          phi = list(
                            prior = "pc",
                            param = c(0.5, 2 / 3)))) +
  f(ID_time,model="rw1", hyper=list(prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01))))+
  f(ID.space.time,model="iid", hyper=list(prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01)))) 

formula_ST_intI = O ~ f(ID_space, model="bym2", graph=shape.adj,
                        hyper=list(prec = list(
                          prior = "pc.prec",
                          param = c(0.5 / 0.31, 0.01)),
                          phi = list(
                            prior = "pc",
                            param = c(0.5, 2 / 3)))) +
  f(ID_time,model="rw1", hyper=list(prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01))))+
  f(ID.space.time,model="iid", hyper=list(prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01)))) + tmin + rebanded_pdsi + water_network


stIntI1.BYM.model = inla(formula=formula_ST_intI, family="poisson", data=RESP_DATA_ST, offset=log(E),
                        control.compute=list(dic=TRUE, waic=TRUE))
stIntI1.BYM.model_nc = inla(formula=formula_ST_intI_nc, family="poisson", data=RESP_DATA_ST, offset=log(E),
                         control.compute=list(dic=TRUE, waic=TRUE))

summary(stIntI1.BYM.model)
#Spatial Relative risks
RR_stIntI1.BYM = c()

for(i in 1:557){
  RR_stIntI1.BYM[i] = inla.emarginal(function(x) exp(x),
                                    stIntI1.BYM.model$marginals.random$ID_space[[i]])
}

#Posterior probabilities (for spatial RR)
RR_stIntI1.BYM_marg = stIntI1.BYM.model$marginals.random$ID_space[1:557]
PP_stIntI1.BYM = lapply(RR_stIntI1.BYM_marg, function(x) {1-inla.pmarginal(0,x)})

#Temporal Relative risks and CI95
RR_stIntI1.RW_RR = c()
RR_stIntI1.RW_lo = c()
RR_stIntI1.RW_hi = c()


for(i in 1:60){
  #Posterior mean
  RR_stIntI1.RW_RR[i] = inla.emarginal(function(x) exp(x),
                                      stIntI1.BYM.model$marginals.random$ID_time[[i]])
  #2.5% quantile
  RR_stIntI1.RW_lo[i] = inla.qmarginal(0.025,inla.tmarginal(function(x) exp(x), stIntI1.BYM.model$marginals.random$ID_time[[i]]))
  #97.5% quantile
  RR_stIntI1.RW_hi[i] = inla.qmarginal(0.975, inla.tmarginal(function(x) exp(x), stIntI1.BYM.model$marginals.random$ID_time[[i]]))
}

RR_stIntI1.RW = data.frame(RR=RR_stIntI1.RW_RR,low=RR_stIntI1.RW_lo,high=RR_stIntI1.RW_hi)
#plotting the residual RR
Temp2 = ggplot(RR_stIntI1.RW, aes(seq(1,60), RR)) + geom_line() + ggtitle("ST model Int I") + geom_ribbon(aes(ymin=low,ymax=high), alpha=0.2) + labs(x="month")

Temp1 | Temp2
#with intercation we see an even higher relative risk than before
resRR_PP_stIntI1 = data.frame(resRR=RR_stIntI1.BYM,
                             PP=unlist(PP_stIntI1.BYM),
                             code=data_agg[,2])
# breakpoints
resRR_PP_stIntI1$resRRcat = cut(resRR_PP_stIntI1$resRR, breaks=c(min(resRR_PP_stIntI1$resRR),
                                                               0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6,
                                                               max(resRR_PP_stIntI1$resRR)),include.lowest = T)

resRR_PP_stIntI1$PPcat = cut(resRR_PP_stIntI1$PP, c(0, 0.2, 0.8, 1.00), include.lowest = TRUE)

map_RR_ST.IntI1 = left_join(shape, resRR_PP_stIntI1, by = c("code" = "code"))
#plotting these maps and comparing
Interaction_1_RR = ggplot() + geom_sf(data = map_RR_ST.IntI1) + aes(fill = resRRcat) +
  theme_bw() + scale_fill_brewer(palette = "PuOr") +
  guides(fill=guide_legend(title="RR")) +  ggtitle("RR ST model Int I") +
  theme(text = element_text(size=15),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

Interaction_1_PP = ggplot() + geom_sf(data = map_RR_ST.IntI1) + aes(fill = PPcat) +
  theme_bw() +
  scale_fill_viridis(
    option = "plasma",
    name = "PP ST model Int I",
    discrete = T,
    direction = -1,
    guide = guide_legend(
      title.position = 'top',
      reverse = T
    )) +  ggtitle("PP ST model Int I") + theme(text = element_text(size=15),
                                               axis.text.x = element_blank(),
                                               axis.text.y = element_blank(), plot.title = element_text(size = 12, face = "bold"))

(Interaction_1_RR|Interaction_1_PP)
summary(stIntI.BYM.model)
#when looking at the temporal interaction we are able to interpret the interactive relationship more
#plotting the space time interaction
RESP_DATA_ST$intI = stIntI.BYM.model$summary.random$ID.space.time$mean

RESP_DATA_ST$intI_cat = cut(RESP_DATA_ST$intI,  breaks=c(-5,-0.05,
                                                         -0.01, 0.01, 0.05, 1,Inf),include.lowest = T)

ggplot() +
  geom_sf(data = RESP_DATA_ST, aes(fill = intI_cat))+ theme_bw() +  scale_fill_brewer(palette = "PuOr") +
  guides(fill=guide_legend(title=NULL)) +
  theme(text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  facet_wrap(~ year, ncol = 3, labeller=labeller(ID.year=c("1"="2015","2"="2016","3"="2017","4"="2018","5"="2019"))) +
  labs("")
ggplot() +
  geom_sf(data = RESP_DATA_ST, aes(fill = intI_cat)) +
  scale_fill_brewer(palette = "PuOr") +
  guides(fill=guide_legend(title=NULL)) +
  theme_bw() +
  theme(text = element_text(size=20),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  facet_wrap(~ year + month, ncol = 12, nrow = 5,
             labeller = labeller(mget(c("year", "month"), envir = as.environment(RESP_DATA_ST)))) +
  labs("")    
#creating table
dat.hyper2 =
  round(
    data.frame(median = stIntI.BYM.model$summary.hyperpar[,4],
               LL = stIntI.BYM.model$summary.hyperpar[,3],
               UL = stIntI.BYM.model$summary.hyperpar[,5]),
    digits = 3)

row.names(dat.hyper2) =
  rownames(stIntI.BYM.model$summary.hyperpar)

knitr::kable(dat.hyper2, caption = "Posterior median and 95% CrI of hyperparameters.") %>%  
  kable_styling(bootstrap_options = "striped", full_width = F, position = "center")

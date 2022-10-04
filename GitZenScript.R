### Loading packages ####
library(MRFcov)
library(tidyverse)
library(data.table)
library(igraph)
library(graphlayouts)
library(ggraph)
library(influential)
library(skimr)
library(car) 
library(sf)
library(grid)
library(ROCit)

# setwd()
source('prettyMRFfull.R')

# load("fulldf.Rda)
# load("Crf025.Rda")
# load("CRFdata025.Rda")
# load(specieslist.Rda)
# load("speciesClassification.Rda")
# fishtherm <- read_csv("fishthermalR.csv")


########### Fit CRF model ##############
myst=Sys.time()
CRFboot025 <- bootstrap_MRF(data = Crf025,
                            n_nodes = nrow(Occ025),
                            family = 'binomial',
                            sample_prop = 1,
                            n_bootstraps = 100,
                            spatial=TRUE,
                            coords = SP,
                            n_cores = 6)
x1 <- Sys.time() - myst

# save(CRFboot025,file = "CRFboot025.Rda")



# Code prep for joining from different scripts
nodenames025 <- names(CRFboot025$mean_key_coefs)

covnames025 <- names(Crf025)[(length(nodenames025)+1):ncol(Crf025)]

names(CRFboot025$indirect_coef_mean) <- covnames025

specplotref <- data.frame(species=nodenames025,fignam=paste0("S",1:length(nodenames025)))

########### Predictions and model fit ###########
rawpred <- predict_MRF(data=Crf025,
                       MRF_mod=CRFboot025)

predscore <- as.data.frame(rawpred$Probability_predictions) %>%
  mutate(featureid=fulldf$featureid) %>%
  pivot_longer(cols=-featureid,names_to = "species",values_to = "predocc") %>%
  left_join(fulldf %>%
              select(featureid,all_of(nodenames025)) %>%
              pivot_longer(cols=-featureid,names_to = "species",values_to = "actocc"),
            by=c("featureid","species"))


rocsumm <- rocit(score=predscore$predocc,
                 class=predscore$actocc,
                 negref = 0)

J <- rocsumm$TPR - rocsumm$FPR
jmax <- which.max(J)
optcut <- rocsumm$Cutoff[jmax]

modsens <- rocsumm$TPR[jmax]
modspcfc <- 1 - rocsumm$FPR[jmax]

newpred <- as_tibble((rawpred$Probability_predictions > optcut) + 0)

# data frame of observed vs expected
preddf025 <- newpred %>%
  mutate(featureid=fulldf$featureid) %>%
  pivot_longer(cols=-featureid,names_to = "species",values_to = "predocc") %>%
  left_join(fulldf %>%
              select(featureid,all_of(nodenames025)) %>%
              pivot_longer(cols=-featureid,names_to = "species",values_to = "actocc"),
            by=c("featureid","species"))

# ROC by species
specROC <- function(specpar){
  rocdf = filter(predscore,species==specpar)
  roclist = rocit(score=rocdf$predocc,
                  class=rocdf$actocc,
                  negref = 0)
  Jfn <- roclist$TPR - roclist$FPR
  jmaxfn <- which.max(Jfn)
  optcutfn <- roclist$Cutoff[jmaxfn]
  sensfn <- roclist$TPR[jmaxfn]
  specfn <- 1-roclist$FPR[jmaxfn]
  AUCfn <- roclist$AUC
  
  return(list(AUC=AUCfn,Cutoff=optcutfn,Sensitivity=sensfn,Specificity=specfn))
}

specROCdf <- tibble(species=nodenames025) %>%
  rowwise() %>%
  mutate(n=sum(Crf025[species]),
         sprocl=list(specROC(species))) %>%
  unnest_wider(sprocl)

fullpred <- rbind(
  data.frame(species="Full model",
             n=nrow(Crf025),
             AUC=rocsumm$AUC,
             Cutoff=optcut,
             Sensitivity=modsens,
             Specificity=modspcfc),
  specROCdf
)

save(fullpred,file="fullpred.Rda")

########### Node labels for figures ##########

left_join(fullpred,specplotref,by="species") %>%
  select(species,fignam,n,AUC,Sensitivity,Specificity) %>%
  rowwise() %>%
  mutate(species=paste0("\\",species)) %>%
  knitr::kable(format='latex',linesep="",digits = 3,
               booktabs=TRUE,escape = FALSE,longtable=TRUE) 


########## Key coefficients ###########

abiotics <- covnames025
keycoefs <- bind_rows(CRFboot025$mean_key_coefs, .id = "column_label") %>%
  rowwise() %>%
  mutate(huckey=str_detect(Variable,pattern = "^H"),
         abiokey=sum(str_detect(Variable,pattern=abiotics)),
         interactkey=str_detect(Variable,"_"),
         negpos=Mean_coef>=0) %>%
  group_by(column_label) %>%
  mutate(coefrank=rank(abs(Mean_coef)*-1)) %>%
  ungroup() %>%
  arrange(column_label,coefrank)

write.csv(keycoefs,"keycoefs.csv")

########## Abiotics ##############
abiokey <- filter(keycoefs,abiokey==1) %>%
  mutate(abiochr=sub('_[^_]*$', '', Variable))

abiocount <- count(abiokey,abiochr,interactkey) %>%
  mutate(interactkey=ifelse(interactkey,"Interaction","Main")) %>%
  pivot_wider(names_from = interactkey, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(total=sum(Main,Interaction)) %>%
  ungroup() %>%
  arrange(total)

abiogg <- count(abiokey,abiochr,interactkey) %>%
  mutate(interactkey=ifelse(interactkey,"Interaction","Main")) %>%
  ggplot(aes(fill=interactkey, x=n, 
             y=factor(abiochr,
                      levels=abiocount$abiochr))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic",size = 17),
        axis.text.x = element_text(size=16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=18,face="bold"),
        legend.text = element_text(size=16),
        legend.position = c(0.8,0.2),
        legend.background = element_rect(fill='white',color = "black")) +
  labs(y=NULL,x="N key coefficients",fill="Effect type",title="(a) Abiotic factors") +
  scale_y_discrete(labels=c("agriculture"="Agriculture",
                            "annprcpmm"="Annual precipitation",
                            "AreaSqKM"="Drainage area",
                            "barrier"="Barrier",
                            "developed"="Developed",
                            "H02"="Region",
                            "pctsand"="Sand",
                            "slope"="Slope",
                            "maxmn14"="Stream temp",
                            "water"="Water"))


abiokey %>%
  group_by(abiochr) %>%
  summarise(Mean=mean(Mean_coef),
            Med=median(Mean_coef),
            PercPos=mean(Mean_coef>0)) %>%
  ungroup() %>%
  mutate(abiochr=paste0("\\",abiochr),
         PercPos=paste0(round(PercPos,4)*100,"\\%")) %>%
  knitr::kable(format='latex',linesep="",digits = 2,
               booktabs=TRUE,escape = FALSE) 

########### Biotics ###########
biotics <- names(CRFboot025$mean_key_coefs)

biokey <- keycoefs %>%
  rowwise() %>%
  mutate(biokey=sum(str_detect(Variable,pattern=biotics))) %>%
  filter(biokey>0) %>%
  rowwise() %>%
  mutate(biochr=ifelse(!interactkey,Variable,
                       str_extract(Variable,'(?<=\\_).*')))

biocount <- count(biokey,biochr,interactkey) %>%
  mutate(interactkey=ifelse(interactkey,"Interaction","Main")) %>%
  pivot_wider(names_from = interactkey, values_from = n, values_fill = 0) %>%
  rowwise() %>%
  mutate(total=sum(Main,Interaction)) %>%
  ungroup() %>%
  arrange(total) %>%
  left_join(select(fullspecTSN,nodename,full=species),by=c("biochr"="nodename")) %>%
  mutate(full=str_to_sentence(full))

biogg <- count(biokey,biochr,interactkey) %>%
  mutate(interactkey=ifelse(interactkey,"Interaction","Main")) %>%
  ggplot(aes(fill=interactkey, x=n, 
             y=factor(biochr,
                      levels=biocount$biochr,
                      labels=biocount$full))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic",size = 17),
        axis.text.x = element_text(size=16),
        axis.title = element_text(size=18),
        plot.title = element_text(size=18,face="bold"),
        legend.text = element_text(size=16),
        legend.position = c(0.8,0.2),
        legend.background = element_rect(fill='white',color = "black")) +
  labs(y=NULL,x="N key coefficients",fill="Effect type",title = "(b) Biotic factors")


# Full key coefficient plot
gridExtra::grid.arrange(abiogg,biogg,ncol=2)


######### Full network plot ############
biotic <- filter(keycoefs,!interactkey & abiokey==0)

# 
keyedges <- select(biotic,from=column_label,
                   to=Variable,
                   weight=Mean_coef)

keyedges <- as_tibble(CRFboot025$direct_coef_means,rownames="species") %>%
  select(species,all_of(nodenames025)) %>%
  pivot_longer(-species,names_to = "spec2",values_to = "weight") %>%
  filter(abs(weight) >= 0.001) %>%
  select(from="species","to"=spec2,weight)


keynodes <- data.frame(species=specplotref$species,
                       figname=specplotref$fignam) %>%
  left_join(select(fullspecTSN,nodename,family),by=c("species"="nodename")) %>%
  left_join(select(Occ025,nodename,nsamp),by=c("species"="nodename")) %>%
  left_join(fishtherm,by=c("species"="nodename"))


g <- graph_from_data_frame(keyedges,directed=FALSE,vertices=keynodes)
g <- simplify(g)
E(g)$color <- ifelse(E(g)$weight > 0,"blue","red4")
#E(g)$width <- E(g)$weight
#E(g)$weight <- c(scale(E(g)$weight))
#V(g)$size <- influential::degree(g)

# Note there is some randomness to this, plot will look different with each run
# but general structure should remain
glayo <- layout_with_fr(g)

ggraph(g,"manual",x=glayo[,1],y=glayo[,2]) +
  geom_edge_link0(aes(edge_width = weight, edge_colour=color)) +
  geom_node_point(size=2) +
  geom_node_label(aes(label = figname,
                      fill = thermal),
                  family = "serif",
                  repel=TRUE,
                  max.overlaps=71,
                  size=10) +
  scale_edge_width_continuous(range = c(0.2, 2),guide="none") +
  scale_edge_colour_manual(values = c("blue","red4"),guide="none") +
  scale_fill_viridis_d(begin = 0.18) +
  guides(fill = guide_legend(title = "Thermal regime",
                             override.aes = list(label=""))) +
  theme_graph() +
  theme(legend.title = element_text(size=24),
        legend.text = element_text(size=24))



######## Species network plots #########


# Stonecat Noturus flavus 
sum(Crf025$notuflav)
mean(Crf025$notuflav)

notuflavkey <- CRFboot025$mean_key_coefs$notuflav

notuflavabio <- filter(keycoefs,column_label=="notuflav") %>%
  select(Variable,Mean_coef,Rel_importance) %>%
  mutate(intkey=str_detect(Variable,"_"),
         Variable=str_replace(Variable,"H02","HUC"),
         Variable=str_replace(Variable,"maxmn14","maxmn")) %>%
  rowwise() %>%
  mutate(v2=ifelse(!intkey,paste0("\\",Variable),
                   paste0("\\",str_split(Variable,"_")[[1]][1],
                          "\\ \\TT\\ \\",
                          str_split(Variable,"_")[[1]][2],
                          "\\"))) %>%
  select(v2,Mean_coef,Rel_importance)

knitr::kable(notuflavabio,format='latex',digits=3,linesep="",booktabs=TRUE,escape = FALSE) 

notuflavdat <- filter(Crf025,notuflav==1)

summary(notuflavdat$maxmn14)

notuflavplots <- invisible(lapply(X = c(min(notuflavdat$maxmn14),
                                        (min(notuflavdat$maxmn14)+max(notuflavdat$maxmn14))/2,
                                        max(notuflavdat$maxmn14)),
                                  function(j){
                                    prettyMRF(data = Crf025,
                                              MRF_mod = CRFboot025,
                                              covariate = "maxmn14",
                                              noi = specplotref$fignam[specplotref$species=="notuflav"],
                                              node_names = specplotref$fignam,
                                              obscovquan = j,
                                              cutoff = 0.0,
                                              edge_cutoff = 0)
                                  }))

notuflavgg <- invisible(
  lapply(notuflavplots,
         function(k,nodnm="S52"){
           g <- k
           glay <- layout_as_star(g,center = nodnm)
           
           ggraph(g,"manual",x=glay[,1],y=glay[,2]) +
             geom_edge_link0(aes(filter=!(node1.name == nodnm | node2.name == nodnm),
                                 edge_colour=color,
                                 edge_width=width),
                             alpha=0.3,
                             show.legend = FALSE) +
             geom_edge_link0(aes(filter=(node1.name == nodnm | node2.name == nodnm),
                                 edge_colour=color,
                                 edge_width=width),
                             show.legend = FALSE) +
             geom_node_label(aes(label=name),size=6) +
             scale_edge_width_continuous(range = c(1, 4),guide="none") +
             scale_edge_colour_manual(values = c("blue","red4"),guide="none") +
             theme_graph() +
             coord_cartesian(xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))
         }
  )
)

patchwork::wrap_plots(notuflavgg)


# Banded darter Etheostoma zonale 
sum(Crf025$ethezona)
mean(Crf025$ethezona)

ethezonakey <- CRFboot025$mean_key_coefs$ethezona

ethezonaabio <- filter(keycoefs,column_label=="ethezona") %>%
  select(Variable,Mean_coef,Rel_importance) %>%
  mutate(intkey=str_detect(Variable,"_"),
         Variable=str_replace(Variable,"H02","HUC"),
         Variable=str_replace(Variable,"maxmn14","maxmn")) %>%
  rowwise() %>%
  mutate(v2=ifelse(!intkey,paste0("\\",Variable),
                   paste0("\\",str_split(Variable,"_")[[1]][1],
                          "\\ \\TT\\ \\",
                          str_split(Variable,"_")[[1]][2],
                          "\\"))) %>%
  select(v2,Mean_coef,Rel_importance)

knitr::kable(ethezonaabio,format='latex',digits=3,linesep="",booktabs=TRUE,escape = FALSE) 

ethezonadat <- filter(Crf025,ethezona==1)

hist(ethezonadat$developed,breaks=20)

ethezonaplots <- invisible(lapply(X = c(min(ethezonadat$developed),
                                        (min(ethezonadat$developed)+max(ethezonadat$developed))/2,
                                        max(ethezonadat$developed)),
                                  function(j){
                                    prettyMRF(data = Crf025,
                                              MRF_mod = CRFboot025,
                                              covariate = "developed",
                                              noi = specplotref$fignam[specplotref$species=="ethezona"],
                                              node_names = specplotref$fignam,
                                              obscovquan = j,
                                              cutoff = 0.0,
                                              edge_cutoff = 0)
                                  }))

ethezonagg <- invisible(
  lapply(ethezonaplots,
         function(k,nodnm="S23"){
           g <- k
           glay <- layout_as_star(g,center = nodnm)
           
           ggraph(g,"manual",x=glay[,1],y=glay[,2]) +
             geom_edge_link0(aes(filter=!(node1.name == nodnm | node2.name == nodnm),
                                 edge_colour=color,
                                 edge_width=width),
                             alpha=0.3,
                             show.legend = FALSE) +
             geom_edge_link0(aes(filter=(node1.name == nodnm | node2.name == nodnm),
                                 edge_colour=color,
                                 edge_width=width),
                             show.legend = FALSE) +
             geom_node_label(aes(label=name),size=6) +
             scale_edge_width_continuous(range = c(1, 4),guide="none") +
             scale_edge_colour_manual(values = c("blue","red4"),guide="none") +
             theme_graph() +
             coord_cartesian(xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))
         }
  )
)

patchwork::wrap_plots(ethezonagg)

# Creek chub Semotilus atromaculatus

sum(Crf025$semoatro)
mean(Crf025$semoatro)

semoatrokey <- CRFboot025$mean_key_coefs$semoatro

semoatroabio <- filter(keycoefs,column_label=="semoatro") %>%
  select(Variable,Mean_coef,Rel_importance) %>%
  mutate(intkey=str_detect(Variable,"_"),
         Variable=str_replace(Variable,"H02","HUC"),
         Variable=str_replace(Variable,"maxmn14","maxmn")) %>%
  rowwise() %>%
  mutate(v2=ifelse(!intkey,paste0("\\",Variable),
                   paste0("\\",str_split(Variable,"_")[[1]][1],
                          "\\ \\TT\\ \\",
                          str_split(Variable,"_")[[1]][2],
                          "\\"))) %>%
  select(v2,Mean_coef,Rel_importance)

knitr::kable(semoatroabio,format='latex',digits=3,linesep="",booktabs=TRUE,escape = FALSE) 

semoatrodat <- filter(Crf025,semoatro==1)

summary(semoatrodat$AreaSqKM)

semoatroplots <- invisible(lapply(X = c(min(semoatrodat$AreaSqKM),
                                        (min(semoatrodat$AreaSqKM)+sort(semoatrodat$AreaSqKM,decreasing = T)[2])/2,
                                        max(semoatrodat$AreaSqKM)),
                                  function(j){
                                    prettyMRF(data = Crf025,
                                              MRF_mod = CRFboot025,
                                              covariate = "AreaSqKM",
                                              noi = specplotref$fignam[specplotref$species=="semoatro"],
                                              node_names = specplotref$fignam,
                                              obscovquan = j,
                                              cutoff = 0.0,
                                              edge_cutoff = 0)
                                  }))

semoatrogg <- invisible(
  lapply(semoatroplots,
         function(k,nodnm="S70"){
           g <- k
           glay <- layout_as_star(g,center = nodnm)
           
           ggraph(g,"manual",x=glay[,1],y=glay[,2]) +
             geom_edge_link0(aes(filter=!(node1.name == nodnm | node2.name == nodnm),
                                 edge_colour=color,
                                 edge_width=width),
                             alpha=0.3,
                             show.legend = FALSE) +
             geom_edge_link0(aes(filter=(node1.name == nodnm | node2.name == nodnm),
                                 edge_colour=color,
                                 edge_width=width),
                             show.legend = FALSE) +
             geom_node_label(aes(label=name),size=5) +
             scale_edge_width_continuous(range = c(1, 6),guide="none") +
             scale_edge_colour_manual(values = c("blue","red4"),guide="none") +
             theme_graph() +
             coord_cartesian(xlim=c(-1.1,1.1), ylim=c(-1.1,1.1))
         }
  )
)

patchwork::wrap_plots(semoatrogg)

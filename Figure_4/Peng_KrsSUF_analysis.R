
# ------------------------------------
# Peng et al. 2023
# SUF - Krs script

# Created by Adrien Heymans, feb 2023.
# Use Hydraulic Viper pipeline : https://github.com/HydraulicViper/HydraulicViper/tree/main

# CPlantBox parameter from Felix Bauer, nov 2022.
# ------------------------------------

# Lib
library(plyr)
library(Matrix)
library(tidyverse)
library(data.table)
library(readxl)
'%!in%' <- compose('!', '%in%')

setwd("~/GitHub/Peng2023/")
source("./HydraulicViper/getSUF.R") # MARSHAL (Meunier et al. 2020: https://doi.org/10.1093/insilicoplants/diz012)

# CPlantBox parameter
data_in = read.csv("./www/Felixdesign.csv")

# Load data from Hydraulic Viper
DATA = Krs_d = Kall = RS = NULL
for(i in 1:5){
  
  # Get Macro hydraulic properties
  out = read.csv(paste0("./data/data",i,"/out2_FELIX.csv"))
  
  Krs_d = rbind(Krs_d, out%>%mutate(rep = i))
  # Get Root system volume - proxy for carbon investement
  Cc = read.csv(paste0("./data/data",i,"/Carbon2_FELIX.csv"))
  
  # Get Root system files
  rs = NULL
  ls <- list.files(paste0("./data/data",i,"/"))
  out_files <- ls[grepl("roots", ls)]
  pb <- txtProgressBar(min = 1, max = length(out_files), style = 3)
  for(j in out_files){
    tmp <- try(read.csv(paste0("./data/data",i,"/",j)),silent = T)
    setTxtProgressBar(pb, which(out_files == j))
    rs = rbind(rs, tmp%>%mutate(sam = j))
  }
  rs  = rs %>%mutate(sam = as.numeric(str_remove(str_sub(sam,7),".csv")))
  RS = rbind(RS, rs%>%mutate(rep =i))
  
  # Get Conductivities files
  Condu = read.csv(paste0("./data/data",i,"/conductivities2_FELIX.csv"))
  Kall = rbind(Kall, Condu%>%mutate(rep = i))
  end = out%>%filter(Time == 9)%>%
    transmute(id = sam, Krs = Krs)
  Cc = Cc%>%filter(age == 8.9)%>%
    transmute(id = sam, Cc= Cc)
  data = left_join(data_in, end, by = "id")
  data = left_join(data, Cc, by = "id")
  DATA = rbind(DATA, data%>%mutate(rep = i))
}

pl =DATA%>%
  ggplot()+
  geom_point(aes(Krs, round(max_b), colour = seminal_a*10), alpha = 0.5, size = 2)+
  ylab("Seminal root number [-]")+
  labs(colour = "Seminal root\n radius [mm]")+
  xlab("Krs [cm4 hPa-1 day-1]")+
  scale_y_continuous(breaks=seq(0,12,2))+
  theme_classic()+
  geom_smooth(aes(Krs, round(max_b)), method = "lm", formula = "y ~ x")+
  viridis::scale_color_viridis()

ggsave(pl, filename = "./img/Felix_Peng_Krs_Sem.svg")


# Create a dummy soil for MARSHAL
soil <- data.frame(id=1:101,
                   z = sort(seq(-100,0,1), decreasing = TRUE),
                   value = 1,
                   psi = rep(-73.5021, 101))

# Run MARSHAL to get SUF 
SUF_data = NULL
for (i in unique(RS$sam)){
  tmp_root= RS%>%filter(sam == i)
  root = hydro =  NULL
  for(k in unique(tmp_root$rep)){
    root = tmp_root%>%filter(rep == k)%>%select(-rep)
    tmp_conduct = Kall%>%filter(sam == i, rep == k)
    hydro = try(getSUF(table_data = root, table_cond = tmp_conduct,
                       table_soil = soil, hetero = FALSE,
                       Psi_collar = -15000,
                       soil_param = NULL), silent = T)
    if(class(hydro) == "try-error"){}else{
      # Get parental information for lateral roots
      root$parent = NA
      for(h in unique(root$branchID)){
        root_type = unique(root$type[root$branchID == h])
        if (root_type == 2){
          alt_rs = root%>%filter(type != 2)
          x1 = root$x1[root$branchID == h][1]
          y1 = root$y1[root$branchID == h][1]
          z1 = root$z1[root$branchID == h][1]
          alt_rs$euc = sqrt((alt_rs$x1-x1)^2+ (alt_rs$y1-y1)^2+ (alt_rs$z1-z1)^2)
          parent = alt_rs$type[alt_rs$euc == min(alt_rs$euc)]
          if(length(parent) != 1){
            parent = unique(parent[1])
          }
          root$parent[root$branchID == h] = parent
        }
      }
      root$type[root$parent == 1] = 3
      root$type[root$parent == 4] = 6

      suf_depth = root %>%
        mutate(z = round((z1+z2)/2),
               suf1 = as.vector(hydro$suf1))%>% 
        dplyr::group_by(z, type )%>%
        dplyr::summarise(suf = sum(suf1))%>%
        ungroup()%>%
        mutate(id = i, rep = k)
      SUF_data = rbind(SUF_data, suf_depth)
    }
  }
}
# Group suf per type of root
SUFinfo = SUF_data %>% 
  dplyr::group_by(id,rep,type)%>%
  dplyr::summarise(suf = sum(suf))%>%
  ungroup()
SUFall = left_join(SUFinfo,data_in, by = "id")

# select suf data for the seminal roots
SUFSem = SUFall%>%
  mutate(max_b = ifelse(is.na(max_b), 0, round(max_b)),
         type = ifelse(type %in% c(1,4), "3Seminal", ifelse(type == 2, "2Laterals", "1Nodal")))%>%
  dplyr::group_by(id,rep,type)%>%
  dplyr::summarise(suf = sum(suf), max_b = mean(max_b))%>%
  ungroup()%>%
  filter(suf < 0.99,
         type == "3Seminal")
# regression on suf vs. sem number
fit = lm(suf~max_b, data = SUFSem )
summary(fit)

# suf vs. Seminal root
SUFSem%>%
  ggplot()+
  geom_point(aes(max_b, suf,group = type), size = 5, alpha = 0.3)+
  geom_smooth(aes(max_b, suf,group = type), alpha = 0.5, method = "lm", formula = "y ~ x")+
  xlab("Seminal root number (#)")+
  scale_x_continuous()+
  theme_classic()+
  viridis::scale_colour_viridis(discrete = T)
ggsave("./img/sem_vs_suf.svg", width=10, height=10)

# Stack all suf per seminal root number
g = SUFall%>%
  mutate(N_sem = round(max_b),
         type = ifelse(type %in% c(1,4), "3Seminal", ifelse(type == 2, "2Laterals", "1Nodal")))%>%
  dplyr::group_by(N_sem, type)%>%
  dplyr::summarise(msuf = median(suf))%>% # Get meadian value
  ungroup()
G = NULL
for(i in unique(g$N_sem)){
  tmp = g%>%filter(N_sem == i)
  ss = sum(tmp$msuf)
  tmp$msuf = tmp$msuf/rep(ss, nrow(tmp))
  if(length(tmp$msuf[tmp$type == "1Nodal"])==0){
    tmp = rbind(tmp, tibble(N_sem = tmp$N_sem[1], type = "1Nodal", msuf = 0.0))
  }
  G = rbind(G, tmp)
}

# Create figure for the stacked suf per root type
G%>%
  ggplot(aes(N_sem, msuf, fill = type))+
  geom_area()+
  theme_classic()+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("Seminal root number (#)")+
  ylab("SUF")+
  labs(fill = "root type")

ggsave("./img/suf.svg", width=10, height=10)
ggsave("./img/suf.png", width=10, height=10)

# write suf data.
write.csv(SUF_data, "./data/SUF_Peng.csv")

# No nodal root in data
# Stack all suf per seminal root number 
Id_nodal = unique(SUFall$id[SUFall$type == 5])
g = SUFall%>%
  filter(id %!in% Id_nodal)%>%
  mutate(max_b = ifelse(is.na(max_b), 0, max_b),
         N_sem = round(max_b))%>%
  dplyr::group_by(N_sem, type)%>%
  dplyr::summarise(msuf = median(suf))%>%
  ungroup()

G = NULL
for(i in unique(g$N_sem)){
  tmp = g%>%filter(N_sem == i)
  ss = sum(tmp$msuf)
  tmp$msuf = tmp$msuf/rep(ss, nrow(tmp))
  if(length(tmp$msuf[tmp$type == 4])==0){
    tmp = rbind(tmp, tibble(N_sem = tmp$N_sem[1], type = 4, msuf = 0.0))
    tmp = rbind(tmp, tibble(N_sem = tmp$N_sem[1], type = 6, msuf = 0.0))
  }
  G = rbind(G, tmp)
}
# Create figure for the stacked suf per root type
G%>%
  ggplot(aes(N_sem, msuf, fill = factor(type)))+
  geom_area()+
  theme_classic()+
  viridis::scale_fill_viridis(discrete = T)+
  xlab("Seminal root number (#)")+
  ylab("SUF")+
  labs(fill = "root type")


ggsave("./img/suf_parental.svg", width=10, height=10)



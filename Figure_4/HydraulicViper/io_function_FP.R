# Author: Adrien Heymans
# input - output function to run hydraulic viper version 1.01
# Updated: 16-09-2022

# Get conductivities based on cell arragement and cell hydraulic prop from rf & svm algorithm
GM <- function(MatSobol, sam, shcut = TRUE, all_roots){

  load("./R/GRANAR/rf/svm_model_kr1.RData")
  model_kr1svm <- fit.svm
  
  load("./R/GRANAR/rf/rf_model_kr1.RData")
  model_kr1rf <- fit.rf


  load("./R/GRANAR/rf/svm_model_kr2.RData")
  model_kr2svm <- fit.svm
  load("./R/GRANAR/rf/rf_model_kr2.RData")
  model_kr2rf <- fit.rf


  load("./R/GRANAR/rf/svm_model_kr3.RData")
  model_kr3svm <- fit.svm 
  load("./R/GRANAR/rf/rf_model_kr3.RData")
  model_kr3rf <- fit.rf
  
  v_stele = MatSobol$stele[sam]/100
  v_xylem = MatSobol$xylem[sam]/100
  trans_time_apo =   MatSobol$trans_time_apo[sam]
  trans_time_xyl = MatSobol$trans_time_xyl[sam]
  
  
  CT_root <- all_roots%>%
    dplyr::group_by(K_type)%>%
    dplyr::summarise(radius = mean(radius)*10)%>% # cm to mm
    ungroup()%>%
    mutate(log_RXA = log(pi*radius^2),
           var_stele = v_stele,
           var_xylem = v_xylem,
           RXA = exp(log_RXA),
           log_TSA = -1.421671 + 1.144070 * log_RXA, #Anatomy_proc.Rmd
           TSA = exp(log_TSA)+var_stele*exp(log_TSA),
           log_TSA = log(TSA, exp(1)),
           r_stele = sqrt(TSA/pi),
           log_nX = 2.9116240 + 0.4459025 * log_TSA, #Anatomy_proc.Rmd
           nX = exp(log_nX),
           ratio = (2+0.07456*r_stele*1000)/nX + var_xylem*(2+0.07456*r_stele*1000)/nX, #Anatomy_proc.Rmd
           mod_XVA = 1.395515 - 0.243006 * log_TSA, #Anatomy_proc.Rmd
           MXA = exp(-mod_XVA^2)+var_xylem*exp(-mod_XVA^2),
           log_CW = log(radius-r_stele, exp(1)),
           CF = exp(3.1091221+0.6718735*log_CW),
           OneC = exp(log_CW)/CF,
           oneC = OneC,
           nPX = nX*ratio,
           PXA_1 = 1000^2*(sqrt(radius/35)/10)^2,
           k_protxyl_s = PXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_unM = k_protxyl_s*nPX*200/1E4, # kx when only the proto xylem have their cell wall lignified 
           LMXA = MXA - nPX*PXA_1/1000^2,
           LMXA_1 = LMXA*1000^2/nX,
           k_Mxyl_s = LMXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM, # kx when all xylem elements have their cell wall lignified 
           km = MatSobol$km[sam], 
           kw = MatSobol$kw[sam], 
           kAQP = MatSobol$aquaporine[sam],
           kpl = MatSobol$plasmo[sam], 
           thickness = 1.5,
           aerenchyma = MatSobol$aer[sam],
           TCA = RXA-TSA,
           AA = aerenchyma*TCA,
           vol = pi*radius^2-AA
    )
  
  
  to_analysis <- CT_root%>%
    select(radius, var_stele, var_xylem, aerenchyma, kw, km, kAQP, kpl, thickness, kx_unM, kx_M)%>%
    mutate(sampl_id = CT_root$K_type)
  
  prediction1 <- (predict(model_kr1svm,to_analysis%>%mutate(aerenchyma = 0))+ predict(model_kr1rf,to_analysis%>%mutate(aerenchyma = 0)))/2 
  prediction2 <- (predict(model_kr2svm,to_analysis)+predict(model_kr2rf,to_analysis))/2
  prediction3 <- (predict(model_kr3svm,to_analysis)+predict(model_kr3rf,to_analysis))/2
  data <- to_analysis%>%
    mutate(kr_1 = abs(prediction1),
           kr_2 = abs(prediction2),
           kr_3 = abs(prediction3))

  
  data$K_type <- data$sampl_id
  
  x <- read_xml("CPlantBox/modelparameter/rootsystem/FP.xml")
  elon <- tibble(type = 1:5, r = rep(0,5), la = rep(0,5))
  for(i in 2:6){
    tmp <- xml_children(x)[i]
    temp_r <- xml_find_all(tmp, ".//parameter")[12]
    r <- xml_attr(temp_r, "value")
    elon$r[i-1] = r
    temp_la <- xml_find_all(tmp, ".//parameter")[8]
    la <- xml_attr(temp_la, "value")
    elon$la[i-1] = la
  }
  
  elon$la <- as.numeric(elon$la)
  elon$r <- as.numeric(elon$r)
  
  conductivities <- NULL
  for(i in CT_root$K_type){
    
    la <- elon$la[all_roots$type[all_roots$K_type == i][1]]
    r <- elon$r[all_roots$type[all_roots$K_type == i][1]]
    
    x_kr = c(0, la/r, la/r+trans_time_apo, la/r+2*trans_time_apo, la/r+3*trans_time_apo, la/r+6*trans_time_apo)
    x_kx = sort(c(0, x_kr[3] + MatSobol$gap[sam],  x_kr[3] + MatSobol$gap[sam] + trans_time_xyl, x_kr[3]+MatSobol$gap[sam] + 5*trans_time_xyl))
    
    tmp_kr1 = data$kr_1[data$K_type == i]
    tmp_kr2 = data$kr_2[data$K_type == i]
    tmp_kr3 = data$kr_3[data$K_type == i]
    if(tmp_kr3 > tmp_kr2){tmp_kr3 = tmp_kr2}
    tmp_kxunM = data$kx_unM[data$K_type == i]
    tmp_kxM = data$kx_M[data$K_type == i]
    
    y = c(tmp_kr1, tmp_kr1,
          tmp_kr2, tmp_kr2,
          tmp_kr3, tmp_kr3/10,
          tmp_kxunM, tmp_kxunM, 
          tmp_kxM, tmp_kxM)
    
    tmp <- tibble(order_id = rep(i,10), order = rep(paste0("type_", i),10), type = c(rep("kr",6), rep("kx", 4)), x = c(x_kr,x_kx), y = y)
    
    conductivities <- rbind(conductivities, tmp)
  }
  
  tmp0 <- tibble(order_id = rep(0,4), order = paste0("type_", 0), type = c(rep("kr",2), rep("kx",2)), x = c(0, 10, 0, 10), 
                 y = c(data$kr_3[data$K_type == 1], data$kr_3[data$K_type == 1], data$kx_M[data$K_type == 1], data$kx_M[data$K_type == 1]))
  conductivities <- rbind(conductivities, tmp0)
  conductivities <- cbind(tibble(id = 1:nrow(conductivities)), conductivities)
  conductivities <- left_join(conductivities, CT_root %>%
                                                 mutate(order = paste0("type_", K_type))%>%
                                                 select(order, vol), by = "order")
  return(conductivities)
  
}
# function to give a conductivity type which is a function of the diameter
K_type <- function(all_roots){
  
  radius_latC <- all_roots$radius[all_roots$type == 2]
  radius_latA <- all_roots$radius[all_roots$type == 3]
  radius_B <- all_roots$radius[all_roots$type == 4]
  radius_SBR <- all_roots$radius[all_roots$type == 5]

  all_roots$K_type = NA
  if(length(radius_SBR)> 0){
    breaks_SBRtype = seq(min(radius_SBR), max(radius_SBR), 
                       by=(max(radius_SBR)-min(radius_SBR))/10)
  all_roots$K_type[all_roots$type == 5 & all_roots$radius < breaks_SBRtype[2]] = 13
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[2] & all_roots$radius < breaks_SBRtype[3]] = 14
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[3] & all_roots$radius < breaks_SBRtype[4]] = 15
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[4] & all_roots$radius < breaks_SBRtype[5]] = 16
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[5] & all_roots$radius < breaks_SBRtype[6]] = 17
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[6] & all_roots$radius < breaks_SBRtype[7]] = 18
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[7] & all_roots$radius < breaks_SBRtype[8]] = 19
  all_roots$K_type[all_roots$type == 5 & all_roots$radius >= breaks_SBRtype[8] & all_roots$radius < breaks_SBRtype[9]] = 20
  all_roots$K_type[all_roots$type == 5 & all_roots$radius > breaks_SBRtype[9] ] = 21
  }
  
  
  
  all_roots$K_type[all_roots$type == 1] = 1
  all_roots$K_type[all_roots$type == 2 & all_roots$radius < quantile(radius_latC)[2]] = 2
  all_roots$K_type[all_roots$type == 2 & all_roots$radius >= quantile(radius_latC)[2] & all_roots$radius < quantile(radius_latC)[3]] = 3
  all_roots$K_type[all_roots$type == 2 & all_roots$radius >= quantile(radius_latC)[3] & all_roots$radius < quantile(radius_latC)[4]] = 4
  all_roots$K_type[all_roots$type == 2 & all_roots$radius >= quantile(radius_latC)[4]] = 5
  
  if(length(radius_latA)> 0){
  all_roots$K_type[all_roots$type == 3 & all_roots$radius < quantile(radius_latA)[2]] = 6
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[2] & all_roots$radius < quantile(radius_latA)[3]] = 7
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[3] & all_roots$radius < quantile(radius_latA)[4]] = 8
  all_roots$K_type[all_roots$type == 3 & all_roots$radius >= quantile(radius_latA)[4]] = 9
  }
  if(length(radius_B)> 0){
  all_roots$K_type[all_roots$type == 4 & all_roots$radius < mean(radius_B)] = 10
  all_roots$K_type[all_roots$type == 4 & all_roots$radius == mean(radius_B)] = 11
  all_roots$K_type[all_roots$type == 4 & all_roots$radius > mean(radius_B)] = 12
  }
  
 
  
  return(all_roots)
}

# Create CplantBox param
CPlant_param <- function(MatSobol, sam, age_max = 30){
    
  
  a <- tibble(type = 1:5, radius = c(ifelse(is.na(MatSobol$tap_a[sam]), 0, MatSobol$tap_a[sam]),
                                     ifelse(is.na(MatSobol$lat_a[sam]), 0, MatSobol$lat_a[sam]), 
                                     ifelse(is.na(MatSobol$lat_a[sam]), 0, MatSobol$lat_a[sam]), 
                                     ifelse(is.na(MatSobol$seminal_a[sam]), 0, MatSobol$seminal_a[sam]), 
                                     ifelse(is.na(MatSobol$sb_a[sam]), 0, MatSobol$sb_a[sam])))%>%
    mutate(radius = radius/10) # cm instead of mm
  r <- tibble(type = 1:5, r = c(ifelse(is.na(MatSobol$tap_r[sam]),0,MatSobol$tap_r[sam]),
                                ifelse(is.na(MatSobol$lat_r[sam]),0,MatSobol$lat_r[sam]),
                                ifelse(is.na(MatSobol$lat_r[sam]),0,MatSobol$lat_r[sam]), 
                                ifelse(is.na(MatSobol$seminal_r[sam]),0,MatSobol$seminal_r[sam]),
                                ifelse(is.na(MatSobol$sb_r[sam]),0,MatSobol$sb_r[sam])))
  
  x <- read_xml("./CPlantBox/modelparameter/rootsystem/B73.xml")
  
  pparam <- xml_children(x)[1]
  ptemp <- xml_find_all(pparam, ".//parameter")
  # xml_attr(ptemp[[2]], "value") <- MatSobol$delayB[sam]
  delayRC = 5
  xml_attr(ptemp[[3]], "value") <- delayRC # MatSobol$delayRC[sam]
  xml_attr(ptemp[[4]], "value") <- ifelse(is.na(MatSobol$delay_s[sam]),1000,MatSobol$delay_s[sam])
  xml_attr(ptemp[[5]], "value") <- ifelse(is.na(MatSobol$first_b[sam]),1000,MatSobol$first_b[sam]) # MatSobol$firstB[sam]
  xml_attr(ptemp[[6]], "value") <- ifelse(is.na(MatSobol$first_s[sam]),1000,MatSobol$first_s[sam])
  xml_attr(ptemp[[7]], "value") <- ifelse(is.na(round(MatSobol$max_b[sam])),0,round(MatSobol$max_b[sam]))
  xml_attr(ptemp[[8]], "value") <- ifelse(is.na(round(MatSobol$n_s[sam])),0,round(MatSobol$n_s[sam]))
  
  xml_attr(ptemp[[13]], "value") <- age_max
  
  for(i in 2:6){
    root <- xml_children(x)[i]
    temp <- xml_find_all(root, ".//parameter")
    if(xml_attr(temp[[3]], "name") != "a"){
      print("encoding miss placement in param.xml")
    }
    xml_attr(temp[[1]], "value") <- 2 # growth function 
    xml_attr(temp[[3]], "value") <- a$radius[i-1]
    xml_attr(temp[[3]], "dev") <- a$radius[i-1]*0.1
    xml_attr(temp[[12]], "value") <- r$r[i-1]
    xml_attr(temp[[12]], "dev") <- r$r[i-1]*0.1
    xml_attr(temp[[8]], "value") <- 46.6*r$r[i-1]/20+2
    xml_attr(temp[[8]], "dev") <- (46.6*r$r[i-1]/20+2)*0.1
  }
  
  Tap <- xml_children(x)[2]
  Ttemp <- xml_find_all(Tap, ".//parameter")
  
  xml_attr(Ttemp[[7]], "value") <- 0.5
  xml_attr(Ttemp[[8]], "value")  = ifelse(is.na(MatSobol$tap_la[sam]),xml_attr(Ttemp[[8]], "value"),MatSobol$tap_la[sam])
  xml_attr(Ttemp[[8]], "dev") = ifelse(is.na(MatSobol$tap_la_sd[sam]),xml_attr(Ttemp[[8]], "dev"),MatSobol$tap_la_sd[sam])
  xml_attr(Ttemp[[10]], "value") <- 1000 # ifelse(is.na(MatSobol$tap_lmax[sam]),1000,MatSobol$tap_lmax[sam])
  xml_attr(Ttemp[[11]], "value") <- ifelse(is.na(MatSobol$tap_ln[sam]),1000,MatSobol$tap_ln[sam])
  xml_attr(Ttemp[[11]], "dev") <- ifelse(is.na(MatSobol$tap_ln_sd[sam]),0.1,MatSobol$tap_ln_sd[sam])
  xml_attr(Ttemp[[17]], "percentage") <- 1 # MatSobol$Plc_T[sam]
  xml_attr(Ttemp[[18]], "percentage") <- 0 # 1 - MatSobol$Plc_T[sam]

  LC <- xml_children(x)[3]
  LCtemp <- xml_find_all(LC, ".//parameter")
  xml_attr(LCtemp[[14]], "value") <- ifelse(is.na(MatSobol$lat_theta[sam]),pi/2,MatSobol$lat_theta[sam])
  xml_attr(LCtemp[[14]], "dev") <- ifelse(is.na(MatSobol$lat_theta_sd[sam]),pi/20,MatSobol$lat_theta_sd[sam])
  
  LA <- xml_children(x)[4]
  LAtemp <- xml_find_all(LA, ".//parameter")
  #xml_attr(LAtemp[[11]], "value") <- ifelse(is.na(MatSobol$lat_ln[sam]),1000,MatSobol$lat_ln[sam])
  #xml_attr(LAtemp[[11]], "dev") <- ifelse(is.na(MatSobol$lat_ln_sd[sam]),0.1,MatSobol$lat_ln_sd[sam])
  #xml_attr(LAtemp[[8]], "value")  = ifelse(is.na(MatSobol$lat_la[sam]),xml_attr(LAtemp[[8]], "value"),MatSobol$lat_la[sam])
  #xml_attr(LAtemp[[8]], "dev") = ifelse(is.na(MatSobol$lat_la_sd[sam]),xml_attr(LAtemp[[8]], "dev"),MatSobol$lat_la_sd[sam])
  xml_attr(LAtemp[[10]], "value") <- ifelse(is.na(MatSobol$lat_lmax[sam]),1000,MatSobol$lat_lmax[sam])
  xml_attr(LAtemp[[10]], "dev") <- ifelse(is.na(MatSobol$lat_lmax_sd[sam]),1,MatSobol$lat_lmax_sd[sam])
  xml_attr(LAtemp[[14]], "value") <- ifelse(is.na(MatSobol$lat_theta[sam]),pi/2,MatSobol$lat_theta[sam])
  xml_attr(LAtemp[[14]], "dev") <- ifelse(is.na(MatSobol$lat_theta_sd[sam]),pi/20,MatSobol$lat_theta_sd[sam])
  
  Basal <- xml_children(x)[5]
  Btemp <- xml_find_all(Basal, ".//parameter")
  xml_attr(Btemp[[7]], "value") <- 0.5
  xml_attr(Btemp[[8]], "value")  = ifelse(is.na(MatSobol$seminal_la[sam]),xml_attr(Btemp[[8]], "value"),MatSobol$seminal_la[sam])
  xml_attr(Btemp[[8]], "dev") = ifelse(is.na(MatSobol$seminal_la_sd[sam]),xml_attr(Btemp[[8]], "dev") , MatSobol$seminal_la_sd[sam])
  xml_attr(Btemp[[11]], "value") <- ifelse(is.na(MatSobol$seminal_ln[sam]),1000,MatSobol$seminal_ln[sam])
  xml_attr(Btemp[[11]], "dev") <- ifelse(is.na(MatSobol$seminal_ln_sd[sam]),1,MatSobol$seminal_ln_sd[sam])
  xml_attr(Btemp[[10]], "value") <- ifelse(is.na(MatSobol$seminal_lmax[sam]),1000,MatSobol$seminal_lmax[sam])
  xml_attr(Btemp[[10]], "dev") <- ifelse(is.na(MatSobol$seminal_lmax_sd[sam]),1,MatSobol$seminal_lmax_sd[sam])
 # xml_attr(Btemp[[14]], "value") <- ifelse(is.na(MatSobol$seminal_theta[sam]),0.75,MatSobol$seminal_theta[sam])
 # xml_attr(Btemp[[14]], "dev") <- ifelse(is.na(MatSobol$seminal_theta_sd[sam]),0.075,MatSobol$seminal_theta_sd[sam])
  xml_attr(Btemp[[17]], "percentage") <- 1 #MatSobol$Plc_B[sam]
  xml_attr(Btemp[[18]], "percentage") <- 0 #1 - MatSobol$Plc_B[sam]
  
  SBR <-  xml_children(x)[6]
  SBtemp <- xml_find_all(SBR, ".//parameter")
  xml_attr(SBtemp[[7]], "value") <- 0.5
  xml_attr(SBtemp[[8]], "value")  = ifelse(is.na(MatSobol$sb_la[sam]),xml_attr(SBtemp[[8]], "value"),MatSobol$sb_la[sam])
  xml_attr(SBtemp[[8]], "dev") = ifelse(is.na(MatSobol$sb_la_sd[sam]),xml_attr(SBtemp[[8]], "dev"),MatSobol$sb_la_sd[sam])
  xml_attr(SBtemp[[11]], "value") <- ifelse(is.na(MatSobol$sb_ln[sam]),1000,MatSobol$sb_ln[sam])
  xml_attr(SBtemp[[11]], "dev") <- ifelse(is.na(MatSobol$sb_ln_sd[sam]),1,MatSobol$sb_ln_sd[sam])
  # xml_attr(SBtemp[[14]], "value") <- ifelse(is.na(MatSobol$sb_theta[sam]),pi/2,MatSobol$sb_theta[sam])
  xml_attr(SBtemp[[17]], "percentage") <- 1 # MatSobol$Plc_SBR[sam]
  xml_attr(SBtemp[[18]], "percentage") <- 0 #1 - MatSobol$Plc_SBR[sam]
  
  write_xml(x, file = paste0("./CPlantBox/modelparameter/rootsystem/FP.xml"))
  make_CPscript(age_max+3)
}

# # Make a CplantBox script
make_CPscript <- function(age_max){

  CPscript <- paste0('"""small example in a container"""
print("Importing libraries")
import sys
sys.path.append("./CPlantBox/")
import plantbox as pb
print("libraries loaded")

rootsystem = pb.RootSystem()
# Open plant and root parameter from a file
path = "./CPlantBox/modelparameter/rootsystem/"
output = "',age_max,'_rootsystem"
name = "FP"

print("Read parameter xml file")
rootsystem.readParameters(path + name + ".xml")

# Create and set geometry
# Creates a soil container
print("Create soil container")
rhizotron = pb.SDF_PlantBox(900, 900, 900)

# Pick 1, or 2
print("Set geometry")
rootsystem.setGeometry(rhizotron)  # soilcore, or rhizotron

# Initialize
print("Initialize")
rootsystem.initialize()

# Simulate
print("Run -- CPlantBox -- ")
rootsystem.simulate(',age_max,')  # days
print("analyse")
ana = pb.SegmentAnalyser(rootsystem)

print("write rootsytem")
ana.write("{}.txt".format(str(output)))
  ')
write(CPscript, "./FPscript.py")
}

# Run the RootSystem model of CPlantBox
CRootBox <- function(sam, MatSobol, age_max = 30){
  
  all_roots <- NULL
  
  if(file.exists(paste0("./",age_max+3,"_rootsystem.txt"))){
    file.remove(paste0("./",age_max+3,"_rootsystem.txt"))
  }
  
  CPlant_param(MatSobol, sam, age_max) # change root radius in B37.xml file

  # make_CPscript(age_max+3)
  message("Run CplantBox")
  system("python ./FPscript.py")

  fc <- file.copy(from = paste0("./",age_max+3,"_rootsystem.txt"),
                  to = paste0("./CPlantBox/RS/RS_",sam,".txt"),
                  overwrite = TRUE)
  
  all_roots <- data.table::fread(paste0("./",age_max+3,"_rootsystem.txt"), header = TRUE)

  all_roots <- all_roots%>%arrange(time)
  all_roots <- all_roots%>%dplyr::mutate(length = sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2))
  all_roots$node2ID <- 1:nrow(all_roots)
  all_roots$node1ID[all_roots$branchID == 1][1] <- 0
  all_roots$node1ID[all_roots$branchID == 1][-1] <- which(all_roots$branchID == 1)[-length(which(all_roots$branchID == 1))] # tap root ordination

  for(i in unique(all_roots$branchID)[-1]){
     all_roots$node1ID[all_roots$branchID == i][-1] <- which(all_roots$branchID == i)[-length(which(all_roots$branchID == 1))]
     if(all_roots$type[all_roots$branchID == i][1] %in% c(4,5)){ # connection with the collar
       all_roots$node1ID[all_roots$branchID == i][1] <- 0
     }
     if(all_roots$type[all_roots$branchID == i][1] %in% c(2,3)){ # connection with the parental root
       x1_child <- all_roots$x1[all_roots$branchID == i][1]
       y1_child <- all_roots$y1[all_roots$branchID == i][1]
       z1_child <- all_roots$z1[all_roots$branchID == i][1]
       
       tmp_time <- all_roots$time[all_roots$branchID == i][1]
       
       nearest <- all_roots%>%filter(branchID != i)%>%
                      mutate(euc = sqrt((x1-x1_child)^2+ (y1 - y1_child)^2 + (z1 - z1_child)^2))
       nearest <- nearest[nearest$euc == min(nearest$euc), ]
       all_roots$node1ID[all_roots$branchID == i][1] <- nearest$node2ID[1] # oldest segments
     }
   }
  oups <- which(all_roots$node1ID == all_roots$node2ID)
  if(length(oups) > 0){
    for(o in oups){
      self_seg_age <- all_roots$time[all_roots$node2ID == o][1]
      self_seg_id <- all_roots$branchID[all_roots$node2ID == o][1]
      
      nearest <- all_roots%>%filter(branchID == self_seg_id, time < self_seg_age)
      all_roots$node1ID[all_roots$node2ID == o][1] <- nearest$node2ID[nearest$time == max(nearest$time)]
    }
  }
  # increasing shoot born root diameter as the numbers of nodes increase
  # first node have 100% of the shoot born radius defined in param
  # after 30 days --> SBR have a radius equals to 300% of the shoot born radius defined in param
  first_SBR <- min(age_max - all_roots$age[all_roots$type == 5])
  
  all_roots$radius[all_roots$type == 5] <- all_roots$radius[all_roots$type == 5]+all_roots$radius[all_roots$type == 5]*((age_max - all_roots$age[all_roots$type == 5] - first_SBR)/22)^1.4 
  # 
  pl = all_roots%>%ggplot()+geom_point(aes(age, radius, colour = type))
  print(pl)
  message("end of RSA formating")
  all_roots <- all_roots[all_roots$time <= age_max,]
  
  return(all_roots) 
}

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))


soil_initial <- function(soil_type = stype, field_capacity = 0.5){
  if(soil_type == "sandy loam"){
      soil <- data.frame(id=1:101,
                         z = sort(seq(-100,0,1), decreasing = TRUE),
                         value = 1,
                         psi = rep(-73.5021, 101))
  }else if (soil_type == "loam") {
      soil <- data.frame(id=1:101,
                         z = sort(seq(-100,0,1), decreasing = TRUE),
                         value = 1,
                         psi = rep(-132.869, 101))

  }else if (soil_type == "silty clay") {
      soil <- data.frame(id=1:101,
                         z = sort(seq(-100,0,1), decreasing = TRUE),
                         value = 1,
                         psi = rep(-354.424, 101))
  }else{
    warning("This soil type is not referenced yet")
  }

  soil$psi <- -200 # soil$psi*(1/field_capacity)

  return(soil)
  
}

Ccost <- function(MatSobol, sam, all_roots){
  
  v_stele = MatSobol$stele[sam]/100
  v_xylem = MatSobol$xylem[sam]/100
  trans_time_apo = MatSobol$trans_time_apo[sam]
  
  #-----------------------------
  CT_root <- all_roots%>%
    dplyr::group_by(K_type)%>%
    dplyr::summarise(radius = mean(radius)*10)%>% # cm to mm
    ungroup()%>%
    mutate(log_RXA = log(pi*radius^2),
           var_stele = v_stele,
           var_xylem = v_xylem,
           RXA = exp(log_RXA),
           log_TSA = -1.421671 + 1.144070 * log_RXA, #Anatomy_proc.Rmd
           TSA = exp(log_TSA)+var_stele*exp(log_TSA),
           log_TSA = log(TSA, exp(1)),
           r_stele = sqrt(TSA/pi),
           log_nX = 2.9116240 + 0.4459025 * log_TSA, #Anatomy_proc.Rmd
           nX = exp(log_nX),
           ratio = (2+0.07456*r_stele*1000)/nX, #Anatomy_proc.Rmd
           mod_XVA = 1.395515 - 0.243006 * log_TSA, #Anatomy_proc.Rmd
           MXA = exp(-mod_XVA^2)+var_xylem*exp(-mod_XVA^2),
           log_CW = log(radius-r_stele, exp(1)),
           CF = exp(3.1091221+0.6718735*log_CW),
           OneC = exp(log_CW)/CF,
           oneC = OneC,
           nPX = nX*ratio,
           PXA_1 = 1000^2*(OneC/2.2)^2,
           k_protxyl_s = PXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_unM = k_protxyl_s*nPX*200/1E4, # kx when only the proto xylem have their cell wall lignified 
           LMXA = MXA - nPX*PXA_1/1000^2,
           LMXA_1 = LMXA*1000^2/nX,
           k_Mxyl_s = LMXA_1^2/(8*pi*200*1E-5/3600/24)*1E-12,
           kx_M = k_Mxyl_s*nX*200/1E4 + kx_unM, # kx when all xylem elements have their cell wall lignified 
           km = MatSobol$km[sam], 
           kw = MatSobol$kw[sam], 
           kAQP = MatSobol$aquaporine[sam],
           kpl = MatSobol$plasmo[sam], 
           thickness = 1.5,
           aerenchyma = MatSobol$aer[sam],
           TCA = RXA-TSA,
           AA = aerenchyma*TCA,
           vol = pi*radius^2-AA
    )
  #----------------------------
  
  x <- read_xml("CPlantBox/modelparameter/rootsystem/B73.xml")
  elon <- tibble(type = 1:5, r = rep(0,5), la = rep(0,5))
  for(i in 2:6){
    tmp <- xml_children(x)[i]
    temp_r <- xml_find_all(tmp, ".//parameter")[12]
    r <- xml_attr(temp_r, "value")
    elon$r[i-1] = r
    temp_la <- xml_find_all(tmp, ".//parameter")[8]
    la <- xml_attr(temp_la, "value")
    elon$la[i-1] = la
  }
  
  elon$la <- as.numeric(elon$la)
  elon$r <- as.numeric(elon$r)
  #-------------------------
  C_cost <- NULL
  for(i in CT_root$K_type){
    
    la <- elon$la[all_roots$type[all_roots$K_type == i][1]]
    r <- elon$r[all_roots$type[all_roots$K_type == i][1]]
    
    x = c(0, 2*la, 2*la+r*trans_time_apo, 2*la+2*r*trans_time_apo)
    y = c(CT_root$RXA[CT_root$K_type == i], CT_root$RXA[CT_root$K_type == i],
          CT_root$vol[CT_root$K_type == i], CT_root$vol[CT_root$K_type == i])
    
    tmp <- tibble(x = x, y = y, K_type = i)
    C_cost <- rbind(C_cost, tmp)
  }
  # ---------------------------
  
  order_uni=unique(CT_root$K_type)
  
  Carbon <- NULL
  age_max = max(all_roots$time)
  for(t in seq(0.2,age_max,0.1)){
    table_data <- all_roots%>%
      filter(time <= t)
    Nseg = nrow(table_data)
    order <- table_data$K_type
    seg_age <- table_data$time     # segment age
    seg_age <- max(table_data$time) - seg_age
    
    Cc=matrix(0,Nseg,1) # radial conductivity of the segments
    for ( i in 1:length(order_uni)) {
      
      pos = is.element(order,order_uni[i])
      od <- order_uni[i]
      
      x = C_cost$x[C_cost$K_type == od]
      y = C_cost$y[C_cost$K_type == od]
      x <- c(x, 5000)
      y <- c(y, y[length(y)])
      
      xout = data.frame(seg_age[pos])
      temp=data.frame(approx(x,y,xout[,1]))
      Cc[pos]=temp[,2]
      
    }
    
    table_data$Cc <- Cc[,1]
    table_data$Cc <- table_data$Cc*table_data$length*10 # mm3
    
    tmpC <- tibble(Cc = sum(table_data$Cc), age = t, sam = sam)
    Carbon <- rbind(Carbon, tmpC)
  }
  return(Carbon)
  
}

add_hydraulics <- function(temp_roots, hydraulics){
  
  # Merge output of MARSHAL on specific root segment
  temp_roots$suf <- as.vector(hydraulics$suf)
  temp_roots$suf_eq <- as.vector(hydraulics$suf_eq)
  temp_roots$suf1 <- as.vector(hydraulics$suf1)
  temp_roots$kx <- as.vector(hydraulics$kx)
  temp_roots$kr <- as.vector(hydraulics$kr)
  temp_roots$jr <- as.vector(hydraulics$jr)
  temp_roots$jr_eq <- as.vector(hydraulics$jr_eq)
  temp_roots$jxl_eq <- as.vector(hydraulics$jxl_eq)
  temp_roots$psi_eq <- as.vector(hydraulics$psi_eq)
  temp_roots$psi <- as.vector(hydraulics$psi)
  temp_roots$jxl <- as.vector(hydraulics$jxl)
  temp_roots$psi_soil <- as.vector(hydraulics$psi_soil)
  return(temp_roots)
}

HM_RWUM_lite <- function(all_roots, sam, conductivities, stype){
    age_max <- round(max(all_roots$time))
    time_sequence = seq(0.4,age_max,0.2)
    soil <- soil_initial(soil_type = stype)

    out = tibble(Time = time_sequence, Krs = 0)
    for(t in time_sequence){

      mp_roots = all_roots%>%filter(time <= t)
      hydraulics <- getSUF(mp_roots, conductivities, soil, hetero =FALSE, 
                         Psi_collar = -15000)
      mp_roots <- add_hydraulics(mp_roots, hydraulics)
      out$Krs[which(out$Time == t)] = hydraulics$krs
    }
  return(out)
}




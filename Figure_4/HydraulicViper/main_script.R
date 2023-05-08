

# Felix - Hydraulic Viper over 9 days

DIR_path <- getwd()
setwd(dir = DIR_path)
JOB_NUMBER = "FELIX"

# library and io_function
source("./R/io_function_FP.R")
source("./R/getSUF.R")
source("./R/dependencies.R")

MatSobol<- read.csv(paste0("./www/Felixdesign.csv"))
MatSobol$km <- 3E-5 # hydraulic conductivity of the membrane without aquaporins. 
MatSobol$kw <- 2.4E-4
MatSobol$aquaporine <- 4.3E-4
MatSobol$plasmo <- 5.3E-12
MatSobol$aer <- 0
MatSobol$trans_time_apo = 2
MatSobol$trans_time_xyl = 2
MatSobol$gap = 1
MatSobol$stele = 0
MatSobol$xylem = 0

soil_list <- c("loam")

# allocate table name
HYDRO = CC = MA = OUT = SOIL = NULL
age_max = 9

for(sam in 1:nrow(MatSobol)){ 
    # Run CRootBox ------------------------------
    all_roots <- try(CRootBox(sam, MatSobol, age_max = age_max), silent = T)
    if(class(all_roots) == "try-error"){next()}
    all_roots <- try(K_type(all_roots),silent = T) # gather roots in cluster by their root diameter
    conductivities <- try(GM(MatSobol, sam, shcut = T, all_roots),silent = T)
    if(class(conductivities) == "try-error"){next()}
    carbon <- Ccost(MatSobol, sam, all_roots)
    CC <- rbind(CC, carbon)
    
    pl = all_roots %>%
      ggplot()+
      geom_segment(aes(x = x1, xend = x2, y = z1, yend = z2))+
      coord_fixed()+
      theme_classic()
    print(pl)

    n_try  <- 0
    HYDRO <- rbind(HYDRO, conductivities%>%mutate(sam =  MatSobol$X[sam]))

    for(stype in soil_list){
        marshal <- try(HM_RWUM_lite(all_roots, sam, conductivities,stype),silent = T)
        if(class(marshal) == "try-error"){next()}
       # OUT
        tmp = marshal
        tmp <- tmp %>% mutate(sam = MatSobol$id[sam], soil_type = stype, loc = w)
        OUT = rbind(OUT, tmp)

    }
    write.csv(all_roots, paste0("./res",JOB_NUMBER,"/roots_",MatSobol$id[sam],".csv"))
    write.csv(CC, paste0("./res",JOB_NUMBER,"/Carbon2_",JOB_NUMBER,".csv"))
    write.csv(HYDRO, paste0("./res",JOB_NUMBER,"/conductivities2_",JOB_NUMBER,".csv"))
    write.csv(OUT, paste0("./res",JOB_NUMBER,"/out2_",JOB_NUMBER,".csv"))
}
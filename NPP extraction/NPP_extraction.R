rm(list = ls())
#--------------------------------------------------------------------------------#
# Extracting forest coverage, NPP and cVeg from ISIMIP2b biome simulations       #
#--------------------------------------------------------------------------------#

#Please read the readme for further explanation! 


#set directories
path_isimi2b <- "/p/projects/isimip/isimip/ISIMIP2b/OutputData/biomes/"
path_out <- "/home/maxto/ISIMIP2/ISIMIP2b_Forest/"
path_fig <- paste0(path_out,"ISIMIP2b_forest_",Sys.Date(),"/")
if(!dir.exists(path_fig)){dir.create(path_fig)}

#-----------------#
#### libraries ####
#-----------------#

#.libPaths("/p/projects/open/R/library")
.libPaths("/home/jinchang/R/x86_64-pc-linux-gnu-library/3.6")
library(stats)
library(ncdf4)


#------------------#
#### parameters ####
#------------------#

# set scenarios, model names, and pft names
year_hist <- seq(1861,2005,by=1)
year_future <- seq(2006,2099,by=1)
scen <- c("historical","future")
soc <- c("histsoc","2005soc")
period <- c("1861_2005","2006_2099")
gcm <- c("gfdl-esm2m","hadgem2-es","ipsl-cm5a-lr","miroc5")
#rcp <- c("rcp26","rcp60","rcp85")

model <- c("CARAIB","LPJ-GUESS","LPJmL","ORCHIDEE") #"CLM45","VISIT","DLEM","JULES-ES-55","VEGAS"
model_sub <- tolower(model)

lon <- seq(-179.75,179.75,by=0.5)
lat <- seq(89.75,-89.75,by=-0.5)



## parameters for netcdf-files
fillvalue <- 1e-20; threshold <- 1e8 
londim <- ncdim_def("lon","degrees_east",as.double(lon), longname="Longitude")
latdim <- ncdim_def("lat","degrees_north",as.double(lat), longname="Latitude")
vegetdim <- ncdim_def("model","-",1,unlim=FALSE,
                      create_dimvar=TRUE,longname="Model_Name")#

#-----------------#
#### functions ####
#-----------------#

read_ncvar <- function(x,ncfname) {
  ncin <- nc_open(ncfname)
  var <- x
  ncvar <- ncvar_get(ncin,var)
  nc_close(ncin)
  return(ncvar)
}  


#-----------------------------------#
####           PFTs              ####
#-----------------------------------#

pft_caraib <- c("c3hh","c3dh","c4h","brsuas","brsutecds","brsutewms","brevtecds","brevtewms","brevxs","sds",
                "trs","ndevtecdt","ndevteclt","ndevtedtt","ndevtedttht","ndevstdit","ndsutecdt","ndsustswt","brevdtt","brevdttht",
                "brevstdit","brsutecdt","brsuteclt","brsutewmt","brrgtrt","brevtrt","maizetec","maizetrc","groundnutc","oilrapeseedc",
                "soybeanc","sunflowerc","othertec","othertrc","pulsestec","pulsestrc","ricec","sugarcanec","cerealtec","roottec",
                "cerealtrc","roottrc","c3p","c4p")


# CLM4.5 25 PFTs not include water
#pft_clm45 <- c("not-vegetated",
#               'needleleaf-evergreen-tree-temperate', 'needleleaf-evergreen-tree-boreal', 'needleleaf-deciduous-tree-boreal',
               # 'broadleaf-evergreen-tree-tropical', 'broadleaf-evergreen-tree-temperate', 'broadleaf-deciduous-tree-tropical',
               # 'broadleaf-deciduous-tree-temperate', 'broadleaf-deciduous-tree-boreal', 'broadleaf-evergreen-shrub-temperate',
               # 'broadleaf-deciduous-shrub-temperate', 'broadleaf-deciduous-shrub-boreal', 'c3-arctic-grass',
               # 'c3-non-arctic-grass', 'c4-grass', 'c3-crop-rainfed',
               # 'c3-crop-irrigated', 'corn-rainfed', 'corn-irrigated',
               # 'spring-cereal-temperate-rainfed', 'spring-cereal-temperate-irrigated', 'winter-cereal-temperate-rainfed',
               # 'winter-cereal-temperate-irrigated', 'soybean-rainfed', 'soybean-irrigated')


# LPJ-GUESS 13 PFTs not include baresoil and water
pft_lpj_guess <- c("BNE","BINE","BNS","TeBS","IBS","TeBE","TrBE","TrIBE","TrBR","C3G",
                   "C4G","C3G_agr","C4G_agr")
pft_lpj_guess <- tolower(pft_lpj_guess)


# LPJmL 44 PFTs include baresoil and water
pft_lpjml <- c("bare","biograss-irrigated","biograss-rainfed","biotree-irrigated","biotree-rainfed",
               "boreal-broadleaved-summergreen-tree","boreal-needleleaved-evergreen-tree","c3-perennial-grass","c4-perennial-grass","cas-irrigated",
               "cas-rainfed","lakesrivers","mai-irrigated","mai-rainfed","mgr-irrigated",
               "mgr-rainfed","mil-irrigated","mil-rainfed","nut-irrigated","nut-rainfed",
               "others-irrigated","others-rainfed","pea-irrigated","pea-rainfed","rap-irrigated",
               "rap-rainfed","reservoir","ric-irrigated","ric-rainfed","sgb-irrigated",
               "sgb-rainfed","soy-irrigated","soy-rainfed","sug-irrigated","sug-rainfed",
               "sun-irrigated","sun-rainfed","temperate-broadleaved-evergreen-tree","temperate-broadleaved-summergreen-tree","temperate-needleleaved-evergreen-tree",
               "tropical-broadleaved-evergreen-tree","tropical-broadleaved-raingreen-tree","whe-irrigated","whe-rainfed")

pft_grass_LPJmL <- c("c4-perennial-grass", "c3-perennial-grass", "bare-ground") # possible to extend if necessary

pft_lpjml_var <- c("bare","biograss_irrigated","biograss_rainfed","biotree_irrigated","biotree_rainfed",
                   "boreal_broadleaved_summergreen_tree","boreal_needleleaved_evergreen_tree","c3_perennial_grass","c4_perennial_grass","cas_irrigated",
                   "cas_rainfed","lakesrivers","mai_irrigated","mai_rainfed","mgr_irrigated",
                   "mgr_rainfed","mil_irrigated","mil_rainfed","nut_irrigated","nut_rainfed",
                   "others_irrigated","others_rainfed","pea_irrigated","pea_rainfed","rap_irrigated",
                   "rap_rainfed","reservoir","ric_irrigated","ric_rainfed","sgb_irrigated",
                   "sgb_rainfed","soy_irrigated","soy_rainfed","sug_irrigated","sug_rainfed",
                   "sun_irrigated","sun_rainfed","temperate_broadleaved_evergreen_tree","temperate_broadleaved_summergreen_tree","temperate_needleleaved_evergreen_tree",
                   "tropical_broadleaved_evergreen_tree","tropical_broadleaved_raingreen_tree","whe_irrigated","whe_rainfed")

# ORCHIDEE 17 PFTs include baresoil but not water
pft_orchidee <- c("bare","trbrev","trbrrg","tendev","tebrev","tebrsu","bondev","bobrsu","bondsu","c3gra",
                  "c4gra","c3win","c3sum","c4mai","c4oth","c3pas","c4pas")


# VISIT only provide grid total
pft_visit <- c("","","","","","","","","","",
               "","","","","","","","","","",
               "","","","","")
pft_tree_visit <- rep(0,length(pft_visit))

# DLEM only provide grid total
pft_dlem <- c("","","","","","","","","","",
              "","","","","","","","","","",
              "","","","","")
pft_tree_dlem <- rep(0,length(pft_dlem))

# JULES only provide grid total
pft_jules_es_55 <- c("","","","","","","","","","",
                     "","","","","","","","","","",
                     "","","","","")
pft_tree_jules_es_55 <- rep(0,length(pft_jules_es_55))

# VEGAS only provide grid total
pft_vegas <- c("","","","","","","","","","",
               "","","","","","","","","","",
               "","","","","")
pft_tree_vegas <- rep(0,length(pft_vegas))


#------------------------------------#
#### Calculation of annual NPP values ####
#------------------------------------#

#--------------------# 
#### NPP softwood #### 
#--------------------#

#set softwood PFTs for models 

pft_tree_caraib <- rep(0,length(pft_caraib))
pft_tree_caraib[c(12:18)] <- 1  

#pft_tree_clm45 <- rep(0,length(pft_clm45))
#pft_tree_clm45[c(2:4)] <- 1  

pft_tree_lpj_guess <- rep(0,length(pft_lpj_guess))
pft_tree_lpj_guess[c(1:3)] <- 1  

pft_tree_lpjml <- rep(0,length(pft_lpjml))
pft_tree_lpjml[c(7,40)] <- 1  

pft_tree_orchidee <- rep(0,length(pft_orchidee))
pft_tree_orchidee[c(4,7,9)] <- 1  

pft <- list(pft_caraib,pft_lpj_guess,pft_lpjml,pft_orchidee)
pft_tree <- list(pft_tree_caraib,pft_tree_lpj_guess,pft_tree_lpjml,pft_tree_orchidee)

# read CLM45 fixed pft maps (2005soc)
#ncfname <- paste0(path_clm45_pft,"clm45_2005soc_pft_global_static.nc")
#clm45_pft_fixed <- read_ncvar("pft_pft",ncfname)
#clm45_pft_fixed <- clm45_pft_fixed/100. # from 0-100 to 0-1

#---------------------------#
# historical (1861-2005)----#
#---------------------------#

scen <- "historical"
soc <- "histsoc"
year <- year_hist
period <- "1861_2005"

for (i in 1:length(model)) { 
  for (j in 1:length(gcm)) {
    
          ### setting of PFTs ###
          tmp_pft <- pft[[i]]
          tmp_pft_tree <- pft_tree[[i]]
          tmp_pft <- tmp_pft[tmp_pft_tree == 1]
          
          ### create the array to store values ### (lon,lat,years)
          tmp_pft_annual <- array(NA,c(720,360,length(year)))
          tmp_npp_annual <- array(NA,c(720,360,length(year)))
          
#for each softwood PFT the area and npp values are extracted from the ISIMIP folder 
            for (p in 1:length(tmp_pft)) {
              
# read area 
                ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                                  tolower(model[i]),"_",gcm[j],"_ewembi_",scen,"_",
                                  soc,"_co2_pft-",tmp_pft[p],"_global_annual_",period,".nc4")
                if (!file.exists(ncfname)) {
                  print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                } else {
                  print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                  tmp_pft_annual <- read_ncvar(paste0("pft-",tmp_pft[p]),ncfname)
                  }
          
          
# read npp
                ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                                tolower(model[i]),"_",gcm[j],"_ewembi_",scen,"_",
                                soc,"_co2_npp-",tmp_pft[p],"_global_annual_",period,".nc4")
                #print(ncfname)
                if (!file.exists(ncfname)) {
                print(paste("NO: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                } else {
                  print(paste("YES: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                
                  tmp_npp_annual <- read_ncvar(paste0("npp-",tmp_pft[p]),ncfname)
                }
                
              
# Now the annual NPP values of a grid cell is calculated in two differnt ways, 
# one for LPJmL and one for the other models (adjustable -- just change the names in if commands)
                
                
          if (model[i] %in% c("CARAIB", "ORCHIDEE", "LPJ-GUESS")){
          
          #npp*area
          tmp0 <- tmp_pft_annual*tmp_npp_annual 
          
          #npp1*area1 + npp2*area2 ...
          if(p == 1) {tmp1 <- tmp0} else {
          tmp1 <- tmp1 + tmp0
          }
          
          # sum of pft area
          if(p == 1) {tmp2 <- tmp_pft_annual} else {
            tmp2 <- tmp2 + tmp_pft_annual
          } #final calculation further down 
          
          
          

          } else if (model[i] %in% "LPJmL") {
          
          # npp sum over pfts
          if(p == 1) {npp_forest <- tmp_npp_annual} else {
            npp_forest <- npp_forest + tmp_npp_annual
          }
          # area sum over pfts
          if(p == 1) {pft_forest <- tmp_pft_annual} else {
            pft_forest <- pft_forest + tmp_pft_annual
          }
          
          }# end LPJmL
          }#end pfts

          
# SPECIAL FOR LPJmL: read grass pft area for LPJmL:
          if(model[i] == "LPJmL") {
            
            for(u in 1:length(pft_grass_LPJmL)){
            
            ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                              tolower(model[i]),"_",gcm[j],"_ewembi_",scen,"_",
                              soc,"_co2_pft-",pft_grass_LPJmL[u],"_global_annual_",period,".nc4")
            if (!file.exists(ncfname)) {
              print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
            } else {
              print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
              
              tmp_pft_grass <- read_ncvar(paste0("pft-",pft_grass_LPJmL[p]),ncfname)
            }
            
            # sum over grass area 
            if(u == 1) {pft_grass <- tmp_pft_grass} else {
              pft_grass <- pft_grass + tmp_pft_grass
            }
            
            } # end pft grass
            
            
# Final calculation for LPJmL:
          tmp_annual_hist_small <- npp_forest*sum(pft_forest, pft_grass)/pft_forest
          }
          
# Final calculation for CARAIB, LPJ-GUESS and ORCHIDEE:
          else if (model[i] %in% c("CARAIB", "ORCHIDEE", "LPJ-GUESS")) {tmp_annual_hist_small <- tmp1/tmp2}

          rm(tmp0,tmp1,tmp2,tmp_pft_annual,tmp_npp_annual)
        
        
        #-----------------------------#
        # create the netcdf file -----#
        #-----------------------------#
        
        ### define dimensions ###
        tunit <- paste("years since ",year[1],"-01-01 00:00:00",sep="")
        timedim <- ncdim_def("time_counter",tunit,1:length(year),unlim=TRUE,
                             create_dimvar=TRUE,calendar="noleap")
        v1.def <- ncvar_def("npp","-",list(londim,latdim,vegetdim,timedim),
                            fillvalue,"npp",prec="single")
        dlname <- "npp"
        
        
        ### create netcdf ###
        ncfnameout <- paste0(path_fig,"ISIMIP2b_npp_softwood_hist_",model[i],"_",gcm[j],"_",scen,".nc4")
        ncout <- nc_create(ncfnameout,list(v1.def),force_v4=FALSE)
        ncvar_put(ncout,v1.def,tmp_annual_hist_small)
        ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
        ncatt_put(ncout,"lat","axis","Y")
        ncatt_put(ncout,dlname,"missing_value",fillvalue,prec="float")
        ncatt_put(ncout,0,"contact","changjf@zju.edu.cn")
        ncatt_put(ncout,0,"source",paste0("Source: ISIMIP2b biome simulations; Model name", model[i]))
        nc_close(ncout)
   
  } #end gcm
} # end model
      
#--------------#
# future ------#
#--------------#
        
scen <- "future"
soc <- "2005soc"
year <- year_future
period <- "2006_2099"
        
for (i in 1:length(model)) { 
  
  ### choose the softwood PFTs for the model ###
  tmp_pft <- pft[[i]]
  tmp_pft_tree <- pft_tree[[i]]
  tmp_pft <- tmp_pft[tmp_pft_tree == 1]
  
  ### create the array to store values ### (lon,lat,years)
  tmp_pft_annual <- array(NA,c(720,360,length(year)))
  tmp_npp_annual <- array(NA,c(720,360,length(year)))
  
  ### select rcps ###
  if (model[i] == "CARAIB") {rcp <- c("rcp26", "rcp60")} else { 
    rcp <- c("rcp26", "rcp60", "rcp85")}
  
   for (j in 1:length(gcm)) {
        for (m in 1:length(rcp)) {
          
# for each softwood PFT the area and npp values are extracted from the ISIMIP folder 
              for (p in 1:length(tmp_pft)) {
                
# read area
                  ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                                    tolower(model[i]),"_",gcm[j],"_ewembi_",rcp[m],"_",
                                    soc,"_co2_pft-",tmp_pft[p],"_global_annual_",period,".nc4")
                  if (!file.exists(ncfname)) {
                    print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                  } else {
                    print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                    if (model[i] == "LPJmL" & rcp[m] == "rcp85") { # LPJmL has wield naming for rcp85
                      tmp_pft_annual <- read_ncvar(paste0("pft_",pft_lpjml_var[p]),ncfname)
                    } else {
                      tmp_pft_annual <- read_ncvar(paste0("pft-",tmp_pft[p]),ncfname)
                    }
                  }
                

# read npp 
                ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                                  tolower(model[i]),"_",gcm[j],"_ewembi_",rcp[m],"_",
                                  soc,"_co2_npp-",tmp_pft[p],"_global_annual_",period,".nc4")
                if (!file.exists(ncfname)) {
                  print(paste("NO: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                } else {
                  print(paste("YES: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                  if (model[i] == "LPJmL" & rcp[m] == "rcp85") { # LPJmL has wield naming for rcp85
                    tmp_npp_annual <- read_ncvar(paste0("npp_",pft_lpjml_var[p]),ncfname)
                  } else {
                    tmp_npp_annual <- read_ncvar(paste0("npp-",tmp_pft[p]),ncfname)
                  }
                }
              
# Now the annual NPP values of a grid cell is calculated in two differnt ways, 
# one for LPJmL and one for the other models (adjustable -- just change the names in if commands)            
                
                if (model[i] %in% c("CARAIB", "ORCHIDEE", "LPJ-GUESS")){
                  
                  # npp*area of current PFT
                  tmp0 <- tmp_pft_annual*tmp_npp_annual
                  
                  # npp1*area1 + npp2*area2 ...
                  if(p == 1) {tmp1 <- tmp0} else {
                    tmp1 <- tmp1 + tmp0
                  }
                  
                  # sum of pft area
                  if(p == 1) {tmp2 <- tmp_pft_annual} else {
                    tmp2 <- tmp2 + tmp_pft_annual
                  }
                  
                  
                  
                  
                } else if (model[i] %in% "LPJmL") {
                  
                  # npp sum over pfts
                  if(p == 1) {npp_forest <- tmp_npp_annual} else {
                    npp_forest <- npp_forest + tmp_npp_annual
                  }
                  # area sum over pfts
                  if(p == 1) {pft_forest <- tmp_pft_annual} else {
                    pft_forest <- pft_forest + tmp_pft_annual
                  }
                  
                }# end LPJmL
              } # end pfts
            

# SPECIAL FOR LPJmL: read grass and bareground area for LPJmL:
          if(model[i] == "LPJmL") {
            
            for(u in 1:length(pft_grass_LPJmL)){
              
              ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                                tolower(model[i]),"_",gcm[j],"_ewembi_",rcp[m],scen,"_",
                                soc,"_co2_pft-",pft_grass_LPJmL[u],"_global_annual_",period,".nc4")
              #print(ncfname)
              if (!file.exists(ncfname)) {
                print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
              } else {
                print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
                
                tmp_pft_grass <- read_ncvar(paste0("pft-",pft_grass_LPJmL[p]),ncfname)
              }
              
              # sum over grass area 
              if(u == 1) {pft_grass <- tmp_pft_grass} else {
                pft_grass <- pft_grass + tmp_pft_grass
              }
              
            } # end pft grass
            
            
# Final calculation for LPJmL:
            tmp_annual_hist_small <- npp_forest*sum(pft_forest, pft_grass)/pft_forest
          }
          
# Final calculation for CARAIB, LPJ-GUESS and ORCHIDEE:
          else if (model[i] != "LPJmL") {tmp_annual_hist_small <- tmp1/tmp2}
          
          rm(tmp0,tmp1,tmp2,tmp_pft_annual,tmp_npp_annual)
          
            
            #-----------------------#
            # write the netcdf file #
            #-----------------------#
            
            tunit <- paste("years since ",year[1],"-01-01 00:00:00",sep="")
            timedim <- ncdim_def("time_counter",tunit,1:length(year),unlim=TRUE,
                                 create_dimvar=TRUE,calendar="noleap")
            v1.def <- ncvar_def("npp","-",list(londim,latdim,vegetdim,timedim),
                                fillvalue,"npp",prec="single")
            dlname <- "npp"
            
            ncfnameout <- paste0(path_fig,"ISIMIP2b_npp_softwood_future_",model[i],"_",gcm[j],"_",rcp[m],".nc4")
            ncout <- nc_create(ncfnameout,list(v1.def),force_v4=FALSE)
            ncvar_put(ncout,v1.def,tmp_annual_fut_small)
            ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
            ncatt_put(ncout,"lat","axis","Y")
            ncatt_put(ncout,dlname,"missing_value",fillvalue,prec="float")
            ncatt_put(ncout,0,"contact","changjf@zju.edu.cn")
            ncatt_put(ncout,0,"source",paste0("Source: ISIMIP2b biome simulations; Model name", model[i]))
            nc_close(ncout)
            

          } # end rcp
   }#end gcm
} #end model
      



#--------------------#
#### NPP hardwood ####
#--------------------#

pft_tree_caraib <- rep(0,length(pft_caraib))
pft_tree_caraib[c(19:26)] <- 1  

pft_tree_clm45 <- rep(0,length(pft_clm45))
pft_tree_clm45[c(5:9)] <- 1 

pft_tree_lpj_guess <- rep(0,length(pft_lpj_guess))
pft_tree_lpj_guess[c(4:9)] <- 1  

pft_tree_lpjml <- rep(0,length(pft_lpjml))
pft_tree_lpjml[c(6,38,39,41,42)] <- 1  

pft_tree_orchidee <- rep(0,length(pft_orchidee))
pft_tree_orchidee[c(2,3,5,6,8)] <- 1 

pft <- list(pft_caraib,pft_clm45,pft_lpj_guess,pft_lpjml,pft_orchidee)
pft_tree <- list(pft_tree_caraib,pft_tree_clm45,pft_tree_lpj_guess,pft_tree_lpjml,pft_tree_orchidee)


#---------------------------#
# historical (1861-2005)----#
#---------------------------#

scen <- "historical"
soc <- "histsoc"
year <- year_hist
period <- "1861_2005"

for (i in 1:length(model)) { 
  for (j in 1:length(gcm)) {
    
    ### setting of PFTs ###
    tmp_pft <- pft[[i]]
    tmp_pft_tree <- pft_tree[[i]]
    tmp_pft <- tmp_pft[tmp_pft_tree == 1]
    
    ### create the array to store values ### (lon,lat,years)
    tmp_pft_annual <- array(NA,c(720,360,length(year)))
    tmp_npp_annual <- array(NA,c(720,360,length(year)))
    
# for each hardwood PFT the area and npp values are extracted from the ISIMIP folder 
    for (p in 1:length(tmp_pft)) {
      
# read area 
      ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                        tolower(model[i]),"_",gcm[j],"_ewembi_",scen,"_",
                        soc,"_co2_pft-",tmp_pft[p],"_global_annual_",period,".nc4")
      if (!file.exists(ncfname)) {
        print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
      } else {
        print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
        tmp_pft_annual <- read_ncvar(paste0("pft-",tmp_pft[p]),ncfname)
      }
      
      
# read npp
      ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                        tolower(model[i]),"_",gcm[j],"_ewembi_",scen,"_",
                        soc,"_co2_npp-",tmp_pft[p],"_global_annual_",period,".nc4")
      if (!file.exists(ncfname)) {
        print(paste("NO: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
      } else {
        print(paste("YES: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
        
        tmp_npp_annual <- read_ncvar(paste0("npp-",tmp_pft[p]),ncfname)
      }
      
      
# Now the annual NPP values of a grid cell is calculated in two differnt ways, 
# one for LPJmL and one for the other models (adjustable -- just change the names in if commands)
      
      
      if (model[i] %in% c("CARAIB", "ORCHIDEE", "LPJ-GUESS")){
        
        #npp*area
        tmp0 <- tmp_pft_annual*tmp_npp_annual 
        
        #npp1*area1 + npp2*area2 ...
        if(p == 1) {tmp1 <- tmp0} else {
          tmp1 <- tmp1 + tmp0
        }
        
        # sum of pft area
        if(p == 1) {tmp2 <- tmp_pft_annual} else {
          tmp2 <- tmp2 + tmp_pft_annual
        } #final calculation further down 
        
        
        
        
      } else if (model[i] %in% "LPJmL") {
        
        # npp sum over pfts
        if(p == 1) {npp_forest <- tmp_npp_annual} else {
          npp_forest <- npp_forest + tmp_npp_annual
        }
        # area sum over pfts
        if(p == 1) {pft_forest <- tmp_pft_annual} else {
          pft_forest <- pft_forest + tmp_pft_annual
        }
        
      }# end LPJmL
    }#end pfts
    
    
# SPECIAL FOR LPJmL: read grass pft area for LPJmL:
    if(model[i] == "LPJmL") {
      
      for(u in 1:length(pft_grass_LPJmL)){
        
        ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                          tolower(model[i]),"_",gcm[j],"_ewembi_",scen,"_",
                          soc,"_co2_pft-",pft_grass_LPJmL[u],"_global_annual_",period,".nc4")
        if (!file.exists(ncfname)) {
          print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
        } else {
          print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
          
          tmp_pft_grass <- read_ncvar(paste0("pft-",pft_grass_LPJmL[p]),ncfname)
        }
        
        # sum over grass area 
        if(u == 1) {pft_grass <- tmp_pft_grass} else {
          pft_grass <- pft_grass + tmp_pft_grass
        }
        
      } # end pft grass
      
      
# Final calculation for LPJmL:
      tmp_annual_hist_small <- npp_forest*sum(pft_forest, pft_grass)/pft_forest
    }
    
# Final calculation for CARAIB, LPJ-GUESS and ORCHIDEE:
    else if (model[i] %in% c("CARAIB", "ORCHIDEE", "LPJ-GUESS")) {tmp_annual_hist_small <- tmp1/tmp2}
    
    rm(tmp0,tmp1,tmp2,tmp_pft_annual,tmp_npp_annual)
    
    
    #-----------------------------#
    # create the netcdf file -----#
    #-----------------------------#
    
    ### define dimensions ###
    tunit <- paste("years since ",year[1],"-01-01 00:00:00",sep="")
    timedim <- ncdim_def("time_counter",tunit,1:length(year),unlim=TRUE,
                         create_dimvar=TRUE,calendar="noleap")
    v1.def <- ncvar_def("npp","-",list(londim,latdim,vegetdim,timedim),
                        fillvalue,"npp",prec="single")
    dlname <- "npp"
    
    
    ### create netcdf ###
    ncfnameout <- paste0(path_fig,"ISIMIP2b_npp_hardwood_hist_",model[i],"_",gcm[j],"_",scen,".nc4")
    ncout <- nc_create(ncfnameout,list(v1.def),force_v4=FALSE)
    ncvar_put(ncout,v1.def,tmp_annual_hist_small)
    ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
    ncatt_put(ncout,"lat","axis","Y")
    ncatt_put(ncout,dlname,"missing_value",fillvalue,prec="float")
    ncatt_put(ncout,0,"contact","changjf@zju.edu.cn")
    ncatt_put(ncout,0,"source",paste0("Source: ISIMIP2b biome simulations; Model name", model[i]))
    nc_close(ncout)
    
  } # end gcm
} # end model

#--------------#
# future ------#
#--------------#

scen <- "future"
soc <- "2005soc"
year <- year_future
period <- "2006_2099"

for (i in 1:length(model)) { 
  
  ### choose the hardwood PFTs for the model ###
  tmp_pft <- pft[[i]]
  tmp_pft_tree <- pft_tree[[i]]
  tmp_pft <- tmp_pft[tmp_pft_tree == 1]
  
  ### create the array to store values ### (lon,lat,years)
  tmp_pft_annual <- array(NA,c(720,360,length(year)))
  tmp_npp_annual <- array(NA,c(720,360,length(year)))
  
  ### select rcps ###
  if (model[i] == "CARAIB") {rcp <- c("rcp26", "rcp60")} else { 
    rcp <- c("rcp26", "rcp60", "rcp85")}
  
  for (j in 1:length(gcm)) {
    for (m in 1:length(rcp)) {
      
# for each hardwood PFT the area and npp values are extracted from the ISIMIP folder 
      for (p in 1:length(tmp_pft)) {
        
# read area
        ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                          tolower(model[i]),"_",gcm[j],"_ewembi_",rcp[m],"_",
                          soc,"_co2_pft-",tmp_pft[p],"_global_annual_",period,".nc4")
        if (!file.exists(ncfname)) {
          print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
        } else {
          print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
          if (model[i] == "LPJmL" & rcp[m] == "rcp85") { # LPJmL has wield naming for rcp85
            tmp_pft_annual <- read_ncvar(paste0("pft_",pft_lpjml_var[p]),ncfname)
          } else {
            tmp_pft_annual <- read_ncvar(paste0("pft-",tmp_pft[p]),ncfname)
          }
        }
        
        
# read npp 
        ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                          tolower(model[i]),"_",gcm[j],"_ewembi_",rcp[m],"_",
                          soc,"_co2_npp-",tmp_pft[p],"_global_annual_",period,".nc4")
        if (!file.exists(ncfname)) {
          print(paste("NO: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
        } else {
          print(paste("YES: npp_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
          if (model[i] == "LPJmL" & rcp[m] == "rcp85") { # LPJmL has wield naming for rcp85
            tmp_npp_annual <- read_ncvar(paste0("npp_",pft_lpjml_var[p]),ncfname)
          } else {
            tmp_npp_annual <- read_ncvar(paste0("npp-",tmp_pft[p]),ncfname)
          }
        }
        
# Now the annual NPP values of a grid cell is calculated in two differnt ways, 
# one for LPJmL and one for the other models (adjustable -- just change the names in if commands)            
        
        if (model[i] %in% c("CARAIB", "ORCHIDEE", "LPJ-GUESS")){
          
          # npp*area of current PFT
          tmp0 <- tmp_pft_annual*tmp_npp_annual
          
          # npp1*area1 + npp2*area2 ...
          if(p == 1) {tmp1 <- tmp0} else {
            tmp1 <- tmp1 + tmp0
          }
          
          # sum of pft area
          if(p == 1) {tmp2 <- tmp_pft_annual} else {
            tmp2 <- tmp2 + tmp_pft_annual
          }
          
          
          
          
        } else if (model[i] %in% "LPJmL") {
          
          # npp sum over pfts
          if(p == 1) {npp_forest <- tmp_npp_annual} else {
            npp_forest <- npp_forest + tmp_npp_annual
          }
          # area sum over pfts
          if(p == 1) {pft_forest <- tmp_pft_annual} else {
            pft_forest <- pft_forest + tmp_pft_annual
          }
          
        }# end LPJmL
      } # end pfts
      
      
# SPECIAL FOR LPJmL: read grass and bareground area for LPJmL:
      if(model[i] == "LPJmL") {
        
        for(u in 1:length(pft_grass_LPJmL)){
          
          ncfname <- paste0(path_isimi2b,model[i],"/",gcm[j],"/",scen,"/",
                            tolower(model[i]),"_",gcm[j],"_ewembi_",rcp[m],scen,"_",
                            soc,"_co2_pft-",pft_grass_LPJmL[u],"_global_annual_",period,".nc4")
          if (!file.exists(ncfname)) {
            print(paste("NO: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
          } else {
            print(paste("YES: pft_annual",strsplit(ncfname,path_isimi2b)[[1]][2]))
            
            tmp_pft_grass <- read_ncvar(paste0("pft-",pft_grass_LPJmL[p]),ncfname)
          }
          
          # sum over grass area 
          if(u == 1) {pft_grass <- tmp_pft_grass} else {
            pft_grass <- pft_grass + tmp_pft_grass
          }
          
        } # end pft grass
        
        
# Final calculation for LPJmL:
        tmp_annual_hist_small <- npp_forest*sum(pft_forest, pft_grass)/pft_forest
      }
      
# Final calculation for CARAIB, LPJ-GUESS and ORCHIDEE:
      else if (model[i] != "LPJmL") {tmp_annual_hist_small <- tmp1/tmp2}
      
      rm(tmp0,tmp1,tmp2,tmp_pft_annual,tmp_npp_annual)
      
      
      #-----------------------#
      # write the netcdf file #
      #-----------------------#
      
      tunit <- paste("years since ",year[1],"-01-01 00:00:00",sep="")
      timedim <- ncdim_def("time_counter",tunit,1:length(year),unlim=TRUE,
                           create_dimvar=TRUE,calendar="noleap")
      v1.def <- ncvar_def("npp","-",list(londim,latdim,vegetdim,timedim),
                          fillvalue,"npp",prec="single")
      dlname <- "npp"
      
      ncfnameout <- paste0(path_fig,"ISIMIP2b_npp_hardwood_future_",model[i],"_",gcm[j],"_",rcp[m],".nc4")
      ncout <- nc_create(ncfnameout,list(v1.def),force_v4=FALSE)
      ncvar_put(ncout,v1.def,tmp_annual_fut_small)
      ncatt_put(ncout,"lon","axis","X") #,verbose=FALSE) #,definemode=FALSE)
      ncatt_put(ncout,"lat","axis","Y")
      ncatt_put(ncout,dlname,"missing_value",fillvalue,prec="float")
      ncatt_put(ncout,0,"contact","changjf@zju.edu.cn")
      ncatt_put(ncout,0,"source",paste0("Source: ISIMIP2b biome simulations; Model name", model[i]))
      nc_close(ncout)
      
      
    } # end rcp
  }#end gcm
} #end model




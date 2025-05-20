

# Read ecolinc raster to use as template
r0 <- rast("C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/Manuscripts/dss_paper/sabah_case_study_2024/resist_ecolinc_baseline_0m.tif")
res(r0) # 90x90
dim(r0) # 528 395 1
crs(r0) # Asia_South_Albers_Equal_Area_Conic

# Create raster with values of 1 and save
r1 <- rast(r0, vals=1)
r1[is.na(r1)] <- -9999
writeRaster(r1, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/res_1.tif", overwrite=T)
writeRaster(r1, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/res_1.asc",
            filetype="AAIGrid", NAflag=-9999, overwrite=T)

# Create raster with values of 2 and save
r2 <- rast(r0, vals=2)
r2[is.na(r2)] <- -9999
writeRaster(r2, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/res_2.tif", overwrite=T)
writeRaster(r2, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/res_2.asc",
            filetype="AAIGrid", NAflag=-9999, overwrite=T)

# Create random uniform raster
r3 <- rast(r0, vals=runif(ncell(r1), min=1, max=5))
writeRaster(r3, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/res_3.tif", overwrite=T)
writeRaster(r3, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/res_3.asc",
            filetype="AAIGrid", NAflag=-9999, overwrite=T)

# Create two points 20km apart and one in the bottom right and save
x1 <- c(-920000,-935000,-935000)
x2 <- c(2390000,2420000,2400000)
coordinates <- data.frame(
  ids = c(1,2,3),
  lon = x1,
  lat = x2
)
points <- vect(coordinates, geom = c("lon", "lat"), crs = crs(r0))
writeVector(points, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/pts.shp", overwrite=T)

# Reformat for csv
coordinates <- data.frame(
  X = x1,
  Y = x2
)
write.csv(coordinates, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/pts.xy", quote = T, row.names = F)

# Random points
set.seed(1)
# Create random points
x1 <- runif(20, min=xmin(ext(r1)), max=xmax(ext(r1)))
x2 <- runif(20, min=ymin(ext(r1)), max=ymax(ext(r1)))

coordinates <- data.frame(
  ids = c(1:length(x1)),
  lon = x1,
  lat = x2
)
points <- vect(coordinates, geom = c("lon", "lat"), crs = crs(r0))
writeVector(points, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/pts20.shp", overwrite=T)

# Reformat for csv
coordinates <- data.frame(
  X = x1,
  Y = x2
)
write.csv(coordinates, "C:/Users/pj276/OneDrive - Northern Arizona University/All/Projects/USFSIP_Connectivity/SoftwareDevelopment/cola_unicor_comparison/pts20.xy", quote = F, row.names = F)




input_shp <- "C:/Users/pj276/cola_unicor_comparison/pts.shp"
input_tif <- "C:/Users/pj276/cola_unicor_comparison/res_2.tif"
out_tif <- "C:/Users/pj276/cola_unicor_comparison/res_2_crk_console.tif"

kernels <- crk_py( inshp = input_shp, intif = input_tif, outtif = out_tif,
                   param4 =  25000, 
                   param5 = 'linear',
                   param6 = 1)

if(file.exists(kernels$file)){
  terra::plot(terra::rast(kernels$file))  
}
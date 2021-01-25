library("sp")
library("raster")
library("knitr")
library("elevatr")
library("httr")
library(wesanderson)

NOT_CRAN <- identical(tolower(Sys.getenv("NOT_CRAN")), "true")
knitr::opts_chunk$set(purl = NOT_CRAN, 
                      eval = NOT_CRAN,
                      fig.width = 7, 
                      fig.height = 7, 
                      tidy = TRUE,
                      dpi=600)
                      
distri <- read.delim("../Genome_Passport_firstVersion_Li.txt", header=T)
distri <- distri[, c(5,4,3)]
distri2 <- distri[, 1:2]


prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
distri_sp <- SpatialPoints(distri2, proj4string = CRS(prj_dd))
elevation_df <- get_elev_raster(distri,prj=prj_dd, z = 5)

ext <- extent(-120, -70, -15, 38)
s.crop <- crop(elevation_df, ext)

##for ggplot###
s.points <- rasterToPoints(s.crop)
s.df <- as.data.frame(s.points)

names(s.df) <- c("longitude", "latitude", "elevation")
names(distri)[3] <- "populations"
s.df$elevation <- ifelse(s.df$elevation<0, -100, s.df$elevation)

cols <- c("AN" = "purple", "GH" = wes_palette("Royal2")[3], "MH" = wes_palette("Darjeeling1")[1], "ML" = wes_palette("Royal1")[1], "SL"="black", "US"=wes_palette("Darjeeling1")[3])


distri$population <- ifelse(distri$populations=="Andes", "AN", distri$population)
distri$population <- ifelse(distri$populations=="Gua_Highland", "GH", distri$population)
distri$population <- ifelse(distri$populations=="Mex_Highland", "MH", distri$population)
distri$population <- ifelse(distri$populations=="Mex_Lowland", "ML", distri$population)
distri$population <- ifelse(distri$populations=="SA_Lowland", "SL", distri$population)
distri$population <- ifelse(distri$populations=="SW_US", "US", distri$population)


p1 <- ggplot(data=s.df, aes(y=latitude, x=longitude)) +
geom_raster(aes(fill=elevation)) +
geom_point(data=distri, aes(x=Long, y=Lat, color=population), size=2) +
coord_equal() +
scale_fill_gradient(limits=c(0, 5630), low="blue", high="brown", breaks=c(0, 2000, 4000)) +
theme(legend.position = c(0.18, 0.33)) +
ggtitle("")

p1 + scale_colour_manual(values = cols)

ggsave("Fig1elevation.pdf")
#####

SW_US <- distri[1:6, 1:2]
MexHigh <- distri[c(7,8,9,11,12,13), 1:2]
MexLow <- distri[c(10,14,15,19,20), 1:2]
Andes <- distri[c(21:23, 25, 31), 1:2]
SA_Low <- distri[c(24, 26:30), 1:2]
GuaHigh <- distri[16:18, 1:2]

plot(s.crop)
#plot(distri_sp, add = TRUE)

#wes_palette("Darjeeling1")[c(3,1)], wes_palette("Royal2")[3], wes_palette("FantasticFox1")[2]
pdf("distributionMapRaster.pdf", height=7, width=7)
plot(s.crop)
points(MexHigh$Long, MexHigh$Lat, col=wes_palette("Darjeeling1")[1], pch=20, cex=2)
points(MexLow$Long, MexLow$Lat, col=wes_palette("Royal1")[1], pch=20, cex=2)
points(GuaHigh$Long, GuaHigh$Lat, col=wes_palette("Royal2")[3], pch=20, cex=2)
points(SW_US$Long, SW_US$Lat, col=wes_palette("Darjeeling1")[3], pch=20, cex=2)
points(Andes$Long, Andes$Lat, col="purple", pch=20, cex=2)
points(SA_Low$Long, SA_Low$Lat, col="black", pch=20, cex=2)
legend(-112, 2, bg="white", col=c(wes_palette("Darjeeling1")[c(3,1)], wes_palette("Royal1")[1], wes_palette("Royal2")[3], "black", "purple"), pch=19, legend=c("US", "MH", "ML", "GH", "SL", "AN"), title="populations")
dev.off()

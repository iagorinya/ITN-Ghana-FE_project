##----------------------------------------------------------------------------------
## Script to obtain administrative boundaries for Ghana and save as csv file
##----------------------------------------------------------------------------------
library(tidyverse)
library(rdhs)
library(raster)
library(malariaAtlas)

## Custom functions
f_get_admin_df <- function(admin_sp) {
  admin_sp$area_sqkm <- area(admin_sp) / 1000000
  cents <- coordinates(admin_sp)
  cents_sp <- SpatialPointsDataFrame(coords = cents, data = admin_sp@data, proj4string = CRS(proj_str))
  cents_name <- over(cents_sp, admin_sp)
  admin_df <- as.data.frame(cents_sp)

  admin_df <- admin_df %>%
    dplyr::select(starts_with("NAME"), "GID_0", "coords.x1", "coords.x2", "area_sqkm") %>%
    dplyr::rename("cent_long" = coords.x1, "cent_lat" = coords.x2)

  return(admin_df)
}

proj_str <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
admin0_sp <- getData("GADM", country = "GHA", level = 0)
admin1_sp <- getData("GADM", country = "GHA", level = 1)
#admin2_sp <- getData("GADM", country = "GHA", level = 2)

admin1_sp.f <- as.MAPshp(admin1_sp)
#admin2_df <- f_get_admin_df(admin2_sp)
admin1_df <- f_get_admin_df(admin1_sp) %>% mutate(Region=NAME_1)

###TODO add ecocone column, and overlay province map with ecozone boundaries


## write csv file
fwrite(admin1_df, file.path('GHA_admin_dat.csv'))

### plot administrative units
pmap <- ggplot(warnings = FALSE) +
    geom_polygon(
      data = admin1_sp.f,
      aes(x = long, y = lat, group = group), col = "black", fill = NA
    )  #+theme_map()+  # requires ggthemes R package
  labs(title='Ghana province map')

ggsave(file.path( 'figures','GHA_region_map.png'), plot = pmap,  width = 8, height = 8, device = "png")

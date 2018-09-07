library(SDMTools)
library(raster)
library(rgdal) # loads sp package
library(rgeos)
library(maptools)
library(abind)
library(pixmap)
library(offsetsim)

bin_raster <- function(current_feature_raster, current_filename, shape_data_files, agg_factor){
  current_feature_raster[is.na(current_feature_raster)] = 0
  
  if (current_filename %in% shape_data_files){
    current_feature_raster = aggregate(current_feature_raster, fact = agg_factor, fun = modal)
  } else {
    current_feature_raster = aggregate(current_feature_raster, fact = agg_factor, fun = mean)
  }
  return(current_feature_raster)
}


convert_asc_to_raster <- function(output_data_folder, data_folder, asc_data_filenames, shape_data_files, agg_factor){
  
  if (!file.exists(output_data_folder)){
    dir.create(output_data_folder)
  }
  
  for (file_ind in seq_along(asc_data_filenames)){
    
    current_filename = asc_data_filenames[file_ind]
    current_feature_raster = raster(paste0(data_folder, current_filename))
    
    if (agg_factor > 1){
      current_feature_raster <- bin_raster(current_feature_raster, current_filename, shape_data_files, agg_factor)
    }
    
    current_filename = paste0(output_data_folder, gsub('.asc', '', current_filename), '.tif')
    writeRaster(current_feature_raster, current_filename, overwrite = TRUE)
    
    removeTmpFiles(h = 0)
    print(file_ind)
  }
  
}

simulation_inputs_folder = paste0(path.expand('~'), '/offset_data/uruguay/simulation_inputs/')
data_folder = '~/offset_data/uruguay/uruguay_data/'
output_data_folder = '~/offset_data/uruguay/uruguay_data/raster_tiff/species_features/'

asc_data_filenames <- list.files(path = paste0(data_folder, 'uruguay_raw_data/'), pattern = '.asc', all.files = FALSE, 
                                 full.names = FALSE, recursive = FALSE, ignore.case = FALSE, 
                                 include.dirs = FALSE, no.. = FALSE)
shape_data_files = vector()
data_characteristics <- as.list(read.csv(file=paste0(data_folder, 'group_defs.csv'), header=TRUE, sep=","))
names(data_characteristics) = c('group_index', 'group', 'filename')
group_characteristics = as.list(table(data_characteristics$group))
total_group_names = c(unique(as.vector(data_characteristics$group)), 'misc')

feature_type = 'species' # suitable_for_crops_&_forestry, 
group_names_to_use = feature_type #c('amphibians', 'birds', 'plants', 'mammals', 'fish', 'ecoregions', 'landscape_units', 'Spp_VU_CC')
convert_asc_layers = TRUE

build_site_characteristics = FALSE
overwrite_simulation_inputs = FALSE
save_current_layers = FALSE
agg_factor = 1


if (build_site_characteristics == TRUE){
  LGA_raster = load_rasters("~/offset_data/uruguay/uruguay_data/uruguay_raw_data/parcelas_uy.asc", features_to_use = 'all')
  LGA_array = raster_to_array(LGA_raster)
  site_characteristics <- build_site_characteristics(LGA_array)
  objects_to_save$site_characteristics <- site_characteristics
} else {
  site_characteristics = readRDS(paste0(simulation_inputs_folder, 'site_characteristics.rds'))
}


if (feature_type == 'ecosystem_services'){
  ecoservice_group_types = c('agua', 'amort', 'calidad', 'clima', 'enferm', 'genet') #(drinking water, , climactic regulation, , genetic resources)
  features_to_use = which(!is.na(match(data_characteristics$group, feature_type)))
  datalist_filenames = data_characteristics$filename[features_to_use]
  feature_rasters = load_rasters(paste0(data_folder, 'uruguay_raw_data/', datalist_filenames), 'all')
  feature_layers = lapply(seq_along(features_to_use), function(i) raster_to_array(subset(feature_rasters, i)))
  
  for (group_ind in seq_along(ecoservice_group_types)){
    current_reduced_group = Reduce('+', feature_layers[grep(pattern = ecoservice_group_types[group_ind], x = datalist_filenames)])
    current_feature_raster = raster(current_reduced_group)
    current_file_name = paste0(simulation_inputs_folder, 'ecosystem_service_feature_', ecoservice_group_types[group_ind], '.tif')
    writeRaster(current_feature_raster, current_file_name, overwrite = TRUE)
    print(paste(ecoservice_group_types[group_ind], 'done'))
  }
  
} else {
  convert_asc_to_raster(output_data_folder, paste0(data_folder, 'uruguay_raw_data/'), asc_data_filenames, shape_data_files, agg_factor)
}

euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}

outer(mat1,mat2, Vectorize(euclidean_distance))

# group_inds_to_use = match(group_names_to_use, names(group_characteristics))
# 
# current_group_characteristics = group_characteristics[group_inds_to_use]
# 
# data_length = sum(unlist(current_group_characteristics))
# data_list <- vector('list', length(data_length))
# data_list_names = vector()
# 
#   for (file_ind in seq_along(asc_data_filenames)){
#     
#     match_ind = which(data_characteristics$filename == asc_data_filenames[file_ind])
#     
#     if (length(match_ind) > 0){
#       current_group = data_characteristics$group[match_ind]
#     } else {
#       current_group = 'misc'
#     }
#     
#     if (current_group %in% group_names_to_use){
#       
#       #current_data = readRDS(paste0(rds_data_folder, current_file))
#       data_list_names = append(data_list_names, asc_data_filenames[file_ind])
#       #data_list[[data_ind]] = current_data
#       #data_ind = data_ind + 1
#     }
#     
#   }

names(data_list) = data_list_names



dev_weight_layer = readRDS(paste0(rds_data_folder, 'agut_clase1.rds'))
scale_fac = sum(as.vector(dev_weight_layer))

objects_to_save$dev_weights = lapply(seq_along(site_characteristics$land_parcels), 
                                     function(i) sum(dev_weight_layer[site_characteristics$land_parcels[[i]] ])/scale_fac)

offset_weight_layer = readRDS(paste0(rds_data_folder, 'agut_clase3.rds'))
scale_fac = sum(as.vector(offset_weight_layer))

objects_to_save$offset_weights = lapply(seq_along(site_characteristics$land_parcels), 
                                     function(i) sum(offset_weight_layer[site_characteristics$land_parcels[[i]] ])/scale_fac)

if (overwrite_simulation_inputs == TRUE){
  save_simulation_inputs(objects_to_save, paste0(simulation_inputs_folder))
  saveRDS(site_characteristics, paste0(simulation_inputs_folder, 'site_characteristics.rds'))
}

if (save_current_layers){
  data_dir <- paste0(simulation_inputs_folder, '/agg_factor_', agg_factor, '/', feature_type, '/')
  save_simulation_inputs(objects_to_save, data_dir)
}

print('all simulation objects saved')

# data_dir = paste0(simulation_inputs_folder, 'agg_factor_', agg_factor, '/', feature_type, '/')
# object_to_image = objects_to_save$landscape_ecology
# folder_to_output = paste0(simulation_inputs_folder, '/pdf_layers/')
# 
# if (write_layers_to_pdf){
#   graphics.off()
# #   output_pdf_filename = paste0(folder_to_output, 'agg_factor_', agg_factor, '_', feature_type, '.pdf')
# #   pdf(output_pdf_filename, width = 8.3, height = 11.7)
#   jpeg(paste0(folder_to_output, 'test.jpg'))
#   setup_sub_plots(nx = 3, ny = 3, x_space = 0.5, y_space = 0)
#   for (data_ind in seq_along(data_list)){
#     
#     image(data_list[[data_ind]], main = names(data_list)[data_ind], axes = FALSE)
#     
#     print(paste(data_ind, 'done'))
#   }
#   graphics.off()
# } 
# 


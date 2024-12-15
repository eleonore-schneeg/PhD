########## Adding Expanded OnPlaque objects to spe_filtered
## Expanded by 10
setwd("/home/ems2817/mtg_dna/deepcell/matched_objects_10")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe_filtered$OnPlaqueObject_10=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj_10=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj_10) == (spe_filtered$OnPlaqueObject_10=='OutsidePlaque'))


## Expanded by 20
setwd("/home/ems2817/mtg_dna/deepcell/matched_objects_20")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe_filtered$OnPlaqueObject_20=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj_20=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj_20) == (spe_filtered$OnPlaqueObject_20=='OutsidePlaque'))


## Expanded by 30
setwd("/home/ems2817/mtg_dna/deepcell/matched_objects_30")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe_filtered$OnPlaqueObject_30=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj_30=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj_30) == (spe_filtered$OnPlaqueObject_30=='OutsidePlaque'))

## Expanded by 40
setwd("/home/ems2817/mtg_dna/deepcell/matched_objects_40")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe_filtered$OnPlaqueObject_40=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj_40=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj_40) == (spe_filtered$OnPlaqueObject_40=='OutsidePlaque'))

## Expanded by 50
setwd("/home/ems2817/mtg_dna/deepcell/matched_objects_50")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe_filtered$OnPlaqueObject_50=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj_50=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj_50) == (spe_filtered$OnPlaqueObject_50=='OutsidePlaque'))

## Expanded by 100
setwd("/home/ems2817/mtg_dna/deepcell/matched_objects_100")
file_names=list.files( recursive = TRUE, pattern = "*.csv")

read_csv_column <- function(x){ data=data.frame(read_csv(x))
if(nrow(data)>0){
  data$name =sub(x = x, pattern = ".csv", replacement = "")
}
return(data)}
all_df <- do.call(rbind,lapply(file_names,read_csv_column))

all_df$DataFrame=paste0(all_df$name,'_',all_df$masks)

spe_filtered$OnPlaqueObject_100=ifelse(rownames(colData(spe_filtered)) %in% all_df$DataFrame,'OnPlaque','OutsidePlaque')
spe_filtered$PlaqueObj_100=all_df$masks_plaques[match(rownames(colData(spe_filtered)),all_df$DataFrame)]

#sanity check
all(is.na(spe_filtered$PlaqueObj_100) == (spe_filtered$OnPlaqueObject_100=='OutsidePlaque'))



# Create a vector of distances
distances <- c("1", "10", "20", "30", "40", "50", "100")

# Initialize the distance column with an empty string
spe_filtered$distance <- ""

# Iterate over the "OnPlaqueObject" columns and update the distance column
for (col in c("OnPlaqueObject", "OnPlaqueObject_10", "OnPlaqueObject_20", "OnPlaqueObject_30", "OnPlaqueObject_40", "OnPlaqueObject_50", "OnPlaqueObject_100")) {
  spe_filtered$distance <- ifelse(spe_filtered$distance == "" & spe_filtered[[col]] == "OnPlaque", distances[match(col, c("OnPlaqueObject", "OnPlaqueObject_10", "OnPlaqueObject_20", "OnPlaqueObject_30", "OnPlaqueObject_40", "OnPlaqueObject_50", "OnPlaqueObject_100"))], spe_filtered$distance)
}

# If no distance is assigned, set it to "NA"
spe_filtered$distance <- ifelse(spe_filtered$distance == "", "NA", spe_filtered$distance)

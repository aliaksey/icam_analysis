##load data and filter images that are out of focus and have two high min intensity
source("icam_load_data_filter_in_focus.R")
##filtering cells based on area aamd perimeter, compactness and solidity
source("icam_filter_cells.R")
##correct for intensities in different repeats
source("icam_correct_for_intensity_in_different_repeats.R")
##find hits based on treshhold intensity value
source("icam_analysis_frequences.R")
##find hits based on intensity value
source("icam_analysis_intensities.R")
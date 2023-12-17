# Header ----------------------------------------------------------------
# Project: PhanCurve
# File name: phanerozoic_diversity_curve.R
# Last updated: 2023-11-14
# Author: Lewis A. Jones
# Email: LewisA.Jones@outlook.com
# Repository: https://github.com/LewisAJones/PhanCurve

# Load libraries --------------------------------------------------------
library(palaeoverse)
library(ggplot2)
library(stringr)
library(deeptime)
library(RCurl)
library(tidyverse)

# Parameters ------------------------------------------------------------
# Prevent timeout for large API request
RCurl::curlSetOpt(timeout = 3000)
# PBDB download parameters
groups = c("Bivalvia,Brachiopoda,Cephalopoda,Gastropoda,Trilobita")
taxonomic_resolution = c("genus")
preservation = c("regular")
modifiers = c("genus_certain")
taxon_status = c("valid")
identification = c("latest")
environment = c("marine")
study_period = c("Fortunian,Holocene")
show = c("genus,pres,strat,coll,coords,paleoloc,loc,class")
# Bin size
bin_size = 1
# Download data ---------------------------------------------------------
# Skip if file already exists
if (any(list.files("./data") != "raw_pbdb_data_.RDS")) {
  # Setup API link
  API <- paste0("https://paleobiodb.org/data1.2/occs/list.csv?",
                "base_name=", groups,
                "&taxon_reso=", taxonomic_resolution, 
                "&ident=",identification,
                "&taxon_status=", taxon_status,
                "&idqual=", modifiers,
                "&pres=", preservation,
                "&interval=", study_period,
                "&envtype=", environment,
                "&show=", show)
  # Get data
  occdf <- RCurl::getURL(url = API, ssl.verifypeer = FALSE)
  # Read data
  occdf <- read.csv(textConnection(occdf))
  # Save raw data
  saveRDS(occdf, "./data/raw_pbdb_data_.RDS")
} else {
  occdf <- readRDS("./data/raw_pbdb_data_.RDS")
}

# Fix taxonomy ----------------------------------------------------------
# Filter unconstrained taxa
occdf <- occdf[which(occdf$class != "NO_CLASS_SPECIFIED"), ]
occdf <- occdf[which(occdf$class != "NO_CLASS_SPECIFIED"), ]
occdf <- occdf[which(occdf$class != "NO_FAMILY_SPECIFIED"), ]
occdf <- occdf[which(occdf$class != "NO_GENUS_SPECIFIED"), ]
# Use genus-level identifications
accepted_name <- str_split_fixed(string = occdf$accepted_name, 
                                 pattern = " ", 
                                 n = 2)
occdf$genus <- accepted_name[, 1]
occdf$rank <- "genus"
# Collapse Brachiopoda classes (for plotting later)
brachs <- c("Chileata", "Craniata", "Kutorginata", "Lingulata",
            "Obolellata", "Paterinata", "Rhynchonellata", "Strophomenata")
occdf$group <- occdf$class
occdf$group[which(occdf$group %in% brachs)] <- c("Brachiopoda")

# Calculate range through -----------------------------------------------
# Generate dataframe for storing counts
bins <- data.frame(bin = seq(from = 1, to = (541 / bin_size)),
                   max_ma = seq(from = 541, to = (0 + bin_size), by = -bin_size),
                   mid_ma = seq(from = (541 - (bin_size / 2)), to = (0 + (bin_size / 2)), by = -bin_size),
                   min_ma = seq(from = (541 - bin_size), to = 0, by = -bin_size),
                   total_counts = NA,
                   bivalvia_counts = NA,
                   brachiopoda_counts = NA,
                   cephalopoda_counts = NA,
                   gastropoda_counts = NA,
                   trilobita_counts = NA)
# Create unique taxon (to avoid genera with the same name but belonging to 
# different families)
occdf$family_genus <- paste0(occdf$family, "_", occdf$genus)
# Calculate ranges for all taxa
all <- tax_range_time(
  occdf = occdf, name = "family_genus",
  min_ma = "min_ma", max_ma = "max_ma",
  plot = TRUE)
# Calculate ranges for bivalvia
bivalvia <- tax_range_time(
  occdf = occdf[which(occdf$group == "Bivalvia"), ],
  name = "family_genus",
  min_ma = "min_ma", max_ma = "max_ma",
  plot = TRUE)
# Calculate ranges for brachiopoda
brachiopoda <- tax_range_time(
  occdf = occdf[which(occdf$group == "Brachiopoda"), ],
  name = "family_genus",
  min_ma = "min_ma", max_ma = "max_ma",
  plot = TRUE)
# Calculate ranges for cephalopoda
cephalopoda <- tax_range_time(
  occdf = occdf[which(occdf$group == "Cephalopoda"), ],
  name = "family_genus",
  min_ma = "min_ma", max_ma = "max_ma",
  plot = TRUE)
# Calculate ranges for gastropoda
gastropoda <- tax_range_time(
  occdf = occdf[which(occdf$group == "Gastropoda"), ],
  name = "family_genus",
  min_ma = "min_ma", max_ma = "max_ma",
  plot = TRUE)
# Calculate ranges for trilobita
trilobita <- tax_range_time(
  occdf = occdf[which(occdf$group == "Trilobita"), ],
  name = "family_genus",
  min_ma = "min_ma", max_ma = "max_ma",
  plot = TRUE)

# Count taxa
for (i in 1:nrow(bins)) {
  # All taxa
  bins$total_counts[i] <- length(
    which(all$min_ma <= bins$max_ma[i] & all$max_ma >= bins$min_ma[i]))
  # Bivalvia
  bins$bivalvia_counts[i] <- length(
    which(bivalvia$min_ma <= bins$max_ma[i] & bivalvia$max_ma >= bins$min_ma[i]))  
  # Brachiopoda
  bins$brachiopoda_counts[i] <- length(
    which(brachiopoda$min_ma <= bins$max_ma[i] & brachiopoda$max_ma >= bins$min_ma[i]))  
  # Cephalopoda
  bins$cephalopoda_counts[i] <- length(
    which(cephalopoda$min_ma <= bins$max_ma[i] & cephalopoda$max_ma >= bins$min_ma[i]))  
  # Gastropoda
  bins$gastropoda_counts[i] <- length(
    which(gastropoda$min_ma <= bins$max_ma[i] & gastropoda$max_ma >= bins$min_ma[i]))  
  # Trilobita
  bins$trilobita_counts[i] <- length(
    which(trilobita$min_ma <= bins$max_ma[i] & trilobita$max_ma >= bins$min_ma[i]))
}

bins %>% 
  gather(key = "group", value = "counts", 6:10) -> counts

ggplot(data = counts, aes(x = mid_ma, y = counts)) +
  geom_area(aes(fill = group)) +
  geom_vline(aes(xintercept = 66), alpha = 1) + # Cretaceous
  geom_vline(aes(xintercept = 201.400), alpha = 1) + # end-Triassic
  geom_vline(aes(xintercept = 251.902), alpha = 1) + # end-Permian
  geom_vline(aes(xintercept = 358.900), alpha = 1) + # Late Devonian
  geom_vline(aes(xintercept = 443.800), alpha = 1) + # end-Ordovician
  scale_x_reverse() +
  scale_fill_viridis_d(labels = c("bivalvia_counts" = "Bivalvia", 
                                 "brachiopoda_counts" = "Brachiopoda", 
                                 "cephalopoda_counts" = "Cephalopoda", 
                                 "gastropoda_counts" = "Gastropoda", 
                                 "trilobita_counts" = "Trilobita")) +
  labs(title = "Phanerozoic diversity curve",
       x = "Time (Ma)", y = "Range-through genus richness") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  coord_geo(height = unit(1, "line"))

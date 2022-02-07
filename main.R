library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# library(purrr)

#Sriramteja Veerisetti
# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.

library(tidyverse)
project_metadata <- read.csv("data/proj_metadata.csv")
project_metadata[1:6, 1:5]

#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data.csv')
#' 

read_expression_table <- function(filename){

expression_mat <- readr::read_delim("data/example_intensity_data.csv") %>%
  tidyr::pivot_longer(starts_with("GSM"), names_to = "subject_id", values_to = "value") %>%
  tidyr::pivot_wider(names_from = "probe", values_from = "value") %>%
 
  return ()

}

#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"

period_to_underscore <- function(str) {
 str_replace_all(str, pattern = "\\.", replacement = "_") %>%
  return()
  
}

# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  
  walnut <- readr::read_csv("data/proj_metadata.csv")
  walnut <- rename(walnut, "age" = Age.at.diagnosis)
  walnut <- rename(walnut, "subtype" = SixSubtypesClassification)
  walnut <- rename(walnut, "Batch" = normalizationcombatbatch)
  walnut <- dplyr::select(walnut, Sex, age, TNM.Stage, Tumor.Location, geo_accession, subtype, Batch) %>%

    return ()
}


#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)

stage_as_factor <- function(data) {
  
  data <- readr::read_csv("data/proj_metadata.csv")
  data %>%
  select(everything()) %>%
  mutate (Stage = "Stage") %>%
  mutate(New_Stage = paste(Stage, TNM.Stage)) %>%
  select(New_Stage) %>%

   return ()
}

#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")

mean_age_by_sex <- function(data, sex) {
  
  proj_metadata <- readr::read_csv("data/proj_metadata.csv")
  dplyr::group_by(proj_metadata, 
      Sex
      )%>% dplyr::summarize(mean_age_by_sex = mean(age.at.diagnosis))%>%
  
    return ()
}


#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage.
#' @export
#'
#' @examples age_by_stage(data)

age_by_stage <- function(data) {
  
  data <- readr::read_csv("data/proj_metadata.csv") %>%
    
  temp_data <-
    data %>%
    select(everything()) %>%
    mutate (Stage = "Stage") %>%
    mutate(New_Stage = paste(Stage, TNM.Stage)) %>%
    select(New_Stage, everything()) 
  
  dplyr::group_by(temp_data, 
      New_Stage
        )%>% dplyr::summarize(mean_by_stage = mean(age.at.diagnosis))
  
  return ()
}

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
  
  data <- readr::read_csv("data/proj_metadata.csv")
  
  candy_data <-
    data %>%
    select(everything()) %>%
    mutate (Stage = "Stage") %>%
    mutate(New_Stage = paste(Stage, TNM.Stage)) %>%
    select(New_Stage, everything()) 
    
  dplyr::group_by(candy_data, New_Stage, SixSubtypesClassification)%>%
  dplyr::summarize(C3 = n()) %>%
    pivot_wider(names_from = SixSubtypesClassification, values_from = C3) %>%

  return ()
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.

summarize_expression <- function(exprs) {
  
  
  nikhil <- readr::read_delim("data/example_intensity_data.csv")
  
  install.packages("matrixStats")
  library(matrixStats)
  
  combination <- readr::read_delim("data/example_intensity_data.csv") 
    Df <- as.data.frame(combination)
    Df$row_var = rowVars(as.matrix(Df[,c(2:36)]))
    Df$row_means = rowMeans(as.matrix(Df[,c(2:36)]))
   
    Df2 <- Df[, c("probe", "row_var", "row_means")]
    
    head(Df2) %>%
    
  
  return ()
}

## code to prepare `DATASET` dataset goes here
contact_all <- list()
contact_home <- list()
for(country_id in 1:4){
  if(country_id == 1){
    age_dist_country <- c(0.481, 0.203, 0.316) # Uganda
  }else if(country_id == 2){
    age_dist_country <- c(0.292, 0.193, 0.515) # South africa
  }else if(country_id == 3){
    age_dist_country <- c(0.419, 0.195, 0.386) # kenya
  }else if(country_id == 4){
    age_dist_country <- c(0.440, 0.190, 0.370) # nigeria
  }


  if(country_id == 1){
    # note: Uganda data is compute from survey data, which require a separate prcessing process; other contries were from the Prem et al PLOS Comp. Bio 2017 paper
    age_mix <- read.csv("data-raw/Uganda_age_mixture_population_level.csv")
    contact_all[[country_id]] <- age_mix$total_age
    contact_home[[country_id]] <- age_mix$home_age
  }else{
    if(country_id == 2){
      AGE_matrix <- read.csv("data-raw/south_africa_age.csv", header = F)
      HOME_age_matrix <- read.csv("data-raw/home_south_africa_age.csv", header = F)
    }else if(country_id == 3){
      AGE_matrix <- read.csv("data-raw/kenya_age.csv", header = F) # kenya
      HOME_age_matrix <- read.csv("data-raw/home_kenya_age.csv", header = F) # kenya
    }
    else if(country_id == 4){
      AGE_matrix <- read.csv("data-raw/Nigeria_age.csv", header = F) # nigeria
      HOME_age_matrix <- read.csv("data-raw/home_Nigeria_age.csv", header = F) # Nigeria
    }
    # the age mixing matrix need to collpase into 3*3 since we only have age distribution for 3 groups.
    ##  ##  ##  ##
    ## this step is important: the contact matrix is not symetrical: row and columns represent different things
    ##  ##  ##  ##
    # the columns represent contact, we take colSums
    # the rows are the participant, we take the average number of contact for the participant: rowMeans
    age <- cbind(rowMeans(AGE_matrix[,1:3]),rowMeans(AGE_matrix[,4:5]),rowMeans(AGE_matrix[,6:16]))
    AGE_matrix <- rbind(colSums(age[1:3,]),colSums(age[4:5,]),colSums(age[6:16,]))
    # weight each column by the age distribution: we have to change the matrix to reflect the age distribution in the population
    AGE_matrix <- (as.matrix(rep(1,length(age_dist_country))) %*% t(age_dist_country)) * AGE_matrix

    # matrix in the end should be the number of contact for each age group pair
    AGE_matrix <- (AGE_matrix + t(AGE_matrix))/2
    contact_all[[country_id]] <- as.matrix(AGE_matrix)[which(upper.tri(AGE_matrix,diag = T))]

    age <- cbind(rowMeans(HOME_age_matrix[,1:3]),rowMeans(HOME_age_matrix[,4:5]),rowMeans(HOME_age_matrix[,6:16]))
    HOME_age_matrix <- rbind(colSums(age[1:3,]),colSums(age[4:5,]),colSums(age[6:16,]))
    HOME_age_matrix <- (as.matrix(rep(1,length(age_dist_country))) %*% t(age_dist_country)) * HOME_age_matrix
    # matrix in the end should be a triangle one
    HOME_age_matrix <- (HOME_age_matrix + t(HOME_age_matrix))/2
    contact_home[[country_id]] <- as.matrix(HOME_age_matrix)[which(upper.tri(HOME_age_matrix,diag = T))]
  }
}

usethis::use_data(contact_all, contact_home, overwrite = TRUE)

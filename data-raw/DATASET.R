## code to prepare `DATASET` dataset goes here
SimulatedPrototypes <- read.csv(file='vignettes/Simulated Portotypes.csv', row.names = 1)
colnames(SimulatedPrototypes) <- paste0('Pty', 1:4)

usethis::use_data(SimulatedPrototypes, overwrite = TRUE)

SimulatedMemberships <- read.csv(file='vignettes/Simulated Memberships.csv', row.names = 1)
rownames(SimulatedMemberships) <- paste0('Pty', 1:4)
colnames(SimulatedMemberships) <- 1:ncol(SimulatedMemberships)

usethis::use_data(SimulatedMemberships, overwrite = TRUE)


SimulatedDataset <- read.csv(file='vignettes/Simulated Dataset.csv', row.names = 1)
colnames(SimulatedDataset) <- paste0('Var', 1:ncol(SimulatedDataset))

usethis::use_data(SimulatedDataset, overwrite = TRUE)

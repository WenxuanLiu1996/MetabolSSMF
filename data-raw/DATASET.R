## code to prepare `DATASET` dataset goes here
SimulatedPrototypes <- read.csv(file='vignettes/Simulated Prototypes.csv', row.names = 1)
rownames(SimulatedPrototypes) <- paste0('Pty', 1:4)
colnames(SimulatedPrototypes) <- 1:ncol(SimulatedPrototypes)

usethis::use_data(SimulatedPrototypes, overwrite = TRUE)

SimulatedMemberships <- read.csv(file='vignettes/Simulated Memberships.csv', row.names = 1)
colnames(SimulatedMemberships) <- paste0('Pty', 1:4)
rownames(SimulatedMemberships) <- 1:nrow(SimulatedMemberships)

usethis::use_data(SimulatedMemberships, overwrite = TRUE)


SimulatedDataset <- read.csv(file='vignettes/Simulated Dataset.csv', row.names = 1)
colnames(SimulatedDataset) <- paste0('Var', 1:ncol(SimulatedDataset))

usethis::use_data(SimulatedDataset, overwrite = TRUE)

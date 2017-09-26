# ==============================================================================
# Code to generate matrices for benchmark testing of the econullnetr package. 
#  The accompanies the supplementary material of:
#  Vaughan, I.P., Gotelli, N.J., Memmott, J., Pearson, C.E., Woodward, G. &
#  Symondson, W.O.C. (2017). econullnetr: an R package using null models to 
#  analyse the structure of ecological networks and identify resource 
#  selection. Methods in Ecology and Evolution, in press.
# --------------------------------------------------------------------------
# 100 bipartite communities are created with 1-20 consumer species and 
#  5-50 resource species. Two data sets are created for each community: 
#    i) the resources consumed by each individual consumer (each species 
#       comprises 100 individuals) 
#    ii) the abundance of the different resources. 
# --------------------------------------------------------------------------
# Different versions of the consumer data are generated to reflect different 
#  types of resource preferences. Each one comprises a list of the 100 matrices:
#   i) perfect.generalists. No preferences: resource selection is simply a 
#              function of the relative abundance of the different resources.
#   ii) perfect.specialists.  All individuals of a consumer only consume a single,
#              preferred resource
#   iii) mix.1. Every consumer species shows a different degree of specialisation,
#              ranging between perfect specialists and perfect generalists .
#   iv) add.1. Add in a single specialist, selected at random, to a matrix
#              of otherwise perfect generalists. Three versions: perfect, 
#              'strong' and 'moderate' specialist (add.1.spec, add.1.str.spec
#              and add.1.mod.spec).
#   v)  Add.2. Add two specialists, selected at random, to a matrix of
#              otherwise perfect generalists. Three versions (add.2.spec, 
#              add.2.spec and add.2.mod.spec.
# --------------------------------------------------------------------------
# Three types of data that could be collected in the field are simulated:
#  i) Each individual consumer interacts with one resource
#  ii) Individual consumers can interact with >1 resource, but the number of times
#      each interaction occurs is not recorded (i.e. just a list of resource
#      species with which an interaction occurred.
#  iii) Individual consumers can interact with >1 resource and the number of times
#      each interaction occurs is recorded (e.g. number of individual prey eaten
#      by an individual predator = 'count data').
# Type ii is derived from count data (type iii) by converting counts >1 to 1.
# ==============================================================================



# ==============================================================================
# Stage 1: create a series of objects used in generating the benchmark matrices:
#          Matrix.specifications: dimensions of the benchmarking matrices
#          sim.res.abund: list of resource abundance matrices 
#          pref.res: list of resource preferences for each consumer species
#          indiv.n.res: number of resources consumed by each individual consumer
#          res.pref: list of matrices containing resource preferences of every 
#                    consumer species in the 100 matrices
# ==============================================================================

# --------------------------------------------------------------------------
# Specifications for the 100 matrices:
N.matrices <- 100 
N.ind.cons <- 100 # Number of individuals per consumer species
Matrix.specifications <- matrix(ncol = 4, nrow = N.matrices)
colnames(Matrix.specifications) <- c("Resource.sp", "Consumer.sp", 
                                     "Total.consumers", "Mean.obs.per.cell")

for (i in 1:N.matrices){
  Matrix.specifications[i, 1] <- round(runif(1, 5, 50), digits = 0)
  Matrix.specifications[i, 2] <- round(runif(1, 1, 20), digits = 0)
  Matrix.specifications[i, 3] <- Matrix.specifications[i, 2] * N.ind.cons
  Matrix.specifications[i, 4] <- Matrix.specifications[i, 3] /
    (Matrix.specifications[i, 1] * 
       Matrix.specifications[i, 2]) 
}
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Generate resource abundance matrices (sim.res.abund). Fill 100 matrices with 
#   abundance values drawn from the lognormal distribution, with the standard 
#   deviation (= shape parameter) set to 1.5 (= median from NCEAS pollination 
#   webs: Dormann et al (2009) The Open Ecology Journal, 2, 7-24. 
sim.res.abund <- list()
for (i in 1:N.matrices){
  Resource.sp <- rep("NA", Matrix.specifications[i, 1]) 
  for (j in 1:Matrix.specifications[i, 1]) {
    Resource.sp[j] <- paste('Res.', j, sep = "")
  }
  Resources <- cbind.data.frame("None", t(rlnorm(Matrix.specifications[i, 1], 
                                                 meanlog = 0, sdlog = 1.5)))
  colnames(Resources) <- c("Site", Resource.sp)
  sim.res.abund[[i]] <- Resources
}
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# For count data, set the number of resources consumed by each individual 
#   consumer:
#   1. Determine the mean number of individual resources for each consumer 
#      species. Draw from a uniform random distribution, with values in the 
#      range 1-3.
#   2. Set the number of resources consumed by each individual consumer, by 
#      drawing from a truncated Poisson distribution (aster package), using the 
#      mean value for the species set in stage 1.
# --------------------------------------------------------------------------
library(aster)
for (i in 1:N.matrices) {
  inter.summ <- data.frame(Matrix = paste("Mat.", i, sep = ""),
                           Consumer.sp = paste("C", seq(1, 
                                                        Matrix.specifications[i, 2]), sep = ""),
                           Mean = runif(Matrix.specifications[i, 2], min = 1,
                                        max = 3))
  ifelse(i == 1, Sp.N.res <- inter.summ, 
         Sp.N.res <- rbind(Sp.N.res, inter.summ))
}

for (i in 1:N.matrices) {
  for (j in 1:Matrix.specifications[i, 2]) {
    ind.choice <- data.frame(Matrix = paste("Mat.", i, sep = ""),
                             Consumer.sp = paste("C", j, sep = ""),
                             N.res = rep(0, 100))
    mean.val <- Sp.N.res[(Sp.N.res$Matrix == as.character(ind.choice$Matrix[1]) &
                            Sp.N.res$Consumer.sp == 
                            as.character(ind.choice$Consumer.sp[1])), 3]
    ind.choice$N.res <- rktp(100, 0, mean.val)
    ifelse(j == 1, indiv.res <- ind.choice, 
           indiv.res <- rbind(indiv.res, ind.choice))
  }
  ifelse(i == 1, indiv.n.res <- indiv.res, 
         indiv.n.res <- rbind(indiv.n.res, indiv.res))
}
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Create resource preferences amongst the consumer species. Resources are
#   randomly re-ordered to produce a ranking independent of resource abundance.
res.pref <- list()
for (i in 1:length(sim.res.abund)) {
  for (j in 1:Matrix.specifications[i, 2]) {
    aa <- matrix(runif(Matrix.specifications[i, 1] * 
                         Matrix.specifications[i, 2]), 
                 nrow = Matrix.specifications[i, 1], 
                 ncol = Matrix.specifications[i, 2])
    bb <- apply(aa, 2, rank)
    mat.1 <- data.frame(Resource = paste("Res.", seq(1, 
                                                     length(sim.res.abund[[i]]) - 1), sep = ""), bb)
    colnames(mat.1)[2:ncol(mat.1)] <- paste("C", seq(1, ncol(bb)), sep = "")
  }
  res.pref[[i]] <- mat.1
}
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Set the degree of specialisation of each taxon in communities with variable
#  degree of specialisation (= 'Mix.1'). The relative preference for different
#  resources is set using a beta-binomial distribution from the VGAM package. 
#  Varying the first shape parameter between zero and one represents 
#  patterns of relative preference ranging from perfect specialists (only one 
#  resource selected) to perfect generalists. The value for each taxon
#  is drawn from a uniform random distribution.
for (i in 1:N.matrices) {
  pref.shapes <- data.frame(Mat.num = rep(i, Matrix.specifications[i, 2]),
                            Cons.sp = seq(1, Matrix.specifications[i, 2]),
                            Shape = runif(Matrix.specifications[i, 2]))
  ifelse(i == 1, Mix.1.shapes <- pref.shapes, 
         Mix.1.shapes <- rbind(Mix.1.shapes, pref.shapes))
}
# --------------------------------------------------------------------------


# --------------------------------------------------------------------------
# Set the 1-2 species that will be made into specialists and seeded into
#  matrices of perfect generalists. 
for (i in 1:N.matrices) {
  if(Matrix.specifications[i, 2] > 1) {
    specs <- sample(paste("C", seq(1, Matrix.specifications[i, 2]), sep = ""), 
                    2, replace = FALSE)
  } else {
    specs <- c(sample(paste("C", seq(1, Matrix.specifications[i, 2]), sep = ""), 
                      1), NA)
  }
  specs <- data.frame(Matrix = paste("Mat.", i, sep = ""), Choice1 = specs[1],
                      Choice2 = specs[2])
  ifelse(i == 1, Specialist.spp <- specs, 
         Specialist.spp <- rbind(Specialist.spp, specs))
}
# --------------------------------------------------------------------------




# ==============================================================================
# Stage 2: generate benchmark matrices for data type i i.e. one resource
#          per individual consumer
# --------------------------------------------------------------------------
# Create nine lists of 100 matrices:
#   1. perfect.generalist: no resource choice (probability of selecting a 
#            resource is proportional to its abundance)
#   2. perfect.specialist: all consumers are perfect specialists (individuals 
#            of each consumer sp only interact with one preferred resource sp
#   3. mix.1 = each consumer species differs in its relative generalism v. 
#            specialism
#   4-6. add.1.spec/add.1.str.spec/add.1.mod.spec: make one consumer species a
#            specialist in each community. Either a perfect specialist 
#            (add.1.spec), a 'strong' specialist (beta-binomial shape coef = 0.2, 
#            add.1.str.spec) or a 'moderate' specialist (shape coef = 0.5, 
#            add.1.mod.spec).
#   7-9. add.2.spec/add.2.str.spec/add.2.mod.spec: as 4-6, but adding two
#            randomly-selected specialist species to each matrix.
# ==============================================================================


# --------------------------------------------------------------------------
library(VGAM)

# Create the nine lists
perfect.generalist <- list()
perfect.specialist <- list()
mix.1 <- list() 
add.1.spec <- list()
add.1.str.spec <- list()
add.1.mod.spec <- list()
add.2.spec <- list()
add.2.str.spec <- list()
add.2.mod.spec <- list()

for (i in 1:N.matrices) {
  
  # -------------------------------
  # Generate empty consumer matrices
  for (j in 1:Matrix.specifications[i, 2]){ 
    gen <- data.frame(Consumer.sp = rep(paste('C', j, sep = ""), N.ind.cons), 
                      Cons.id = seq(1, N.ind.cons))
    ifelse(j == 1, gen.mat <- gen, gen.mat <- rbind(gen.mat, gen))
  }
  gen.mat <- cbind.data.frame(gen.mat, matrix(0, ncol = 
                                                Matrix.specifications[i, 1]))
  colnames(gen.mat)[-c(1, 2)] <- paste("Res.", colnames(gen.mat)[-c(1, 2)], 
                                       sep = "")
  spec.mat <- gen.mat
  mix.1.mat <- gen.mat
  seed.1.mat <- gen.mat
  seed.1.ss <- gen.mat
  seed.1.mod <- gen.mat
  seed.2.mat <- gen.mat
  seed.2.ss <- gen.mat
  seed.2.mod <- gen.mat
  # -------------------------------
  
  # -------------------------------
  # Fill the perfect generalist matrix. 
  for (j in 1:Matrix.specifications[i, 2]) {
    for (k in 1:N.ind.cons) {               
      res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                    prob = as.numeric(sim.res.abund[[i]][-1]), 
                    replace = FALSE)
      gen.mat[((j - 1) * N.ind.cons) + k, res] <- 1
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill the perfect specialist matrix. 
  for (j in 1:Matrix.specifications[i, 2]) {  
    pref <- as.character(res.pref[[i]]$Resource[res.pref[[i]][, j + 1] == 1])
    spec.mat[spec.mat$Consumer.sp == colnames(res.pref[[i]][j + 1]), pref] <- 1
  }
  # -------------------------------
  
  # -------------------------------
  # Fill mix.1 matrix. 
  for (j in 1:Matrix.specifications[i, 2]) {  
    # Create preference profile from beta-binomial distribution
    shape.param <- Mix.1.shapes[Mix.1.shapes$Mat.num == i & 
                                  Mix.1.shapes$Cons.sp == j, 3]
    pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                               shape1 = shape.param, shape2 = 1) / 
      max(dbetabinom.ab(0, size = 50, shape1 = 
                          shape.param, shape2 = 1))
    res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                             Ord = seq(1, nrow(res.pref[[i]])),
                             Cons.pref = res.pref[[i]][, j + 1])
    res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
    res.choice$Cons.pref <- pref.prof
    res.choice <- res.choice[with(res.choice, order(Ord)), ]
    weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
    ifelse(j == 1, wr <- weighted.res, wr <- rbind(wr, weighted.res))
    for (k in 1:N.ind.cons) {
      res <- sample(colnames(weighted.res), 1, prob = 
                      as.numeric(weighted.res), replace = FALSE)
      mix.1.mat[((j - 1) * N.ind.cons) + k, res] <- 1
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Select the 1 or 2 consumer species that will be made into specialists and 
  #   seeded into matrices that contain perfect generalists
  sel.1 <- as.numeric(gsub("C", "", Specialist.spp[i, 2]))
  ifelse(is.na(Specialist.spp[i, 3]), sel.2 <- NA, 
         {sel.2 <- as.numeric(gsub("C", "", Specialist.spp[i, 3]))})
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 1 perfect specialist seeded in. 
  for (j in 1:Matrix.specifications[i, 2]) {  
    for (k in 1:N.ind.cons) {
      ifelse(j == sel.1, {
        pref <- as.character(res.pref[[i]]$Resource[res.pref[[i]][, j + 1] == 1])
        seed.1.mat[seed.1.mat$Consumer.sp == colnames(res.pref[[i]][j + 1]), 
                   pref] <-1
      }, {
        res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                      prob = as.numeric(sim.res.abund[[i]][-1]), 
                      replace = FALSE)
        seed.1.mat[((j - 1) * N.ind.cons) + k, res] <- 1
      })
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 2 perfect specialists seeded in. 
  # Select the specialist.
  if(Matrix.specifications[i, 2] >= 2) { 
    for (j in 1:Matrix.specifications[i, 2]) {  
      for (k in 1:N.ind.cons) {
        ifelse(j == sel.1 | j == sel.2, {
          pref <- as.character(res.pref[[i]]$Resource[res.pref[[i]][, j + 1] == 1])
          seed.2.mat[seed.2.mat$Consumer.sp == colnames(res.pref[[i]][j + 1]), 
                     pref] <-1
        }, {
          res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                        prob = as.numeric(sim.res.abund[[i]][-1]))
          seed.2.mat[((j - 1) * N.ind.cons) + k, res] <- 1
        })
      }
    }
  } else {seed.2.mat[, 3:ncol(seed.2.mat)] <- NA}
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 1 moderate specialist seeded in. 
  for (j in 1:Matrix.specifications[i, 2]) {
    if(j == sel.1) {
      pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                 shape1 = 0.5, shape2 = 1) / 
        max(dbetabinom.ab(0, size = 50, shape1 = 
                            0.5, shape2 = 1))
      res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                               Ord = seq(1, nrow(res.pref[[i]])),
                               Cons.pref = res.pref[[i]][, sel.1 + 1])
      res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
      res.choice$Cons.pref <- pref.prof
      res.choice <- res.choice[with(res.choice, order(Ord)), ]
      weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
    }
    
    for (k in 1:N.ind.cons) {
      ifelse(j == sel.1, {
        pref <- as.character(sample(colnames(weighted.res), 1, 
                                    prob = as.numeric(weighted.res)))
        seed.1.mod[((j - 1) * N.ind.cons) + k, pref] <- 1
      }, {
        res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                      prob = as.numeric(sim.res.abund[[i]][-1]), 
                      replace = FALSE)
        seed.1.mod[((j - 1) * N.ind.cons) + k, res] <- 1
      })
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 2 moderate specialists seeded in. 
  # Select the specialist.
  if(Matrix.specifications[i, 2] >= 2) {       
    for (j in 1:Matrix.specifications[i, 2]) {  
      if(j == sel.1 | j == sel.2) {
        pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                   shape1 = 0.5, shape2 = 1) / 
          max(dbetabinom.ab(0, size = 50, shape1 = 
                              0.5, shape2 = 1))
        res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                                 Ord = seq(1, nrow(res.pref[[i]])),
                                 Cons.pref = res.pref[[i]][, j + 1])
        
        res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
        res.choice$Cons.pref <- pref.prof 
        res.choice <- res.choice[with(res.choice, order(Ord)), ]
        weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
      }
      
      for (k in 1:N.ind.cons) {
        ifelse(j == sel.1 | j == sel.2, {
          pref <- as.character(sample(colnames(weighted.res), 1, 
                                      prob = as.numeric(weighted.res), replace = FALSE))
          seed.2.mod[((j - 1) * N.ind.cons) + k, pref] <- 1
        }, {
          res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                        prob = as.numeric(sim.res.abund[[i]][-1]), 
                        replace = FALSE)
          seed.2.mod[((j - 1) * N.ind.cons) + k, res] <- 1
        })
      }
    }
  } else {seed.2.mod[, 3:ncol(seed.2.mod)] <- NA}
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 1 strong specialist seeded in. 
  for (j in 1:Matrix.specifications[i, 2]) {  
    if(j == sel.1) {
      pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                 shape1 = 0.2, shape2 = 1) / 
        max(dbetabinom.ab(0, size = 50, shape1 = 
                            0.2, shape2 = 1))
      res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                               Ord = seq(1, nrow(res.pref[[i]])),
                               Cons.pref = res.pref[[i]][, sel.1 + 1])
      res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
      res.choice$Cons.pref <- pref.prof
      res.choice <- res.choice[with(res.choice, order(Ord)), ]
      weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
    }
    
    for (k in 1:N.ind.cons) {
      ifelse(j == sel.1, {
        pref <- as.character(sample(colnames(weighted.res), 1, 
                                    prob = as.numeric(weighted.res)))
        seed.1.ss[((j - 1) * N.ind.cons) + k, pref] <- 1
      }, {
        res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                      prob = as.numeric(sim.res.abund[[i]][-1]), 
                      replace = FALSE)
        seed.1.ss[((j - 1) * N.ind.cons) + k, res] <- 1
      })
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 2 strong specialists seeded in. 
  # Select the specialist.
  if(Matrix.specifications[i, 2] >= 2) {       
    for (j in 1:Matrix.specifications[i, 2]) {  
      if(j == sel.1 | j == sel.2) {
        pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                   shape1 = 0.2, shape2 = 1) / 
          max(dbetabinom.ab(0, size = 50, shape1 = 
                              0.2, shape2 = 1))
        res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                                 Ord = seq(1, nrow(res.pref[[i]])),
                                 Cons.pref = res.pref[[i]][, j + 1])
        
        res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
        res.choice$Cons.pref <- pref.prof 
        res.choice <- res.choice[with(res.choice, order(Ord)), ]
        weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
      }
      
      for (k in 1:N.ind.cons) {
        ifelse(j == sel.1 | j == sel.2, {
          pref <- as.character(sample(colnames(weighted.res), 1, 
                                      prob = as.numeric(weighted.res), replace = FALSE))
          seed.2.ss[((j - 1) * N.ind.cons) + k, pref] <- 1
        }, {
          res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                        prob = as.numeric(sim.res.abund[[i]][-1]), 
                        replace = FALSE)
          seed.2.ss[((j - 1) * N.ind.cons) + k, res] <- 1
        })
      }
    }
  } else {seed.2.ss[, 3:ncol(seed.2.mod)] <- NA}
  # -------------------------------
  
  # -------------------------------
  perfect.generalist[[i]] <- gen.mat
  perfect.specialist[[i]] <- spec.mat
  mix.1[[i]] <- mix.1.mat
  add.1.spec[[i]] <- seed.1.mat
  add.1.str.spec[[i]] <- seed.1.ss
  add.1.mod.spec[[i]] <- seed.1.mod
  add.2.spec[[i]] <- seed.2.mat
  add.2.str.spec[[i]] <- seed.2.ss
  add.2.mod.spec[[i]] <- seed.2.mod
  
  ## Progress count
  Sys.sleep(0.02); print(i); flush.console() 
}
# --------------------------------------------------------------------------





# ==============================================================================
# Stage 3: generate benchmark matrices for data type iii i.e. counts of the 
#          number of resources consumed per individual consumer
# --------------------------------------------------------------------------
# Equivalent code to 'Stage 2' above, producing lists of matrices with the 
#  same names
# ==============================================================================


# Create list objects to hold the matrices
perfect.generalist <- list()
perfect.specialist <- list()
mix.1 <- list() 
add.1.spec <- list()
add.1.str.spec <- list()
add.1.mod.spec <- list()
add.2.spec <- list()
add.2.str.spec <- list()
add.2.mod.spec <- list()


for (i in 1:N.matrices) {
  # -------------------------------
  for (j in 1:Matrix.specifications[i, 2]){ 
    gen <- data.frame(Consumer.sp = rep(paste('C', j, sep = ""), N.ind.cons), 
                      Cons.id = seq(1, N.ind.cons))
    ifelse(j == 1, gen.mat <- gen, gen.mat <- rbind(gen.mat, gen))
  }
  gen.mat <- cbind.data.frame(gen.mat, matrix(0, ncol = 
                                                Matrix.specifications[i, 1]))
  colnames(gen.mat)[-c(1, 2)] <- paste("Res.", colnames(gen.mat)[-c(1, 2)], 
                                       sep = "")
  spec.mat <- gen.mat
  mix.1.mat <- gen.mat
  seed.1.mat <- gen.mat
  seed.2.mat <- gen.mat
  seed.1.mod <- gen.mat
  seed.2.mod <- gen.mat
  seed.1.ss <- gen.mat
  seed.2.ss <- gen.mat
  # -------------------------------
  
  # -------------------------------
  # Set the number of resource items (not necessarily resource spp, as species 
  #   can be selected >1) for each individual consumer
  n.res <- eval(parse(text = paste('indiv.n.res[indiv.n.res$Matrix == "Mat.', 
                                   i, '", ]', sep = '')))
  # -------------------------------
  
  # -------------------------------
  # Fill the generalist matrix. 
  for (k in 1:nrow(n.res)) {
    for (l in 1:n.res[k, 3]) {
      res.choice <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                           prob = as.numeric(sim.res.abund[[i]][, -1]))
      gen.mat[k, res.choice] <- gen.mat[k, res.choice] + 1
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill the perfect specialist matrix. 
  for (j in 1:Matrix.specifications[i, 2]) {  
    pref <- as.character(res.pref[[i]]$Resource[res.pref[[i]][, j + 1] == 1])
    spec.mat[spec.mat$Consumer.sp == colnames(res.pref[[i]][j + 1]), pref] <- 
      n.res[n.res$Matrix == paste("Mat.", i, sep = "") & 
              n.res$Consumer.sp == paste("C", j, sep = ""), 3]
  }
  # -------------------------------
  
  # -------------------------------
  # Fill mix.1 matrix. 
  shape.mat <- eval(parse(text = paste("Mix.1.shapes[Mix.1.shapes$Mat.num == ",
                                       i, ", ]", sep = "")))
  for (j in 1:Matrix.specifications[i, 2]) {   
    # Create preference profile from beta-binomial distribution
    pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                               shape1 = shape.mat[j, 3], shape2 = 1) / 
      max(dbetabinom.ab(0, size = 50, shape1 = 
                          shape.mat[j, 3], shape2 = 1))
    res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                             Ord = seq(1, nrow(res.pref[[i]])),
                             Cons.pref = res.pref[[i]][, j + 1])
    res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
    res.choice$Cons.pref <- pref.prof
    res.choice <- res.choice[with(res.choice, order(Ord)), ]
    weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
    
    for (k in 1:N.ind.cons) {
      for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
        res.choice <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                             prob = as.numeric(weighted.res))
        mix.1.mat[((j - 1) * N.ind.cons) + k, res.choice] <- 
          mix.1.mat[((j - 1) * N.ind.cons) + k, res.choice] + 1
      }
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Select the 1 or 2 consumer species that will be made into specialists and 
  #   seeded into matrices that contain perfect generalists
  sel.1 <- as.numeric(gsub("C", "", Specialist.spp[i, 2]))
  ifelse(is.na(Specialist.spp[i, 3]), sel.2 <- NA, 
         {sel.2 <- as.numeric(gsub("C", "", Specialist.spp[i, 3]))})
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 1 perfect specialist seeded in. 
  for (j in 1:Matrix.specifications[i, 2]) { 
    if(j == sel.1) {
      pref <- as.character(res.pref[[i]]$Resource[res.pref[[i]][, j + 1] == 1])
      seed.1.mat[seed.1.mat$Consumer.sp == colnames(res.pref[[i]][j + 1]), 
                 pref] <- n.res[n.res$Matrix == paste("Mat.", i, sep = "") & 
                                  n.res$Consumer.sp == paste("C", j, sep = ""), 3]
    } else {
      for (k in 1:N.ind.cons) {
        for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
          res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                        prob = as.numeric(sim.res.abund[[i]][-1]), 
                        replace = TRUE)
          seed.1.mat[((j - 1) * N.ind.cons) + k, res] <- 
            seed.1.mat[((j - 1) * N.ind.cons) + k, res] + 1
        }
      }
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 2 perfect specialists seeded in. 
  # Select the specialist.
  if(Matrix.specifications[i, 2] >= 2) {       
    for (j in 1:Matrix.specifications[i, 2]) { 
      if(j == sel.1 | j == sel.2) {
        pref <- as.character(res.pref[[i]]$Resource[res.pref[[i]][, j + 1] == 1])
        seed.2.mat[seed.2.mat$Consumer.sp == colnames(res.pref[[i]][j + 1]), 
                   pref] <- n.res[n.res$Matrix == paste("Mat.", i, sep = "") & 
                                    n.res$Consumer.sp == paste("C", j, sep = ""), 3]
      } else {
        for (k in 1:N.ind.cons) {
          for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
            res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                          prob = as.numeric(sim.res.abund[[i]][-1]), 
                          replace = TRUE)
            seed.2.mat[((j - 1) * N.ind.cons) + k, res] <- 
              seed.2.mat[((j - 1) * N.ind.cons) + k, res] + 1
          }
        }
      }
    }
  } else {seed.2.mat[, 3:ncol(seed.2.mat)] <- NA}
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 1 moderate specialist seeded in. 
  for (j in 1:Matrix.specifications[i, 2]) {  
    if(j == sel.1) {
      pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                 shape1 = 0.5, shape2 = 1) / 
        max(dbetabinom.ab(0, size = 50, shape1 = 
                            0.5, shape2 = 1))
      res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                               Ord = seq(1, nrow(res.pref[[i]])),
                               Cons.pref = res.pref[[i]][, sel.1 + 1])
      res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
      res.choice$Cons.pref <- pref.prof
      res.choice <- res.choice[with(res.choice, order(Ord)), ]
      weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
      
      for (k in 1:N.ind.cons) {
        for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
          res <- sample(colnames(weighted.res), 1,
                        prob = as.numeric(weighted.res))
          seed.1.mod[((j - 1) * N.ind.cons) + k, res] <- 
            seed.1.mod[((j - 1) * N.ind.cons) + k, res] + 1
        }
      }
    } else {
      for (k in 1:N.ind.cons) {
        for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
          res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                        prob = as.numeric(sim.res.abund[[i]][-1]))
          seed.1.mod[((j - 1) * N.ind.cons) + k, res] <- 
            seed.1.mod[((j - 1) * N.ind.cons) + k, res] + 1
        }
      }
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 1 strong specialist seeded in. 
  for (j in 1:Matrix.specifications[i, 2]) {
    if(j == sel.1) {
      pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                 shape1 = 0.2, shape2 = 1) / 
        max(dbetabinom.ab(0, size = 50, shape1 = 
                            0.2, shape2 = 1))
      res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                               Ord = seq(1, nrow(res.pref[[i]])),
                               Cons.pref = res.pref[[i]][, sel.1 + 1])
      res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
      res.choice$Cons.pref <- pref.prof
      res.choice <- res.choice[with(res.choice, order(Ord)), ]
      weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
      
      for (k in 1:N.ind.cons) {
        for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
          res <- sample(colnames(weighted.res), 1,
                        prob = as.numeric(weighted.res))
          seed.1.ss[((j - 1) * N.ind.cons) + k, res] <- 
            seed.1.ss[((j - 1) * N.ind.cons) + k, res] + 1
        }
      }
    } else {
      for (k in 1:N.ind.cons) {
        for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
          res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                        prob = as.numeric(sim.res.abund[[i]][-1]))
          seed.1.ss[((j - 1) * N.ind.cons) + k, res] <- 
            seed.1.ss[((j - 1) * N.ind.cons) + k, res] + 1
        }
      }
    }
  }
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 2 moderate specialists seeded in. 
  if(Matrix.specifications[i, 2] >= 2) {    
    for (j in 1:Matrix.specifications[i, 2]) {  
      if(j == sel.1 | j == sel.2) {
        pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                   shape1 = 0.5, shape2 = 1) / 
          max(dbetabinom.ab(0, size = 50, shape1 = 
                              0.5, shape2 = 1))
        res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                                 Ord = seq(1, nrow(res.pref[[i]])),
                                 Cons.pref = res.pref[[i]][, j + 1])
        res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
        res.choice$Cons.pref <- pref.prof
        res.choice <- res.choice[with(res.choice, order(Ord)), ]
        weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
        
        for (k in 1:N.ind.cons) {
          for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
            res <- sample(colnames(weighted.res), 1,
                          prob = as.numeric(weighted.res))
            seed.2.mod[((j - 1) * N.ind.cons) + k, res] <- 
              seed.2.mod[((j - 1) * N.ind.cons) + k, res] + 1
          }
        }
      } else {
        for (k in 1:N.ind.cons) {
          for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
            res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                          prob = as.numeric(sim.res.abund[[i]][-1]))
            seed.2.mod[((j - 1) * N.ind.cons) + k, res] <- 
              seed.2.mod[((j - 1) * N.ind.cons) + k, res] + 1
          }
        }
      }
    } 
  } else {seed.2.mod[, 3:ncol(seed.2.mod)] <- NA}
  # -------------------------------
  
  # -------------------------------
  # Fill matrix with 2 strong specialists seeded in. 
  if(Matrix.specifications[i, 2] >= 2) {
    for (j in 1:Matrix.specifications[i, 2]) {
      if(j == sel.1 | j == sel.2) {
        pref.prof <- dbetabinom.ab(seq(0, nrow(res.pref[[i]]) - 1), size = 50, 
                                   shape1 = 0.2, shape2 = 1) / 
          max(dbetabinom.ab(0, size = 50, shape1 = 
                              0.2, shape2 = 1))
        res.choice <- data.frame(Res = res.pref[[i]]$Resource, 
                                 Ord = seq(1, nrow(res.pref[[i]])),
                                 Cons.pref = res.pref[[i]][, j + 1])
        res.choice <- res.choice[with(res.choice, order(Cons.pref)), ]
        res.choice$Cons.pref <- pref.prof
        res.choice <- res.choice[with(res.choice, order(Ord)), ]
        weighted.res <- sim.res.abund[[i]][-1] * res.choice$Cons.pref
        
        for (k in 1:N.ind.cons) {
          for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
            res <- sample(colnames(weighted.res), 1,
                          prob = as.numeric(weighted.res))
            seed.2.ss[((j - 1) * N.ind.cons) + k, res] <- 
              seed.2.ss[((j - 1) * N.ind.cons) + k, res] + 1
          }
        }
      } else {
        for (k in 1:N.ind.cons) {
          for (l in 1:n.res[((j - 1) * N.ind.cons) + k, 3]) {
            res <- sample(colnames(sim.res.abund[[i]][-1]), 1,
                          prob = as.numeric(sim.res.abund[[i]][-1]))
            seed.2.ss[((j - 1) * N.ind.cons) + k, res] <- 
              seed.2.ss[((j - 1) * N.ind.cons) + k, res] + 1
          }
        }
      }
    } 
  } else {seed.2.ss[, 3:ncol(seed.2.ss)] <- NA}
  # -------------------------------
  
  # -------------------------------
  perfect.generalist[[i]] <- gen.mat
  perfect.specialist[[i]] <- spec.mat
  mix.1[[i]] <- mix.1.mat
  add.1.spec [[i]] <- seed.1.mat
  add.2.spec [[i]] <- seed.2.mat
  add.1.mod.spec[[i]] <- seed.1.mod
  add.2.mod.spec[[i]] <- seed.2.mod
  add.1.str.spec[[i]] <- seed.1.ss
  add.2.str.spec[[i]] <- seed.2.ss
  
  ## Progress count
  Sys.sleep(0.02); print(i); flush.console() 
  
}
# --------------------------------------------------------------------------





# ==============================================================================
# Stage 4: generate benchmark matrices for data type ii (i.e. lists of 
#          resource species with which each individual consumer interacts 
#          (i.e. nominal data)) from the count (type iii) data.
# --------------------------------------------------------------------------
# Derived directly from the count data in Stage 3, so that code needs to be 
#  run first.
# ==============================================================================


perfect.generalist.PA <- list()
perfect.specialist.PA <- list()
mix.1.PA <- list()
add.1.spec.PA <- list()
add.2.spec.PA <- list()
add.1.mod.spec.PA <- list()
add.2.mod.spec.PA <- list()
add.1.str.spec.PA <- list()
add.2.str.spec.PA <- list()


for (h in 1:N.matrices) {
  perfect.generalist.PA[[h]] <- perfect.generalist[[h]]
  aa <- perfect.generalist.PA[[h]][, 3:ncol(perfect.generalist.PA[[h]])]
  aa[aa > 1] <- 1
  perfect.generalist.PA[[h]][, 3:ncol(perfect.generalist.PA[[h]])] <- aa
  
  perfect.specialist.PA[[h]] <- perfect.specialist[[h]]
  aa <- perfect.specialist.PA[[h]][, 3:ncol(perfect.specialist.PA[[h]])]
  aa[aa > 1] <- 1
  perfect.specialist.PA[[h]][, 3:ncol(perfect.specialist.PA[[h]])] <- aa
  
  mix.1.PA[[h]] <- mix.1[[h]]
  aa <- mix.1.PA[[h]][, 3:ncol(mix.1.PA[[h]])]
  aa[aa > 1] <- 1
  mix.1.PA[[h]][, 3:ncol(mix.1.PA[[h]])] <- aa
  
  add.1.spec.PA[[h]] <- add.1.spec[[h]]
  aa <- add.1.spec.PA[[h]][, 3:ncol(add.1.spec.PA[[h]])]
  aa[aa > 1] <- 1
  add.1.spec.PA[[h]][, 3:ncol(add.1.spec.PA[[h]])] <- aa
  
  add.2.spec.PA[[h]] <- add.2.spec[[h]]
  aa <- add.2.spec.PA[[h]][, 3:ncol(add.2.spec.PA[[h]])]
  aa[aa > 1] <- 1
  add.2.spec.PA[[h]][, 3:ncol(add.2.spec.PA[[h]])] <- aa
  
  add.1.mod.spec.PA[[h]] <- add.1.mod.spec[[h]]
  aa <- add.1.mod.spec.PA[[h]][, 3:ncol(add.1.mod.spec.PA[[h]])]
  aa[aa > 1] <- 1
  add.1.mod.spec.PA[[h]][, 3:ncol(add.1.mod.spec.PA[[h]])] <- aa
  
  add.2.mod.spec.PA[[h]] <- add.2.mod.spec[[h]]
  aa <- add.2.mod.spec.PA[[h]][, 3:ncol(add.2.mod.spec.PA[[h]])]
  aa[aa > 1] <- 1
  add.2.mod.spec.PA[[h]][, 3:ncol(add.2.mod.spec.PA[[h]])] <- aa
  
  add.1.str.spec.PA[[h]] <- add.1.str.spec[[h]]
  aa <- add.1.str.spec.PA[[h]][, 3:ncol(add.1.str.spec.PA[[h]])]
  aa[aa > 1] <- 1
  add.1.str.spec.PA[[h]][, 3:ncol(add.1.str.spec.PA[[h]])] <- aa
  
  add.2.str.spec.PA[[h]] <- add.2.str.spec[[h]]
  aa <- add.2.str.spec.PA[[h]][, 3:ncol(add.2.str.spec.PA[[h]])]
  aa[aa > 1] <- 1
  add.2.str.spec.PA[[h]][, 3:ncol(add.2.str.spec.PA[[h]])] <- aa
  
}
# --------------------------------------------------------------------------



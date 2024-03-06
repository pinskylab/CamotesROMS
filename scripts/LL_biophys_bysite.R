# calculates the log-likelihoods of dispersal observations conditional on an oceanographic dispersal matrix
# from Bode et al. 2019 PLOS Bio
# like LL_biophys.R, but returns values separately for each destination site
# the one argument is a list with the following components:
#       sampled_reefs_vec: indices of the sites that were sampled
#       pop_size_vec: a vector of population size at each sampling site
#       BioPhysMat: a matrix of simulated larval dispersal from oceanographic currents, produced in part by format_biophys_normalized_matrix()
#       prop_samp_vec: proportion of anemones at each site that were sampled
#       unassigned_vec: vector of sampled juveniles without assigned parents
#       Assignments: a matrix of parentage assignments of larva, produced in part by format_genetic_parentage_matrix()


LL_biophys_bysite <- function(data){
    sampled_reefs_vec <- data$sampled_reefs_vec
    pop_size_vec <- data$pop_size_vec
    BioPhysMat <- data$BioPhysMat
    prop_samp_vec <- data$prop_samp_vec 
    unassigned_vec <- data$unassigned_vec
    Assignments <- data$Assignments
    
    #term 1 parenthesis- probability that sampled recruit is unassigned
    #multiply the source normalized connectivity matrix by a vector of adult abundance at each source to generate an expected connectivity matrix of the number of recruits
    ExpectedRecruits = (BioPhysMat)*(matrix(pop_size_vec, nrow = nrow(BioPhysMat), ncol=nrow(BioPhysMat)))
    UnassignedRecruits = matrix(0, nrow=nrow(prop_samp_vec), ncol=nrow(prop_samp_vec)) # rows are sources, columns are
    
    for(i in 1:nrow(prop_samp_vec)){
        prop_samp_source = prop_samp_vec[i] #value of prop_samp at time point t at site i (source), This_SS_A
        for(j in 1:nrow(prop_samp_vec)){
            RecruitsFromAssignedReefs = ExpectedRecruits[i,j]
            UnassignedRecruits[i,j] = RecruitsFromAssignedReefs*(1-prop_samp_source)^2
        }
    }
    #replace the -Inf values with 0? -Inf happens when the prop_samp for a source is 0, because in the equation above those routes have a zero value. ln(0) = -Inf in R
    UnassignedRecruitsLn <- log(UnassignedRecruits)
    UnassignedRecruitsLn[is.infinite(UnassignedRecruitsLn)] <- 0
    
    #term 2 parenthesis- probability that sampled recruit is assigned
    #multiply the source normalized connectivity matrix by a vector of adult abundance at each source to generate an expected connectivity matrix of the number of recruits
    ExpectedRecruits = (BioPhysMat)*(matrix(pop_size_vec, nrow = nrow(BioPhysMat), ncol=nrow(BioPhysMat)))
    AssignedRecruits = matrix(0, nrow=nrow(prop_samp_vec),ncol=nrow(prop_samp_vec)) #should there be an extra row for unassigned- NO
    
    for(i in 1:nrow(prop_samp_vec)){
        prop_samp_source = prop_samp_vec[i] #value of prop_samp at time point t at site i (source), This_SS_A
        for(j in 1:nrow(prop_samp_vec)){
            RecruitsFromAssignedReefs = ExpectedRecruits[i,j]
            AssignedRecruits[i,j] = RecruitsFromAssignedReefs*(prop_samp_source^2 + 2*prop_samp_source*(1 - prop_samp_source))
        }
    }
    #replace the -Inf values with 0? -Inf happens when the prop_samp for a source is 0, because in the equation above those routes have a zero value. ln(0) = -Inf in R
    AssignedRecruitsLn <- log(AssignedRecruits)
    AssignedRecruitsLn[is.infinite(AssignedRecruitsLn)] <- 0
    
    #to complete term 1, 
    ProbUnassignedFromParentage = (UnassignedRecruitsLn)*(matrix(unassigned_vec, nrow = nrow(UnassignedRecruitsLn), ncol=nrow(UnassignedRecruitsLn))) #this is bxjt in Eqn.S3.4
    term1 <- as.matrix(colSums(ProbUnassignedFromParentage)) # sum across source sites. this completes term1
    
    #to complete term 2, 
    ProbAssignedFromParentage = (AssignedRecruitsLn)*(matrix(Assignments, nrow = nrow(AssignedRecruitsLn), ncol=nrow(AssignedRecruitsLn))) #this is bijt in Eqn.S3.4
    term2 <- as.matrix(colSums(ProbAssignedFromParentage)) # sum across source sites. this completes term1
    
    #calculate  log likelihood. leave separate for each destination site.
    LLbysite <- term1 + term2
    return(LLbysite)
    
}
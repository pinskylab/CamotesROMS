PredProp <- function (k, theta, Assignments, Sampled_reefs, Distances, Reef_sizes, 
    Adult_sample_proportions) 
{
    d_pred <- function(k, theta, Distances) {
        d <- Distances
        disp <- exp(k) * theta * exp(-(exp(k) * d)^theta)/gamma(1/theta)
        return(disp)
    }
    NumReefs = nrow(Distances)
    NumSampledReefs = length(Sampled_reefs)
    Proportions = d_pred(k, theta, Distances)
    Settlers = Proportions * (matrix(Reef_sizes, nrow = nrow(Distances), 
        ncol = nrow(Distances)))
    AllSettlers = colSums(Settlers)
    AssignedSettlers = matrix(0, nrow = NumSampledReefs + 1, 
        ncol = NumSampledReefs)
    for (i in 1:NumSampledReefs) {
        This_SS_A = Adult_sample_proportions[i]
        for (j in 1:NumSampledReefs) {
            SettlersFromAssignedReefs = Settlers[Sampled_reefs[i], 
                Sampled_reefs[j]]
            AssignedSettlers[i, j] = SettlersFromAssignedReefs * 
                (This_SS_A^2 + 2 * This_SS_A * (1 - This_SS_A))
            AssignedSettlers[NumSampledReefs + 1, j] = AssignedSettlers[NumSampledReefs + 
                1, j] + SettlersFromAssignedReefs * (1 - This_SS_A)^2
        }
    }
    Unsampled = as.matrix(setdiff(1:NumReefs, Sampled_reefs))
    for (j in 1:length(Sampled_reefs)) {
        AssignedSettlers[NumSampledReefs + 1, j] = sum(Settlers[Unsampled, 
            Sampled_reefs[j]]) + AssignedSettlers[NumSampledReefs + 
            1, j]
    }
    PredictedProportions = AssignedSettlers/(matrix(rep(t(colSums(AssignedSettlers)), 
        NumSampledReefs + 1), ncol = ncol(AssignedSettlers), 
        byrow = TRUE))
    colSums(PredictedProportions)
    PredictedProportions[PredictedProportions == 0] = 1e-12
    PredictedProportions[PredictedProportions <= 1e-12] = 1e-12
    PredictedProportions = PredictedProportions/(matrix(rep(t(colSums(PredictedProportions)), 
        nrow(PredictedProportions)), ncol = ncol(PredictedProportions), 
        byrow = TRUE))
    log_like = 0
    for (j in 1:NumSampledReefs) {
        ObsVector = Assignments[, j]
        ProbVector = PredictedProportions[, j]
        x <- 1e+12
        y <- sum(ObsVector * log(ProbVector))
        log_like = log_like + ifelse(is.finite(y), y, x)
    }
    log_like = -log_like
    return(PredictedProportions)
}

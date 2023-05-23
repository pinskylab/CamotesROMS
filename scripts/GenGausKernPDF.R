predicted_disp <- function(d, k, theta){ 
    
    disp <- exp(k)*theta*exp(-(exp(k)*d)^theta)/gamma(1/theta)
    return(disp)

}

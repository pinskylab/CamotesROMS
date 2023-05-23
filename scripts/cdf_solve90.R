cdf_solve90 <- function (d, theta = theta_eval, k = k_eval) 
{
    return(cdf(d, theta, k) - 0.45)$x
}
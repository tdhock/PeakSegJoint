#include <math.h>
#include <stdio.h>

double OptimalPoissonLoss(double cumsum_value, double mean_value){
  if(cumsum_value == 0){
    return 0.0;
  }
  double log_mean = log(mean_value);
  double mean_term = 1-log_mean;
  //printf("mean_value=%f log_mean=%f mean_term=%f cumsum_value=%f\n", mean_value, log_mean, mean_term, cumsum_value);
  return cumsum_value * mean_term;
}


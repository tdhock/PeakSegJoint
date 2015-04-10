#include <math.h>

double OptimalPoissonLoss(int cumsum_value, double mean_value){
  if(cumsum_value == 0){
    return 0.0;
  }
  return cumsum_value * (1-log(mean_value));
}


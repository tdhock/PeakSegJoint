#include <math.h>
#include "OptimalPoissonLoss.h"

/*
  solver for the DPA recursion: 
  L_{s,t} = min_{t' < t} L_{s-1, t'} + c_(t', t]
  where L is the optimal loss in s segments up to t
  and c is the optimal loss of segment (t', t]

  - PrevSegs_loss_vec is the L_{s-1} vector,
  - sample_cumsum_mat is the matrix of cumsums
    which is used to compute the optimal loss c.
  - first_possible_index = s.
  - last_possible_index = t.
  - LastSeg_FirstIndex = t'.
 */
void get_best_FirstIndex(
  double *PrevSegs_loss_vec, // n_data
  int n_data,
  int *sample_cumsum_mat, // n_data x n_samples
  int n_samples,
  int first_possible_index,
  int last_possible_index,
  int bases_per_bin,
  int *best_FirstIndex, //output
  double *best_loss //output
  ){
  int sample_i;
  int* cumsum_vec, cumsum_value, cumsum_prev_end, cumsum_this_end;
  int LastSeg_FirstIndex;
  double candidate_loss, mean_value, bases_value,
    prev_loss, this_loss, loss_value;
  *best_loss = INFINITY;
  for(LastSeg_FirstIndex = first_possible_index; 
      LastSeg_FirstIndex <= last_possible_index;
      LastSeg_FirstIndex++){
    // start with previous segment optimal loss.
    this_loss = 0.0;
    for(sample_i=0; sample_i < n_samples; sample_i++){
      cumsum_vec = sample_cumsum_mat + n_data * sample_i;
      /*  
	  the optimal loss of segment (t', t]
	  can be computed in O(1) time given the vector of cumsums.
      */
      cumsum_prev_end = cumsum_vec[LastSeg_FirstIndex-1];
      cumsum_this_end = cumsum_vec[last_possible_index];
      cumsum_value = cumsum_this_end - cumsum_prev_end;
      bases_value = bases_per_bin*(last_possible_index-LastSeg_FirstIndex+1);
      mean_value = ((double)cumsum_value)/bases_value;
      loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
      /* printf("%dsample%d prev=%d this=%d bases=%f mean=%f loss=%f\n",  */
      /* 	     LastSeg_FirstIndex, sample_i, */
      /* 	     //cumsum_value,  */
      /* 	     cumsum_prev_end, cumsum_this_end, */
      /* 	     bases_value, mean_value, loss_value); */
      this_loss += loss_value;
    }
    prev_loss = PrevSegs_loss_vec[LastSeg_FirstIndex-1];
    /* printf("%d prev=%f this=%f\n", LastSeg_FirstIndex,  */
    /* 	   prev_loss, this_loss); */
    candidate_loss = prev_loss + this_loss;
    if(candidate_loss < *best_loss){
      //printf("LastSeg_FirstIndex=%d\n", LastSeg_FirstIndex);
      *best_loss = candidate_loss;
      *best_FirstIndex = LastSeg_FirstIndex;
    }
  }
}

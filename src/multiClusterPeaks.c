/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiClusterPeaks.h"
#include <stdlib.h>

int multiClusterPeaks
(int *peakStart, int *peakEnd, int *sample, int n_peaks,
 int *cluster){//out
  int peak_i, n_samples=0, candidate_samples;
  for(peak_i=0; peak_i<n_peaks; peak_i++){
    if(sample < 0){
      return ERROR_NEGATIVE_SAMPLE_ID;
    }
    candidate_samples = sample[peak_i] + 1;
    if(n_samples < candidate_samples){
      n_samples = candidate_samples;
    }
  }
  int *start_vec = (int*)malloc(sizeof(int)*n_samples);
  int *end_vec = (int*)malloc(sizeof(int)*n_samples);
  // start_vec and end_vec contain the candidates for the current peak
  // cluster, and are initialized to -1 to indicate unused.
  int sample_i;
  for(sample_i=0; sample_i<n_samples; sample_i++){
    start_vec[sample_i] = -1;
    end_vec[sample_i] = -1;
  }
  for(peak_i=0; peak_i<n_peaks; peak_i++){
    //keep adding peaks to start_vec and end_vec
    
    //end a cluster when next start is after the first end of this
    //cluster -- this means that two peaks on the same sample can
    //never be in the same cluster.
  }  
  free(start_vec);
  free(end_vec);
  return 0;
}

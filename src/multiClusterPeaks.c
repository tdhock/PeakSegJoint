/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiClusterPeaks.h"
#include <stdlib.h>

int multiClusterPeaks
(int *peakStart, int *peakEnd, int *cluster, int n_peaks){//out
  int peak_i;
  int cluster_i = 0;
  int cluster_first_end = 0;
  for(peak_i=0; peak_i<n_peaks; peak_i++){
    //end a cluster when next start is after the first end of this
    //cluster -- this means that two peaks on the same sample can
    //never be in the same cluster.
    if(cluster_first_end < peakStart[peak_i]){
      cluster_i++;
      cluster_first_end = peakEnd[peak_i];
    }else{
      //this is the same cluster, so update cluster_first_end only if
      //necessary.
      if(peakEnd[peak_i] < cluster_first_end){
	cluster_first_end = peakEnd[peak_i];
      }
    }
    cluster[peak_i]=cluster_i;
  }  
  return 0;
}

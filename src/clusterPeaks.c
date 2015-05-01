#include "clusterPeaks.h"

int clusterPeaks(int *peakStart, int *peakEnd,
		 int *cluster,
		 int peaks){
  int peak_i=1, cluster_i=0, 
    clusterStart=peakStart[0], 
    clusterEnd=peakEnd[0];
  cluster[0] = 0;
  while(peak_i < peaks){
    if(peakStart[peak_i] < clusterStart){
      return ERROR_PEAKSTART_DECREASING;
    }
    if(peakStart[peak_i] < clusterEnd){
      if(clusterEnd < peakEnd[peak_i]){
	clusterEnd = peakEnd[peak_i];
      }
    }else{
      clusterStart = peakStart[peak_i];
      clusterEnd = peakEnd[peak_i];
      cluster_i++;
    }
    cluster[peak_i] = cluster_i;
    peak_i++;
  }
  return 0;
}

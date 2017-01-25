#include "multiClusterPeaks.h"
#include <R.h>

void multiClusterPeaks_interface
(int *peakStart, int *peakEnd, int *sample,
 int *cluster, 
 int *peaks) {
  int status;
  status = multiClusterPeaks
    (peakStart, peakEnd, sample, 
     *peaks,
     cluster);
  if(status != 0){
    error("unrecognized error code");
  }
}


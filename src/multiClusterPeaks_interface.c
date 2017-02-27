#include "multiClusterPeaks.h"
#include <R.h>

void multiClusterPeaks_interface
(int *peakStart, int *peakEnd, int *cluster, 
 int *peaks) {
  int status;
  status = multiClusterPeaks
    (peakStart, peakEnd, cluster,
     *peaks);
  if(status != 0){
    error("unrecognized error code");
  }
}


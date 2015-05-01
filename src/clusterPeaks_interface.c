#include "clusterPeaks.h"
#include <R.h>

void clusterPeaks_interface
(int *peakStart, int *peakEnd,
 int *cluster,
 int *peaks) {
  int status;
  status = clusterPeaks(peakStart, peakEnd, 
			cluster,
			*peaks);
  if(status == ERROR_PEAKSTART_DECREASING){
    error("peakStart decreasing");
  }
  if(status != 0){
    error("unrecognized error code");
  }
}


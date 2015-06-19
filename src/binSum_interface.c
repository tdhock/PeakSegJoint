#include "binSum.h"
#include <R.h>

void binSum_interface(
  int *profile_chromStart,
  int *profile_chromEnd,
  int *profile_coverage,
  int *n_profiles,
  int *bin_total,
  int *bin_size,
  int *n_bins,
  int *bin_chromStart){
  int status;
  status = binSum(profile_chromStart, 
		  profile_chromEnd,
		  profile_coverage,
		  *n_profiles,
		  bin_total,
		  *bin_size,
		  *n_bins,
		  *bin_chromStart, 0);
  if(status == ERROR_CHROMSTART_NOT_LESS_THAN_CHROMEND){
    error("chromStart not less than chromEnd");
  }
  if(status == ERROR_CHROMSTART_BEFORE_PREVIOUS_CHROMEND){
    error("chromStart before previous chromEnd");
  }
  if(status == ERROR_CHROMSTART_CHROMEND_MISMATCH){
    error("chromStart[i] != chromEnd[i-1]");
  }
  if(status != 0){
    error("error code %d", status);
  }
}

/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "binSum.h"
#include <stdlib.h>
#include <stdio.h>

int oneBin
(int *profile_chromStart, 
 int *profile_chromEnd, 
 int *profile_coverage, 
 int n_profiles,
 int *bin_total,
 int bin_chromStart,
 int bin_chromEnd){
  return binSum(
    profile_chromStart,
    profile_chromEnd,
    profile_coverage,
    n_profiles,
    bin_total,
    bin_chromEnd - bin_chromStart,
    1,
    bin_chromStart,
    EMPTY_AS_ZERO);
}

int binSum
(int *profile_chromStart, 
 int *profile_chromEnd, 
 int *profile_coverage, 
 int n_profiles,
 int *bin_total, 
 int bin_size,
 int n_bins, 
 int bin_chromStart,
 int status_for_empty_bin){
  int profile_i, bin_i;
  // check that chromEnd < chromStart for all profile data.
  for(profile_i = 0; profile_i < n_profiles; profile_i++){
    if(profile_chromEnd[profile_i] <= profile_chromStart[profile_i]){
      return ERROR_CHROMSTART_NOT_LESS_THAN_CHROMEND;
    }
  }
  // check that chromEnd[i-1] <= chromStart[i] for all i>0.
  for(profile_i = 1; profile_i < n_profiles; profile_i++){
    if(profile_chromStart[profile_i] < profile_chromEnd[profile_i-1]){
      return ERROR_CHROMSTART_BEFORE_PREVIOUS_CHROMEND;
    }
  }
  /* printf("bin_chromStart=%d bin_size=%d n_bins=%d status=%d\n", */
  /* 	 bin_chromStart, bin_size, n_bins, status_for_empty_bin); */
  int *bin_touched = (int*) malloc(sizeof(int) * n_bins);
  for(bin_i = 0; bin_i < n_bins; bin_i++){
    bin_total[bin_i] = 0;
    bin_touched[bin_i] = 0;
  }
  // bin_chromStart gives the base before the first position that we
  // want to count, for example 1000 means we want to start counting at
  // 1001. so we should ignore profile entries of (0, 10], (0, 1000], 
  // (999, 1000], but start counting (0, 1001], (1000, 1001], (1000, 1002].
  bin_i = 0;
  profile_i = 0;
  while(profile_i < n_profiles && profile_chromEnd[profile_i] <= bin_chromStart){
    profile_i ++;
  }
  int count_until, bases, bin_add, profile_add;
  int begin_count_after;
  int first_chromStart = profile_chromStart[profile_i];
  if(bin_chromStart < first_chromStart){
    begin_count_after = first_chromStart;
  }else{
    begin_count_after = bin_chromStart;
  }
  int bin_end = bin_chromStart + bin_size;
  while(bin_i < n_bins && profile_i < n_profiles){
    // at this point there are two cases.
    if(bin_end <= profile_chromEnd[profile_i]){
      // 1. the profile segment continues to the end of this bin,
      //    so add profile_coverage * bin_size to this bin total.
      // -profile----]
      //             (-----------
      //          bin]
      //         bin]
      //  bin]
      count_until = bin_end;
      if(bin_end == profile_chromEnd[profile_i]){
	profile_add = 1; // done adding from this profile segment.
      }else{
	profile_add = 0; // not done adding from this profile segment.
      }
      bin_add = 1;
    }else{      
      // 2. the profile segment ends before this bin ends.
      // -profile----]
      //             (-----------
      //           bin]
      //             bin]
      count_until = profile_chromEnd[profile_i];
      profile_add = 1; // done adding from this profile segment.
      bin_add = 0; // not done adding to this bin total.
    }
    bases = count_until - begin_count_after;
    bin_total[bin_i] += profile_coverage[profile_i] * bases;
    bin_touched[bin_i] = 1;
    // setup next iteration.
    begin_count_after = count_until;
    profile_i += profile_add;
    if(bin_add){
      bin_i++;
      bin_end += bin_size;
    }
  }
  // If EMPTY_AS_ZERO flag, return now (untouched totals are zero).
  if(status_for_empty_bin == EMPTY_AS_ZERO){
    free(bin_touched);
    return 0;
  }
  // If there was no data at all that overlapped a bin (not even
  // zeros), then set it to -1 to mark that.
  int status = 0;
  for(bin_i=0; bin_i < n_bins; bin_i++){
    if(bin_touched[bin_i] == 0){
      bin_total[bin_i] = -1;
      status = status_for_empty_bin;
    }
  }
  free(bin_touched);
  return status;
}

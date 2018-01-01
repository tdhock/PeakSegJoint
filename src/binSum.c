/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "binSum.h"
#include <stdlib.h>
#include <stdio.h>

int
binSumLR
(int *data_start_end,
 int *chromStart, int *chromEnd,
 int *coverage, int n_entries,
 int *left_bin_vec, int *right_bin_vec,
 int left_chromStart, int right_chromStart,
 int bases_per_bin, int n_bins){
  int bin_chromEnd, bin_chromStart;
  int extra_chromStart, extra_chromEnd, extra_bases, extra_coverage;
  int status;
  /* printf("left bin_size=%d bins=%d start=%d\n",  */
  /* 	 bases_per_bin, n_bins, left_chromStart); */
  status = binSum(chromStart, chromEnd,
		  coverage, n_entries,
		  left_bin_vec,
		  bases_per_bin,
		  n_bins,
		  left_chromStart,
		  EMPTY_AS_ZERO);
  if(status != 0){
    return status;
  }
  /* printf("right bin_size=%d bins=%d start=%d\n",  */
  /* 	 bases_per_bin, n_bins, right_chromStart); */
  status = binSum(chromStart, chromEnd,
		  coverage, n_entries,
		  right_bin_vec,
		  bases_per_bin,
		  n_bins,
		  right_chromStart,
		  EMPTY_AS_ZERO);
  if(status != 0){
    return status;
  }
  for(int bin_i=0; bin_i < n_bins; bin_i++){
    //left bin.
    bin_chromStart = left_chromStart + bases_per_bin * bin_i;
    bin_chromEnd = bin_chromStart + bases_per_bin;
    if(data_start_end[0] < bin_chromEnd){
      if(data_start_end[0] <= bin_chromStart){
	//   (     data   ]
	//   (bin]
	//       (bin]
	// bin is completely data, leave it alone!
      }else{
	//    (  data  ]
	//   (bin]
	//    - extra
	// (bin]
	//  --- extra
	// bin has some data, so subtract the extra.
	extra_chromStart = bin_chromStart;
	extra_chromEnd = data_start_end[0];
	extra_bases = extra_chromEnd - extra_chromStart;
	//printf("left start=%d bases=%d\n", extra_chromStart, extra_bases);
	status = binSum(chromStart, chromEnd,
			coverage, n_entries,
			&extra_coverage,
			extra_bases,
			1,
			extra_chromStart,
			EMPTY_AS_ZERO);
	if(status != 0){
	  return status;
	}
	left_bin_vec[bin_i] -= extra_coverage;
      }
    }else{
      //      (  data  ]
      //  (bin]
      // bin does not overlap data, so set it to zero.
      left_bin_vec[bin_i] = 0;
    }
    //right bin.
    bin_chromStart = right_chromStart + bases_per_bin * bin_i;
    bin_chromEnd = bin_chromStart + bases_per_bin;
    if(bin_chromStart < data_start_end[1]){
      if(bin_chromEnd <= data_start_end[1]){
	// (  data  ]
	//      (bin]    
	//  (bin]
	// bin is completely data, leave it alone!
      }else{
	// (  data  ]
	//        (bin]
	//           -- extra
	extra_chromStart = data_start_end[1];
	extra_chromEnd = bin_chromEnd;
	extra_bases = extra_chromEnd - extra_chromStart;
	//printf("right start=%d bases=%d\n", extra_chromStart, extra_bases);
	status = binSum(chromStart, chromEnd,
			coverage, n_entries,
			&extra_coverage,
			extra_bases,
			1,
			extra_chromStart,
			EMPTY_AS_ZERO);
	if(status != 0){
	  return status;
	}
	right_bin_vec[bin_i] -= extra_coverage;
      }
    }else{
      // (  data  ]
      //          (bin]
      //              (bin]
      // bin does not overlap data, so set it to zero.
      right_bin_vec[bin_i] = 0;
    }
  }
  return 0;
}

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
  int bin_end = bin_chromStart + bin_size;
  while(bin_i < n_bins && profile_i < n_profiles){
    //printf("bin_i=%d profile_i=%d\n", bin_i, profile_i);
    // at this point there are two cases.
    if(bin_end - bin_size < profile_chromStart[profile_i]){
      begin_count_after = profile_chromStart[profile_i];
    }else{
      begin_count_after = bin_end - bin_size;
    }
    if(bin_end <= profile_chromEnd[profile_i]){
      // 1. the profile segment ends after the end of this bin,
      //    so add profile_coverage * bin_size to this bin total.
      // -profile----]
      //             (-----------
      //          bin]
      //         bin]
      //  bin]
      //        (profile] or maybe it is not overlapping!
      count_until = bin_end;
      if(bin_end < profile_chromStart[profile_i]){
	//no overlap!
	begin_count_after = count_until;
      }
      if(bin_end == profile_chromEnd[profile_i]){
	profile_add = 1; // done adding from this profile segment.
      }else{
	profile_add = 0; // not done adding from this profile segment.
      }
      bin_add = 1;
    }else{      
      // 2. the profile segment ends before this bin ends.
      //   -profile----]
      //               (-----------
      //              or gap (------
      //             (bin]bin]bin] next bin with no profile.
      //           ( bin ] bin ] next bin with some profile.
      count_until = profile_chromEnd[profile_i];
      profile_add = 1; // done adding from this profile segment.
      bin_add = 0;
    }
    bases = count_until - begin_count_after;
    /* printf("begin_count_after=%d count_until=%d coverage=%d\n", */
    /* 	   begin_count_after, count_until, */
    /* 	   profile_coverage[profile_i]); */
    bin_total[bin_i] += profile_coverage[profile_i] * bases;
    bin_touched[bin_i] = 1;
    profile_i += profile_add;
    // setup next iteration.
    if(bin_add){
      bin_i++;
      bin_end += bin_size;
    }
  }
  /* printf("binSum %d data", n_profiles); */
  /* for(bin_i=0; bin_i < n_bins; bin_i++){ */
  /*   printf(" %d", bin_total[bin_i]); */
  /* } */
  /* printf("\n"); */

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

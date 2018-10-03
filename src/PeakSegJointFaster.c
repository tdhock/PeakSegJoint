/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegJointFaster.h"
#include "OptimalPoissonLoss.h"
#include "binSum.h"
#include <stdlib.h>
#include <math.h>
//#include <R.h>

int PeakSegJointFaster(
  struct ProfileList *profile_list,
  int bin_factor,
  double *mean_mat,
  double *flat_loss_vec,
  double *peak_loss_vec,
  int *peak_start_end,
  int *data_start_end
  ){
  int n_samples = profile_list->n_profiles;
  if(n_samples == 0){
    return ERROR_FASTER_NO_COVERAGE_DATA;
  }
  int chromStart, chromEnd, unfilled_chromStart, unfilled_chromEnd;
  double seg1_mean=-1, seg2_mean, seg3_mean;
  struct Profile *profile, *samples = profile_list->profile_vec;
  profile = samples;
  unfilled_chromEnd = get_max_chromEnd(profile);
  unfilled_chromStart = get_min_chromStart(profile);
  for(int sample_i=1; sample_i < n_samples; sample_i++){
    profile = samples + sample_i;
    chromStart = get_min_chromStart(profile);
    if(chromStart < unfilled_chromStart){
      unfilled_chromStart = chromStart;
    }
    chromEnd = get_max_chromEnd(profile);
    if(unfilled_chromEnd < chromEnd){
      unfilled_chromEnd = chromEnd;
    }
  }
  data_start_end[0] = unfilled_chromStart;
  data_start_end[1] = unfilled_chromEnd;
  int unfilled_bases = unfilled_chromEnd - unfilled_chromStart;
  double data_bases = (double) unfilled_bases;
  double bin_bases;
  if(unfilled_bases/bin_factor < 4){
    /*
      4 is smallest the number of data points for which the 3-segment
      optimization problem is not trivial.

      If we don't have at least this many data points for the first
      bin step, than we stop with an error.
    */
    return ERROR_FASTER_BIN_FACTOR_TOO_LARGE;
  }
  int bases_per_bin = 1;
  while(unfilled_bases/bases_per_bin/bin_factor >= 4){
    bases_per_bin *= bin_factor;
  }
  int n_bins = unfilled_bases / bases_per_bin;
  /*
    MaxBinSize() from line 1 of the JointZoom algorithm of the
    PeakSegJoint paper returns the value of the C variable
    bases_per_bin.

    Little b in the text of the section that describes the JointZoom
    algorithm is the C variable n_bins.
  */
  if(unfilled_bases % bases_per_bin != 0){
    n_bins ++ ;
  }
  if(n_bins < bin_factor*2){
    n_bins = bin_factor*2;
  }
  int extra_bases = n_bins  * bases_per_bin - unfilled_bases;
  int extra_before = extra_bases/2;
  int extra_after = extra_bases - extra_before;
  //Rprintf("bin_factor=%d n_bins=%d bases_per_bin=%d extra_before=%d extra_after=%d\n", bin_factor, n_bins, bases_per_bin, extra_before, extra_after);
  //int extra_count;
  int seg1_chromStart = unfilled_chromStart - extra_before;
  int seg3_chromEnd = unfilled_chromEnd + extra_after;

  int *sample_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  double *last_cumsum_vec = (double*) malloc(n_samples*sizeof(double));
  int *sample_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));

  double *seg1_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg2_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg3_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  
  double *seg1_loss_vec = (double*)malloc(sizeof(double)*n_samples);
  double *candidate_loss_vec = (double*)malloc(sizeof(double)*n_samples);

  int *left_bin_vec = (int*) malloc(n_bins * sizeof(int));
  int *right_bin_vec = (int*) malloc(n_bins * sizeof(int));
  int *left_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *right_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *left_initial_cumsum_vec = (int*) malloc(n_samples * sizeof(int));
  int *right_initial_cumsum_vec = (int*) malloc(n_samples * sizeof(int));
  
  int *left_cumsum_vec, *right_cumsum_vec;
  int left_cumsum_value, right_cumsum_value;
  int peakStart, peakEnd;
  int best_seg1, best_seg2;

  int *count_vec, *cumsum_vec, cumsum_value;
  int status;
  for(int sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples + sample_i;
    count_vec = sample_count_mat + n_bins*sample_i;
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    count_vec,
		    bases_per_bin, n_bins, seg1_chromStart, 
		    EMPTY_AS_ZERO);
    if(status != 0){
      free(sample_count_mat);
      free(last_cumsum_vec);
      free(sample_cumsum_mat);
  
      free(seg1_mean_vec);
      free(seg2_mean_vec);
      free(seg3_mean_vec);

      free(seg1_loss_vec);
      free(candidate_loss_vec);

      free(left_bin_vec);
      free(right_bin_vec);
      free(left_cumsum_mat);
      free(right_cumsum_mat);
      free(left_initial_cumsum_vec);
      free(right_initial_cumsum_vec);
      return status;
    }
  }//for sample_i
  int bin_i, offset;
  for(int sample_i=0; sample_i < n_samples; sample_i++){
    cumsum_value = 0;
    offset = n_bins * sample_i;
    count_vec = sample_count_mat + offset;
    cumsum_vec = sample_cumsum_mat + offset;
    //printf("[sample%02d] ", sample_i);
    for(bin_i=0; bin_i < n_bins; bin_i++){
      cumsum_value += count_vec[bin_i];
      //printf("%d ", cumsum_value);
      cumsum_vec[bin_i] = cumsum_value;
    }
    last_cumsum_vec[sample_i] = cumsum_value;
    seg1_mean = cumsum_value / data_bases;
    flat_loss_vec[sample_i] = OptimalPoissonLoss(cumsum_value, seg1_mean);
  }
  /*
    The for loops below implement the GridSearch() function mentioned
    on line 2 of the JointZoom algorithm in the PeakSegJoint paper.
  */
  double best_loss = INFINITY, feasible_loss;
  for(int seg1_LastIndex=0; seg1_LastIndex < n_bins-2; seg1_LastIndex++){
    peakStart = seg1_chromStart + (seg1_LastIndex+1)*bases_per_bin;
    if(unfilled_chromStart < peakStart){
      for(int sample_i=0; sample_i < n_samples; sample_i++){
	cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
	cumsum_value = cumsum_vec[seg1_LastIndex];
	bin_bases = (seg1_LastIndex+1)*bases_per_bin;
	data_bases = bin_bases - (double)extra_before;
	/* printf("sample_i=%d extra_before=%d bin_bases=%f data_bases=%f\n", */
	/* 	     sample_i, extra_before, bin_bases, data_bases); */
	seg1_mean = cumsum_value/data_bases;
	seg1_mean_vec[sample_i] = seg1_mean;
	seg1_loss_vec[sample_i] = OptimalPoissonLoss(cumsum_value, seg1_mean);
      }
      for(int seg2_LastIndex=seg1_LastIndex+1; 
	  seg2_LastIndex < n_bins-1; 
	  seg2_LastIndex++){
	feasible_loss=0.0;
	peakEnd = seg1_chromStart + (seg2_LastIndex+1)*bases_per_bin;
	if(peakEnd < unfilled_chromEnd){
	  for(int sample_i=0; sample_i < n_samples; sample_i++){
	    candidate_loss_vec[sample_i] = seg1_loss_vec[sample_i];
	    cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
	    //segment 2.
	    cumsum_value = cumsum_vec[seg2_LastIndex]-cumsum_vec[seg1_LastIndex];
	    data_bases = (seg2_LastIndex-seg1_LastIndex)*bases_per_bin;
	    seg2_mean = cumsum_value/data_bases;
	    seg2_mean_vec[sample_i] = seg2_mean;
	    candidate_loss_vec[sample_i] += OptimalPoissonLoss(cumsum_value, seg2_mean);
	    //segment 3.
	    cumsum_value = cumsum_vec[n_bins-1]-cumsum_vec[seg2_LastIndex];
	    bin_bases = (n_bins-1-seg2_LastIndex)*bases_per_bin;
	    data_bases = bin_bases - extra_after;
	    seg3_mean = cumsum_value/data_bases;
	    seg3_mean_vec[sample_i] = seg3_mean;
	    candidate_loss_vec[sample_i] += OptimalPoissonLoss(cumsum_value, seg3_mean);
	    if(seg1_mean < seg2_mean && seg2_mean > seg3_mean){
	      feasible_loss += candidate_loss_vec[sample_i];
	    }else{
	      feasible_loss += flat_loss_vec[sample_i];
	    }
	  }//sample_i
	  if(feasible_loss < best_loss){
	    best_loss = feasible_loss;
	    peak_start_end[0] = peakStart;
	    peak_start_end[1] = peakEnd;
	    for(int sample_i=0; sample_i < n_samples; sample_i++){
	      peak_loss_vec[sample_i] = candidate_loss_vec[sample_i];
	      mean_mat[sample_i] = seg1_mean_vec[sample_i];
	      mean_mat[sample_i+n_samples] = seg2_mean_vec[sample_i];
	      mean_mat[sample_i+n_samples*2] = seg3_mean_vec[sample_i];
	      // Also save the left/right cumsums, which are needed to
	      // perform the step2 optimization.
	      cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
	      if(seg1_LastIndex == 0){
		left_initial_cumsum_vec[sample_i] = 0;
	      }else{
		left_initial_cumsum_vec[sample_i] = cumsum_vec[seg1_LastIndex-1];
	      }
	      right_initial_cumsum_vec[sample_i] = cumsum_vec[seg2_LastIndex-1];
	    }//for(sample_i
	  }//if(feasible_loss < best_loss
	}//if(peakEnd < unfilled_chromEnd
      }//for(seg2_LastIndex
    }//if(unfilled_chromStart < peakStart
  }//for(seg1_LastIndex
  //return status;//old end of Step1.
  /* When performing the minimization over peakStart/End locations, it
   * is possible that at any given bases_per_bin value, there is no
   * better solution than what we found for the previous bases_per_bin
   * value. In that case, we begin the search anew at a lower
   * resolution, but we need to copy the cumsums from the following
   * index of left_right_vec: */
  int left_chromStart, right_chromStart;
  /*
    The while loop below corresponds to line 3 of the JointZoom
    algorithm from the PeakSegJoint paper.
  */
  while(1 < bases_per_bin){
    best_seg1 = -1; // indicates no min found.
    left_chromStart = peak_start_end[0] - bases_per_bin;
    right_chromStart = peak_start_end[1] - bases_per_bin;
    /*
      Below in the C code we decrease the value of bases_per_bin,
      as in line 4 of the JointZoom algorithm in the PeakSegJoint
      paper.
    */
    bases_per_bin /= bin_factor;
    //printf("bases_per_bin=%d left cumsum before:\n", bases_per_bin);
    for(int sample_i=0; sample_i < n_samples; sample_i++){
      profile = profile_list->profile_vec + sample_i;
      //printf("binSumLR sample_i=%d\n", sample_i);
      status = binSumLR(data_start_end,
			profile->chromStart, profile->chromEnd,
			profile->coverage, profile->n_entries,
			left_bin_vec, right_bin_vec,
			left_chromStart, right_chromStart,
			bases_per_bin, n_bins);
      if(status != 0){
	//printf("binSumLR bad status\n");
	free(sample_count_mat);
	free(last_cumsum_vec);
	free(sample_cumsum_mat);
  
	free(seg1_mean_vec);
	free(seg2_mean_vec);
	free(seg3_mean_vec);

	free(seg1_loss_vec);
	free(candidate_loss_vec);

	free(left_bin_vec);
	free(right_bin_vec);
	free(left_cumsum_mat);
	free(right_cumsum_mat);
	free(left_initial_cumsum_vec);
	free(right_initial_cumsum_vec);
	return status;
      }
      left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
      left_cumsum_value = left_initial_cumsum_vec[sample_i];
      right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
      right_cumsum_value = right_initial_cumsum_vec[sample_i];
      //printf("%d ", left_cumsum_value);
      //printf("%d ", right_cumsum_value);
      for(int bin_i=0; bin_i < n_bins; bin_i++){
	left_cumsum_value += left_bin_vec[bin_i];
	left_cumsum_vec[bin_i] = left_cumsum_value;
	right_cumsum_value += right_bin_vec[bin_i];
	right_cumsum_vec[bin_i] = right_cumsum_value;
      }
    }//for(sample_i
    /* 
       cumsum matrices have been computed, so now use them to
       compute the loss and feasibility of all models.

       The for loops below correspond to SearchNearPeak() on line
       5 of the JointZoom algorithm in the PeakSegJoint paper.
    */
    for(int seg1_LastIndex=0; seg1_LastIndex < n_bins; seg1_LastIndex++){
      peakStart = left_chromStart + (seg1_LastIndex+1)*bases_per_bin;
      if(unfilled_chromStart < peakStart){
	//printf("[seg1last=%d] seg1 cumsum bases ", seg1_LastIndex);
	for(int sample_i=0; sample_i < n_samples; sample_i++){
	  left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	  cumsum_value = left_cumsum_vec[seg1_LastIndex];
	  bin_bases = peakStart - seg1_chromStart;
	  data_bases = bin_bases - extra_before;
	  seg1_mean = cumsum_value/data_bases;
	  //printf("%d %f ", cumsum_value, bases_value);
	  seg1_mean_vec[sample_i] = seg1_mean;
	  seg1_loss_vec[sample_i] = OptimalPoissonLoss(cumsum_value, seg1_mean);
	}
	//printf("\n");
	for(int seg2_LastIndex=0; seg2_LastIndex < n_bins; seg2_LastIndex++){
	  peakEnd = right_chromStart + (seg2_LastIndex+1)*bases_per_bin;
	  if(peakStart < peakEnd && peakEnd < unfilled_chromEnd){
	    feasible_loss = 0.0;
	    //printf("[seg2last=%d]\n", seg2_LastIndex);
	    for(int sample_i=0; sample_i < n_samples; sample_i++){
	      left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	      right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	      //segment 1.
	      candidate_loss_vec[sample_i] = seg1_loss_vec[sample_i];
	      //segment 2.
	      cumsum_value = 
		right_cumsum_vec[seg2_LastIndex]-
		left_cumsum_vec[seg1_LastIndex];
	      data_bases = peakEnd-peakStart;
	      seg2_mean = cumsum_value/data_bases;
	      seg2_mean_vec[sample_i] = seg2_mean;
	      candidate_loss_vec[sample_i] += OptimalPoissonLoss(cumsum_value, seg2_mean);
	      //segment 3.
	      cumsum_value = 
		last_cumsum_vec[sample_i]-
		right_cumsum_vec[seg2_LastIndex];
	      bin_bases = seg3_chromEnd - peakEnd;
	      data_bases = bin_bases - extra_after;
	      seg3_mean = cumsum_value/data_bases;
	      seg3_mean_vec[sample_i] = seg3_mean;
	      candidate_loss_vec[sample_i] += OptimalPoissonLoss(cumsum_value, seg3_mean);
	      if(seg1_mean < seg2_mean && seg2_mean > seg3_mean){
		feasible_loss += candidate_loss_vec[sample_i];
	      }else{
		feasible_loss += flat_loss_vec[sample_i];
	      }
	    }//for(sample_i
	    if(feasible_loss < best_loss){
	      best_loss = feasible_loss;
	      peak_start_end[0] = peakStart;
	      peak_start_end[1] = peakEnd;
	      best_seg1 = seg1_LastIndex;
	      best_seg2 = seg2_LastIndex;
	      for(int sample_i=0; sample_i < n_samples; sample_i++){
		peak_loss_vec[sample_i] = candidate_loss_vec[sample_i];
		mean_mat[sample_i] = seg1_mean_vec[sample_i];
		mean_mat[sample_i+n_samples] = seg2_mean_vec[sample_i];
		mean_mat[sample_i+n_samples*2] = seg3_mean_vec[sample_i];
	      }
	    }//if(feasible_loss<best_loss
	  }//if(peakStart<peakEnd
	}//for(seg2_LastIndex
      }//if(unfilled_chromStart < peakStart
    }//for(seg1_LastIndex
    for(int sample_i=0; sample_i < n_samples; sample_i++){
      left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
      right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
      if(best_seg1 == -1){
	/* the bin_factor is the number of new bins that get put into old bins.
	   [ ] old bin   
	   _ _ new bins
	   bin_factor=2 [ _ _ ]   peakStart [ _ _ ]
	   bin_factor=3 [ _ _ _ ] peakStart [ _ _ _ ]
     
	   We always start the search on the new level from the old bin just
	   before the peakStart/End. Therefore if no new min is found, the
	   best choice is the old one, seg1_LastIndex=bin_factor-1 */
	best_seg1 = bin_factor-1;
	best_seg2 = bin_factor-1;
	//Rprintf("no new min found! best_seg1=%d best_seg2=%d\n", best_seg1, best_seg2);
      }
      if(best_seg1 != 0){
	left_initial_cumsum_vec[sample_i] = left_cumsum_vec[best_seg1-1];
      }
      if(best_seg2 != 0){
	right_initial_cumsum_vec[sample_i] = right_cumsum_vec[best_seg2-1];
      }
    }
  }//while(1 < bases_per_bin)
  //printf("free at end\n");
  free(sample_count_mat);
  free(last_cumsum_vec);
  free(sample_cumsum_mat);
  
  free(seg1_mean_vec);
  free(seg2_mean_vec);
  free(seg3_mean_vec);

  free(seg1_loss_vec);
  free(candidate_loss_vec);

  free(left_bin_vec);
  free(right_bin_vec);
  free(left_cumsum_mat);
  free(right_cumsum_mat);
  free(left_initial_cumsum_vec);
  free(right_initial_cumsum_vec);
  
  return 0;
}

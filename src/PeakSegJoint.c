/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegJoint.h"
#include "DPA.h"
#include "OptimalPoissonLoss.h"
#include "binSum.h"
#include <stdio.h>
#include <stdlib.h>

int LossIndex_compare(const void *a, const void *b){
  const struct LossIndex *A = a, *B = b;
  //printf("compare %f %f\n", A->loss, B->loss);
  return (int)(A->loss - B->loss);
}

int PeakSegJointHeuristicStep1(
  struct ProfileList *profile_list,
  int bin_factor,
  struct PeakSegJointModelList *model_list
  ){
  int n_samples = profile_list->n_profiles;
  int sample_i, coverage_i, 
    chromStart, chromEnd, unfilled_chromStart, unfilled_chromEnd;
  struct Profile *profile, *samples = profile_list->profile_vec;
  struct PeakSegJointModel *model;
  profile = samples;
  unfilled_chromEnd = get_max_chromEnd(profile);
  unfilled_chromStart = get_min_chromStart(profile);
  for(sample_i=1; sample_i < n_samples; sample_i++){
    profile = samples + sample_i;
    chromStart = get_min_chromStart(profile);
    if(unfilled_chromStart < chromStart){
      unfilled_chromStart = chromStart;
    }
    chromEnd = get_max_chromEnd(profile);
    if(chromEnd < unfilled_chromEnd){
      unfilled_chromEnd = chromEnd;
    }
  }
  model_list->data_start_end[0] = unfilled_chromStart;
  model_list->data_start_end[1] = unfilled_chromEnd;
  //printf("data_start_end=[%d,%d]\n", unfilled_chromStart, unfilled_chromEnd);
  int unfilled_bases = unfilled_chromEnd - unfilled_chromStart;
  double bases_value, seg1_loss_value;
  if(unfilled_bases/bin_factor < 4){
    /*
      4 is smallest the number of data points for which the 3-segment
      optimization problem is not trivial.

      If we don't have at least this many data points for the first
      bin step, than we stop with an error.
    */
    return ERROR_BIN_FACTOR_TOO_LARGE;
  }
  int bases_per_bin = 1;
  while(unfilled_bases/bases_per_bin/bin_factor >= 4){
    bases_per_bin *= bin_factor;
  }
  int n_bins = unfilled_bases / bases_per_bin;
  if(unfilled_bases % bases_per_bin != 0){
    n_bins ++ ;
  }
  model_list->n_bins[0] = n_bins;
  model_list->bases_per_bin[0] = bases_per_bin;
  model_list->bin_factor[0] = bin_factor;
  int extra_bases = n_bins  * bases_per_bin - unfilled_bases;
  int extra_before = extra_bases/2;
  int extra_after = extra_bases - extra_before;
  int seg1_chromStart = unfilled_chromStart - extra_before;
  int seg3_chromEnd = unfilled_chromEnd + extra_after;
  int bases = seg3_chromEnd - seg1_chromStart;
  model_list->seg_start_end[0] = seg1_chromStart;
  model_list->seg_start_end[1] = seg3_chromEnd;
  //printf("seg_start_end=[%d,%d]\n", seg1_chromStart, seg3_chromEnd);
  // sample_*_mat variables are matrices n_bins x n_samples (in
  // contrast to model_*_mat which are n_bins x n_segments=3).
  int *sample_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *count_vec, *cumsum_vec, cumsum_value;
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples + sample_i;
    count_vec = sample_count_mat + n_bins*sample_i;
    //printf("initial sample_i=%d\n", sample_i);
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    count_vec,
		    bases_per_bin, n_bins, seg1_chromStart, 
		    ERROR_EMPTY_BIN);
    if(status != 0){
      free(sample_count_mat);
      return status;
    }
    /* Profiles may not have the same first chromStart and last
     * chromEnd, so subtract any values that are before or after the
     * maximum region that overlaps all profiles.
     */

    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    &extra_before,
		    unfilled_chromStart - seg1_chromStart,
		    1,
		    seg1_chromStart,
		    EMPTY_AS_ZERO);
    if(status != 0){
      free(sample_count_mat);
      return status;
    }
    count_vec[0] -= extra_before;
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    &extra_after,
		    seg3_chromEnd - unfilled_chromEnd,
		    1,
		    unfilled_chromEnd,
		    EMPTY_AS_ZERO);
    if(status != 0){
      free(sample_count_mat);
      return status;
    }
    count_vec[n_bins - 1] -= extra_after;
  }//for sample_i
  int bin_i, offset;
  double mean_value, loss_value;
  double flat_loss_total = 0.0;
  int *sample_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  struct LossIndex *diff_index_vec = 
    (struct LossIndex *)malloc(sizeof(struct LossIndex)*n_samples);
  int n_feasible;
  for(sample_i=0; sample_i < n_samples; sample_i++){
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
    model_list->last_cumsum_vec[sample_i] = cumsum_value;
    //printf("\n");
    bases_value = (double)bases;
    mean_value = cumsum_value / bases_value;
    //printf("sample_i=%d cumsum=%d bases=%d\n", sample_i, cumsum_value, bases);
    model_list->sample_mean_vec[sample_i] = mean_value;
    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
    model_list->flat_loss_vec[sample_i] = loss_value;
    flat_loss_total += loss_value;
  }
  model_list->model_vec[0].loss[0] = flat_loss_total;

  int seg1_LastIndex, seg2_LastIndex;
  int model_i, diff_i, n_peaks;
  double *seg1_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg2_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg3_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *peak_loss_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg1_loss_vec = (double*)malloc(sizeof(double)*n_samples);
  for(seg1_LastIndex=0; seg1_LastIndex < n_bins-2; seg1_LastIndex++){
    for(sample_i=0; sample_i < n_samples; sample_i++){
      cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
      cumsum_value = cumsum_vec[seg1_LastIndex];
      bases_value = (seg1_LastIndex+1)*bases_per_bin;
      mean_value = cumsum_value/bases_value;
      seg1_mean_vec[sample_i] = mean_value;
      loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
      seg1_loss_vec[sample_i] = loss_value;
    }
    for(seg2_LastIndex=seg1_LastIndex+1; 
	seg2_LastIndex < n_bins-1; 
	seg2_LastIndex++){
      n_feasible=0;
      for(sample_i=0; sample_i < n_samples; sample_i++){
	peak_loss_vec[sample_i] = seg1_loss_vec[sample_i];
	cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
	//segment 2.
	cumsum_value = cumsum_vec[seg2_LastIndex]-cumsum_vec[seg1_LastIndex];
	bases_value = (seg2_LastIndex-seg1_LastIndex)*bases_per_bin;
	mean_value = cumsum_value/bases_value;
	seg2_mean_vec[sample_i] = mean_value;
	loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	peak_loss_vec[sample_i] += loss_value;
	//segment 3.
	cumsum_value = cumsum_vec[n_bins-1]-cumsum_vec[seg2_LastIndex];
	bases_value = (n_bins-1-seg2_LastIndex)*bases_per_bin;
	mean_value = cumsum_value/bases_value;
	seg3_mean_vec[sample_i] = mean_value;
	loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	peak_loss_vec[sample_i] += loss_value;
	//if feasible, add to list of loss differences.
	if(seg1_mean_vec[sample_i] < seg2_mean_vec[sample_i] &&
	   seg3_mean_vec[sample_i] < seg2_mean_vec[sample_i]){
	  diff_index_vec[n_feasible].sample_i = sample_i;
	  diff_index_vec[n_feasible].loss = 
	    peak_loss_vec[sample_i]-model_list->flat_loss_vec[sample_i];
	  n_feasible++;
	}
      }//sample_i
      if(0 < n_feasible){
	/* printf("[0,%d][%d,%d][%d,%d] %d feasible\n", */
	/* 	     seg1_LastIndex, */
	/* 	     seg1_LastIndex+1, */
	/* 	     seg2_LastIndex, */
	/* 	     seg2_LastIndex+1, */
	/* 	     n_bins-1, */
	/* 	     n_feasible); */
	/* printf("before sort"); */
	/* for(sample_i=0; sample_i < n_feasible; sample_i++){ */
	/* 	printf(" %f", diff_index_vec[sample_i].loss); */
	/* } */
	/* printf("\n"); */
	qsort(diff_index_vec, n_feasible, sizeof(struct LossIndex), 
	      LossIndex_compare);
	/* printf("after sort"); */
	/* for(sample_i=0; sample_i < n_feasible; sample_i++){ */
	/* 	printf(" %f", diff_index_vec[sample_i].loss); */
	/* } */
	/* printf("\n"); */
	for(model_i=0; model_i < n_feasible; model_i++){
	  // start from loss of all samples with 1 segment.
	  loss_value = flat_loss_total;
	  n_peaks = model_i + 1;
	  for(diff_i=0; diff_i < n_peaks; diff_i++){
	    sample_i = diff_index_vec[diff_i].sample_i;
	    // subtract the loss of this sample with 1 segment.
	    loss_value -= model_list->flat_loss_vec[sample_i];
	    // add loss from this sample with 3 segments (1 peak).
	    loss_value += peak_loss_vec[sample_i];
	  }
	  model = model_list->model_vec + n_peaks; 
	  if(loss_value < model->loss[0]){
	    model->loss[0] = loss_value;
	    model->peak_start_end[0] = 
	      seg1_chromStart + (seg1_LastIndex+1)*bases_per_bin;
	    model->peak_start_end[1] = 
	      seg1_chromStart + (seg2_LastIndex+1)*bases_per_bin;
	    for(diff_i=0; diff_i < n_peaks; diff_i++){
	      sample_i = diff_index_vec[diff_i].sample_i;
	      cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
	      model->samples_with_peaks_vec[diff_i] = sample_i;
	      if(seg1_LastIndex == 0){
		model->left_cumsum_vec[diff_i] = 0;
	      }else{
		model->left_cumsum_vec[diff_i] = cumsum_vec[seg1_LastIndex-1];
	      }
	      model->right_cumsum_vec[diff_i] = cumsum_vec[seg2_LastIndex-1];
	      model->seg1_mean_vec[diff_i] = seg1_mean_vec[sample_i];
	      model->seg2_mean_vec[diff_i] = seg2_mean_vec[sample_i];
	      model->seg3_mean_vec[diff_i] = seg3_mean_vec[sample_i];
	    }
	  }
	}//model_i
      }//if(n_feasible)
    }//seg2_LastIndex
  }//seg2_FirstIndex
  free(sample_cumsum_mat);
  free(sample_count_mat);
  free(peak_loss_vec);
  free(seg1_loss_vec);
  free(seg1_mean_vec);
  free(seg2_mean_vec);
  free(seg3_mean_vec);
  free(diff_index_vec);
  return status;
}

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

int 
PeakSegJointHeuristicStep2
(struct ProfileList *profile_list,
 struct PeakSegJointModelList *model_list
  ){
  int n_peaks, diff_i, sample_i, bin_i;
  int n_bins = model_list->bin_factor[0] * 2;
  int n_samples = model_list->n_models - 1;
  struct PeakSegJointModel *model;
  struct Profile *profile;
  int bases_per_bin;
  int left_chromStart, right_chromStart;
  //printf("before malloc n_bins=%d n_samples=%d\n", n_bins, n_samples);
  int *left_bin_vec = (int*) malloc(n_bins * sizeof(int));
  int *right_bin_vec = (int*) malloc(n_bins * sizeof(int));
  int *left_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *right_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  double *seg1_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg2_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg3_mean_vec = (double*)malloc(sizeof(double)*n_samples);
  double *seg1_loss_vec = (double*)malloc(sizeof(double)*n_samples);
  //printf("after malloc\n");
  double total_loss, loss_value, bases_value, mean_value;
  int *left_cumsum_vec, *right_cumsum_vec;
  int seg1_LastIndex, seg2_LastIndex;
  int status;
  int left_cumsum_value, right_cumsum_value, cumsum_value;
  int peakStart, peakEnd;
  int min_found;
  for(n_peaks=1; n_peaks < model_list->n_models; n_peaks++){
    model = model_list->model_vec + n_peaks;
    if(model->loss[0] < INFINITY){
      bases_per_bin = model_list->bases_per_bin[0];
      while(1 < bases_per_bin){
	min_found = 0;
	left_chromStart = model->peak_start_end[0] - bases_per_bin;
	right_chromStart = model->peak_start_end[1] - bases_per_bin;
	bases_per_bin /= model_list->bin_factor[0];
	for(diff_i=0; diff_i < n_peaks; diff_i++){
	  sample_i = model->samples_with_peaks_vec[diff_i];
	  profile = profile_list->profile_vec + sample_i;
	  //printf("binSumLR sample_i=%d\n", sample_i);
	  status = binSumLR(model_list->data_start_end,
			    profile->chromStart, profile->chromEnd,
			    profile->coverage, profile->n_entries,
			    left_bin_vec, right_bin_vec,
			    left_chromStart, right_chromStart,
			    bases_per_bin, n_bins);
	  if(status != 0){
	    //printf("binSumLR bad status\n");
	    free(left_bin_vec);
	    free(right_bin_vec);
	    free(left_cumsum_mat);
	    free(right_cumsum_mat);
	    free(seg1_mean_vec);
	    free(seg2_mean_vec);
	    free(seg3_mean_vec);
	    free(seg1_loss_vec);
	    return status;
	  }
	  left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	  left_cumsum_value = model->left_cumsum_vec[diff_i];
	  right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	  right_cumsum_value = model->right_cumsum_vec[diff_i];
	  for(bin_i=0; bin_i < n_bins; bin_i++){
	    left_cumsum_value += left_bin_vec[bin_i];
	    left_cumsum_vec[bin_i] = left_cumsum_value;
	    right_cumsum_value += right_bin_vec[bin_i];
	    right_cumsum_vec[bin_i] = right_cumsum_value;
	  }
	}//for(diff_i

	/* printf("left bases_per_bin=%d n_peaks=%d\n", bases_per_bin, n_peaks); */
	/* for(diff_i=0;diff_i<n_peaks;diff_i++){ */
	/* 	sample_i = model->samples_with_peaks_vec[diff_i]; */
	/* 	left_cumsum_vec = left_cumsum_mat + n_bins*sample_i; */
	/* 	for(bin_i=0;bin_i<n_bins;bin_i++){ */
	/* 	  printf("%d ", left_cumsum_vec[bin_i]); */
	/* 	} */
	/* 	printf("\n"); */
	/* } */
	/* printf("\n"); */

	/* 
	   cumsum matrices have been computed, so now use them to
	   compute the loss and feasibility of all models.
	*/
	for(seg1_LastIndex=0; seg1_LastIndex < n_bins; seg1_LastIndex++){
	  peakStart = left_chromStart + (seg1_LastIndex+1)*bases_per_bin;
	  //printf("[seg1last=%d] seg1 cumsum bases ", seg1_LastIndex);
	  for(diff_i=0; diff_i < n_peaks; diff_i++){
	    sample_i = model->samples_with_peaks_vec[diff_i];
	    left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	    cumsum_value = left_cumsum_vec[seg1_LastIndex];
	    bases_value = peakStart - model_list->seg_start_end[0];
	    mean_value = cumsum_value/bases_value;
	    //printf("%d %f ", cumsum_value, bases_value);
	    seg1_mean_vec[sample_i] = mean_value;
	    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	    seg1_loss_vec[sample_i] = loss_value;
	  }
	  //printf("\n");
	  for(seg2_LastIndex=0; seg2_LastIndex < n_bins; seg2_LastIndex++){
	    peakEnd = right_chromStart + (seg2_LastIndex+1)*bases_per_bin;
	    /* printf("\npeaks=%d[%d,%d]bases_per_bin=%d\n",  */
	    /* 	 n_peaks, peakStart, peakEnd, bases_per_bin); */
	    total_loss = model_list->model_vec[0].loss[0];
	    if(peakEnd <= peakStart){
	      total_loss = INFINITY;
	    }
	    //printf("[seg2last=%d]\n", seg2_LastIndex);
	    for(diff_i=0; diff_i < n_peaks; diff_i++){
	      sample_i = model->samples_with_peaks_vec[diff_i];
	      left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	      right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	      total_loss -= model_list->flat_loss_vec[sample_i];
	      //segment 1.
	      total_loss += seg1_loss_vec[sample_i];
	      //segment 2.
	      cumsum_value = 
		right_cumsum_vec[seg2_LastIndex]-left_cumsum_vec[seg1_LastIndex];
	      bases_value = peakEnd-peakStart;
	      mean_value = cumsum_value/bases_value;
	      /* printf("[sample=%d][seg=2] %d %f", */
	      /* 	   sample_i, */
	      /* 	   cumsum_value, */
	      /* 	   bases_value); */
	      seg2_mean_vec[sample_i] = mean_value;
	      loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	      total_loss += loss_value;
	      //segment 3.
	      cumsum_value = 
		model_list->last_cumsum_vec[sample_i]-
		right_cumsum_vec[seg2_LastIndex];
	      bases_value = model_list->seg_start_end[1] - peakEnd;
	      /* printf("[sample=%d][seg=3] %d %f", */
	      /* 	   sample_i, */
	      /* 	   cumsum_value, */
	      /* 	   bases_value); */
	      mean_value = cumsum_value/bases_value;
	      seg3_mean_vec[sample_i] = mean_value;
	      loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	      total_loss += loss_value;
	      /* printf("[sample=%d] %f %f %f\n", */
	      /* 	   sample_i, */
	      /* 	   seg1_mean_vec[sample_i], */
	      /* 	   seg2_mean_vec[sample_i], */
	      /* 	   seg3_mean_vec[sample_i]); */
	      //if not feasible, loss is infinite.
	      if(seg2_mean_vec[sample_i] <= seg1_mean_vec[sample_i] ||
		 seg2_mean_vec[sample_i] <= seg3_mean_vec[sample_i]){
		total_loss = INFINITY;
	      }
	    }	  
	    //printf("loss=%f\n", total_loss);
	    if(total_loss < model->loss[0]){
	      min_found = 1;
	      model->loss[0] = total_loss;
	      model->peak_start_end[0] = peakStart;
	      model->peak_start_end[1] = peakEnd;
	      /* printf("new best loss=%f [%d,%d]\n", */
	      /* 	   total_loss, seg1_LastIndex, seg2_LastIndex); */
	      for(diff_i=0; diff_i < n_peaks; diff_i++){
		sample_i = model->samples_with_peaks_vec[diff_i];
		if(seg1_LastIndex != 0){
		  left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
		  model->left_cumsum_vec[diff_i] = 
		    left_cumsum_vec[seg1_LastIndex-1];
		}
		if(seg2_LastIndex != 0){
		  right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
		  model->right_cumsum_vec[diff_i] = 
		    right_cumsum_vec[seg2_LastIndex-1];
		}
		model->seg1_mean_vec[diff_i] = seg1_mean_vec[sample_i];
		model->seg2_mean_vec[diff_i] = seg2_mean_vec[sample_i];
		model->seg3_mean_vec[diff_i] = seg3_mean_vec[sample_i];
	      }//diff_i
	    }//total_loss
	  }//seg2_LastIndex
	}//seg1_LastIndex
	/* printf("n_peaks=%d bases_per_bin=%d [%d,%d] loss=%f\n", */
	/* 	     n_peaks, bases_per_bin,  */
	/* 	     model->peak_start_end[0], model->peak_start_end[1], */
	/* 	     model->loss[0]); */
	if(min_found == 0){
	  //printf("no min found\n");
	  for(diff_i=0; diff_i < n_peaks; diff_i++){
	    sample_i = model->samples_with_peaks_vec[diff_i];
	    left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	    model->left_cumsum_vec[diff_i] = left_cumsum_vec[0];
	    right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	    model->right_cumsum_vec[diff_i] = right_cumsum_vec[0];
	  }//diff_i
	}//min_found
      }//while(1 < bases_per_bin)
    }//if(loss < INFINITY
  }//for(n_peaks
  //printf("free at end\n");
  free(left_bin_vec);
  free(right_bin_vec);
  free(left_cumsum_mat);
  free(right_cumsum_mat);
  free(seg1_mean_vec);
  free(seg2_mean_vec);
  free(seg3_mean_vec);
  free(seg1_loss_vec);
  return 0;
}

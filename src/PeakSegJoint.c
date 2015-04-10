/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegJoint.h"
#include "DPA.h"
#include "OptimalPoissonLoss.h"
#include "binSum.h"
#include <stdio.h>
#include <stdlib.h>

struct PeakSegJointModelList * 
malloc_model_list(
  int n_models
  ){
  int model_i;
  struct PeakSegJointModel *model;
  struct PeakSegJointModelList *model_list = 
    malloc(sizeof(struct PeakSegJointModelList));
  model_list->n_models = n_models;
  model_list->seg_start_end = malloc(2*sizeof(int));
  model_list->model_vec = 
    malloc(n_models * sizeof(struct PeakSegJointModel *));
  for(model_i=0; model_i < n_models; model_i++){
    model = malloc(sizeof(struct PeakSegJointModel));
    if(0 < model_i){
      model->loss = malloc(sizeof(double));
      model->peak_start_end = malloc(2*sizeof(int));
      model->samples_with_peaks_vec = malloc(model_i * sizeof(int));
      model->left_cumsum_vec = malloc(model_i * sizeof(int));
      model->right_cumsum_vec = malloc(model_i * sizeof(int));
      model->last_cumsum_vec = malloc(model_i * sizeof(int));
      model->seg1_mean_vec = malloc(model_i * sizeof(double));
      model->seg2_mean_vec = malloc(model_i * sizeof(double));
      model->seg3_mean_vec = malloc(model_i * sizeof(double));
    }
    model->loss[0] = INFINITY;
    model_list->model_vec[model_i] = model;
  }
  return model_list;
}

void free_model_list(struct PeakSegJointModelList *model_list){
  int model_i;
  struct PeakSegJointModel *model;
  for(model_i = 0; model_i < model_list->n_models; model_i++){
    model = model_list->model_vec[model_i];
    if(0 < model_i){
      free(model->peak_start_end);
      free(model->loss);
      free(model->samples_with_peaks_vec);
      free(model->left_cumsum_vec);
      free(model->right_cumsum_vec);
      free(model->last_cumsum_vec);
      free(model->seg1_mean_vec);
      free(model->seg2_mean_vec);
      free(model->seg3_mean_vec);
    }
    free(model);
  }
  free(model_list->model_vec);
  free(model_list->seg_start_end);
  free(model_list);
}

int PeakSegJointHeuristicStep1(
  struct Profile **samples,
  int n_samples,
  int bin_factor,
  struct PeakSegJointModelList *model_list
  ){
  int sample_i, coverage_i, min_chromEnd, max_chromStart, 
    chromStart, chromEnd, unfilled_chromStart, unfilled_chromEnd;
  struct Profile *profile;
  profile = samples[0];
  unfilled_chromEnd = get_max_chromEnd(profile);
  unfilled_chromStart = get_min_chromStart(profile);
  for(sample_i=1; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    chromStart = get_min_chromStart(profile);
    if(unfilled_chromStart < chromStart){
      unfilled_chromStart = chromStart;
    }
    chromEnd = get_max_chromEnd(profile);
    if(chromEnd < unfilled_chromEnd){
      unfilled_chromEnd = chromEnd;
    }
  }
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
  int extra_bases = n_bins  * bases_per_bin - unfilled_bases;
  int extra_before = extra_bases/2;
  int extra_after = extra_bases - extra_before;
  max_chromStart = unfilled_chromStart - extra_before;
  min_chromEnd = unfilled_chromEnd + extra_after;
  model_list->seg_start_end[0] = max_chromStart;
  model_list->seg_start_end[1] = min_chromEnd;
  return 0;
  // sample_*_mat variables are matrices n_bins x n_samples (in
  // contrast to model_*_mat which are n_bins x n_segments=3).
  int *sample_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    sample_count_mat + n_bins*sample_i,
		    bases_per_bin, n_bins, max_chromStart, 
		    ERROR_EMPTY_BIN);
  }//for sample_i
  if(status != 0){
    /* printf("first sample_i=%d bases_per_bin=%d n_bins=%d\n",  */
    /* 	     sample_i, bases_per_bin, n_bins); */
    free(sample_count_mat);
    return status;
  }
  int bin_i, offset;
  int *count_vec, *cumsum_vec, cumsum_value;
  double *loss_vec, mean_value, loss_value;
  int *sample_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
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
    //printf("\n");
  }
  double *seg1_loss_vec = (double*) malloc(n_bins * sizeof(double));
  for(bin_i=0; bin_i < n_bins; bin_i++){
    seg1_loss_vec[bin_i] = 0;
    for(sample_i=0; sample_i < n_samples; sample_i++){
      cumsum_value = sample_cumsum_mat[n_bins*sample_i+bin_i];
      /*
      printf("bin=%d sample=%d cumsum=%d\n", 
	     bin_i, sample_i, cumsum_value);
      */
      bases_value = (bin_i+1) * bases_per_bin;
      mean_value = ((double) cumsum_value) / bases_value;
      loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
      seg1_loss_vec[bin_i] += loss_value;
    }
    //printf("seg1_loss[%d]=%f\n", bin_i, seg1_loss_vec[bin_i]);
  }
  /* 
     Second step of DPA: compute optimal loss in 2 segments up to data
     point t, for all data points = bins.
   */
  double *seg12_loss_vec = (double*) malloc(n_bins * sizeof(double));
  int *seg2_first_vec = (int*) malloc(n_bins * sizeof(int));
  int seg2_FirstIndex, seg2_LastIndex, best_FirstIndex;
  double seg1_loss, min_loss, seg2_loss, candidate_loss, seg3_loss;
  for(seg2_LastIndex=1; 
      seg2_LastIndex <= n_bins-2; 
      seg2_LastIndex++){
    //printf("seg2_LastIndex=%d\n", seg2_LastIndex);
    get_best_FirstIndex(
      seg1_loss_vec,
      n_bins,
      sample_cumsum_mat,
      n_samples,
      1, // first_possible_index
      seg2_LastIndex, // last_possible_index
      bases_per_bin,
      seg2_first_vec + seg2_LastIndex,
      seg12_loss_vec + seg2_LastIndex);
  }
  /*
    For the best segmentation in 3 segments up to n_bins-1, the first
    index of the 3rd segment is computed using the following code. No
    need for a for loop on seg3_LastIndex, since we only want to know
    the optimal model which ends at the last data point = bin (we are
    not continuing the DPA past 3 segments).
  */
  int seg3_FirstIndex;
  double best_loss;
  get_best_FirstIndex(
    seg12_loss_vec,
    n_bins,
    sample_cumsum_mat,
    n_samples,
    2, // first_possible_index
    n_bins-1, // last_possible_index
    bases_per_bin,
    &seg3_FirstIndex,
    &best_loss);
  seg2_LastIndex = seg3_FirstIndex-1;
  seg2_FirstIndex = seg2_first_vec[seg2_LastIndex];
  int seg1_LastIndex = seg2_FirstIndex-1;

  free(seg1_loss_vec);
  free(seg12_loss_vec);
  free(seg2_first_vec);

  /* printf("[0,%d] %d[%d,%d]%d [%d,%d]\n", */
  /* 	 seg1_LastIndex, */
  /* 	 optimal_start_end[0], */
  /* 	 seg2_FirstIndex, seg2_LastIndex, */
  /* 	 optimal_start_end[1], */
  /* 	 seg3_FirstIndex, n_bins-1); */

  int last_chromEnd = n_bins * bases_per_bin + max_chromStart;
  int *left_cumsum_vec = (int*) malloc(n_samples * sizeof(int));
  int *right_cumsum_vec = (int*) malloc(n_samples * sizeof(int));
  int *last_cumsum_vec = (int*) malloc(n_samples * sizeof(int));
  for(sample_i=0; sample_i < n_samples; sample_i++){
    if(seg1_LastIndex > 0){
      left_cumsum_vec[sample_i] = 
	sample_cumsum_mat[n_bins*sample_i+seg1_LastIndex-1];
    }else{ // there is no data before the first data point.
      left_cumsum_vec[sample_i] = 0;
    }
    right_cumsum_vec[sample_i] = 
      sample_cumsum_mat[n_bins*sample_i+seg2_LastIndex-1];
    last_cumsum_vec[sample_i] = sample_cumsum_mat[n_bins*(sample_i+1)-1];
  }
  free(sample_cumsum_mat);
  free(sample_count_mat);

  //cleanup!
  free(left_cumsum_vec);
  free(right_cumsum_vec);
  free(last_cumsum_vec);

  return status;
}


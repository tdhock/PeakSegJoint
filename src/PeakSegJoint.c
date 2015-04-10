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

int LossIndex_compare(const void *a, const void *b){
  const struct LossIndex *A = a, *B = b;
  //printf("compare %f %f\n", A->loss, B->loss);
  return (int)(A->loss - B->loss);
}

int PeakSegJointHeuristicStep1(
  struct Profile **samples,
  int n_samples,
  int bin_factor,
  struct PeakSegJointModelList *model_list
  ){
  int sample_i, coverage_i, 
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
  int seg1_chromStart = unfilled_chromStart - extra_before;
  int seg3_chromEnd = unfilled_chromEnd + extra_after;
  int bases = seg3_chromEnd - seg1_chromStart;
  model_list->seg_start_end[0] = seg1_chromStart;
  model_list->seg_start_end[1] = seg3_chromEnd;
  // sample_*_mat variables are matrices n_bins x n_samples (in
  // contrast to model_*_mat which are n_bins x n_segments=3).
  int *sample_count_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  int *count_vec, *cumsum_vec, cumsum_value;
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    count_vec = sample_count_mat + n_bins*sample_i;
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
  int *sample_cumsum_mat = (int*) malloc(n_bins * n_samples * sizeof(int));
  double *flat_loss_vec = malloc(sizeof(double)*n_samples);
  struct LossIndex *diff_index_vec = 
    malloc(sizeof(struct LossIndex)*n_samples);
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
    //printf("\n");
    bases_value = (double)bases;
    mean_value = cumsum_value / bases_value;
    //printf("sample_i=%d cumsum=%d bases=%d\n", sample_i, cumsum_value, bases);
    model_list->sample_mean_vec[sample_i] = mean_value;
    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
    flat_loss_vec[sample_i] = loss_value;
  }

  int seg1_LastIndex, seg2_LastIndex;
  double *seg1_mean_vec = malloc(sizeof(double)*n_samples);
  double *seg2_mean_vec = malloc(sizeof(double)*n_samples);
  double *seg3_mean_vec = malloc(sizeof(double)*n_samples);
  double *peak_loss_vec = malloc(sizeof(double)*n_samples);
  double *seg1_loss_vec = malloc(sizeof(double)*n_samples);
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
	    peak_loss_vec[sample_i]-flat_loss_vec[sample_i];
	  /* printf("%d peak=%f flat=%f diff=%f\n",  */
	  /* 	 sample_i, */
	  /* 	 peak_loss_vec[sample_i],  */
	  /* 	 flat_loss_vec[sample_i], */
	  /* 	 diff_index_vec[sample_i].loss); */
	  n_feasible++;
	}
      }//sample_i
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
    }//seg2_LastIndex
  }//seg2_FirstIndex
  free(sample_cumsum_mat);
  free(sample_count_mat);
  free(flat_loss_vec);
  free(peak_loss_vec);
  free(seg1_loss_vec);
  free(seg1_mean_vec);
  free(seg2_mean_vec);
  free(seg3_mean_vec);
  free(diff_index_vec);
  return status;
}


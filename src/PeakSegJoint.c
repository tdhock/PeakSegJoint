#include "PeakSegJoint.h"
#include "OptimalPoissonLoss.h"
#include "binSum.h"
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <limits.h>

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
  if(n_samples == 0){
    return ERROR_NO_COVERAGE_DATA;
  }
  int chromStart, chromEnd, unfilled_chromStart, unfilled_chromEnd;
  struct Profile *profile, *samples = profile_list->profile_vec;
  struct PeakSegJointModel *model;
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
  model_list->data_start_end[0] = unfilled_chromStart;
  model_list->data_start_end[1] = unfilled_chromEnd;
  //printf("data_start_end=[%d,%d]\n", unfilled_chromStart, unfilled_chromEnd);
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
    return ERROR_BIN_FACTOR_TOO_LARGE;
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
  model_list->n_bins[0] = n_bins;
  model_list->bases_per_bin[0] = bases_per_bin;
  model_list->bin_factor[0] = bin_factor;
  double extra_bases = n_bins  * bases_per_bin - unfilled_bases;//int overflow
  double extra_before = extra_bases/2;
  double extra_after = extra_bases - extra_before;
  //int extra_count;
  int seg1_chromStart = unfilled_chromStart - extra_before;
  double seg3_chromEnd_dbl = unfilled_chromEnd + extra_after;
  int seg3_chromEnd = (seg3_chromEnd_dbl > INT_MAX) ? INT_MAX : (unfilled_chromEnd + extra_after);
  model_list->bin_start_end[0] = seg1_chromStart;
  model_list->bin_start_end[1] = seg3_chromEnd;
  //printf("bin_start_end=[%d,%d]\n", seg1_chromStart, seg3_chromEnd);
  // sample_*_mat variables are matrices n_bins x n_samples (in
  // contrast to model_*_mat which are n_bins x n_segments=3).
  double *sample_count_mat = Calloc(n_bins * n_samples, double);
  double *count_vec, *cumsum_vec, cumsum_value;
  int status;
  for(int sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples + sample_i;
    count_vec = sample_count_mat + n_bins*sample_i;
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    count_vec,
		    bases_per_bin, n_bins, seg1_chromStart, 
		    EMPTY_AS_ZERO);
    /* printf("initial sample_i=%d start=%d\n", sample_i, seg1_chromStart); */
    /* for(int bin_i=0; bin_i < n_bins; bin_i++){ */
    /*   printf("%d ", count_vec[bin_i]); */
    /* } */
    /* printf("\n"); */
    if(status != 0){
      Free(sample_count_mat);
      return status;
    }
    /* Profiles may not have the same first chromStart and last
     * chromEnd, so assume any values outside the observed range are
     * zeros.
     */

    /* The old code below would be useful for the case where we would
     * want to subtract away those data:

    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    &extra_count,
		    unfilled_chromStart - seg1_chromStart,
		    1,
		    seg1_chromStart,
		    EMPTY_AS_ZERO);
    if(status != 0){
      Free(sample_count_mat);
      return status;
    }
    count_vec[0] -= extra_count;
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    &extra_count,
		    seg3_chromEnd - unfilled_chromEnd,
		    1,
		    unfilled_chromEnd,
		    EMPTY_AS_ZERO);
    if(status != 0){
      Free(sample_count_mat);
      return status;
    }
    count_vec[n_bins - 1] -= extra_count;
    */
  }//for sample_i
  int bin_i, offset;
  double mean_value, loss_value;
  double flat_loss_total = 0.0;
  double *sample_cumsum_mat = Calloc(n_bins * n_samples, double);
  struct LossIndex *diff_index_vec = Calloc(n_samples,struct LossIndex);
  int n_feasible;
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
    model_list->last_cumsum_vec[sample_i] = cumsum_value;
    //printf("\n");
    mean_value = cumsum_value / data_bases;
    /* printf("sample_i=%d cumsum=%d bases=%f\n",  */
    /* 	   sample_i, cumsum_value, data_bases); */
    model_list->sample_mean_vec[sample_i] = mean_value;
    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
    model_list->flat_loss_vec[sample_i] = loss_value;
    flat_loss_total += loss_value;
  }
  model_list->model_vec[0].loss[0] = flat_loss_total;

  int n_peaks;
  int sample_i;
  double *seg1_mean_vec = model_list->mean_mat;
  double *seg2_mean_vec = model_list->mean_mat + n_samples;
  double *seg3_mean_vec = model_list->mean_mat + n_samples*2;
  double *peak_loss_vec = Calloc(n_samples,double);
  double *seg1_loss_vec = Calloc(n_samples,double);
  /*
    The for loops below implement the GridSearch() function mentioned
    on line 2 of the JointZoom algorithm in the PeakSegJoint paper.
  */
  for(int seg1_LastIndex=0; seg1_LastIndex < n_bins-2; seg1_LastIndex++){
    for(int sample_i=0; sample_i < n_samples; sample_i++){
      cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
      cumsum_value = cumsum_vec[seg1_LastIndex];
      bin_bases = (seg1_LastIndex+1)*bases_per_bin;
      data_bases = bin_bases - (double)extra_before;
      /* printf("sample_i=%d extra_before=%d bin_bases=%f data_bases=%f\n", */
      /* 	     sample_i, extra_before, bin_bases, data_bases); */
      mean_value = cumsum_value/data_bases;
      seg1_mean_vec[sample_i] = mean_value;
      loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
      seg1_loss_vec[sample_i] = loss_value;
    }
    for(int seg2_LastIndex=seg1_LastIndex+1; 
	seg2_LastIndex < n_bins-1; 
	seg2_LastIndex++){
      n_feasible=0;
      for(sample_i=0; sample_i < n_samples; sample_i++){
	peak_loss_vec[sample_i] = seg1_loss_vec[sample_i];
	cumsum_vec = sample_cumsum_mat + n_bins*sample_i;
	//segment 2.
	cumsum_value = cumsum_vec[seg2_LastIndex]-cumsum_vec[seg1_LastIndex];
	data_bases = (seg2_LastIndex-seg1_LastIndex)*bases_per_bin;
	mean_value = cumsum_value/data_bases;
	seg2_mean_vec[sample_i] = mean_value;
	loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	peak_loss_vec[sample_i] += loss_value;
	//segment 3.
	cumsum_value = cumsum_vec[n_bins-1]-cumsum_vec[seg2_LastIndex];
	bin_bases = (n_bins-1-seg2_LastIndex)*bases_per_bin;
	data_bases = bin_bases - extra_after;
	mean_value = cumsum_value/data_bases;
	seg3_mean_vec[sample_i] = mean_value;
	loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	peak_loss_vec[sample_i] += loss_value;
	/* printf("sample_i=%d means %f %f %f\n", */
	/*        sample_i, */
	/*        seg1_mean_vec[sample_i], */
	/*        seg2_mean_vec[sample_i], */
	/*        seg3_mean_vec[sample_i]); */
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
	for(int model_i=0; model_i < n_feasible; model_i++){
	  // start from loss of all samples with 1 segment.
	  loss_value = flat_loss_total;
	  n_peaks = model_i + 1;
	  for(int diff_i=0; diff_i < n_peaks; diff_i++){
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
	    for(int diff_i=0; diff_i < n_peaks; diff_i++){
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
  Free(sample_cumsum_mat);
  Free(sample_count_mat);
  Free(peak_loss_vec);
  Free(seg1_loss_vec);
  Free(diff_index_vec);
  return status;
}

int 
PeakSegJointHeuristicStep2
(struct ProfileList *profile_list,
 struct PeakSegJointModelList *model_list
  ){
  int n_bins = model_list->bin_factor[0] * 2;
  int n_samples = model_list->n_models - 1;
  struct PeakSegJointModel *model;
  struct Profile *profile;
  int bases_per_bin;
  int left_chromStart, right_chromStart;
  double *left_bin_vec = Calloc(n_bins, double);
  double *right_bin_vec = Calloc(n_bins, double);
  double *left_cumsum_mat = Calloc(n_bins * n_samples, double);
  double *right_cumsum_mat = Calloc(n_bins * n_samples, double);
  double *seg1_mean_vec = model_list->mean_mat;
  double *seg2_mean_vec = model_list->mean_mat + n_samples;
  double *seg3_mean_vec = model_list->mean_mat + n_samples*2;
  double *seg1_loss_vec = Calloc(n_samples,double);
  double total_loss, loss_value, mean_value;
  double bin_bases, data_bases;
  int extra_before = model_list->data_start_end[0] - 
    model_list->bin_start_end[0];
  int extra_after = model_list->bin_start_end[1] -
    model_list->data_start_end[1];
  double *left_cumsum_vec, *right_cumsum_vec;
  int status;
  double left_cumsum_value, right_cumsum_value, cumsum_value;
  int peakStart, peakEnd;
  int best_seg1, best_seg2=-1, sample_i;
  /* When performing the minimization over peakStart/End locations, it
   * is possible that at any given bases_per_bin value, there is no
   * better solution than what we found for the previous bases_per_bin
   * value. In that case, we begin the search anew at a lower
   * resolution, but we need to copy the cumsums from the following
   * index of left_right_vec: */
  int no_min_index = model_list->bin_factor[0] - 2;
  for(int n_peaks=1; n_peaks < model_list->n_models; n_peaks++){
    model = model_list->model_vec + n_peaks;
    if(model->loss[0] < INFINITY){
      bases_per_bin = model_list->bases_per_bin[0];
      /*
	The while loop below corresponds to line 3 of the JointZoom
	algorithm from the PeakSegJoint paper.
      */
      while(1 < bases_per_bin){
	best_seg1 = -1; // indicates no min found.
	left_chromStart = model->peak_start_end[0] - bases_per_bin;
	right_chromStart = model->peak_start_end[1] - bases_per_bin;
	/*
	  Below in the C code we decrease the value of bases_per_bin,
	  as in line 4 of the JointZoom algorithm in the PeakSegJoint
	  paper.
	*/
	bases_per_bin /= model_list->bin_factor[0];
	//printf("bases_per_bin=%d left cumsum before:\n", bases_per_bin);
	for(int diff_i=0; diff_i < n_peaks; diff_i++){
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
	    Free(left_bin_vec);
	    Free(right_bin_vec);
	    Free(left_cumsum_mat);
	    Free(right_cumsum_mat);
	    Free(seg1_loss_vec);
	    return status;
	  }
	  left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	  left_cumsum_value = model->left_cumsum_vec[diff_i];
	  right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	  right_cumsum_value = model->right_cumsum_vec[diff_i];
	  //printf("%d ", left_cumsum_value);
	  //printf("%d ", right_cumsum_value);
	  for(int bin_i=0; bin_i < n_bins; bin_i++){
	    left_cumsum_value += left_bin_vec[bin_i];
	    left_cumsum_vec[bin_i] = left_cumsum_value;
	    right_cumsum_value += right_bin_vec[bin_i];
	    right_cumsum_vec[bin_i] = right_cumsum_value;
	  }
	}//for(diff_i
	//printf("\n");

	/* printf("left bases_per_bin=%d n_peaks=%d\n", bases_per_bin, n_peaks); */
	/* for(int diff_i=0;diff_i<n_peaks;diff_i++){ */
	/*   sample_i = model->samples_with_peaks_vec[diff_i]; */
	/*   left_cumsum_vec = left_cumsum_mat + n_bins*sample_i; */
	/*   for(int bin_i=0;bin_i<n_bins;bin_i++){ */
	/*     printf("%d ", left_cumsum_vec[bin_i]); */
	/*   } */
	/*   printf("\n"); */
	/* } */
	/* printf("\n"); */

	/* printf("right bases_per_bin=%d n_peaks=%d\n", bases_per_bin, n_peaks); */
	/* for(int diff_i=0;diff_i<n_peaks;diff_i++){ */
	/*   sample_i = model->samples_with_peaks_vec[diff_i]; */
	/*   right_cumsum_vec = right_cumsum_mat + n_bins*sample_i; */
	/*   for(int bin_i=0;bin_i<n_bins;bin_i++){ */
	/*     printf("%d ", right_cumsum_vec[bin_i]); */
	/*   } */
	/*   printf("\n"); */
	/* } */

	/* 
	   cumsum matrices have been computed, so now use them to
	   compute the loss and feasibility of all models.

	   The for loops below correspond to SearchNearPeak() on line
	   5 of the JointZoom algorithm in the PeakSegJoint paper.
	*/
	for(int seg1_LastIndex=0; seg1_LastIndex < n_bins; seg1_LastIndex++){
	  peakStart = left_chromStart + (seg1_LastIndex+1)*bases_per_bin;
	  //printf("[seg1last=%d] seg1 cumsum bases ", seg1_LastIndex);
	  for(int diff_i=0; diff_i < n_peaks; diff_i++){
	    sample_i = model->samples_with_peaks_vec[diff_i];
	    left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	    cumsum_value = left_cumsum_vec[seg1_LastIndex];
	    bin_bases = peakStart - model_list->bin_start_end[0];
	    data_bases = bin_bases - extra_before;
	    mean_value = cumsum_value/data_bases;
	    //printf("%d %f ", cumsum_value, bases_value); 
	    seg1_mean_vec[sample_i] = mean_value;
	    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	    seg1_loss_vec[sample_i] = loss_value;
	  }
	  //printf("\n");
	  for(int seg2_LastIndex=0; seg2_LastIndex < n_bins; seg2_LastIndex++){
	    peakEnd = right_chromStart + (seg2_LastIndex+1)*bases_per_bin;
	    /* printf("\npeaks=%d[%d,%d]bases_per_bin=%d\n",  */
	    /* 	 n_peaks, peakStart, peakEnd, bases_per_bin); */
	    total_loss = model_list->model_vec[0].loss[0];
	    if(peakEnd <= peakStart){
	      total_loss = INFINITY;
	    }
	    //printf("[seg2last=%d]\n", seg2_LastIndex);
	    for(int diff_i=0; diff_i < n_peaks; diff_i++){
	      sample_i = model->samples_with_peaks_vec[diff_i];
	      left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	      right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	      total_loss -= model_list->flat_loss_vec[sample_i];
	      //segment 1.
	      total_loss += seg1_loss_vec[sample_i];
	      //segment 2.
	      cumsum_value = 
		right_cumsum_vec[seg2_LastIndex]-
		left_cumsum_vec[seg1_LastIndex];
	      data_bases = peakEnd-peakStart;
	      mean_value = cumsum_value/data_bases;
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
	      bin_bases = model_list->bin_start_end[1] - peakEnd;
	      data_bases = bin_bases - extra_after;
	      mean_value = cumsum_value/data_bases;
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
		 seg2_mean_vec[sample_i] <= seg3_mean_vec[sample_i] ||
		 peakStart <= model_list->data_start_end[0] ||
		 model_list->data_start_end[1] <= peakEnd){
		total_loss = INFINITY;
	      }
	    }	  
	    //printf("loss=%f\n", total_loss);
	    if(total_loss < model->loss[0]){
	      model->loss[0] = total_loss;
	      model->peak_start_end[0] = peakStart;
	      model->peak_start_end[1] = peakEnd;
	      /* printf("new best loss=%f [%d,%d]\n", */
	      /* 	   total_loss, seg1_LastIndex, seg2_LastIndex); */
	      best_seg1 = seg1_LastIndex;
	      best_seg2 = seg2_LastIndex;
	      for(int diff_i=0; diff_i < n_peaks; diff_i++){
		sample_i = model->samples_with_peaks_vec[diff_i];
		model->seg1_mean_vec[diff_i] = seg1_mean_vec[sample_i];
		model->seg2_mean_vec[diff_i] = seg2_mean_vec[sample_i];
		model->seg3_mean_vec[diff_i] = seg3_mean_vec[sample_i];
	      }
	    }//total_loss
	  }//seg2_LastIndex
	}//seg1_LastIndex
	/* printf("n_peaks=%d bases_per_bin=%d [%d,%d] loss=%f\n", */
	/*        n_peaks, bases_per_bin, */
	/*        model->peak_start_end[0], model->peak_start_end[1], */
	/*        model->loss[0]); */
	if(best_seg1 == -1){
	  //printf("no min found\n");
	  for(int diff_i=0; diff_i < n_peaks; diff_i++){
	    sample_i = model->samples_with_peaks_vec[diff_i];
	    left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	    model->left_cumsum_vec[diff_i] = 
	      left_cumsum_vec[no_min_index];
	    right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	    model->right_cumsum_vec[diff_i] = 
	      right_cumsum_vec[no_min_index];
	  }//diff_i
	}else{
	  for(int diff_i=0; diff_i < n_peaks; diff_i++){
	    sample_i = model->samples_with_peaks_vec[diff_i];
	    if(best_seg1 != 0){
	      left_cumsum_vec = left_cumsum_mat + n_bins*sample_i;
	      model->left_cumsum_vec[diff_i] = 
		left_cumsum_vec[best_seg1-1];
	    }
	    if(best_seg2 != 0){
	      right_cumsum_vec = right_cumsum_mat + n_bins*sample_i;
	      model->right_cumsum_vec[diff_i] = 
		right_cumsum_vec[best_seg2-1];
	    }
	  }//diff_i
	}
	//printf("\n");
      }//while(1 < bases_per_bin)
    }//if(loss < INFINITY
  }//for(n_peaks
  //printf("Free at end\n");
  Free(left_bin_vec);
  Free(right_bin_vec);
  Free(left_cumsum_mat);
  Free(right_cumsum_mat);
  Free(seg1_loss_vec);
  return 0;
}

int PeakSegJointHeuristicStep3
(struct ProfileList *profile_list,
 struct PeakSegJointModelList *model_list
  ){
  struct PeakSegJointModel *model, *prev_model;
  int n_samples = model_list->n_models - 1;
  double flat_loss_total = model_list->model_vec[0].loss[0];
  double *seg1_mean_vec = model_list->mean_mat;
  double *seg2_mean_vec = model_list->mean_mat + n_samples;
  double *seg3_mean_vec = model_list->mean_mat + n_samples*2;
  struct LossIndex *diff_index_vec = Calloc(n_samples,struct LossIndex);
  int n_feasible, peakStart, peakEnd, status;
  double total;
  int dataStart = model_list->data_start_end[0];
  int dataEnd = model_list->data_start_end[1];
  double data_bases, mean_value, loss_value;
  struct Profile *profile;
  for(int n_peaks=2; n_peaks < model_list->n_models; n_peaks++){
    n_feasible=0;
    model = model_list->model_vec + n_peaks;
    prev_model = model_list->model_vec + n_peaks - 1;
    if(prev_model->loss[0] < INFINITY){
      peakStart = prev_model->peak_start_end[0];
      peakEnd = prev_model->peak_start_end[1];
      for(int sample_i=0; sample_i < n_samples; sample_i++){
	profile = profile_list->profile_vec + sample_i;
	//segment 1.
	status = oneBin(profile->chromStart, profile->chromEnd,
			profile->coverage, profile->n_entries,
			&total, dataStart, peakStart);
	if(status != 0){
	  Free(diff_index_vec);
	  return status;
	}
	data_bases = peakStart - dataStart;
	mean_value = total/data_bases;
	seg1_mean_vec[sample_i] = mean_value;
	loss_value = OptimalPoissonLoss(total, mean_value);
	//segment 2.
	status = oneBin(profile->chromStart, profile->chromEnd,
			profile->coverage, profile->n_entries,
			&total, peakStart, peakEnd);
	if(status != 0){
	  Free(diff_index_vec);
	  return status;
	}
	data_bases = peakEnd - peakStart;
	mean_value = total/data_bases;
	seg2_mean_vec[sample_i] = mean_value;
	loss_value += OptimalPoissonLoss(total, mean_value);
	//segment 3.
	status = oneBin(profile->chromStart, profile->chromEnd,
			profile->coverage, profile->n_entries,
			&total, peakEnd, dataEnd);
	if(status != 0){
	  Free(diff_index_vec);
	  return status;
	}
	data_bases = dataEnd - peakEnd;
	mean_value = total/data_bases;
	seg3_mean_vec[sample_i] = mean_value;
	loss_value += OptimalPoissonLoss(total, mean_value);
	//if feasible, add to list of loss differences.
	/* printf("n_peaks=%d %f %f %f", */
	/* 	     n_peaks, seg1_mean_vec[sample_i], */
	/* 	     seg2_mean_vec[sample_i], */
	/* 	     seg3_mean_vec[sample_i]); */
	double sample_loss_diff =
	  loss_value-model_list->flat_loss_vec[sample_i];
	double sample_sign;
	if(seg1_mean_vec[sample_i] < seg2_mean_vec[sample_i] &&
	   seg3_mean_vec[sample_i] < seg2_mean_vec[sample_i]){
	  //printf(" FEASIBLE");
	  sample_sign = -1.0;
	  diff_index_vec[n_feasible].sample_i = sample_i;
	  diff_index_vec[n_feasible].loss = sample_loss_diff;
	  n_feasible++;
	}else{
	  sample_sign = 1.0;
	}
	model_list->loss_change_vec[sample_i] = sample_loss_diff*sample_sign;
	//printf("\n");
      }//sample_i
      if(n_peaks <= n_feasible){
	qsort(diff_index_vec, n_feasible, sizeof(struct LossIndex), 
	      LossIndex_compare);
	loss_value = flat_loss_total;
	for(int diff_i=0; diff_i < n_peaks; diff_i++){
	  // add loss difference.
	  loss_value += diff_index_vec[diff_i].loss;
	}
	/* printf("n_peaks=%d loss_value=%f model loss=%f",  */
	/* 	     n_peaks, loss_value, model->loss[0]); */
	if(loss_value < model->loss[0]){
	  //printf(" NEW OPTIMUM!");
	  model->loss[0] = loss_value;
	  model->peak_start_end[0] = peakStart;
	  model->peak_start_end[1] = peakEnd;
	  for(int diff_i=0; diff_i < n_peaks; diff_i++){
	    int sample_i = diff_index_vec[diff_i].sample_i;
	    model->samples_with_peaks_vec[diff_i] = sample_i;
	    model->seg1_mean_vec[diff_i] = seg1_mean_vec[sample_i];
	    model->seg2_mean_vec[diff_i] = seg2_mean_vec[sample_i];
	    model->seg3_mean_vec[diff_i] = seg3_mean_vec[sample_i];
	  }
	}//if(loss_value < model->loss[0]
	//printf("\n");
      }//if(n_feasible)
    }//if(prev_model->loss[0] < INFINITY
  }//for(n_peaks
  Free(diff_index_vec);
  return 0;
}

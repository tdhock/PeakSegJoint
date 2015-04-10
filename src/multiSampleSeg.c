/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "multiSampleSeg.h"
#include "binSum.h"
#include "DPA.h"
#include <stdio.h>
#include <stdlib.h>

int multiSampleSegZoom
(int bases_per_bin, // initial bin size in bases.
 int bin_factor, // divide bases_per_bin by this.
 int *optimal_start_end, // initial guess for peakStart/peakEnd.
 int n_samples,
 struct Profile **samples, //n_samples
 int *left_cumsum_vec, int *right_cumsum_vec, //n_samples
 int *last_cumsum_vec, //n_samples
 int first_chromStart,
 int last_chromEnd
  ){
  int sample_i, bin_i;
  int status;
  int n_bins_zoom = bin_factor * 2; 
  int n_cumsum_zoom = n_bins_zoom + 1;
  int left_cumsum_value, right_cumsum_value;
  int *left_count_mat = (int*) malloc(n_bins_zoom * n_samples * sizeof(int));
  int *right_count_mat = (int*) malloc(n_bins_zoom * n_samples * sizeof(int));
  int *left_cumsum_mat = (int*) malloc(
    n_cumsum_zoom * n_samples * sizeof(int));
  int *right_cumsum_mat = (int*) malloc(
    n_cumsum_zoom * n_samples * sizeof(int));
  double min_loss;
  int seg2_FirstIndex, seg1_LastIndex, seg2_LastIndex, seg3_FirstIndex;
  int cumsum_value, bases_value;
  double seg1_loss, seg2_loss, seg3_loss, candidate_loss,
    loss_value, mean_value;
  struct Profile *profile;
  //Now we zoom in, and search on the left and right bins.
  while(bases_per_bin > 1){
    int left_chromStart = optimal_start_end[0] - bases_per_bin;
    int right_chromStart = optimal_start_end[1] - bases_per_bin;
    // Important to determine chromStart before!
    bases_per_bin /= bin_factor;
    /* printf("bases/bin (zoom)=%d [%d,%d]\n", */
    /* 	   bases_per_bin, */
    /* 	   peakStart, */
    /* 	   peakEnd); */

    /* printf("bases_per_bin=%d n_bins=%d left=%d right=%d\n", */
    /* 	   bases_per_bin, n_bins_zoom, */
    /* 	   left_chromStart, right_chromStart); */

    // one bin before and after estimated start/end:
    for(sample_i=0; sample_i < n_samples; sample_i++){
      profile = samples[sample_i];
      status = binSum(profile->chromStart, profile->chromEnd,
		      profile->coverage, profile->n_entries,
		      left_count_mat + n_bins_zoom*sample_i,
		      bases_per_bin, n_bins_zoom, 
		      left_chromStart, 
		      ERROR_EMPTY_BIN);
      if(status != 0){
	//printf("bases/bin=%d left sample_i=%d\n", bases_per_bin, sample_i);
	free(left_count_mat);
	free(right_count_mat);
	free(left_cumsum_mat);
	free(right_cumsum_mat);
	return status;
      }
      status = binSum(profile->chromStart, profile->chromEnd,
		      profile->coverage, profile->n_entries,
		      right_count_mat + n_bins_zoom*sample_i,
		      bases_per_bin, n_bins_zoom, 
		      right_chromStart, 
		      EMPTY_AS_ZERO);
      if(status != 0){
	//printf("bases/bin=%d right sample_i=%d\n", bases_per_bin, sample_i);
	free(left_count_mat);
	free(right_count_mat);
	free(left_cumsum_mat);
	free(right_cumsum_mat);
	return status;
      }
      left_cumsum_value = left_cumsum_vec[sample_i];
      left_cumsum_mat[n_cumsum_zoom*sample_i] = left_cumsum_value;
      right_cumsum_value = right_cumsum_vec[sample_i];
      right_cumsum_mat[n_cumsum_zoom*sample_i] = right_cumsum_value;
      for(bin_i=0; bin_i < n_bins_zoom; bin_i++){
	left_cumsum_value += left_count_mat[n_bins_zoom*sample_i + bin_i];
	left_cumsum_mat[n_cumsum_zoom*sample_i+bin_i+1] = left_cumsum_value;
	right_cumsum_value += right_count_mat[n_bins_zoom*sample_i + bin_i];
	right_cumsum_mat[n_cumsum_zoom*sample_i+bin_i+1] = right_cumsum_value;
      }
    }//for sample_i

    //debug: print left and right cumsum matrices! OK.
    /* for(bin_i=0; bin_i < n_cumsum_zoom; bin_i++){ */
    /*   for(sample_i=0; sample_i < n_samples; sample_i++){ */
    /* 	printf("%d ", left_cumsum_mat[n_cumsum_zoom*sample_i + bin_i]); */
    /*   } */
    /*   printf("\n"); */
    /* } */
    /* printf("\n"); */
    /* for(bin_i=0; bin_i < n_cumsum_zoom; bin_i++){ */
    /*   for(sample_i=0; sample_i < n_samples; sample_i++){ */
    /* 	printf("%d ", right_cumsum_mat[n_cumsum_zoom*sample_i + bin_i]); */
    /*   } */
    /*   printf("\n"); */
    /* } */
  
    int best_seg2_FirstIndex, best_seg3_FirstIndex, 
      seg2_chromEnd, seg1_chromEnd, after_cumsum,
      last_bin_cumsum, cumsum_seg2_end;
    min_loss = INFINITY;
    /* 
       the first possible seg2_FirstIndex value is in fact 2.
       - not 0 since that is the cumsum before,
       - and not 1, since if that was optimal then we 
         would not know the cumsum before in the next step.
	 for example consider bin_factor=2:
	 [     |     |     |     |     |     |     |     ] original bins
	             [                 ] best peak
	       [1 |2 |3 |4 ]     [1 |2 |3 |4 ] possible starts/ends.
	       [     ] seg2_FirstIndex=2 is fine
	                [     ] seg2_FirstIndex=4 is fine
	    [     ] seg2_FirstIndex=1 NOT OK -- unknown cumsum before.
	    so in this case 3 possible starts (2-4),
	    and 4 possible ends (1-4).
    */
    for(seg2_FirstIndex=2; seg2_FirstIndex <= n_bins_zoom; seg2_FirstIndex++){
      seg1_LastIndex = seg2_FirstIndex-1;
      seg1_chromEnd = left_chromStart + seg1_LastIndex*bases_per_bin;
      seg1_loss = 0.0;
      for(sample_i=0; sample_i < n_samples; sample_i++){
	cumsum_value = left_cumsum_mat[n_cumsum_zoom*sample_i+seg1_LastIndex];
	bases_value = seg1_chromEnd - first_chromStart;
	mean_value = ((double)cumsum_value)/bases_value;
	loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	/* printf("sample%dseg1[%d] cumsum=%d bases=%f mean=%f loss=%f\n", */
	/*        sample_i, seg1_LastIndex, */
	/*        cumsum_value, bases_value, mean_value, loss_value); */
	seg1_loss += loss_value;
      }
      for(seg3_FirstIndex=1; seg3_FirstIndex <= n_bins_zoom; seg3_FirstIndex++){
	seg2_LastIndex = seg3_FirstIndex -1;
	seg2_chromEnd = right_chromStart + seg3_FirstIndex*bases_per_bin;
	if(seg1_chromEnd < seg2_chromEnd &&  //seg2 size positive?
	   seg2_chromEnd < last_chromEnd){ //seg3 size positive
	  //printf("%d,%d\n", seg2_FirstIndex, seg3_FirstIndex);
	  seg2_loss = 0.0;
	  seg3_loss = 0.0;
	  for(sample_i=0; sample_i < n_samples; sample_i++){
	    // first compute seg2_loss.
	    cumsum_seg2_end = right_cumsum_mat[
	      n_cumsum_zoom*sample_i+seg3_FirstIndex];
	    cumsum_value = cumsum_seg2_end-
	      left_cumsum_mat[n_cumsum_zoom*sample_i+seg1_LastIndex];
	    bases_value = seg2_chromEnd - seg1_chromEnd;	  
	    mean_value = ((double)cumsum_value)/bases_value;
	    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	    seg2_loss += loss_value;
	    /* printf("sample%dseg2[%d,%d] cumsum=%d bases=%f mean=%f loss=%f\n", */
	    /* 	   sample_i, seg2_FirstIndex, seg2_LastIndex, */
	    /* 	   cumsum_value, bases_value, mean_value, loss_value); */
	    // then compute seg3_loss.
	    last_bin_cumsum = last_cumsum_vec[sample_i];
	    //after_cumsum = -
	    cumsum_value = last_bin_cumsum - cumsum_seg2_end;
	    bases_value = last_chromEnd - seg2_chromEnd;
	    mean_value = ((double)cumsum_value)/bases_value;
	    loss_value = OptimalPoissonLoss(cumsum_value, mean_value);
	    /* printf("sample%dseg3[%d] cumsum=%d bases=%f mean=%f loss=%f\n", */
	    /* 	   sample_i, seg3_FirstIndex, */
	    /* 	   cumsum_value, bases_value, mean_value, loss_value); */
	    seg3_loss += loss_value;
	  }
	  candidate_loss = seg1_loss + seg2_loss + seg3_loss;
	  /* printf("loss %d %d loss=%f\n", */
	  /* 	 seg2_FirstIndex, seg3_FirstIndex, */
	  /* 	 candidate_loss); */
	  if(candidate_loss < min_loss){
	    min_loss = candidate_loss;
	    best_seg2_FirstIndex = seg2_FirstIndex;
	    best_seg3_FirstIndex = seg3_FirstIndex;
	  }
	}
      }//seg3_FirstIndex
    }//seg2_FirstIndex
    /* printf("best %d %d loss=%f final\n", */
    /* 	   best_seg2_FirstIndex, best_seg3_FirstIndex, */
    /* 	   min_loss); */
    for(sample_i=0; sample_i < n_samples; sample_i++){
      left_cumsum_vec[sample_i] = left_cumsum_mat[
	n_cumsum_zoom*sample_i+best_seg2_FirstIndex-2];
      right_cumsum_vec[sample_i] = right_cumsum_mat[
	n_cumsum_zoom*sample_i+best_seg3_FirstIndex-1];
    }
    optimal_start_end[0] = 
      left_chromStart + bases_per_bin * (best_seg2_FirstIndex-1);
    optimal_start_end[1] = 
      right_chromStart + bases_per_bin * (best_seg3_FirstIndex);
  }
  free(left_count_mat);
  free(right_count_mat);
  free(left_cumsum_mat);
  free(right_cumsum_mat);
  return 0;
}

/*
  Fast heuristic solver for the multiple samples with 1 shared peak
  problem. The sub-optimality parameter bin_factor is used to
  down-sample the profiles, which results in a faster segmentation
  problem.
 */
int
multiSampleSegHeuristic(
  struct Profile **samples,
  int n_samples,
  int bin_factor,
  int *optimal_start_end // array of length 2.
  ) {
  int sample_i, coverage_i, min_chromEnd, max_chromStart,
    chromStart, chromEnd;
  struct Profile *profile;
  profile = samples[0];
  min_chromEnd = get_max_chromEnd(profile);
  max_chromStart = get_min_chromStart(profile);
  for(sample_i=1; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    chromStart = get_min_chromStart(profile);
    if(max_chromStart < chromStart){
      max_chromStart = chromStart;
    }
    chromEnd = get_max_chromEnd(profile);
    if(chromEnd < min_chromEnd){
      min_chromEnd = chromEnd;
    }
    //printf("sample_i=%d profile_size=%d\n", sample_i, profile->n_entries);
    //for(coverage_i=0; coverage_i < profile->n_entries; coverage_i++){
    //profile->chromStart[0];
    //}
  }
  /* printf("max_chromStart=%d min_chromEnd=%d\n", */
  /* 	 max_chromStart, min_chromEnd); */
  int bases = min_chromEnd - max_chromStart;
  double bases_value, seg1_loss_value;
  if(bases/bin_factor < 4){
    /*
      4 is smallest the number of data points for which the 3-segment
      optimization problem is not trivial.

      If we don't have at least this many data points for the first
      bin step, than we stop with an error.
    */
    return ERROR_BIN_FACTOR_TOO_LARGE;
  }
  int bases_per_bin = 1;
  while(bases/bases_per_bin/bin_factor >= 4){
    bases_per_bin *= bin_factor;
  }
  int n_bins = bases / bases_per_bin;
  if(bases % bases_per_bin != 0){
    n_bins ++ ;
  }
  
  //printf("n_bins=%d bases_per_bin=%d\n", n_bins, bases_per_bin);

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
  /*
    First step of DPA: compute optimal loss for 1 segment up to data
    point t, for all data points = bins.
  */
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

  optimal_start_end[0] = max_chromStart + bases_per_bin * seg2_FirstIndex;
  optimal_start_end[1] = max_chromStart + bases_per_bin * seg3_FirstIndex;
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

  status = multiSampleSegZoom(
    bases_per_bin,
    bin_factor,
    optimal_start_end,
    n_samples,
    samples,
    left_cumsum_vec, right_cumsum_vec,
    last_cumsum_vec,
    max_chromStart,
    last_chromEnd);

  //cleanup!
  free(left_cumsum_vec);
  free(right_cumsum_vec);
  free(last_cumsum_vec);

  return status;
}

/*
  Implements base-pair level DPA = Dynamic Programming Algorithm (not
  constrained), in order to recover the most likely Poisson model with
  the same 3 segments (but different mean values) on each of n_samples
  profiles.
 */
int
multiSampleSegOptimal(
  struct Profile **samples,
  int n_samples,
  int *optimal_start_end // array of length 2.
  ) {
  int sample_i, coverage_i, min_chromEnd, max_chromStart,
    chromStart, chromEnd;
  struct Profile *profile;
  profile = samples[0];
  min_chromEnd = get_max_chromEnd(profile);
  max_chromStart = get_min_chromStart(profile);
  for(sample_i=1; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    chromStart = get_min_chromStart(profile);
    if(max_chromStart < chromStart){
      max_chromStart = chromStart;
    }
    chromEnd = get_max_chromEnd(profile);
    if(chromEnd < min_chromEnd){
      min_chromEnd = chromEnd;
    }
  }
  int n_bases = min_chromEnd - max_chromStart;

  // sample_*_mat variables are matrices n_bases x n_samples (in
  // contrast to model_*_mat which are n_bases x n_segments=3).
  int *sample_count_mat = (int*) malloc(n_bases * n_samples * sizeof(int));
  int *sample_cumsum_mat = (int*) malloc(n_bases * n_samples * sizeof(int));
  double *sample_mean1_mat = (double*) malloc(
    n_bases * n_samples * sizeof(double));
  double *sample_loss1_mat = (double*) malloc(
    n_bases * n_samples * sizeof(double));
  int status;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    profile = samples[sample_i];
    status = binSum(profile->chromStart, profile->chromEnd,
		    profile->coverage, profile->n_entries,
		    sample_count_mat + n_bases*sample_i,
		    1, n_bases, max_chromStart, 
		    ERROR_EMPTY_BIN);
    if(status != 0){
      return status;
    }
  }//for sample_i
  int base_i, offset;
  int *count_vec, *cumsum_vec, cumsum_value;
  double *mean_vec, *loss_vec, mean_value, loss_value;
  for(sample_i=0; sample_i < n_samples; sample_i++){
    cumsum_value = 0;
    offset = n_bases * sample_i;
    count_vec = sample_count_mat + offset;
    cumsum_vec = sample_cumsum_mat + offset;
    mean_vec = sample_mean1_mat + offset;
    loss_vec = sample_loss1_mat + offset;
    for(base_i=0; base_i < n_bases; base_i++){
      cumsum_value += count_vec[base_i];
      cumsum_vec[base_i] = cumsum_value;
      mean_value = ((double) cumsum_value) / ((double)base_i+1);
      mean_vec[base_i] = mean_value;
      loss_vec[base_i] = OptimalPoissonLoss(cumsum_value, mean_value);
    }
  }
  int maxSegments = 3;
  double *model_loss_mat = (double*) malloc(
    n_bases * maxSegments * sizeof(double));
  int *model_first_mat = (int*) malloc(
    n_bases * maxSegments * sizeof(int));
  for(base_i=0; base_i < n_bases; base_i++){
    model_loss_mat[base_i] = 0;
    for(sample_i=0; sample_i < n_samples; sample_i++){
      model_loss_mat[base_i] += sample_loss1_mat[base_i + n_bases*sample_i];
    }
  }
  int segment_i, LastSeg_FirstIndex, LastSeg_LastIndex, best_FirstIndex;
  double prev_loss, min_loss, LastSeg_loss, candidate_loss;
  for(segment_i=1; segment_i < maxSegments; segment_i++){
    for(LastSeg_LastIndex=segment_i; 
	LastSeg_LastIndex < n_bases; 
	LastSeg_LastIndex++){
      min_loss = INFINITY;
      for(LastSeg_FirstIndex=segment_i; 
	  LastSeg_FirstIndex <= LastSeg_LastIndex;
	  LastSeg_FirstIndex++){
	prev_loss = model_loss_mat[
	  LastSeg_FirstIndex-1 + n_bases*(segment_i-1)];
	LastSeg_loss = 0.0;
	for(sample_i=0; sample_i < n_samples; sample_i++){
	  cumsum_vec = sample_cumsum_mat + n_bases * sample_i;
	  cumsum_value = cumsum_vec[LastSeg_LastIndex] - 
	    cumsum_vec[LastSeg_FirstIndex-1];
	  mean_value = ((double)cumsum_value)/
	    ((double)LastSeg_LastIndex-LastSeg_FirstIndex+1);
	  LastSeg_loss += OptimalPoissonLoss(cumsum_value, mean_value);
	}
	candidate_loss = LastSeg_loss + prev_loss;
	if(candidate_loss < min_loss){
	  min_loss = candidate_loss;
	  best_FirstIndex = LastSeg_FirstIndex;
	}
      }
      model_first_mat[LastSeg_LastIndex + n_bases*segment_i] = best_FirstIndex;
      model_loss_mat[LastSeg_LastIndex + n_bases*segment_i] = min_loss;
    }
  }
  // For the best segmentation in 3 segments up to n_bases-1, the first
  // index of the 3rd segment is
  int seg3_FirstIndex = model_first_mat[n_bases*3-1];
  int seg2_LastIndex = seg3_FirstIndex-1;
  int seg2_FirstIndex = model_first_mat[n_bases*1+seg2_LastIndex];
  int peakStart = max_chromStart + seg2_FirstIndex;
  int peakEnd = max_chromStart + seg3_FirstIndex;

  free(model_loss_mat);
  free(model_first_mat);
  free(sample_count_mat);
  free(sample_cumsum_mat);
  free(sample_mean1_mat);
  free(sample_loss1_mat);
  optimal_start_end[0] = peakStart;
  optimal_start_end[1] = peakEnd;
  return 0;
}


/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "PeakSegJoint.h"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>

// TODO: for each malloc, check for out of memory error.
void * malloc_or_error(size_t size, const char * msg){
  void * ptr = malloc(size);
  if(ptr == 0){
    //memory_error();
  }
  return ptr;
}

SEXP
PeakSegJointHeuristicStep1_interface(
  SEXP profile_list,
  SEXP bin_factor
  ){
  int n_profiles = length(profile_list);
  int profile_i;
  SEXP df, chromStart, chromEnd, coverage;
  // malloc Profiles for input data.
  struct Profile *profile;
  struct Profile *samples = malloc(n_profiles * sizeof(struct Profile));
  for(profile_i=0; profile_i < n_profiles; profile_i++){
    df = VECTOR_ELT(profile_list, profile_i);
    profile = samples + profile_i;
    chromStart = VECTOR_ELT(df, 0);
    chromEnd = VECTOR_ELT(df, 1);
    coverage = VECTOR_ELT(df, 2);
    profile->chromStart = INTEGER(chromStart);
    profile->chromEnd = INTEGER(chromEnd);
    profile->coverage = INTEGER(coverage);
    profile->n_entries = length(chromStart);
  }
  // allocVector for outputs.
  int model_i;
  struct PeakSegJointModel *model;
  struct PeakSegJointModelList model_list;
  model_list.n_models = n_profiles + 1;
  model_list.model_vec = 
    malloc(model_list.n_models * sizeof(struct PeakSegJointModel));
  SEXP model_list_names, model_list_sexp, model_vec_sexp,
    seg_start_end_sexp, sample_mean_sexp, last_cumsum_sexp;
  PROTECT(model_list_sexp = allocVector(VECSXP, 4));

  PROTECT(model_list_names = allocVector(STRSXP, 4));
  SET_STRING_ELT(model_list_names,0,mkChar("models"));
  SET_STRING_ELT(model_list_names,1,mkChar("seg_start_end"));
  SET_STRING_ELT(model_list_names,2,mkChar("sample_mean_vec"));
  SET_STRING_ELT(model_list_names,3,mkChar("last_cumsum_vec"));
  namesgets(model_list_sexp, model_list_names);
  UNPROTECT(1);

  PROTECT(model_vec_sexp = allocVector(VECSXP, model_list.n_models));
  PROTECT(seg_start_end_sexp = allocVector(INTSXP, 2));
  PROTECT(sample_mean_sexp = allocVector(REALSXP, n_profiles));
  PROTECT(last_cumsum_sexp = allocVector(INTSXP, n_profiles));
  SET_VECTOR_ELT(model_list_sexp,0,model_vec_sexp);
  SET_VECTOR_ELT(model_list_sexp,1,seg_start_end_sexp);
  SET_VECTOR_ELT(model_list_sexp,2,sample_mean_sexp);
  SET_VECTOR_ELT(model_list_sexp, 3, last_cumsum_sexp);
  model_list.seg_start_end = INTEGER(seg_start_end_sexp);
  model_list.sample_mean_vec = REAL(sample_mean_sexp);
  model_list.last_cumsum_vec = INTEGER(last_cumsum_sexp);
  UNPROTECT(4);

  SEXP loss_sexp, peak_start_end_sexp, samples_with_peaks_sexp,
    left_cumsum_sexp, right_cumsum_sexp, 
    seg1_mean_sexp, seg2_mean_sexp, seg3_mean_sexp,
    model_sexp, model_names;
  PROTECT(model_names = allocVector(STRSXP, 8));
  SET_STRING_ELT(model_names,0,mkChar("loss"));
  SET_STRING_ELT(model_names,1,mkChar("peak_start_end"));
  SET_STRING_ELT(model_names,2,mkChar("samples_with_peaks_vec"));
  SET_STRING_ELT(model_names,3,mkChar("left_cumsum_vec"));
  SET_STRING_ELT(model_names,4,mkChar("right_cumsum_vec"));
  SET_STRING_ELT(model_names,5,mkChar("seg1_mean_vec"));
  SET_STRING_ELT(model_names,6,mkChar("seg2_mean_vec"));
  SET_STRING_ELT(model_names,7,mkChar("seg3_mean_vec"));
  
  for(model_i=0; model_i < model_list.n_models; model_i++){
    model = model_list.model_vec + model_i;

    PROTECT(model_sexp = allocVector(VECSXP, 8));
    SET_VECTOR_ELT(model_vec_sexp, model_i, model_sexp);
    namesgets(model_sexp, model_names);
    UNPROTECT(1);

    PROTECT(loss_sexp = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(model_sexp, 0, loss_sexp);
    model->loss = REAL(loss_sexp);
    model->loss[0] = INFINITY;
    UNPROTECT(1);

    if(0 < model_i){
      PROTECT(peak_start_end_sexp = allocVector(INTSXP, 2));
      SET_VECTOR_ELT(model_sexp, 1, peak_start_end_sexp);
      model->peak_start_end = INTEGER(peak_start_end_sexp);
      UNPROTECT(1);

      PROTECT(samples_with_peaks_sexp = allocVector(INTSXP, model_i));
      SET_VECTOR_ELT(model_sexp, 2, samples_with_peaks_sexp);
      model->samples_with_peaks_vec = INTEGER(samples_with_peaks_sexp);
      UNPROTECT(1);

      PROTECT(left_cumsum_sexp = allocVector(INTSXP, model_i));
      SET_VECTOR_ELT(model_sexp, 3, left_cumsum_sexp);
      model->left_cumsum_vec = INTEGER(left_cumsum_sexp);
      UNPROTECT(1);

      PROTECT(right_cumsum_sexp = allocVector(INTSXP, model_i));
      SET_VECTOR_ELT(model_sexp, 4, right_cumsum_sexp);
      model->right_cumsum_vec = INTEGER(right_cumsum_sexp);
      UNPROTECT(1);

      PROTECT(seg1_mean_sexp = allocVector(REALSXP, model_i));
      SET_VECTOR_ELT(model_sexp, 5, seg1_mean_sexp);
      model->seg1_mean_vec = REAL(seg1_mean_sexp);
      UNPROTECT(1);

      PROTECT(seg2_mean_sexp = allocVector(REALSXP, model_i));
      SET_VECTOR_ELT(model_sexp, 6, seg2_mean_sexp);
      model->seg2_mean_vec = REAL(seg2_mean_sexp);
      UNPROTECT(1);

      PROTECT(seg3_mean_sexp = allocVector(REALSXP, model_i));
      SET_VECTOR_ELT(model_sexp, 7, seg3_mean_sexp);
      model->seg3_mean_vec = REAL(seg3_mean_sexp);
      UNPROTECT(1);
    }// if(0 < model_i)
  }
  int status;
  status = PeakSegJointHeuristicStep1(
    samples, n_profiles, INTEGER(bin_factor)[0], &model_list);
  // free inputs.
  free(samples);
  free(model_list.model_vec);
  UNPROTECT(2); //model_list_sexp and model_names.
  // TODO: if known status codes...
  if(status != 0){
    error("unrecognized error code %d", status);
  }
  return model_list_sexp;
}


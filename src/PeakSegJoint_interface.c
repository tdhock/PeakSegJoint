/* -*- compile-command: "R CMD INSTALL .." -*- */

#include "binSum.h"
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

void 
Ralloc_profile_list
(SEXP profile_list_sexp, 
 struct ProfileList *profile_list
  ){
  struct Profile *profile;
  int profile_i;
  SEXP df, chromStart, chromEnd, coverage;
  profile_list->n_profiles = length(profile_list_sexp);
  profile_list->profile_vec = 
    malloc(profile_list->n_profiles * sizeof(struct Profile));
  for(profile_i=0; profile_i < profile_list->n_profiles; profile_i++){
    df = VECTOR_ELT(profile_list_sexp, profile_i);
    profile = profile_list->profile_vec + profile_i;
    chromStart = VECTOR_ELT(df, 0);
    chromEnd = VECTOR_ELT(df, 1);
    coverage = VECTOR_ELT(df, 2);
    profile->chromStart = INTEGER(chromStart);
    profile->chromEnd = INTEGER(chromEnd);
    profile->coverage = INTEGER(coverage);
    profile->n_entries = length(chromStart);
  }
}  

void
free_profile_list
(struct ProfileList *profile_list){
  free(profile_list->profile_vec);
}

SEXP allocPeakSegJointModelList(){
  return allocVector(VECSXP, 9);
}

void
Ralloc_model_struct
(SEXP model_list_sexp,
 struct PeakSegJointModelList *model_list
  ){
  SEXP     
    model_list_names, flat_loss_sexp,
    model_vec_sexp, data_start_end_sexp,
    n_bins_sexp, bases_per_bin_sexp, bin_factor_sexp,
    bin_start_end_sexp, sample_mean_sexp, last_cumsum_sexp;
  int model_i;
  int n_profiles = model_list->n_models - 1;
  struct PeakSegJointModel *model;

  PROTECT(model_list_names = allocVector(STRSXP, 9));
  SET_STRING_ELT(model_list_names,0,mkChar("models"));
  SET_STRING_ELT(model_list_names,1,mkChar("bin_start_end"));
  SET_STRING_ELT(model_list_names,2,mkChar("sample_mean_vec"));
  SET_STRING_ELT(model_list_names,3,mkChar("last_cumsum_vec"));
  SET_STRING_ELT(model_list_names,4,mkChar("flat_loss_vec"));
  SET_STRING_ELT(model_list_names,5,mkChar("n_bins"));
  SET_STRING_ELT(model_list_names,6,mkChar("bases_per_bin"));
  SET_STRING_ELT(model_list_names,7,mkChar("bin_factor"));
  SET_STRING_ELT(model_list_names,8,mkChar("data_start_end"));
  namesgets(model_list_sexp, model_list_names);
  UNPROTECT(1);

  PROTECT(model_vec_sexp = allocVector(VECSXP, model_list->n_models));
  PROTECT(bin_start_end_sexp = allocVector(INTSXP, 2));
  PROTECT(sample_mean_sexp = allocVector(REALSXP, n_profiles));
  PROTECT(last_cumsum_sexp = allocVector(INTSXP, n_profiles));
  PROTECT(flat_loss_sexp = allocVector(REALSXP, n_profiles));
  PROTECT(n_bins_sexp = allocVector(INTSXP, 1));
  PROTECT(bases_per_bin_sexp = allocVector(INTSXP, 1));
  PROTECT(bin_factor_sexp = allocVector(INTSXP, 1));
  PROTECT(data_start_end_sexp = allocVector(INTSXP, 2));
  SET_VECTOR_ELT(model_list_sexp,0,model_vec_sexp);
  SET_VECTOR_ELT(model_list_sexp,1,bin_start_end_sexp);
  SET_VECTOR_ELT(model_list_sexp,2,sample_mean_sexp);
  SET_VECTOR_ELT(model_list_sexp, 3, last_cumsum_sexp);
  SET_VECTOR_ELT(model_list_sexp, 4, flat_loss_sexp);
  SET_VECTOR_ELT(model_list_sexp, 5, n_bins_sexp);
  SET_VECTOR_ELT(model_list_sexp, 6, bases_per_bin_sexp);
  SET_VECTOR_ELT(model_list_sexp, 7, bin_factor_sexp);
  SET_VECTOR_ELT(model_list_sexp, 8, data_start_end_sexp);
  model_list->bin_start_end = INTEGER(bin_start_end_sexp);
  model_list->sample_mean_vec = REAL(sample_mean_sexp);
  model_list->last_cumsum_vec = INTEGER(last_cumsum_sexp);
  model_list->flat_loss_vec = REAL(flat_loss_sexp);
  model_list->n_bins = INTEGER(n_bins_sexp);
  model_list->bases_per_bin = INTEGER(bases_per_bin_sexp);
  model_list->bin_factor = INTEGER(bin_factor_sexp);
  model_list->data_start_end = INTEGER(data_start_end_sexp);
  UNPROTECT(9);

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
  
  for(model_i=0; model_i < model_list->n_models; model_i++){
    model = model_list->model_vec + model_i;

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
  }//for(model_i)
  UNPROTECT(1);//model_names.
}

struct PeakSegJointModelList *
malloc_PeakSegJointModelList(int n_models){
  struct PeakSegJointModelList *model_list = 
    malloc(sizeof(struct PeakSegJointModelList));
  model_list->n_models = n_models;
  model_list->model_vec = 
    malloc(n_models * sizeof(struct PeakSegJointModel));
  return model_list;
}

void
free_PeakSegJointModelList(struct PeakSegJointModelList *model_list){
  free(model_list->model_vec);
  free(model_list);
}

SEXP
PeakSegJointHeuristicStep1_interface(
  SEXP profile_list_sexp,
  SEXP bin_factor
  ){
  int n_profiles = length(profile_list_sexp);
  // malloc Profiles for input data.
  struct ProfileList profile_list;
  Ralloc_profile_list(profile_list_sexp, &profile_list);
  // allocVector for outputs.
  struct PeakSegJointModelList *model_list = 
    malloc_PeakSegJointModelList(n_profiles+1);
  SEXP model_list_sexp;
  PROTECT(model_list_sexp = allocPeakSegJointModelList());

  /* 
     This calls allocVector and assigns the resulting pointers to the
     elements of model_list.
  */
  Ralloc_model_struct(model_list_sexp, model_list);

  int status;
  status = PeakSegJointHeuristicStep1(
    &profile_list, INTEGER(bin_factor)[0], model_list);
  // free inputs.
  free_profile_list(&profile_list);
  free_PeakSegJointModelList(model_list);
  UNPROTECT(1); //model_list_sexp.
  // TODO: if known status codes...
  if(status != 0){
    error("unrecognized error code %d", status);
  }
  return model_list_sexp;
}

SEXP
PeakSegJointHeuristicStep2_interface(
  SEXP profile_list_sexp,
  SEXP bin_factor
  ){
  int n_profiles = length(profile_list_sexp);
  // malloc Profiles for input data.
  struct ProfileList profile_list;
  Ralloc_profile_list(profile_list_sexp, &profile_list);
  // allocVector for outputs.
  struct PeakSegJointModelList *model_list = 
    malloc_PeakSegJointModelList(n_profiles+1);
  SEXP model_list_sexp;
  PROTECT(model_list_sexp = allocPeakSegJointModelList());

  /* 
     This calls allocVector and assigns the resulting pointers to the
     elements of model_list.
  */
  Ralloc_model_struct(model_list_sexp, model_list);
  
  //printf("before step1\n");
  int status;
  status = PeakSegJointHeuristicStep1(
    &profile_list, INTEGER(bin_factor)[0], model_list);
  //printf("before step2\n");
  if(status == 0){
    status = PeakSegJointHeuristicStep2(&profile_list, model_list);
  }
  //printf("after step2 status=%d\n", status);
  free_profile_list(&profile_list);
  free_PeakSegJointModelList(model_list);
  UNPROTECT(1); //model_list_sexp.
  if(status != 0){
    if(status == ERROR_CHROMSTART_NOT_LESS_THAN_CHROMEND)
      error("chromStart not less than chromEnd");
    if(status == ERROR_CHROMSTART_CHROMEND_MISMATCH)
      error("chromStart[i] != chromEnd[i-1]");
    if(status == ERROR_BIN_FACTOR_TOO_LARGE)
      error("bin factor too large");
    if(status == ERROR_EMPTY_BIN)
      error("empty bin");
    error("unrecognized error code %d", status);
  }
  return model_list_sexp;
}

SEXP
PeakSegJointHeuristic_interface(
  SEXP profile_list_sexp,
  SEXP bin_factor
  ){
  int n_profiles = length(profile_list_sexp);
  // malloc Profiles for input data.
  struct ProfileList profile_list;
  Ralloc_profile_list(profile_list_sexp, &profile_list);
  // allocVector for outputs.
  struct PeakSegJointModelList *model_list = 
    malloc_PeakSegJointModelList(n_profiles+1);
  SEXP model_list_sexp;
  PROTECT(model_list_sexp = allocPeakSegJointModelList());

  /* 
     This calls allocVector and assigns the resulting pointers to the
     elements of model_list.
  */
  Ralloc_model_struct(model_list_sexp, model_list);
  
  //printf("before step1\n");
  int status;
  status = PeakSegJointHeuristicStep1(
    &profile_list, INTEGER(bin_factor)[0], model_list);
  //printf("before step2\n");
  if(status == 0){
    status = PeakSegJointHeuristicStep2(&profile_list, model_list);
  }
  if(status == 0){
    status = PeakSegJointHeuristicStep3(&profile_list, model_list);
  }
  //printf("after step2 status=%d\n", status);
  free_profile_list(&profile_list);
  free_PeakSegJointModelList(model_list);
  UNPROTECT(1); //model_list_sexp.
  if(status != 0){
    if(status == ERROR_CHROMSTART_NOT_LESS_THAN_CHROMEND)
      error("chromStart not less than chromEnd");
    if(status == ERROR_CHROMSTART_CHROMEND_MISMATCH)
      error("chromStart[i] != chromEnd[i-1]");
    if(status == ERROR_BIN_FACTOR_TOO_LARGE)
      error("bin factor too large");
    if(status == ERROR_EMPTY_BIN)
      error("empty bin");
    error("unrecognized error code %d", status);
  }
  return model_list_sexp;
}


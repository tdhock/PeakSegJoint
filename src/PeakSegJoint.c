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
  return 0;
}

#include <math.h>
#include "profile.h"

#define ERROR_BIN_FACTOR_TOO_LARGE 1
#define ERROR_NO_COVERAGE_DATA 2

struct PeakSegJointModel {
  double *loss; 
  // size 1 but pointer so we can set it as the address of an R vector's array.
  int *peak_start_end; // size 2.
  // size of the following vectors are of size
  // n_samples (number of samples with peaks),
  // which is determined by the position in the list:
  // 0, ..., n_samples.
  // (there are a total of n_samples + 1 models).
  int *samples_with_peaks_vec;
  int *left_cumsum_vec;
  int *right_cumsum_vec;
  double *seg1_mean_vec;
  double *seg2_mean_vec;
  double *seg3_mean_vec;
};

struct PeakSegJointModelList {
  int n_models;
  struct PeakSegJointModel *model_vec; // n_models = n_samples + 1.
  int *bin_start_end; // size 2.
  double *sample_mean_vec; //n_samples;
  double *flat_loss_vec; //n_samples;
  int *last_cumsum_vec; //n_samples;
  int *n_bins; //size 1.
  int *bases_per_bin; //size 1.
  int *bin_factor; //size 1.
  int *data_start_end; //size 2.
  double *mean_mat;//size n_samples x 3.
  double *loss_change_vec;//size n_samples.
};

struct LossIndex {
  int sample_i;
  double loss;
};

struct PeakSegJointModelList * 
malloc_model_list(
  int n_samples
  );

void free_model_list(struct PeakSegJointModelList *);

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
  );

int PeakSegJointOptimal(
  struct Profile **samples,
  int n_samples,
  int *optimal_start_end // array of length 2.
  );

int PeakSegJointHeuristic(
  struct Profile **samples,
  int n_samples,
  int bin_factor,
  struct PeakSegJointModelList *model_list
  );

int PeakSegJointHeuristicStep1(
  struct ProfileList *,
  int bin_factor,
  struct PeakSegJointModelList *model_list
  );

int PeakSegJointHeuristicStep2(
  struct ProfileList *,
  struct PeakSegJointModelList *
  );

int PeakSegJointHeuristicStep3(
  struct ProfileList *,
  struct PeakSegJointModelList *
  );

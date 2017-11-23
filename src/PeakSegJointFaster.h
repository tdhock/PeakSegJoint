#include "profile.h"

#define ERROR_FASTER_NO_COVERAGE_DATA 1
#define ERROR_FASTER_BIN_FACTOR_TOO_LARGE 2

int PeakSegJointFaster(
  struct ProfileList *,
  int bin_factor,
  double *mean_mat,
  double *flat_loss_vec,
  double *peak_loss_vec,
  int *peak_start_end,
  int *data_start_end
  );

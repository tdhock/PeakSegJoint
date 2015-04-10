void get_best_FirstIndex(
  double *PrevSegs_loss_vec, // n_data
  int n_data,
  int *sample_cumsum_mat, // n_data x n_samples
  int n_samples,
  int first_possible_index,
  int last_possible_index,
  int bases_per_bin,
  int *best_FirstIndex, //output
  double *best_loss //output
  );

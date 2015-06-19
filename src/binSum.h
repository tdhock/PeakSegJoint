#define ERROR_CHROMSTART_NOT_LESS_THAN_CHROMEND 12
#define ERROR_CHROMSTART_CHROMEND_MISMATCH 13
#define ERROR_EMPTY_BIN 14
#define EMPTY_AS_ZERO 15
#define ERROR_CHROMSTART_BEFORE_PREVIOUS_CHROMEND 16

int binSum
(int *profile_chromStart, 
 int *profile_chromEnd, 
 int *profile_coverage, 
 int n_profiles,
 int *bin_total, 
 int bin_size,
 int n_bins, 
 int bin_chromStart,
  int status_for_empty_bins);

int oneBin
(int *profile_chromStart, 
 int *profile_chromEnd, 
 int *profile_coverage, 
 int n_profiles,
 int *bin_total,
 int bin_chromStart,
 int bin_chromEnd);

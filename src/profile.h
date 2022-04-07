#define PROFILEH

struct ProfileList {
  struct Profile * profile_vec;
  int n_profiles;
};

struct Profile {
  int *chromStart;
  int *chromEnd;
  int *coverage;
  int n_entries;
};

int get_min_chromStart(struct Profile *profile);
int get_max_chromEnd(struct Profile *profile);

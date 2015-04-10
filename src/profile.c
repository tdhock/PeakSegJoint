#include "profile.h"

int get_min_chromStart(struct Profile *profile){
  return profile->chromStart[0];
}

int get_max_chromEnd(struct Profile *profile){
  return profile->chromEnd[profile->n_entries-1];
}


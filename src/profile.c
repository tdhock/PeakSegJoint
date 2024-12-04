#include "profile.h"

int get_min_chromStart(struct Profile *profile){
  if(profile->n_entries < 1){
    return -1;
  }
  return profile->chromStart[0];
}

int get_max_chromEnd(struct Profile *profile){
  if (profile->n_entries < 1){
    return -1;
  }
  return profile->chromEnd[profile->n_entries-1];
}


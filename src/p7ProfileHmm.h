#ifndef P7_HMM_READER_PROFILE_HMM_H
#define P7_HMM_READER_PROFILE_HMM_H

#include <stdbool.h>
#include <stdint.h>
#include "p7HmmReader.h"


void p7HmmInit(struct P7Hmm *phmm);
void p7HmmDealloc(struct P7Hmm *phmm);
struct P7Hmm *p7HmmListAppendHmm(struct P7HmmList *phmmList);
uint32_t hmmReaderGetAlphabetCardinality(const struct P7Hmm *const currentPhmm);

//alphabet cardinality is also required to note to the user that this should only be called
//after the application knows the alphabet of the profile hmm.
enum P7HmmReturnCode p7HmmAllocateModelData(struct P7Hmm *currentPhmm, const uint32_t alphabetCardinality);


#endif

#ifndef P7_HMM_READER_PROFILE_HMM_H
#define P7_HMM_READER_PROFILE_HMM_H

#include <stdbool.h>
#include <stdint.h>
#include "p7HmmReader.h"


/*
 * Function:  p7HmmListInit
 * --------------------
 * Initializes a given P7HmmList struct, but does not allocate any memory.
 *  Functionally, this sets the phmms array to NULL, and count to zero.
 *  When a profile hmm is to be added to the list, the p7HmmListAppendHmm function
 *  will allocate necessary data and update the count.
 *
 *  Inputs:
 *    phmmList: pointer to list struct to initialize.
 */
void p7HmmListInit(struct P7HmmList *phmmList);

/*
 * Function:  p7HmmInit
 * --------------------
 * Initializes all the member data in the phmm to a default configuration.
 *  after this function concludes, all scalars are set to 0, and all arrays are
 *  set to NULL.
 *
 *  Inputs:
 *    phmm: pointer to the profile Hmm to initialize.
 */
void p7HmmInit(struct P7Hmm *phmm);

/*
 * Function:  p7HmmDealloc
 * --------------------
 * Deallocates all arrays in the given profile hmm, and sets their pointers to NULL.
 *
 *  Inputs:
 *    phmm: pointer to profile hmm struct to deallocate.
 */
void p7HmmDealloc(struct P7Hmm *phmm);

/*
 * Function:  p7HmmListAppendHmm
 * --------------------
 * Grows the phmms array in the given phmmList by one, and (if successful)
 *  updates the count to new correct number of elements in the list.
 *
 *  Inputs:
 *    phmmList: pointer to P7HmmList struct that will contain an additional profile hmm.
 *
 *  Returns:
 *    Pointer to the new phmm that now resides at the end of the phmms list, or
 *      returns NULL if the realloc fails to allocate memory.
 *      If NULL is returned, the existing phmmList should not be affected.
 */
struct P7Hmm *p7HmmListAppendHmm(struct P7HmmList *phmmList);

/*
 * Function:  p7HmmGetAlphabetCardinality
 * --------------------
 * Returns the number of symbols in the alphabet of the given P7Hmm.
 *  This function is only meant to be called after parsing the header far enough
 *  to have encountered the 'ALPH' tag and have the alphabet set in the P7Hmm struct.
 *
 *  Inputs:
 *    currentPhmm: pointer to the phmm that you would like to know the alphabet cardinality of.
 *
 *  Returns:
 *    The total number of symbols in the given alphabet.
 *      Possible alphabets, and their cardinalities are as follows:
 *      amino: 20,    DNA: 4,   RNA: 4,   dice: 6,  coins: 2,
 *      or 0 if the alphabet has not been set (ALPHABET_NOT_SET).
 */
uint32_t p7HmmGetAlphabetCardinality(const struct P7Hmm *const currentPhmm);


/*
 * Function:  p7HmmAllocateModelData
 * --------------------
 * Once the model length and alphabet are known, this function allocates all data
 *  for the profile Hmm that has values per node. For example, this includes
 *  emission probabilities, state transition probabilities, and the ancillary data
 *  on the match emission lines, like reference annotations and model masks, if they
 *  are described to be present in the header. When these ancillary data are set to 'no'
 *  in the header, they will be left NULL.
 *
 *  The alphabet and model length MUST be set in the currentPhmm before this function is called.
 *  Otherwise, the function wouldn't know how much memory to allocate.
 *
 *  Inputs:
 *    currentPhmm: pointer to the profile hmm to allocate data for.
 *
 *  Returns:
 *    p7HmmSuccess on success,
 *    p7HmmFormatError if either the alphabet or modelLength are uninitialized.
 *    p7HmmAllocationFailure if there was a failure to allocate data for the profile hmm.
 */
enum P7HmmReturnCode p7HmmAllocateModelData(struct P7Hmm *currentPhmm);


#endif

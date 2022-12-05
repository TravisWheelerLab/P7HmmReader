#ifndef P7_HMM_READER_HMM_PARSE_H
#define P7_HMM_READER_HMM_PARSE_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>


#define P7_HMM_READER_ALPHABET_AMINO "amino"
#define P7_HMM_READER_ALPHABET_DNA   "DNA"
#define P7_HMM_READER_ALPHABET_RNA   "RNA"
#define P7_HMM_READER_ALPHABET_COINS "coins"
#define P7_HMM_READER_ALPHABET_DICE  "dice"


enum P7HmmReturnCode{
  p7HmmSuccess = 0, p7HmmAllocationFailure = -1, p7HmmFormatError = -2
};

enum P7Alphabet{
  amino, DNA, RNA, coins, dice, ALPHABET_NOT_SET
};

struct P7InitialTransitions{
  float beginToM1;
  float beginToInsert0;
  float beginToDelete1;
  float insert0ToMatch1;
  float insert0ToInsert0;
  //the last 2 transitions, (insert0 to match1 and insert0 to insert0) are guaranteed to be
  //values 0.0 and 1.0 respectively, as per the file spec.
};

struct P7Stats{
  float msvGumbelMu;
  float msvGumbelLambda;
  float viterbiGumbelMu;
  float viterbiGumbelLambda;
  float forwardTau;
  float forwardLambda;
};

struct P7StateTransitions{
  float *matchToMatch;    //Mk -> Mk+1
  float *matchToInsert;   //Mk -> Ik
  float *matchToDelete;   //Mk -> Dk+1
  float *insertToMatch;   //Ik -> Mk+1
  float *insertToInsert;  //Ik -> Ik+1
  float *deleteToMatch;   //Dk -> Mk+1
  float *deleteToDelete;  //Dk -> Dk+1
};

struct P7Model{
  float *compo;
  float *insert0Emissions;
  struct P7InitialTransitions initialTransitions;
  float *matchEmissionScores;
  float *insertEmissionScores;
  struct P7StateTransitions stateTransitions;
  //these remaining fields are the optional fields on the match emission line
  uint32_t *mapAnnotations;
  char *consensusResidues;
  char *referenceAnnotation;
  bool *modelMask;
  char *consensusStructure;
};

struct P7Header{
  bool hasReferenceAnnotation;
  bool hasModelMask;
  bool hasConsensusResidue;
  bool hasConsensusStructure;
  bool hasMapAnnotation;
  uint32_t modelLength;
  uint32_t maxLength;
  uint32_t checksum;
  uint32_t numSequences;
  float effectiveNumSequences;
  char *name;
  char *version;
  char *accessionNumber;
  char *description;
  char *date;
  char *commandLineHistory;
  float gatheringThresholds[2];
  float trustedCutoffs[2];
  float noiseCutoffs[2];
  enum P7Alphabet alphabet;
};

struct P7Hmm{
  struct P7Header header;
  struct P7Stats stats;
  struct P7Model model;
};

struct P7HmmList{
  struct P7Hmm *phmms;
  uint32_t count;
};

/*
 * Function:  readP7Hmm
 * --------------------
 * reads the given fileSrc as a profile Hmm file, like those generated by HMMER3.
 *    This function allocates all necessary data. If the file doesn't appear to
 *    meet the format requirements for the HMMER3 spec, readP7Hmm will print an error
 *    message to sterr with the line and error description.
 *
 *  Inputs:
 *    fileSrc: Location of the hmm file to open.
 *    phmmList: Pointer to a P7HmmList, either dynamically allocated by the user,
 *      or allocated on the stack, but uninitialized.
 *
 *  Returns:
 *    P&P7HmmReturnCode represnting the result of the read. Possible returns are:
 *      p7HmmSuccess on successful file read.
 *      p7HmmFormatError if there appears to be a file formatting error,
 *        and therefore the parser could not read the file correctly.
 *      p7HmmAllocationFailure if the file could not be read sucessfully.
 */
enum P7HmmReturnCode readP7Hmm(const char *const fileSrc, struct P7HmmList *phmmList);

/*
 * Function:  p7HmmListDealloc
 * --------------------
 * Deallocates and cleans up the given phmmList. This function will walk through
 *  all hmms in the list, deallocating all their allocated data, until finally deallocating
 *  the list its self.
 *
 *  Inputs:
 *    phmmList: struct containing the list of profile hmms, generated and allocated
 *    by the readP7Hmm function.
 *
 */
void p7HmmListDealloc(struct P7HmmList *phmmList);


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
 * Function:  p7HmmGetMatchEmissionScore
 * --------------------
 * Gets the match emission score from the given phmm, for the specified nodeIndex and symbol of the alphabet.
 *
 *  Inputs:
 *    phmm: pointer to the phmm to extract the match emission score from
 *    nodeIndex: position in the profile hmm, as specified by the node indices in the hmm file.
 *    symbolIndex: which symbol of the profile hmm's alphabet to get the score for.
 *
 *  Returns:
 *    float value of the match emission score, or NaN if the given nodeIndex or symbolIndex is out of range.
 */
float p7HmmGetMatchEmissionScore(const struct P7Hmm *const phmm, uint32_t nodeIndex, uint32_t symbolIndex);

/*
 * Function:  p7HmmGetInsertEmissionScores
 * --------------------
 * Gets the insert emission score from the given phmm, for the specified nodeIndex and symbol of the alphabet.
 *
 *  Inputs:
 *    phmm: pointer to the phmm to extract the match emission score from
 *    nodeIndex: position in the profile hmm, as specified by the node indices in the hmm file.
 *    symbolIndex: which symbol of the profile hmm's alphabet to get the score for.
 *
 *  Returns:
 *    float value of the insert emission score, or NaN if the given nodeIndex or symbolIndex is out of range.
 */
float p7HmmGetInsertEmissionScores(const struct P7Hmm *const phmm, uint32_t nodeIndex, uint32_t symbolIndex);

#endif

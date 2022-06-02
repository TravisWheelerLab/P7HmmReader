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
  amino, DNA, RNA, coins, dice
};

struct P7InitialTransitions{
  float beginToM1;
  float beginToInsert0;
  float beginToDelete1;
  float insert0ToMatch1;
  float insert0ToInsert0;
  //the last 2 transitions, (insert0 to match1 and insert0 to insert0) are guarnteed to be
  //values 0.0 and 1.0 respectively, as per the file spec.
};

struct P7Stats{
  float msvGumbelMu;
  float msvGumbleLambda;
  float viterbiGumbleMu;
  float viterbiGumbleLambda;
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
  float *comp0;
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
  enum P7Alphabet alphabet;
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


enum P7HmmReturnCode readP7Hmm(const char *const fileSrc, struct P7HmmList **phmmList);
void p7HmmListInit(struct P7HmmList *phmmList);
void p7HmmListDealloc(struct P7HmmList *phmmList);

#endif

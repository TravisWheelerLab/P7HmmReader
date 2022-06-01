#include "p7ProfileHmm.h"
#include <stdlib.h>


void p7HmmListInit(struct P7HmmList *phmmList){
  phmmList->phmms = NULL;
  phmmList->count = 0;
}

//returns NULL on error
struct P7Hmm *p7HmmListAppendHmm(struct P7HmmList *phmmList){
  void *profileHmmListPointer = realloc(phmmList->phmms, phmmList->count + 1);
  if(profileHmmListPointer == NULL){
    return NULL;
  }
  else{
    phmmList->phmms = profileHmmListPointer;
    p7HmmInit(&phmmList->phmms[phmmList->count]);  //initialize the newly allocated phmm.
    struct P7Hmm *newlyAllocatedPhmm = &phmmList->phmms[phmmList->count];
    phmmList->count++;
    return newlyAllocatedPhmm;
  }
  // return false; //fall-through case if the malloc failed
}

void p7HmmInit(struct P7Hmm *phmm){
  phmm->header.hasReferenceAnnotation = false;
  phmm->header.hasModelMask = false;
  phmm->header.hasConsensusResidue = false;
  phmm->header.hasConsensusStructure = false;
  phmm->header.hasMapAnnotation = false;
  phmm->header.alphabet = DNA;
  phmm->header.modelLength = 0;
  phmm->header.maxLength = 0;
  phmm->header.numSequences = 0;
  phmm->header.effectiveNumSequences = 0;
  phmm->header.checksum = 0;
  phmm->header.version = NULL;
  phmm->header.name = NULL;
  phmm->header.accessionNumber = NULL;
  phmm->header.description = NULL;
  phmm->header.date = NULL;
  phmm->header.commandLineHistory = NULL;
  phmm->header.gatheringThresholds[0] = 0.0f;
  phmm->header.gatheringThresholds[1] = 0.0f;
  phmm->header.trustedCutoffs[0] = 0.0f;
  phmm->header.trustedCutoffs[1] = 0.0f;
  phmm->header.noiseCutoffs[0] = 0.0f;
  phmm->header.noiseCutoffs[1] = 0.0f;
  phmm->stats.msvGumbelMu = 0;
  phmm->stats.msvGumbleLambda = 0;
  phmm->stats.viterbiGumbleMu = 0;
  phmm->stats.viterbiGumbleLambda = 0;
  phmm->stats.forwardTau = 0;
  phmm->stats.forwardLambda = 0;
  phmm->model.comp0 = NULL;
  phmm->model.insert0Emissions = NULL;
  phmm->model.initialTransitions.beginToM1 = 0;
  phmm->model.initialTransitions.beginToInsert0 = 0;
  phmm->model.initialTransitions.beginToDelete1 = 0;
  phmm->model.matchEmissionScores = NULL;
  phmm->model.insertEmissionScores = NULL;
  phmm->model.stateTransitions.matchToMatch = NULL;
  phmm->model.stateTransitions.matchToInsert = NULL;
  phmm->model.stateTransitions.matchToDelete = NULL;
  phmm->model.stateTransitions.insertToMatch = NULL;
  phmm->model.stateTransitions.insertToInsert = NULL;
  phmm->model.stateTransitions.deleteToMatch = NULL;
  phmm->model.stateTransitions.deleteToDelete = NULL;
  phmm->model.mapAnnotations = NULL;
  phmm->model.consensusResidues = NULL;
  phmm->model.referenceAnnotation = NULL;
  phmm->model.modelMask = NULL;
  phmm->model.consensusStructure = NULL;
}

void p7HmmDealloc(struct P7Hmm *phmm){
  free(phmm->header.version);
  free(phmm->header.accessionNumber);
  free(phmm->header.description);
  free(phmm->header.date);
  free(phmm->model.comp0);
  free(phmm->model.insert0Emissions);
  free(phmm->model.matchEmissionScores);
  free(phmm->model.insertEmissionScores);
  free(phmm->model.stateTransitions.matchToMatch);
  free(phmm->model.stateTransitions.matchToInsert);
  free(phmm->model.stateTransitions.matchToDelete);
  free(phmm->model.stateTransitions.insertToMatch);
  free(phmm->model.stateTransitions.insertToInsert);
  free(phmm->model.stateTransitions.deleteToMatch);
  free(phmm->model.stateTransitions.deleteToDelete);
  free(phmm->model.mapAnnotations);
  free(phmm->model.consensusResidues);
  free(phmm->model.referenceAnnotation);
  free(phmm->model.modelMask);
  free(phmm->model.consensusStructure);
  phmm->header.version = NULL;
  phmm->header.accessionNumber = NULL;
  phmm->header.description = NULL;
  phmm->header.date = NULL;
  phmm->model.comp0 = NULL;
  phmm->model.insert0Emissions = NULL;
  phmm->model.matchEmissionScores = NULL;
  phmm->model.insertEmissionScores = NULL;
  phmm->model.stateTransitions.matchToMatch = NULL;
  phmm->model.stateTransitions.matchToInsert = NULL;
  phmm->model.stateTransitions.matchToDelete = NULL;
  phmm->model.stateTransitions.insertToMatch = NULL;
  phmm->model.stateTransitions.insertToInsert = NULL;
  phmm->model.stateTransitions.deleteToMatch = NULL;
  phmm->model.stateTransitions.deleteToDelete = NULL;
  phmm->model.mapAnnotations = NULL;
  phmm->model.consensusResidues = NULL;
  phmm->model.referenceAnnotation = NULL;
  phmm->model.modelMask = NULL;
  phmm->model.consensusStructure = NULL;
}

void p7HmmListDealloc(struct P7HmmList *phmmList){
  for(size_t i = 0; i < phmmList->count; i++){
    p7HmmDealloc(&phmmList->phmms[i]);
  }
  free(phmmList->phmms);
  phmmList->phmms = NULL;
  phmmList->count = 0;
}

//returns 0 if the alphabet type is unsupported, although this shouldn't happen in practice
//unless the struct has been tampered with.
uint32_t hmmReaderGetAlphabetCardinality(const struct P7Hmm *const currentPhmm){
  switch(currentPhmm->header.alphabet){
    case(amino): return 20;
    case(DNA):  return 4;
    case(RNA): return 4;
    case(coins): return 2;
    case(dice): return 6;
    default: return 0;
  }
}

//allocates model arrays for the given phmm, based on its header data.
//the application must know the alphabet being used in order to allocate memory correctly,
//so this will likely be done after reading the header.
enum P7HmmReturnCode p7HmmAllocateModelData(struct P7Hmm *currentPhmm, const uint32_t alphabetCardinality){
  currentPhmm->model.insert0Emissions       = malloc(alphabetCardinality * sizeof(float));
  currentPhmm->model.matchEmissionScores    = malloc(alphabetCardinality * sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.insertEmissionScores   = malloc(alphabetCardinality * sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.matchToMatch    = malloc(sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.matchToInsert   = malloc(sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.matchToDelete   = malloc(sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.insertToMatch   = malloc(sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.insertToInsert  = malloc(sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.deleteToMatch   = malloc(sizeof(float) * currentPhmm->header.modelLength);
  currentPhmm->model.stateTransitions.deleteToDelete  = malloc(sizeof(float) * currentPhmm->header.modelLength);

  //bitwise OR the allocated arrays togeter to determine if the allocation suceeded
  bool allAllocationsSuccessful = currentPhmm->model.insert0Emissions != NULL &&
    currentPhmm->model.matchEmissionScores != NULL && currentPhmm->model.insertEmissionScores != NULL &&
    currentPhmm->model.stateTransitions.matchToMatch != NULL && currentPhmm->model.stateTransitions.matchToInsert != NULL &&
    currentPhmm->model.stateTransitions.matchToDelete != NULL && currentPhmm->model.stateTransitions.insertToMatch != NULL &&
    currentPhmm->model.stateTransitions.insertToInsert != NULL && currentPhmm->model.stateTransitions.deleteToMatch != NULL &&
    currentPhmm->model.stateTransitions.deleteToDelete != NULL;

  if(allAllocationsSuccessful){
    return p7HmmSuccess;
  }
  else{
    return p7HmmAllocationFailure;
  }
}

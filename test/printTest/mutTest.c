#include <stdio.h>
#include <string.h>
#include "../../src/p7HmmReader.h"
#include "../../src/p7ProfileHmm.h"
#include "../test.h"

char *mutFileSrc = "mut.hmm";


char printBuffer[2048];

bool floatCompare(float f1, float f2){
  const float threshold = .00001f;
  float difference = f1 - f2;
  if(difference < 0){
    difference *= -1;
  }
  return difference < threshold;
}

int main(int argc, char ** argv){
  printf("\n\tstarting amylase test\n");
  struct P7HmmList phmmList;
  enum P7HmmReturnCode rc = readP7Hmm(mutFileSrc, &phmmList);
  sprintf(printBuffer, "read hmm gave return code %u, expected success (%u)\n", rc, p7HmmSuccess);

  amalyseHmmTest(&phmmList.phmms[0]);
  p7HmmListDealloc(&phmmList);

  printf("\n\tstarting ox test\n");
  rc = readP7Hmm(oxFileSrc, &phmmList);
  testAssertString(rc == p7HmmSuccess, printBuffer);
  testAssertString(phmmList.phmms != NULL, "phmmList Returned Null");
  testAssertString(phmmList.count == 1, "phmmList did not have expected count of 1");
  oxHmmTest(&phmmList.phmms[0]);
  p7HmmListDealloc(&phmmList);

  printf("\n\tstarting t2 test\n");
  rc = readP7Hmm(t2FileSrc, &phmmList);
  testAssertString(rc == p7HmmSuccess, printBuffer);
  testAssertString(phmmList.phmms != NULL, "phmmList Returned Null");
  testAssertString(phmmList.count == 1, "phmmList did not have expected count of 1");
  t2HmmTest(&phmmList.phmms[0]);
  p7HmmListDealloc(&phmmList);

  printf("\n\tstarting thio test\n");
  rc = readP7Hmm(thioFileSrc, &phmmList);
  testAssertString(rc == p7HmmSuccess, printBuffer);
  testAssertString(phmmList.phmms != NULL, "phmmList Returned Null");
  testAssertString(phmmList.count == 1, "phmmList did not have expected count of 1");
  thioHmmTest(&phmmList.phmms[0]);
  p7HmmListDealloc(&phmmList);

  printf("\n\tstarting tae test\n");
  rc = readP7Hmm(taeFileSrc, &phmmList);
  testAssertString(rc == p7HmmSuccess, printBuffer);
  testAssertString(phmmList.phmms != NULL, "phmmList Returned Null");
  testAssertString(phmmList.count == 1, "phmmList did not have expected count of 1");
  taeHmmTest(&phmmList.phmms[0]);
  p7HmmListDealloc(&phmmList);

  printf("\n\tstarting combined test\n");
  rc = readP7Hmm(combinedFileSrc, &phmmList);
  testAssertString(rc == p7HmmSuccess, printBuffer);
  testAssertString(phmmList.phmms != NULL, "phmmList Returned Null");
  testAssertString(phmmList.count == 5, "phmmList did not have expected count of 1");
  combinedHmmTest(&phmmList);
  p7HmmListDealloc(&phmmList);
}


void testHmmHeader(struct P7Hmm *phmm, const char *expectedVersion, const char *expectedName,
  const char *expectedAccession, const char *expectedDescription, const uint32_t expectedLength,
  const enum P7Alphabet expectedAlphabet, const bool expectedRF, const bool expectedModelMask,
  const bool expectedConsResidue, const bool expectedConsStructure, const bool expectedMap,
  const char *expectedDate, const int expectedNseq, const float expectedEffn,
  const int expectedChecksum, const float expectedGA[2], const float expectedTC[2] ,
  const float expectedNC[2], const float expectedMsvStats[2], const float expectedViterbiStats[2],
  const float expectedForwardStats[2]){
    sprintf(printBuffer, "expected version %s, but got version %s.", expectedVersion, phmm->header.version);
    testAssertString(strcmp(expectedVersion, phmm->header.version) == 0, printBuffer);
    sprintf(printBuffer, "expected name %s, but got name %s.", expectedName, phmm->header.name);
    testAssertString(strcmp(expectedName, phmm->header.name) == 0, printBuffer);

    sprintf(printBuffer, "expected accession %s, but got %s.", expectedAccession, phmm->header.accessionNumber);
    testAssertString(strcmp(expectedAccession, phmm->header.accessionNumber) == 0, printBuffer);

    sprintf(printBuffer, "expected description '%s', but got '%s'.", expectedDescription, phmm->header.description);
    testAssertString(strcmp(expectedDescription, phmm->header.description) == 0, printBuffer);

    sprintf(printBuffer, "expected length %u, but got %u.", expectedLength, phmm->header.modelLength);
    testAssertString(expectedLength == phmm->header.modelLength, printBuffer);

    sprintf(printBuffer, "expected alphabet %u, but got %u.", expectedAlphabet, phmm->header.alphabet);
    testAssertString(expectedAlphabet == phmm->header.alphabet, printBuffer);

    sprintf(printBuffer, "expected ref annot %u, but got %u.", expectedRF, phmm->header.hasReferenceAnnotation);
    testAssertString(expectedRF == phmm->header.hasReferenceAnnotation, printBuffer);

    sprintf(printBuffer, "expected model mask %u, but got %u.", expectedModelMask, phmm->header.hasModelMask);
    testAssertString(expectedModelMask == phmm->header.hasModelMask, printBuffer);

    sprintf(printBuffer, "expected consensus residue %u, but got %u.", expectedConsResidue, phmm->header.hasConsensusResidue);
    testAssertString(expectedConsResidue == phmm->header.hasConsensusResidue, printBuffer);

    sprintf(printBuffer, "expected consensus structure %u, but got %u.", expectedConsStructure, phmm->header.hasConsensusStructure);
    testAssertString(expectedConsStructure == phmm->header.hasConsensusStructure, printBuffer);

    sprintf(printBuffer, "expected map annotation %u, but got %u.", expectedMap, phmm->header.hasMapAnnotation);
    testAssertString(expectedMap == phmm->header.hasMapAnnotation, printBuffer);

    sprintf(printBuffer, "expected date %s, but got %s.", expectedDate, phmm->header.date);
    testAssertString(strcmp(expectedDate, phmm->header.date) == 0, printBuffer);

    sprintf(printBuffer, "expected num sequences %u, but got %u.", expectedNseq, phmm->header.numSequences);
    testAssertString(expectedNseq == phmm->header.numSequences, printBuffer);

    sprintf(printBuffer, "expected effective seq num %f, but got %f.", expectedEffn, phmm->header.effectiveNumSequences);
    testAssertString(floatCompare(expectedEffn, phmm->header.effectiveNumSequences), printBuffer);

    sprintf(printBuffer, "expected effective seq num %u, but got %u.", expectedChecksum, phmm->header.checksum);
    testAssertString(expectedChecksum == phmm->header.checksum, printBuffer);

    sprintf(printBuffer, "expected gathering threshold 0 as %f, but got %f.", expectedGA[0], phmm->header.gatheringThresholds[0]);
    testAssertString(floatCompare(expectedGA[0], phmm->header.gatheringThresholds[0]), printBuffer);

    sprintf(printBuffer, "expected gathering threshold 1 as %f, but got %f.", expectedGA[1], phmm->header.gatheringThresholds[1]);
    testAssertString(floatCompare(expectedGA[1], phmm->header.gatheringThresholds[1]), printBuffer);

    sprintf(printBuffer, "expected trusted cutoff 0 as %f, but got %f.", expectedTC[0], phmm->header.trustedCutoffs[0]);
    testAssertString(floatCompare(expectedTC[0], phmm->header.trustedCutoffs[0]), printBuffer);

    sprintf(printBuffer, "expected trusted cutoff 1 as %f, but got %f.", expectedTC[1], phmm->header.trustedCutoffs[1]);
    testAssertString(floatCompare(expectedTC[1], phmm->header.trustedCutoffs[1]), printBuffer);

    sprintf(printBuffer, "expected noise cutoff 0 as %f, but got %f.", expectedNC[0], phmm->header.noiseCutoffs[0]);
    testAssertString(floatCompare(expectedNC[0], phmm->header.noiseCutoffs[0]), printBuffer);

    sprintf(printBuffer, "expected noise cutoff 1 as %f, but got %f.", expectedNC[0], phmm->header.noiseCutoffs[0]);
    testAssertString(floatCompare(expectedNC[1], phmm->header.noiseCutoffs[1]), printBuffer);

    sprintf(printBuffer, "expected msv gumbel mu %f, but got %f.", expectedMsvStats[0], phmm->stats.msvGumbelMu);
    testAssertString(floatCompare(expectedMsvStats[0], phmm->stats.msvGumbelMu), printBuffer);

    sprintf(printBuffer, "expected msv gumbel lambda %f, but got %f.", expectedMsvStats[1], phmm->stats.msvGumbelLambda);
    testAssertString(floatCompare(expectedMsvStats[1], phmm->stats.msvGumbelLambda), printBuffer);

    sprintf(printBuffer, "expected viterbi gumbel mu %f, but got %f.", expectedViterbiStats[0], phmm->stats.viterbiGumbelMu);
    testAssertString(floatCompare(expectedViterbiStats[0], phmm->stats.viterbiGumbelMu), printBuffer);

    sprintf(printBuffer, "expected viterbi gumbel lambda %f, but got %f.", expectedViterbiStats[1], phmm->stats.viterbiGumbelLambda);
    testAssertString(floatCompare(expectedViterbiStats[1], phmm->stats.viterbiGumbelLambda), printBuffer);

    sprintf(printBuffer, "expected forward gumbel tau %f, but got %f.", expectedForwardStats[0], phmm->stats.forwardTau);
    testAssertString(floatCompare(expectedForwardStats[0], phmm->stats.forwardTau), printBuffer);

    sprintf(printBuffer, "expected forward gumbel lambda %f, but got %f.", expectedForwardStats[1], phmm->stats.forwardLambda);
    testAssertString(floatCompare(expectedForwardStats[1], phmm->stats.forwardLambda), printBuffer);
}

void printSanityCheckInfo(struct P7Hmm *phmm){
  printf("---------------------------------\n");
  uint32_t alphabetCardinality = p7HmmGetAlphabetCardinality(phmm);
  if(phmm->model.compo != NULL){
    printf("\tCOMPO");
    for(size_t i = 0; i < alphabetCardinality; i++){
      printf("  %f", phmm->model.compo[i]);
    }
    printf("]\n  ");
  }

  for(size_t i = 0; i < alphabetCardinality; i++){
    printf("  %f",phmm->model.insert0Emissions[i]);
  }
  printf("]\n  ");
  printf("  %f  %f  %f  %f  %f]\n", phmm->model.initialTransitions.beginToM1,
  phmm->model.initialTransitions.beginToInsert0, phmm->model.initialTransitions.beginToDelete1,
  phmm->model.initialTransitions.insert0ToMatch1, phmm->model.initialTransitions.insert0ToInsert0);


  for(int nodeIndex = 0; nodeIndex < 3; nodeIndex++){
    printf("  %u\t", nodeIndex);
    for(int i = 0; i < alphabetCardinality; i++){
      printf("%f", phmm->model.matchEmissionScores[nodeIndex*alphabetCardinality + i]);
    }
    printf("  ");
    if(phmm->header.hasMapAnnotation){
      printf("  %u", phmm->model.mapAnnotations[nodeIndex]);
    }
    else{
      printf(" -");
    }
    if(phmm->header.hasConsensusResidue){
      printf("  %c", phmm->model.consensusResidues[nodeIndex]);
    }
    else{
      printf("  -");
    }
    if(phmm->header.hasReferenceAnnotation){
      printf("  %c", phmm->model.referenceAnnotation[nodeIndex]);
    }
    else{
      printf("  -");
    }
    if(phmm->header.hasModelMask){
      printf("  %c", phmm->model.modelMask[nodeIndex]);
    }
    else{
      printf("  -");
    }
    if(phmm->header.hasConsensusStructure){
      printf("  %c", phmm->model.consensusStructure[nodeIndex]);
    }
    else{
      printf(" -");
    }

    for(int i = 0; i < alphabetCardinality; i++){
      printf("  %f", phmm->model.insertEmissionScores[nodeIndex*alphabetCardinality + i]);
    }
    printf("]\n  %f  %f  %f  %f  %f  %f  %f]\n", phmm->model.stateTransitions.matchToMatch[nodeIndex],
    phmm->model.stateTransitions.matchToInsert[nodeIndex], phmm->model.stateTransitions.matchToDelete[nodeIndex],
    phmm->model.stateTransitions.insertToMatch[nodeIndex], phmm->model.stateTransitions.insertToInsert[nodeIndex],
    phmm->model.stateTransitions.deleteToMatch[nodeIndex], phmm->model.stateTransitions.deleteToDelete[nodeIndex]);
  }
  printf("...\n");
  for(int nodeIndex = phmm->header.modelLength - 3; nodeIndex < phmm->header.modelLength; nodeIndex++){
    printf("  %u", nodeIndex);
    for(int i = 0; i < alphabetCardinality; i++){
      printf(" %f", phmm->model.matchEmissionScores[nodeIndex*alphabetCardinality + i]);
    }
    printf("  ");
    if(phmm->header.hasMapAnnotation){
      printf("  %u", phmm->model.mapAnnotations[nodeIndex]);
    }
    else{
      printf(" -");
    }
    if(phmm->header.hasConsensusResidue){
      printf("  %c", phmm->model.consensusResidues[nodeIndex]);
    }
    else{
      printf("  -");
    }
    if(phmm->header.hasReferenceAnnotation){
      printf("  %c", phmm->model.referenceAnnotation[nodeIndex]);
    }
    else{
      printf("  -");
    }
    if(phmm->header.hasModelMask){
      printf("  %c", phmm->model.modelMask[nodeIndex]);
    }
    else{
      printf("  -");
    }
    if(phmm->header.hasConsensusStructure){
      printf("  %c", phmm->model.consensusStructure[nodeIndex]);
    }
    else{
      printf(" -");
    }

    for(int i = 0; i < alphabetCardinality; i++){
      printf("  %f", phmm->model.insertEmissionScores[nodeIndex*alphabetCardinality + i]);
    }
    printf("\n %f %f %f %f %f %f %f \n", phmm->model.stateTransitions.matchToMatch[nodeIndex],
    phmm->model.stateTransitions.matchToInsert[nodeIndex], phmm->model.stateTransitions.matchToDelete[nodeIndex],
    phmm->model.stateTransitions.insertToMatch[nodeIndex], phmm->model.stateTransitions.insertToInsert[nodeIndex],
    phmm->model.stateTransitions.deleteToMatch[nodeIndex], phmm->model.stateTransitions.deleteToDelete[nodeIndex]);
  }
}

void amalyseHmmTest(struct P7Hmm *phmm){
  const char *expectedVersion = "HMMER3/f [3.1b2 | February 2015]";
  const char *expectedName  = "Alpha-amylase";
  const char *expectedAccession = "PF00128.27";
  const char * expectedDescription = "Alpha amylase, catalytic domain";
  const uint32_t expectedLength = 336;
  const enum P7Alphabet expectedAlphabet = amino;
  const bool expectedRF = false;
  const bool expectedModelMask = false;
  const bool expectedConsResidue = true;
  const bool expectedConsStructure = true;
  const bool expectedMap = true;
  const char *expectedDate = "Wed Oct 20 22:19:58 2021";
  const int expectedNseq = 18;
  const float expectedEffn = 2.008301;
  const int expectedChecksum = 181186707;
  const float expectedGA[2] = {22, 22};
  const float expectedTC[2] = {22, 22};
  const float expectedNC[2] = {21.9, 21.9};
  const float expectedMsvStats[2] = {-11.1495, 0.70041};
  const float expectedViterbiStats[2] = {-11.9808,  0.70041};
  const float expectedForwardStats[2] = {-5.1114, 0.70041};

  testHmmHeader(phmm, expectedVersion, expectedName, expectedAccession,
    expectedDescription, expectedLength, expectedAlphabet, expectedRF,
    expectedModelMask, expectedConsResidue, expectedConsStructure, expectedMap,
    expectedDate, expectedNseq, expectedEffn, expectedChecksum, expectedGA,
    expectedTC, expectedNC, expectedMsvStats, expectedViterbiStats,
    expectedForwardStats);

  printSanityCheckInfo(phmm);
}

void oxHmmTest(struct P7Hmm *phmm){
  const char *expectedVersion = "HMMER3/f [3.1b2 | February 2015]";
  const char *expectedName  = "OxRdtase_C";
  const char *expectedAccession = "PF19858.2";
  const char * expectedDescription = "Bacterial oxidoreductases, C-terminal";
  const uint32_t expectedLength = 163;
  const enum P7Alphabet expectedAlphabet = amino;
  const bool expectedRF = false;
  const bool expectedModelMask = false;
  const bool expectedConsResidue = true;
  const bool expectedConsStructure = false;
  const bool expectedMap = true;
  const char *expectedDate = "Mon Oct 18 06:29:33 2021";
  const int expectedNseq = 26;
  const float expectedEffn = 0.593506;
  const int expectedChecksum = 2542544526;
  const float expectedGA[2] = {26.3, 26.3};
  const float expectedTC[2] = {26.5, 26.4};
  const float expectedNC[2] = {26.1, 26};
  const float expectedMsvStats[2] = {-10.3066,  0.70808};
  const float expectedViterbiStats[2] = {-11.0414,  0.70808};
  const float expectedForwardStats[2] = {-5.0756,  0.70808};

  testHmmHeader(phmm, expectedVersion, expectedName, expectedAccession,
    expectedDescription, expectedLength, expectedAlphabet, expectedRF,
    expectedModelMask, expectedConsResidue, expectedConsStructure, expectedMap,
    expectedDate, expectedNseq, expectedEffn, expectedChecksum, expectedGA,
    expectedTC, expectedNC, expectedMsvStats, expectedViterbiStats,
    expectedForwardStats);

  printSanityCheckInfo(phmm);
}

void t2HmmTest(struct P7Hmm *phmm){
  const char *expectedVersion = "HMMER3/f [3.1b2 | February 2015]";
  const char *expectedName  = "T2SSL";
  const char *expectedAccession = "PF05134.16";
  const char * expectedDescription = "Type II secretion system (T2SS), protein L";
  const uint32_t expectedLength = 233;
  const enum P7Alphabet expectedAlphabet = amino;
  const bool expectedRF = false;
  const bool expectedModelMask = false;
  const bool expectedConsResidue = true;
  const bool expectedConsStructure = true;
  const bool expectedMap = true;
  const char *expectedDate = "Fri Oct 15 11:00:22 2021";
  const int expectedNseq = 9;
  const float expectedEffn = 2.192871;
  const int expectedChecksum = 3510166814;
  const float expectedGA[2] = {25.5, 25.5};
  const float expectedTC[2] = {25.6, 25.5};
  const float expectedNC[2] = {25.2, 25.4};
  const float expectedMsvStats[2] = {-10.4089,  0.70362};
  const float expectedViterbiStats[2] = {-11.1159,  0.70362};
  const float expectedForwardStats[2] = {-5.2362,  0.70362};

  testHmmHeader(phmm, expectedVersion, expectedName, expectedAccession,
    expectedDescription, expectedLength, expectedAlphabet, expectedRF,
    expectedModelMask, expectedConsResidue, expectedConsStructure, expectedMap,
    expectedDate, expectedNseq, expectedEffn, expectedChecksum, expectedGA,
    expectedTC, expectedNC, expectedMsvStats, expectedViterbiStats,
    expectedForwardStats);

  printSanityCheckInfo(phmm);
}

void thioHmmTest(struct P7Hmm *phmm){
  const char *expectedVersion = "HMMER3/f [3.1b2 | February 2015]";
  const char *expectedName  = "Thioredoxin_10";
  const char *expectedAccession = "PF17991.4";
  const char * expectedDescription = "Thioredoxin like C-terminal domain";
  const uint32_t expectedLength = 142;
  const enum P7Alphabet expectedAlphabet = amino;
  const bool expectedRF = false;
  const bool expectedModelMask = false;
  const bool expectedConsResidue = true;
  const bool expectedConsStructure = false;
  const bool expectedMap = true;
  const char *expectedDate = "Sat Oct 16 13:57:11 2021";
  const int expectedNseq = 158;
  const float expectedEffn = 2.415710;
  const int expectedChecksum = 1083335185;
  const float expectedGA[2] = {27.1, 27.1};
  const float expectedTC[2] = {29.9, 28.7};
  const float expectedNC[2] = {25.4, 26.8};
  const float expectedMsvStats[2] = {-10.0039,  0.71034};
  const float expectedViterbiStats[2] = {-10.9949,  0.71034};
  const float expectedForwardStats[2] = {-3.8790,  0.71034};

  testHmmHeader(phmm, expectedVersion, expectedName, expectedAccession,
    expectedDescription, expectedLength, expectedAlphabet, expectedRF,
    expectedModelMask, expectedConsResidue, expectedConsStructure, expectedMap,
    expectedDate, expectedNseq, expectedEffn, expectedChecksum, expectedGA,
    expectedTC, expectedNC, expectedMsvStats, expectedViterbiStats,
    expectedForwardStats);

  printSanityCheckInfo(phmm);

}

void taeHmmTest(struct P7Hmm *phmm){
  const char *expectedVersion = "HMMER3/f [3.1b2 | February 2015]";
  const char *expectedName  = "Tae4";
  const char *expectedAccession = "PF14113.9";
  const char * expectedDescription = "Type VI secretion system (T6SS), amidase effector protein 4";
  const uint32_t expectedLength = 121;
  const enum P7Alphabet expectedAlphabet = amino;
  const bool expectedRF = false;
  const bool expectedModelMask = false;
  const bool expectedConsResidue = true;
  const bool expectedConsStructure = false;
  const bool expectedMap = true;
  const char *expectedDate = "Fri Oct 22 09:24:18 2021";
  const int expectedNseq = 26;
  const float expectedEffn = 2.250244;
  const int expectedChecksum = 1451623049;
  const float expectedGA[2] = {27, 27};
  const float expectedTC[2] = {27.3, 27.6};
  const float expectedNC[2] = {26.8, 26.6};
  const float expectedMsvStats[2] = {-10.0638,  0.71332};
  const float expectedViterbiStats[2] = {-11.1668,  0.71332};
  const float expectedForwardStats[2] = {-4.1337,  0.71332};

  testHmmHeader(phmm, expectedVersion, expectedName, expectedAccession,
    expectedDescription, expectedLength, expectedAlphabet, expectedRF,
    expectedModelMask, expectedConsResidue, expectedConsStructure, expectedMap,
    expectedDate, expectedNseq, expectedEffn, expectedChecksum, expectedGA,
    expectedTC, expectedNC, expectedMsvStats, expectedViterbiStats,
    expectedForwardStats);

  printSanityCheckInfo(phmm);
}


void combinedHmmTest(struct P7HmmList *phmmList){
  char printBuffer[256];
  sprintf(printBuffer, "combined phmm list should contain 5 phmms, but contained %u\n", phmmList->count);
  testAssertString(phmmList->count == 5, printBuffer);
  if(phmmList->count == 5){
    amalyseHmmTest(&phmmList->phmms[0]);
    oxHmmTest(&phmmList->phmms[1]);
    t2HmmTest(&phmmList->phmms[2]);
    taeHmmTest(&phmmList->phmms[3]);
    thioHmmTest(&phmmList->phmms[4]);

  }
}

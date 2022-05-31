#define  _POSIX_C_SOURCE 200809L     //required for the getline function
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "p7HmmReader.h"
#include "p7ProfileHmm.h"
#include "p7HmmReaderLog.h"


#define P7_HEADER_FORMAT_FLAG "HMMER3"
#define P7_HEADER_NAME_FLAG "NAME"
#define P7_HEADER_ACCESSION_FLAG "ACC"
#define P7_HEADER_DESCRIPTION_FLAG "DESC"
#define P7_HEADER_LENGTH_FLAG "LENG"
#define P7_HEADER_MAXL_FLAG "MAXL"
#define P7_HEADER_ALPHABET_FLAG "ALPH"
#define P7_HEADER_REFERENCE_FLAG "RF"
#define P7_HEADER_MASK_FLAG "MM"
#define P7_HEADER_CONSENSUS_RESIDUE_FLAG "CONS"
#define P7_HEADER_CONSENSUS_STRUCTURE_FLAG "CS"
#define P7_HEADER_MAP_FLAG "MAP"
#define P7_HEADER_DATE_FLAG "DATE"
#define P7_HEADER_COMMAND_FLAG "COM"
#define P7_HEADER_NSEQ_FLAG "NSEQ"
#define P7_HEADER_EFFN_FLAG "EFFN"
#define P7_HEADER_CHECKSUM_FLAG "CKSUM"
#define P7_HEADER_GATHERING_FLAG "GA"
#define P7_HEADER_TRUSTED_FLAG "TC"
#define P7_HEADER_NOISE_FLAG "NC"
#define P7_HEADER_STATS_FLAG "STATS"

#define P7_BODY_HMM_MODEL_START_FLAG "HMM"
#define P7_BODY_COMP0_FLAG "COMP0"
#define P7_BODY_END_FLAG "//"

//TODO: debug running through a phmm file



enum HmmReaderParserState{
  parsingHmmIdle, parsingHmmHeader, parsingHmmModelHead, parsingHmmModelBody
};


enum P7HmmReturnCode readP7Hmm(const char *const fileSrc, struct P7HmmList **phmmList){
  size_t lineBufferLength = 1 << 10; //1024
  char *lineBuffer = malloc(lineBufferLength * sizeof(char));
  if(!lineBuffer){
    printAllocationError(fileSrc, 0, "failed to allocate memory for internal line buffer.");
    return p7HmmAllocationFailure;
  }

    //init the profilehmm list
    *phmmList = malloc(sizeof(struct P7HmmList));
    if(phmmList == NULL){
      printAllocationError(fileSrc, 0, "failed to allocate memory for profileHmmList.");
      return p7HmmAllocationFailure;
    }
    p7HmmListInit(*phmmList);

    struct P7Hmm *currentPhmm = NULL;

    FILE *openedFile = fopen(fileSrc, "r");

    enum HmmReaderParserState parserState = parsingHmmHeader;
    uint32_t alphabetCardinality = 0;
    //for counting which node number we're in when we get to the model body
    uint32_t expectedNodeIndex = 1;
    //for printing errors
    size_t lineNumber = 0;
    while(true){
      lineNumber++;
      size_t numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);

      //check to make sure getline didn't have an allocation failure
      if(numCharactersRead == -1){
        printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer.");
        return p7HmmAllocationFailure;
      }
      else if(numCharactersRead <= 2){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
        if(feof(openedFile)){
          if(parserState == parsingHmmIdle){
            return p7HmmSuccess;
          }
          else return p7HmmFormatError;
        }
      }

      //remove the newline if present at the final buffer position
      if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
        lineBuffer[strlen(lineBuffer) - 1] = 0;
      }

      char *firstTokenLocation = strtok(lineBuffer, " ");
      if(firstTokenLocation == NULL){
        //if the token is null, we can assume that this line only contained whitespace, so we should skip this line
        continue;
      }

      switch(parserState){
        case parsingHmmIdle:
        if(strcmp(firstTokenLocation, P7_HEADER_FORMAT_FLAG) == 0){
          //we've found another header, so append a new hmm to the list
          currentPhmm = p7HmmListAppendHmm(*phmmList);
          if(currentPhmm == NULL){
            printAllocationError(fileSrc, lineNumber, "could not allocate memory to grow the P7ProfileHmmList list.");
            return p7HmmAllocationFailure;
          }
          size_t versionLength = strlen(firstTokenLocation);
          currentPhmm->header.version = malloc((versionLength + 1)*sizeof(char)); //+1 for null terminator
          if(currentPhmm->header.version == NULL){
            printAllocationError(fileSrc, lineNumber, "couldn't allocate buffer for format tag.");
            return p7HmmAllocationFailure;
          }
          strcpy(currentPhmm->header.version, firstTokenLocation);

          //now switch modes to parsing the header
          parserState = parsingHmmHeader;
        }

        //if we're idle and we encounter a line that doesn't start with a "HMMER3" version tag, skip the line until we do.
        break;
        case parsingHmmHeader:
          if(strcmp(firstTokenLocation, P7_HEADER_NAME_FLAG) == 0){
            char *nameText = strtok(NULL, " ");
            //set a default name if this field is missing
            if(nameText == NULL){
              nameText = "none_given";
            }
            else{
              currentPhmm->header.name = malloc(strlen(nameText)+1);
              if(currentPhmm->header.name == NULL){
                printAllocationError(fileSrc, lineNumber, "unalble to allocate memory for name.");
                return p7HmmAllocationFailure;
              }
              strcpy(currentPhmm->header.name, nameText);

            }
          }

          if(strcmp(firstTokenLocation, P7_HEADER_ACCESSION_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse accession number tag (ACC).");
              return p7HmmFormatError;
            }
            currentPhmm->header.accessionNumber = malloc(strlen(flagText)+1);
            if(currentPhmm->header.accessionNumber == NULL){
              printAllocationError(fileSrc, lineNumber, "couldn't allocate buffer for accession number tag (ACC).");
              return p7HmmAllocationFailure;
            }
            strcpy(currentPhmm->header.accessionNumber, flagText);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_DESCRIPTION_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse description tag (DESC).");
              return p7HmmFormatError;
            }
            currentPhmm->header.description = malloc(strlen(flagText)+1);
            if(currentPhmm->header.description == NULL){
              printAllocationError(fileSrc, lineNumber, "couldn't allocate buffer for description tag (DESC).");
              return p7HmmAllocationFailure;
            }
            strcpy(currentPhmm->header.description, flagText);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_LENGTH_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse model length tag (LENG).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, "%u", &currentPhmm->header.modelLength);
            if(scanVariablesFilled < 1){
              printFormatError(fileSrc, lineNumber, "expected positive nonzero integer after model length tag (LENG).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_MAXL_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse max length tag (MAXL).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %u", &currentPhmm->header.maxLength);
            if(scanVariablesFilled == EOF || (currentPhmm->header.maxLength == 0)){
              printFormatError(fileSrc, lineNumber, "expected positive nonzero integer after max length tag (MAXL).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_ALPHABET_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse alphabet tag (ALPH).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, P7_HMM_READER_ALPHABET_AMINO) == 0){
              currentPhmm->header.alphabet = amino;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_DNA) == 0){
              currentPhmm->header.alphabet = DNA;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_RNA) == 0){
              currentPhmm->header.alphabet = RNA;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_COINS) == 0){
              currentPhmm->header.alphabet = coins;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_DICE) == 0){
              currentPhmm->header.alphabet = dice;
            }
            else{
              printFormatError(fileSrc, lineNumber, "expected 'amino', 'DNA', 'RNA', 'coins', or 'dice' after alphabet tag (ALPH).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_REFERENCE_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse consensus annotation tag (CONS).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, "yes") == 0){
              currentPhmm->header.hasReferenceAnnotation = true;
            }
            else if(strcmp(flagText, "no") == 0){
              currentPhmm->header.hasReferenceAnnotation = false;
            }
            else{
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after reference annotation tag (RF).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_MASK_FLAG) == 0){
            char *flagText = strtok(NULL ," ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse model mask tag (MM).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, "yes") == 0){
              currentPhmm->header.hasModelMask = true;
            }
            else if(strcmp(flagText, "no") == 0){
              currentPhmm->header.hasModelMask = false;
            }
            else{
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after model mask tag (MM).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_CONSENSUS_RESIDUE_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse consensus residue tag (CONS).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, "yes") == 0){
              currentPhmm->header.hasConsensusResidue = true;
            }
            else if(strcmp(flagText, "no") == 0){
              currentPhmm->header.hasConsensusResidue = false;
            }
            else{
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after consensus residue tag (CONS).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_CONSENSUS_STRUCTURE_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse consensus structure tag (CS).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, "yes") == 0){
              currentPhmm->header.hasConsensusStructure = true;
            }
            else if(strcmp(flagText, "no") == 0){
              currentPhmm->header.hasConsensusStructure = false;
            }
            else{
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after consensus structure tag (CS).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_MAP_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse map annotation tag (MAP).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, "yes") == 0){
              currentPhmm->header.hasMapAnnotation = true;
            }
            else if(strcmp(flagText, "no") == 0){
              currentPhmm->header.hasMapAnnotation = false;
            }
            else{
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after map annotation tag (MAP).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_DATE_FLAG) == 0){
            char *flagText = strtok(NULL, "");  //leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse date tag (DATE).");
            }
            currentPhmm->header.date = malloc(sizeof(flagText) + 1);
            if(currentPhmm->header.date == NULL){
              printAllocationError(fileSrc, lineNumber, "couldn't allocate memory for date buffer.");
              return p7HmmFormatError;
            }
            strcpy(currentPhmm->header.date, flagText);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_COMMAND_FLAG) == 0){
            char *flagText = strtok(NULL, ""); //leaving the delimeter empty goes until the string's null terminator

            size_t currentCmdHistoryLength = currentPhmm->header.commandLineHistory == NULL?
              0:  strlen(currentPhmm->header.commandLineHistory);
            size_t expandedCmdHistoryLength = currentCmdHistoryLength + strlen(flagText) + 2;
            //+2 to the new length is for the null terminator and a separating newline

            //resize the current commandLineHistory to include the
            char *expandedCmdHistory = realloc(currentPhmm->header.commandLineHistory, expandedCmdHistoryLength * sizeof(char));
            if(expandedCmdHistory == NULL){
              printAllocationError(fileSrc, lineNumber, "failed to allocate memory for command line history buffer.");
              return p7HmmAllocationFailure;
            }
            //add a newline separating the current history from the newly added line
            //if there was already a line in the history
            if(currentCmdHistoryLength != 0){
              strcat(currentPhmm->header.commandLineHistory, "\n");
            }
            strcat(currentPhmm->header.commandLineHistory, flagText);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_NSEQ_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse sequence number tag (NSEQ).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %d", &currentPhmm->header.numSequences);
            if(scanVariablesFilled < 1){
              printFormatError(fileSrc, lineNumber, "expected 1 float value after sequence number tag (NSEQ).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_EFFN_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse effective sequence number tag (EFFN).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %f", &currentPhmm->header.effectiveNumSequences);
            if(scanVariablesFilled < 1){
              printFormatError(fileSrc, lineNumber, "expected 1 float value after effective sequence number tag (EFFN).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_CHECKSUM_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse checksum tag (CKSUM).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %u", &currentPhmm->header.checksum);
            if(scanVariablesFilled < 1){
              printFormatError(fileSrc, lineNumber, " unsigned 32-bit int value is required after checksum tag.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_GATHERING_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse GA tag.");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %f %f", &currentPhmm->header.gatheringThresholds[0], &currentPhmm->header.gatheringThresholds[1]);
            if(scanVariablesFilled < 2){
              printFormatError(fileSrc, lineNumber, "expected 2 float values after GA tag.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_TRUSTED_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse TC tag.");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %f %f", &currentPhmm->header.trustedCutoffs[0], &currentPhmm->header.trustedCutoffs[1]);
            if(scanVariablesFilled != 2){
              printFormatError(fileSrc, lineNumber, "expected 2 float values after TC tag.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_NOISE_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse NC flag.");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %f %f", &currentPhmm->header.noiseCutoffs[0], &currentPhmm->header.noiseCutoffs[1]);
            if(scanVariablesFilled != 2){
              printFormatError(fileSrc, lineNumber, "expected 2 float values after NC flag.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_STATS_FLAG)){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse STATS flag.");
              return p7HmmFormatError;
            }
            char scoreDistributionName[16];
            float mu; //first value
            float lambda; //second value
            int scanVariablesFilled = sscanf(flagText, " LOCAL %s %f %f", scoreDistributionName, &mu, &lambda);
            if(scanVariablesFilled != 3){
              printFormatError(fileSrc, lineNumber,
                "expected distribution name and 2 float values after STATS.");
              return p7HmmFormatError;
            }
            if(strcmp(scoreDistributionName, "MSV") == 0){
              currentPhmm->stats.msvGumbelMu = mu;
              currentPhmm->stats.msvGumbleLambda = lambda;
            }
            else if(strcmp(scoreDistributionName, "VITERBI") == 0){
              currentPhmm->stats.viterbiGumbleMu = mu;
              currentPhmm->stats.viterbiGumbleLambda = lambda;
            }

            else if(strcmp(scoreDistributionName, "FORWARD") == 0){
              currentPhmm->stats.fowardTau = mu;
              currentPhmm->stats.forwardLambda = lambda;
            }
            else{
              printFormatError(fileSrc, lineNumber,
                "couldn't parse distribution name. exected distribution name of MSV, VITERBI, or FORWARD.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_BODY_HMM_MODEL_START_FLAG) == 0){
            parserState = parsingHmmModelHead;
            enum P7HmmReturnCode returnCode = p7HmmAllocateModelData(currentPhmm, hmmReaderGetAlphabetCardinality(currentPhmm));
            if(returnCode != p7HmmSuccess){
              printAllocationError(fileSrc, lineNumber, "failed to allocate memory for all buffers for P7 model.");
              return p7HmmAllocationFailure;
            }

            alphabetCardinality = hmmReaderGetAlphabetCardinality(currentPhmm);

            //also consume the next line of labels for the transition characters
            lineNumber++;
            size_t numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);
            if(numCharactersRead < 1){
              printAllocationError(fileSrc, lineNumber, "getline failed to create buffer.");
              return p7HmmAllocationFailure;
            }
          }
          break;

        case parsingHmmModelHead:
          if(strcmp(firstTokenLocation, P7_BODY_COMP0_FLAG) == 0){
            for(uint32_t i = 0; i < alphabetCardinality; i++){
              char *flagText = strtok(NULL, " "); //grab the next float value
              if(flagText == NULL){
                char printBuffer[256];
                sprintf(printBuffer, "Error reading value #%u from comp0 line.", i+1);
                printFormatError(fileSrc, lineNumber, printBuffer);
                return p7HmmFormatError;
              }
              int numValuesScanned = sscanf(flagText, " %f ", &currentPhmm->model.comp0[i]);
              if(numValuesScanned < 1){
                char printBuffer[256];
                sprintf(printBuffer, "Error parsing float value #%u from comp0 line.", i);
                printFormatError(fileSrc, lineNumber, printBuffer);
              }
            }

            //read the insert0 emissions line. If we didn't find the COMP0 tag, this line will still be in the lineBuffer
            lineNumber++;
            numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);
            if(numCharactersRead == -1){
              printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer.");
              return p7HmmAllocationFailure;
            }
            else if(numCharactersRead <= 2 || feof(openedFile)){
              printFormatError(fileSrc, lineNumber,
                "unexpectedly encountered end of file or blank line after hmm tag.");
              return p7HmmFormatError;
            }

            char *floatValuePtr = strtok(lineBuffer, " ");
            for(uint32_t i = 0; i < alphabetCardinality; i++){
              if(floatValuePtr == NULL){
                char printBuffer[256];
                sprintf(printBuffer, "Error reading value #%u from insert0 emissions line.", i+1);
                printFormatError(fileSrc, lineNumber, printBuffer);
                return p7HmmFormatError;
              }
              int numValuesScanned = sscanf(floatValuePtr, " %f ", &currentPhmm->model.insert0Emissions[i]);
              if(numValuesScanned < 1){
                char printBuffer[256];
                sprintf(printBuffer, "Error parsing float value #%u from comp0 line.", i);
                printFormatError(fileSrc, lineNumber, printBuffer);
              }
              floatValuePtr = strtok(NULL, " ");  //load the next value
            }

          //read and parse the the transitions from the begin state and insert state 0
          lineNumber++;
          size_t numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);
          //check to make sure getline didn't have an allocation failure
          if(numCharactersRead == -1){
            printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer.");
            return p7HmmAllocationFailure;
          }
          else if(numCharactersRead <= 2 || feof(openedFile)){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
              printFormatError(fileSrc, lineNumber, "unexpectedly encountered end of file or blank line when expecting .");
              return p7HmmFormatError;
          }
          //remove the newline if present at the final buffer position
          if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
            lineBuffer[strlen(lineBuffer) - 1] = 0;
          }


          const char *const initialTransitionsFormatString = " %f %f %f %f %f";
          int scanVariablesFilled = sscanf(lineBuffer, initialTransitionsFormatString,
            &currentPhmm->model.initialTransitions.beginToM1,
            &currentPhmm->model.initialTransitions.beginToInsert0,
            &currentPhmm->model.initialTransitions.beginToDelete1,
            &currentPhmm->model.initialTransitions.insert0ToMatch1,
            &currentPhmm->model.initialTransitions.insert0ToInsert0);

          if(scanVariablesFilled != 5){
            char errorMessageBuffer[256];
            sprintf(errorMessageBuffer,
              "expected 5 values from initial transitions line, but only got %u (there should be 7, but the last 2 are always 0.0 and *).",
              scanVariablesFilled);
            printFormatError(fileSrc, lineNumber, errorMessageBuffer);
            return p7HmmFormatError;
          }
        }
        break;


        case parsingHmmModelBody:
          if(strcmp(firstTokenLocation, P7_BODY_END_FLAG) == 0){
            //we've encountered an ending profile hmm body tag, so set the parser state and restart
            parserState = parsingHmmIdle;
            continue;
          }

          //parse the node index

          uint32_t nodeIndex = 0;
          int numItemsScanned = sscanf(firstTokenLocation, " %u ", &nodeIndex);
          if(numItemsScanned != 1){
            printFormatError(fileSrc, lineNumber,
              "Error: could not parse node index from match emissions line.");
            return p7HmmFormatError;
          }


          //check to make sure that the node index agrees with what we'd expect
          if(nodeIndex != expectedNodeIndex){
            char errorMessageBuffer[256];
            sprintf(errorMessageBuffer,
              "expected node index value of %u, but received node index value %u.",
              expectedNodeIndex, nodeIndex);
            printFormatError(fileSrc, lineNumber, errorMessageBuffer);
            return p7HmmFormatError;
          }
          expectedNodeIndex++;


          //tokenize the match emissions.
          char *tokenPointer;
          for(size_t matchEmissionIndex = 0; matchEmissionIndex < alphabetCardinality; matchEmissionIndex++){
            tokenPointer = strtok(NULL, " ");
            if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize match emission line.");
              return p7HmmFormatError;
            }
            size_t emissionScoreIndex = ((nodeIndex-1) * alphabetCardinality) + matchEmissionIndex; //-1 is to make the value zero-indexed
            numItemsScanned = sscanf(tokenPointer, " %f", &currentPhmm->model.matchEmissionScores[emissionScoreIndex]);
            if(numItemsScanned != 1){
              char errorMessageBuffer[256];
              sprintf(errorMessageBuffer,
                "Error: could not parse match emission score %zu from match emissions line.",
                matchEmissionIndex);
              printFormatError(fileSrc, lineNumber, errorMessageBuffer);
              return p7HmmFormatError;
            }
          }

          //tokenize the optional character data
          //read map annotation value
          tokenPointer = strtok(NULL, " ");
          if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize map annotation value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] == '-') && (currentPhmm->header.hasMapAnnotation)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file has map annotations, but none given on match line.");
            return p7HmmFormatError;
          }
          else if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasMapAnnotation)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have map annotations, but integer value given on match line.");
            return p7HmmFormatError;
          }
          else if(!(tokenPointer[0] == '-') && (currentPhmm->header.hasMapAnnotation)){
            numItemsScanned = sscanf(tokenPointer, " %u", &currentPhmm->model.mapAnnotations[nodeIndex]);
            if(numItemsScanned != 1){
              printFormatError(fileSrc, lineNumber,
                "Error: could not parse integer value for map annotation value.");
              return p7HmmFormatError;
            }
          }

          //read consensus residue value
          tokenPointer = strtok(NULL, "");
          if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize consensus residue value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] == '-') && (currentPhmm->header.hasConsensusResidue)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file has consensus residues, but none given on match line.");
            return p7HmmFormatError;
          }
          else if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasConsensusResidue)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have consensus residues, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          currentPhmm->model.consensusResidues[nodeIndex] = tokenPointer[0];

          //read reference annotation value
          tokenPointer = strtok(NULL, "");
          if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize reference annotation value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] == '-') && (currentPhmm->header.hasReferenceAnnotation)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file has reference annotation, but none given on match line.");
            return p7HmmFormatError;
          }
          else if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasReferenceAnnotation)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have reference annotation, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          currentPhmm->model.referenceAnnotation[nodeIndex] = tokenPointer[0];

          //read model mask value
          tokenPointer = strtok(NULL, "");
          if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize model mask value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] == '-') && (currentPhmm->header.hasModelMask)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file has reference annotation, but none given on match line.");
            return p7HmmFormatError;
          }
          else if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasModelMask)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have reference annotation, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          currentPhmm->model.modelMask[nodeIndex] = tokenPointer[0] == 'm';


          //read consensus structure value
          tokenPointer = strtok(NULL, "");
          if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize consensus structure value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] == '-') && (currentPhmm->header.hasConsensusStructure)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file has consensus structure, but none given on match line.");
            return p7HmmFormatError;
          }
          else if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasConsensusStructure)){
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have reference annotation, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          currentPhmm->model.consensusStructure[nodeIndex] = tokenPointer[0];


          //get the insert emissions line
          lineNumber++;
          numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);

          //check to make sure getline didn't have an allocation failure
          if(numCharactersRead == -1){
            printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer for insert emissions line.");
            return p7HmmAllocationFailure;
          }
          else if(numCharactersRead <= 2){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
            if(feof(openedFile)){
              if(parserState == parsingHmmIdle){
                return p7HmmSuccess;
              }
              else return p7HmmFormatError;
            }
          }

          //remove the newline if present at the final buffer position
          if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
            lineBuffer[strlen(lineBuffer) - 1] = 0;
          }

          tokenPointer = strtok(lineBuffer, " ");
          if(tokenPointer == NULL){
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize insert emission value.");
              return p7HmmFormatError;
          }
          size_t indexIntoEmissionScores = nodeIndex * alphabetCardinality;
          sscanf(tokenPointer, " %f ", &currentPhmm->model.insertEmissionScores[indexIntoEmissionScores]);
          if(numItemsScanned != 1){
            printFormatError(fileSrc, lineNumber,
              "Error: could not parse float value for insert emissions score.");
            return p7HmmFormatError;
          }

          for(size_t insertEmissionScoreIndex = 1; insertEmissionScoreIndex < alphabetCardinality; insertEmissionScoreIndex++){
            tokenPointer = strtok(NULL, " ");
            if(tokenPointer == NULL){
                printFormatError(fileSrc, lineNumber,
                  "Error: could not tokenize insert emission value.");
                return p7HmmFormatError;
            }
            sscanf(tokenPointer, " %f ", &currentPhmm->model.insertEmissionScores[indexIntoEmissionScores + insertEmissionScoreIndex]);
            if(numItemsScanned != 1){
              printFormatError(fileSrc, lineNumber,
                "Error: could not parse float value for insert emissions score.");
              return p7HmmFormatError;
            }
          }

          //get the insert emissions line
          lineNumber++;
          numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);

          //check to make sure getline didn't have an allocation failure
          if(numCharactersRead == -1){
            printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer for state transitions line.");
            return p7HmmAllocationFailure;
          }
          else if(numCharactersRead <= 2){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
            if(feof(openedFile)){
              if(parserState == parsingHmmIdle){
                return p7HmmSuccess;
              }
              else return p7HmmFormatError;
            }
          }

          //remove the newline if present at the final buffer position
          if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
            lineBuffer[strlen(lineBuffer) - 1] = 0;
          }

          numItemsScanned = sscanf(lineBuffer, " %f %f %f %f %f %f %f",
            &currentPhmm->model.stateTransitions.matchToMatch[nodeIndex],
            &currentPhmm->model.stateTransitions.matchToInsert[nodeIndex],
            &currentPhmm->model.stateTransitions.matchToDelete[nodeIndex],
            &currentPhmm->model.stateTransitions.insertToMatch[nodeIndex],
            &currentPhmm->model.stateTransitions.insertToInsert[nodeIndex],
            &currentPhmm->model.stateTransitions.deleteToMatch[nodeIndex],
            &currentPhmm->model.stateTransitions.deleteToDelete[nodeIndex]);

            if(numItemsScanned != 7){
              char errorMessageBuffer[512];
              sprintf(errorMessageBuffer,
                "Error: expected 7 values on state transitions line, but got %i.",
                numItemsScanned);
              printFormatError(fileSrc, lineNumber, errorMessageBuffer);
              return p7HmmFormatError;

            }
          break;
      }
    }

  return p7HmmSuccess;  //fallthrough condition, should not happen in practice.
}

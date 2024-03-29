#define  _POSIX_C_SOURCE 200809L     //required for the getline function
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
#define P7_BODY_COMPO_FLAG "COMPO"
#define P7_BODY_END_FLAG "//"


enum HmmReaderParserState{
  parsingHmmIdle, parsingHmmHeader, parsingHmmModelHead, parsingHmmModelBody
};


enum P7HmmReturnCode readP7Hmm(const char *const fileSrc, struct P7HmmList *phmmList){
  size_t lineBufferLength = 1 << 10; //1024
  char *lineBuffer = malloc(lineBufferLength * sizeof(char));
  if(!lineBuffer){
    printAllocationError(fileSrc, 0, "failed to allocate memory for internal line buffer.");
    return p7HmmAllocationFailure;
  }


    p7HmmListInit(phmmList);

    struct P7Hmm *currentPhmm = NULL;

    FILE *openedFile = fopen(fileSrc, "r");

    if(openedFile == NULL){
      return p7HmmFileNotFound;
    }

    enum HmmReaderParserState parserState = parsingHmmIdle;
    uint32_t alphabetCardinality = 0;
    //for counting which node number we're in when we get to the model body
    uint32_t expectedNodeIndex = 1;
    //for printing errors
    size_t lineNumber = 0;
    bool completedParsingHmm = false; //used to determine if we're valid when we hit EOF
    while(true){
      lineNumber++;
      ssize_t numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);

      //check to make sure getline didn't have an allocation failure
      if(numCharactersRead == -1){
        if(completedParsingHmm){
          free(lineBuffer);
          return p7HmmSuccess;
        }
        else{
          p7HmmListDealloc(phmmList);
          free(lineBuffer);
          printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer.");
          return p7HmmAllocationFailure;
        }
      }

      //remove the newline if present at the final buffer position
      if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
        lineBuffer[strlen(lineBuffer) - 1] = 0;
      }
      size_t fullLineLength = strlen(lineBuffer);

      char *firstTokenLocation = strtok(lineBuffer, " ");
      if(firstTokenLocation == NULL){
        //if the token is null, we can assume that this line only contained whitespace, so we should skip this line
        continue;
      }
      switch(parserState){
        case parsingHmmIdle:
        //when looking for the format tag, only check to see if it starts with 'HMMER3'
        if(strncmp(firstTokenLocation, P7_HEADER_FORMAT_FLAG, strlen(P7_HEADER_FORMAT_FLAG)) == 0){
          expectedNodeIndex = 1;
          completedParsingHmm = false;
          //we've found another header, so append a new hmm to the list
          currentPhmm = p7HmmListAppendHmm(phmmList);
          if(currentPhmm == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printAllocationError(fileSrc, lineNumber, "could not allocate memory to grow the P7ProfileHmmList list.");
            return p7HmmAllocationFailure;
          }
          currentPhmm->header.version = malloc((fullLineLength + 4)*sizeof(char)); //+4 for null terminator and extra space
          if(currentPhmm->header.version == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printAllocationError(fileSrc, lineNumber, "couldn't allocate buffer for format tag.");
            return p7HmmAllocationFailure;
          }
          currentPhmm->header.version[0] = 0; //null terminate the format since we'll be concatenating to it
          char *tokenLocation = firstTokenLocation;
          while(tokenLocation != NULL){
            strcat(currentPhmm->header.version, tokenLocation);
            strcat(currentPhmm->header.version, " ");
            tokenLocation = strtok(NULL, " ");  //grab the next word of the format tag
          }
          //this adds an extra space at the end, so setting it to a null terminator gets rid of it
          currentPhmm->header.version[strlen(currentPhmm->header.version)-1] = 0;

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
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                printAllocationError(fileSrc, lineNumber, "unalble to allocate memory for name.");
                return p7HmmAllocationFailure;
              }
              strcpy(currentPhmm->header.name, nameText);

            }
          }

          if(strcmp(firstTokenLocation, P7_HEADER_ACCESSION_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse accession number tag (ACC).");
              return p7HmmFormatError;
            }
            currentPhmm->header.accessionNumber = malloc(strlen(flagText)+1);
            if(currentPhmm->header.accessionNumber == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "couldn't allocate buffer for accession number tag (ACC).");
              return p7HmmAllocationFailure;
            }
            strcpy(currentPhmm->header.accessionNumber, flagText);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_DESCRIPTION_FLAG) == 0){
            char *flagText = strtok(NULL, "");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse description tag (DESC).");
              return p7HmmFormatError;
            }
            //fast forward past any spaces to the actual text.
            while(flagText[0] == ' '){
              flagText++;
            }
            currentPhmm->header.description = malloc(strlen(flagText)+1);
            if(currentPhmm->header.description == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "couldn't allocate buffer for description tag (DESC).");
              return p7HmmAllocationFailure;
            }
            strcpy(currentPhmm->header.description, flagText);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_LENGTH_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse model length tag (LENG).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, "%u", &currentPhmm->header.modelLength);
            if(scanVariablesFilled < 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected positive nonzero integer after model length tag (LENG).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_MAXL_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse max length tag (MAXL).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %u", &currentPhmm->header.maxLength);
            if(scanVariablesFilled == EOF || (currentPhmm->header.maxLength == 0)){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected positive nonzero integer after max length tag (MAXL).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_ALPHABET_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse alphabet tag (ALPH).");
              return p7HmmFormatError;
            }
            if(strcmp(flagText, P7_HMM_READER_ALPHABET_AMINO) == 0){
              currentPhmm->header.alphabet = P7HmmReaderAlphabetAmino;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_DNA) == 0){
              currentPhmm->header.alphabet = P7HmmReaderAlphabetDna;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_RNA) == 0){
              currentPhmm->header.alphabet = P7HmmReaderAlphabetRna;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_COINS) == 0){
              currentPhmm->header.alphabet = P7HmmReaderAlphabetCoins;
            }
            else if(strcmp(flagText, P7_HMM_READER_ALPHABET_DICE) == 0){
              currentPhmm->header.alphabet = P7HmmReaderAlphabetDice;
            }
            else{
              currentPhmm->header.alphabet = P7HmmReaderAlphabetNotSet;
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 'amino', 'DNA', 'RNA', 'coins', or 'dice' after alphabet tag (ALPH).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_REFERENCE_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after reference annotation tag (RF).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_MASK_FLAG) == 0){
            char *flagText = strtok(NULL ," ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after model mask tag (MM).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_CONSENSUS_RESIDUE_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after consensus residue tag (CONS).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_CONSENSUS_STRUCTURE_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after consensus structure tag (CS).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_MAP_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 'yes' or 'no' after map annotation tag (MAP).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_DATE_FLAG) == 0){
            char *flagText = strtok(NULL, "");  //leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              printFormatError(fileSrc, lineNumber, "couldn't parse date tag (DATE).");
            }
            //fast forward past any spaces to the actual text.
            while(flagText[0] == ' '){
              flagText++;
            }
            currentPhmm->header.date = malloc(strlen(flagText) + 1);
            if(currentPhmm->header.date == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse sequence number tag (NSEQ).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %d", &currentPhmm->header.numSequences);
            if(scanVariablesFilled < 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 1 float value after sequence number tag (NSEQ).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_EFFN_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse effective sequence number tag (EFFN).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %f", &currentPhmm->header.effectiveNumSequences);
            if(scanVariablesFilled < 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "expected 1 float value after effective sequence number tag (EFFN).");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_CHECKSUM_FLAG) == 0){
            char *flagText = strtok(NULL, " ");
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse checksum tag (CKSUM).");
              return p7HmmFormatError;
            }
            int scanVariablesFilled = sscanf(flagText, " %u", &currentPhmm->header.checksum);
            if(scanVariablesFilled < 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, " unsigned 32-bit int value is required after checksum tag.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_HEADER_GATHERING_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse GA tag.");
              return p7HmmFormatError;
            }
            //set the default values for the gathering thresholds to NAN, these will be overwritten if the GA line contained values
            currentPhmm->header.gatheringThresholds[0] = NAN;
            currentPhmm->header.gatheringThresholds[1] = NAN;
            sscanf(flagText, " %f %f", 
              &currentPhmm->header.gatheringThresholds[0], &currentPhmm->header.gatheringThresholds[1]);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_TRUSTED_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse TC tag.");
              return p7HmmFormatError;
            }
            //set the default values for the trusted cutoffs to NAN, these will be overwritten if the TC line contained values
            currentPhmm->header.trustedCutoffs[0] = NAN;
            currentPhmm->header.trustedCutoffs[1] = NAN;
            sscanf(flagText, " %f %f", 
              &currentPhmm->header.trustedCutoffs[0], &currentPhmm->header.trustedCutoffs[1]);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_NOISE_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse NC flag.");
              return p7HmmFormatError;
            }
            //set the default values for the noise cutoffs to NAN, these will be overwritten if the NC line contained values
            currentPhmm->header.noiseCutoffs[0] = NAN;
            currentPhmm->header.noiseCutoffs[1] = NAN;
            sscanf(flagText, " %f %f", &currentPhmm->header.noiseCutoffs[0], &currentPhmm->header.noiseCutoffs[1]);
          }
          if(strcmp(firstTokenLocation, P7_HEADER_STATS_FLAG) == 0){
            char *flagText = strtok(NULL, "");//leaving the delimeter empty goes until the string's null terminator
            if(flagText == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "couldn't parse STATS flag.");
              return p7HmmFormatError;
            }
            char scoreDistributionName[16];
            float mu; //first value
            float lambda; //second value
            int scanVariablesFilled = sscanf(flagText, " LOCAL %s %f %f", scoreDistributionName, &mu, &lambda);
            if(scanVariablesFilled != 3){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "expected distribution name and 2 float values after STATS.");
              return p7HmmFormatError;
            }
            if(strcmp(scoreDistributionName, "MSV") == 0){
              currentPhmm->stats.msvGumbelMu = mu;
              currentPhmm->stats.msvGumbelLambda = lambda;
            }
            else if(strcmp(scoreDistributionName, "VITERBI") == 0){
              currentPhmm->stats.viterbiGumbelMu = mu;
              currentPhmm->stats.viterbiGumbelLambda = lambda;
            }

            else if(strcmp(scoreDistributionName, "FORWARD") == 0){
              currentPhmm->stats.forwardTau = mu;
              currentPhmm->stats.forwardLambda = lambda;
            }
            else{
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "couldn't parse distribution name. exected distribution name of MSV, VITERBI, or FORWARD.");
              return p7HmmFormatError;
            }
          }
          if(strcmp(firstTokenLocation, P7_BODY_HMM_MODEL_START_FLAG) == 0){
            parserState = parsingHmmModelHead;
            enum P7HmmReturnCode returnCode = p7HmmAllocateModelData(currentPhmm);
            if(returnCode == p7HmmFormatError){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "model alphabet and/or model length was not set.");
              return p7HmmAllocationFailure;
            }
            else if(returnCode == p7HmmAllocationFailure){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "failed to allocate memory for all buffers for P7 model.");
              return p7HmmAllocationFailure;
            }

            alphabetCardinality = p7HmmGetAlphabetCardinality(currentPhmm);
            if(alphabetCardinality == 0){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "model alphabet is required at this point in the file, but was not set.");
              return p7HmmFormatError;
            }

            //also consume the next line of labels for the transition characters
            lineNumber++;
            ssize_t numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);
            if(numCharactersRead < 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "getline failed to create buffer.");
              return p7HmmAllocationFailure;
            }
          }
          break;

        case parsingHmmModelHead:
          if(strcmp(firstTokenLocation, P7_BODY_COMPO_FLAG) == 0){
            currentPhmm->model.compo = malloc(alphabetCardinality * sizeof(float));
            if(currentPhmm->model.compo == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "unable to allocate memory for COMPO array.");
              return p7HmmAllocationFailure;
            }

            for(uint32_t i = 0; i < alphabetCardinality; i++){
              char *flagText = strtok(NULL, " "); //grab the next float value
              if(flagText == NULL){
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                char printBuffer[256];
                sprintf(printBuffer, "Error reading value #%u from compo line.", i+1);
                printFormatError(fileSrc, lineNumber, printBuffer);
                return p7HmmFormatError;
              }
              int numValuesScanned = sscanf(flagText, " %f ", &currentPhmm->model.compo[i]);
              if(numValuesScanned < 1){
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                char printBuffer[256];
                sprintf(printBuffer, "Error parsing float value #%u from compo line.", i);
                printFormatError(fileSrc, lineNumber, printBuffer);
                return p7HmmFormatError;
              }
            }

            //read the insert0 emissions line. If we didn't find the COMPO tag, this line will still be in the lineBuffer
            lineNumber++;
            numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);
            if(numCharactersRead == -1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer.");
              return p7HmmAllocationFailure;
            }
            else if(numCharactersRead <= 2 || feof(openedFile)){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "unexpectedly encountered end of file or blank line after hmm tag.");
              return p7HmmFormatError;
            }

            char *floatValuePtr = strtok(lineBuffer, " ");
            for(uint32_t i = 0; i < alphabetCardinality; i++){
              if(floatValuePtr == NULL){
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                char printBuffer[256];
                sprintf(printBuffer, "Error reading value #%u from insert0 emissions line.", i+1);
                printFormatError(fileSrc, lineNumber, printBuffer);
                return p7HmmFormatError;
              }
              int numValuesScanned = sscanf(floatValuePtr, " %f ", &currentPhmm->model.insert0Emissions[i]);
              if(numValuesScanned < 1){
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                char printBuffer[256];
                sprintf(printBuffer, "Error parsing float value #%u from beginning transition probabilities line.", i);
                printFormatError(fileSrc, lineNumber, printBuffer);
                return p7HmmFormatError;
              }
              floatValuePtr = strtok(NULL, " ");  //load the next value
            }

          //read and parse the the transitions from the begin state and insert state 0
          lineNumber++;
          ssize_t numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);
          //check to make sure getline didn't have an allocation failure
          if(numCharactersRead == -1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer.");
            return p7HmmAllocationFailure;
          }
          else if(numCharactersRead <= 2 || feof(openedFile)){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            char errorMessageBuffer[256];
            sprintf(errorMessageBuffer,
              "expected 5 values from initial transitions line, but only got %u (there should be 7, but the last 2 are always 0.0 and *).",
              scanVariablesFilled);
            printFormatError(fileSrc, lineNumber, errorMessageBuffer);
            return p7HmmFormatError;
          }
          parserState = parsingHmmModelBody;
        }
        break;


        case parsingHmmModelBody:
          if(strcmp(firstTokenLocation, P7_BODY_END_FLAG) == 0){
            //we've encountered an ending profile hmm body tag, so set the parser state and restart
            completedParsingHmm = true;
            parserState = parsingHmmIdle;
            continue;
          }

          //parse the node index

          uint32_t nodeIndex = 0;
          int numItemsScanned = sscanf(firstTokenLocation, " %u ", &nodeIndex);
          if(numItemsScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: could not parse node index from match emissions line.");
            return p7HmmFormatError;
          }


          //check to make sure that the node index agrees with what we'd expect
          if(nodeIndex != expectedNodeIndex){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
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
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize match emission line.");
              return p7HmmFormatError;
            }
            size_t emissionScoreIndex = ((nodeIndex-1) * alphabetCardinality) + matchEmissionIndex; //-1 is to make the value zero-indexed
            numItemsScanned = sscanf(tokenPointer, " %f", &currentPhmm->model.matchEmissionScores[emissionScoreIndex]);
            if(numItemsScanned != 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize map annotation value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasMapAnnotation)){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have map annotations, but integer value given on match line.");
            return p7HmmFormatError;
          }
          else if(currentPhmm->header.hasMapAnnotation){
            numItemsScanned = sscanf(tokenPointer, " %u", &currentPhmm->model.mapAnnotations[nodeIndex - 1]);
            if(numItemsScanned != 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not parse integer value for map annotation value.");
              return p7HmmFormatError;
            }
          }

          //read consensus residue value
          tokenPointer = strtok(NULL, " ");
          if(tokenPointer == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize consensus residue value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] != '-') && (!currentPhmm->header.hasConsensusResidue)){
            printFormatError(fileSrc, lineNumber,
              "Warning: header declared the file does not have consensus residues, but character residue value was given on match line.");
          }
          else if(tokenPointer[0] != '-'){
            currentPhmm->model.consensusResidues[nodeIndex - 1] = tokenPointer[0];
          }

          //read reference annotation value
          tokenPointer = strtok(NULL, " ");
          if(tokenPointer == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize reference annotation value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasReferenceAnnotation)){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have reference annotation, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          else if(tokenPointer[0] != '-'){
            currentPhmm->model.referenceAnnotation[nodeIndex - 1] = tokenPointer[0];
          }

          //read model mask value
          tokenPointer = strtok(NULL, " ");
          if(tokenPointer == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize model mask value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if((tokenPointer[0] != '-') && (!currentPhmm->header.hasModelMask)){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have reference annotation, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          else if(tokenPointer[0] != '-'){
            currentPhmm->model.modelMask[nodeIndex - 1] = tokenPointer[0] == 'm';
          }

          //read consensus structure value
          tokenPointer = strtok(NULL, " ");
          if(tokenPointer == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
              printFormatError(fileSrc, lineNumber,
                "Error: could not tokenize consensus structure value");
              return p7HmmFormatError;
          }
          //check to see if the map annotation value's existance agrees with what we'd expect from hasMapAnnotation
          if(!(tokenPointer[0] == '-') && (!currentPhmm->header.hasConsensusStructure)){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: header declared the file does not have reference annotation, but character residue value was given on match line.");
            return p7HmmFormatError;
          }
          else if(currentPhmm->header.hasConsensusStructure){
            currentPhmm->model.consensusStructure[nodeIndex - 1] = tokenPointer[0];
          }


          //get the insert emissions line
          lineNumber++;
          numCharactersRead = getline(&lineBuffer, &lineBufferLength, openedFile);

          //check to make sure getline didn't have an allocation failure
          if(numCharactersRead == -1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer for insert emissions line.");
            return p7HmmAllocationFailure;
          }
          else if(numCharactersRead <= 2){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
            if(feof(openedFile)){
              free(lineBuffer);
              if(parserState == parsingHmmIdle){
                return p7HmmSuccess;
              }
              else{
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                printFormatError(fileSrc, lineNumber,
                  "file ended unexpectedly when still parsing an Hmm. Is the file missing an expected model termination flag ('//')?");
                return p7HmmFormatError;
              }
            }
          }

          //remove the newline if present at the final buffer position
          if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
            lineBuffer[strlen(lineBuffer) - 1] = 0;
          }

          tokenPointer = strtok(lineBuffer, " ");
          if(tokenPointer == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: could not tokenize insert emission value.");
              return p7HmmFormatError;
          }
          size_t indexIntoEmissionScores = (nodeIndex-1) * alphabetCardinality;
          sscanf(tokenPointer, " %f ", &currentPhmm->model.insertEmissionScores[indexIntoEmissionScores]);
          if(numItemsScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber,
              "Error: could not parse float value for insert emissions score.");
            return p7HmmFormatError;
          }

          for(size_t insertEmissionScoreIndex = 1; insertEmissionScoreIndex < alphabetCardinality; insertEmissionScoreIndex++){
            tokenPointer = strtok(NULL, " ");
            if(tokenPointer == NULL){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
                printFormatError(fileSrc, lineNumber,
                  "Error: could not tokenize insert emission value.");
                return p7HmmFormatError;
            }
            sscanf(tokenPointer, " %f ", &currentPhmm->model.insertEmissionScores[indexIntoEmissionScores + insertEmissionScoreIndex]);
            if(numItemsScanned != 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
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
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printAllocationError(fileSrc, lineNumber, "getline failed to allocate buffer for state transitions line.");
            return p7HmmAllocationFailure;
          }
          else if(numCharactersRead <= 2){  //2 is a reasonable number to distinguish a trailing newline from a line with a tag
            if(feof(openedFile)){
              free(lineBuffer);
              if(parserState == parsingHmmIdle){
                return p7HmmSuccess;
              }
              else{
                p7HmmListDealloc(phmmList);
                free(lineBuffer);
                printFormatError(fileSrc, lineNumber,
                  "unexpectedly encountered end of file while expecting an insert emissions line (2nd of set of 3 lines for each node index).");
                return p7HmmFormatError;
              }
            }
          }

          //remove the newline if present at the final buffer position
          if(lineBuffer[strlen(lineBuffer) - 1] == '\n'){
            lineBuffer[strlen(lineBuffer) - 1] = 0;
          }
          //parse out the state transition scores.
          char *tokenLocation = strtok(lineBuffer, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse match to match state transition score (1st value on state transition line).");
            return p7HmmFormatError;
          }
          int numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.matchToMatch[nodeIndex - 1]);
          if(numScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse match to match state transition score (1st value on state transition line).");
            return p7HmmFormatError;
          }

          tokenLocation = strtok(NULL, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse match to insert state transition score (2nd value on state transition line).");
            return p7HmmFormatError;
          }
          numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.matchToInsert[nodeIndex - 1]);
          if(numScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse match to insert state transition score (2nd value on state transition line).");
            return p7HmmFormatError;
          }

          tokenLocation = strtok(NULL, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse match to delete state transition score(3rd value on state transition line).");
            return p7HmmFormatError;
          }
          if(nodeIndex == currentPhmm->header.modelLength){
            //this is the last node, so this will always be '*', for infinity. since -log(INF) is undefined, we set to NAN
            currentPhmm->model.stateTransitions.matchToDelete[nodeIndex - 1] = NAN;
          }
          else{
            numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.matchToDelete[nodeIndex - 1]);
            if(numScanned != 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "failed to parse match to delete state transition score(3rd value on state transition line).");
              return p7HmmFormatError;
            }
          }

          tokenLocation = strtok(NULL, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse insert to match state transition score (4th value on state transition line).");
            return p7HmmFormatError;
          }
          numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.insertToMatch[nodeIndex - 1]);
          if(numScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse insert to match state transition score (4th value on state transition line).");
            return p7HmmFormatError;
          }

          tokenLocation = strtok(NULL, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse insert to insert state transition score (5th value on state transition line).");
            return p7HmmFormatError;
          }
          numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.insertToInsert[nodeIndex - 1]);
          if(numScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse insert to insert state transition score (5th value on state transition line).");
            return p7HmmFormatError;
          }

          tokenLocation = strtok(NULL, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse delete to match state transition score (6th value on state transition line).");
            return p7HmmFormatError;
          }
          numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.deleteToMatch[nodeIndex - 1]);
          if(numScanned != 1){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse delete to match state transition score (6th value on state transition line).");
            return p7HmmFormatError;
          }

          tokenLocation = strtok(NULL, " ");
          if(tokenLocation == NULL){
            p7HmmListDealloc(phmmList);
            free(lineBuffer);
            printFormatError(fileSrc, lineNumber, "failed to parse delete to delete state transition score (7th value on state transition line).");
            return p7HmmFormatError;
          }
          if(nodeIndex == currentPhmm->header.modelLength){
            //this is the last node, so this will always be '*', for infinity. since -log(INF) is undefined, we set to NAN
            currentPhmm->model.stateTransitions.deleteToDelete[nodeIndex - 1] = NAN;
          }
          else{
            numScanned = sscanf(tokenLocation, "%f", &currentPhmm->model.stateTransitions.deleteToDelete[nodeIndex - 1]);
            if(numScanned != 1){
              p7HmmListDealloc(phmmList);
              free(lineBuffer);
              printFormatError(fileSrc, lineNumber, "failed to parse delete to delete state transition score (7th value on state transition line).");
              return p7HmmFormatError;
            }
          }

        break;
      }
    }

  return p7HmmSuccess;  //fallthrough condition, should not happen in practice.
}

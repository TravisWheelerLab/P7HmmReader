#ifndef P7_HMM_READER_LOG_H
#define P7_HMM_READER_LOG_H

#include "p7HmmReader.h"
#include <stdint.h>

/*
 * Function:  printFormatError
 * --------------------
 * prints an error message detailing a formatting error to stderr.
 *  this is functionally equivalent to printAllocationError, other than the
 *  type of error in the output error message.
 *
 *  Inputs:
 *    fileSrc: Location of the hmm file that the error occurred in.
 *    lineNumber: line number where the formatting error took place.
 *    errormessage: text string describing the nature of the formatting error.
 */
void printFormatError(const char *const fileSrc, size_t lineNumber, char *errorMessage);

/*
 * Function:  printAllocationError
 * --------------------
 * prints an error message detailing an error due to a failed dyanamic allocation to stderr.
 *  this is functionally equivalent to printFormatError, other than the
 *  type of error in the output error message.
 *
 *  Inputs:
 *    fileSrc: Location of the hmm file that the error occurred in.
 *    lineNumber: line number where the formatting error took place.
 *    errormessage: text string describing what data element failed to allocate.
 */
void printAllocationError(const char *const fileSrc, size_t lineNumber, char *errorMessage);

#endif

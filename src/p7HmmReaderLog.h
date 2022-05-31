#ifndef P7_HMM_READER_LOG_H
#define P7_HMM_READER_LOG_H

#include "p7HmmReader.h"
#include <stdint.h>

void printFormatError(const char *const fileSrc, size_t lineNumber, char *errorMessage);
void printAllocationError(const char *const fileSrc, size_t lineNumber, char *errorMessage);

#endif

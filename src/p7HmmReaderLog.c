#include "p7HmmReaderLog.h"
#include <stdio.h>

void printFormatError(const char *const fileSrc, size_t lineNumber, char *errorMessage){
  fprintf(stderr, "\033[0;31mHmmReader format error\033[0m:File %s Line %zu.\n", fileSrc, lineNumber);
  fprintf(stderr, "\t%s\n", errorMessage);
}
void printAllocationError(const char *const fileSrc, size_t lineNumber, char *errorMessage){
  fprintf(stderr, "\033[0;31mHmmReader allocation failure\033[0m:File %s Line %zu.\n", fileSrc, lineNumber);
  fprintf(stderr, "\t%s\n", errorMessage);
}

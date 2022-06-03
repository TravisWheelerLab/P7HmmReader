TEST_NAME = printTest
MAIN_SRC = printTest.c

GCC = gcc
CFLAGS 	= -std=c11 -Wall -mtune=native -O0 -g -fPIC -fsanitize=address
OTHER_SRCS = $(wildcard ../../src/*.c)
SRCS = $(MAIN_SRC)  $(OTHER_SRCS)
TEST_BIN_NAME = $(TEST_NAME).run

printTest: $(SRC)
	$(GCC) $(CFLAGS) $(MAIN_SRC) $(OTHER_SRCS) -o $(TEST_BIN_NAME)

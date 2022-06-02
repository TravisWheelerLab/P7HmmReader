LIB_NAME = libP7HmmReader
MAJOR_VERSION = 1
MINOR_VERSION = 0
VERSION = $(MAJOR_VERSION).$(MINOR_VERSION)

CC 														= gcc
CFLAGS 												= -std=c11 -Wall -mtune=native -O3 -fPIC
LDFLAGS_SHARED_LIB 						= -shared
STATIC_LIB_FILE_EXTENSION 		= .a


PROJECT_DIR = $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
BUILD_DIR = $(PROJECT_DIR)/build
SRC_DIR = $(PROJECT_DIR)/src
BUILD_INCLUDE_DIR	= $(BUILD_DIR)/include
BUILD_LIB_DIR			=	$(BUILD_DIR)/lib

ifeq ($(PREFIX),)
PREFIX := /usr/local
endif


#determine the current operating system
OS_NAME 	= $(shell uname -s)
ifeq ($(OS_NAME), Darwin)
SHARED_LIB_FILENAME 	= $(LIB_NAME).dylib
else
SHARED_LIB_FILENAME 	= $(LIB_NAME).so
endif
STATIC_LIB_FILENAME 	= $(FASTA_VECTOR_LIB_NAME).a
GLOBAL_HEADER_FILENAME 									= p7HmmReader.h


PROJECT_HEADER_SRC 									= $(SRC_DIR)/$(GLOBAL_HEADER_FILENAME)
BUILD_HEADER_SRC 										= $(BUILD_INCLUDE_DIR)/$(GLOBAL_HEADER_FILENAME)
SHARED_LIB_BUILD_SRC 								= $(BUILD_LIB_DIR)/$(SHARED_LIB_FILENAME)
STATIC_LIB_BUILD_SRC 								= $(BUILD_LIB_DIR)/$(STATIC_LIB_FILENAME)


INSTALL_DIR 								= $(DESTDIR)$(PREFIX)
INSTALL_LIB_DIR 						= $(INSTALL_DIR)/lib
INSTALL_INCLUDE_DIR 				= $(INSTALL_DIR)/include
INSTALL_HEADER_SRC 				= $(INSTALL_INCLUDE_DIR)/$(GLOBAL_HEADER_FILENAME)
INSTALL_STATIC_LIB_SRC 		= $(INSTALL_LIB_DIR)/$(STATIC_LIB_FILENAME)
INSTALL_SHARED_LIB_SRC 		= $(INSTALL_LIB_DIR)/$(SHARED_LIB_FILENAME)

BUILD_LIB_TARGETS = $(STATIC_LIB_BUILD_SRC) $(SHARED_LIB_BUILD_SRC) $(BUILD_HEADER_SRC)

SRCS := $(shell find $(SRC_DIR) -name *.c)
OBJS := $(patsubst $(SRC_DIR)/%, $(BUILD_DIR)/%, $(SRCS:.c=.o))


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS)  -c $< -o $@


.PHONY: shared
shared: $(BUILD_INCLUDE_DIR) $(BUILD_LIB_DIR) $(OBJS)
	$(CC) -o $(SHARED_LIB_BUILD_SRC) $(LDFLAGS_SHARED_LIB) $(OBJS)
	cp $(PROJECT_HEADER_SRC) $(BUILD_HEADER_SRC)


.PHONY: static
static: $(BUILD_INCLUDE_DIR) $(BUILD_LIB_DIR) $(OBJS)
	ar rcs $(STATIC_LIB_BUILD_SRC) $(OBJS)
	cp $(PROJECT_HEADER_SRC) $(BUILD_HEADER_SRC)

.PHONY: clean
clean:
	rm -f $(BUILD_LIB_TARGETS) $(OBJS)
	$(RM)

.PHONY: install
install:$(INSTALL_LIB_DIR) $(INSTALL_INCLUDE_DIR)
	cp $(PROJECT_HEADER_SRC) $(INSTALL_HEADER_SRC)
	#copy the library files to the install src if they exist
ifneq ("$(wildcard $(STATIC_LIB_BUILD_SRC))","")
	cp $(STATIC_LIB_BUILD_SRC) $(INSTALL_STATIC_LIB_SRC)
endif
ifneq ("$(wildcard $(SHARED_LIB_BUILD_SRC))","")
	cp $(SHARED_LIB_BUILD_SRC) $(INSTALL_SHARED_LIB_SRC)
endif


.PHONY: uninstall
uninstall:
	rm -f $(INSTALL_HEADER_SRC)
	rm -f $(INSTALL_STATIC_LIB_SRC)
	rm -f $(INSTALL_SHARED_LIB_SRC)


#make the build directories if they do not exist.
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_INCLUDE_DIR): $(BUILD_DIR)
	mkdir -p $(BUILD_INCLUDE_DIR)

$(BUILD_LIB_DIR): $(BUILD_DIR)
	mkdir -p $(BUILD_LIB_DIR)

$(INSTALL_LIB_DIR):
	mkdir -p $(INSTALL_LIB_DIR)

$(INSTALL_INCLUDE_DIR):
	mkdir -p $(INSTALL_INCLUDE_DIR)
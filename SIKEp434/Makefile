# Hey Emacs, this is a -*- makefile -*-
#----------------------------------------------------------------------------
#
# Makefile for ChipWhisperer SimpleSerial-AES Program
#
#----------------------------------------------------------------------------
# On command line:
#
# make all = Make software.
#
# make clean = Clean out built project files.
#
# make coff = Convert ELF to AVR COFF.
#
# make extcoff = Convert ELF to AVR Extended COFF.
#
# make program = Download the hex file to the device, using avrdude.
#                Please customize the avrdude settings below first!
#
# make debug = Start either simulavr or avarice as specified for debugging,
#              with avr-gdb or avr-insight as the front end for debugging.
#
# make filename.s = Just compile filename.c into the assembler code only.
#
# make filename.i = Create a preprocessed source file for use in submitting
#                   bug reports to the GCC project.
#
# To rebuild project do "make clean" then "make all".
#----------------------------------------------------------------------------

# Target file name (without extension).
# This is the name of the compiled .hex file.
TARGET = sikep434

# List C source files here.
# Header files (.h) are automatically pulled in.
SRC += sikep434.c fp_generic.c P434.c fips202.c keccakf1600.c randombytes.c

# -----------------------------------------------------------------------------


FUNC_SEL = TRACE

ifneq ($(TYPE),)
  CFLAGS += -D$(TYPE)
endif


CFLAGS += -O0 -D$(FUNC_SEL) -D$(KEY)
#CFLAGS += --specs=nosys.specs

ifeq ($(CRYPTO_OPTIONS),)
CRYPTO_OPTIONS = AES128C
endif

#Add simpleserial project to build
include ./simpleserial/Makefile.simpleserial

FIRMWAREPATH = .
include $(FIRMWAREPATH)/Makefile.inc


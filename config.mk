# config.mk -- shared build configuration for miniDSP
#
# Included by the root Makefile, tests/Makefile, and examples/Makefile.

CC = gcc
CFLAGS = -std=c23 -Wall -Wextra -pedantic

# Homebrew on Apple Silicon installs to /opt/homebrew; on Intel to /usr/local.
# If Homebrew isn't present, these lines are skipped and system paths are used.
BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
ifneq ($(BREW_PREFIX),)
  CFLAGS += -I$(BREW_PREFIX)/include
  LDFLAGS += -L$(BREW_PREFIX)/lib
endif

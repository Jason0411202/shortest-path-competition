CXX      ?= c++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra

# macOS workaround: some Command Line Tools installs are missing headers in
# /Library/Developer/CommandLineTools/usr/include/c++/v1 (only a handful of
# files instead of the full ~189). Fall back to the SDK's c++/v1 when the
# default location is broken. No-op on Linux.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    ifeq ($(wildcard /Library/Developer/CommandLineTools/usr/include/c++/v1/cstdint),)
        SDK_CXX_INC := $(shell xcrun --sdk macosx --show-sdk-path 2>/dev/null)/usr/include/c++/v1
        ifneq ($(wildcard $(SDK_CXX_INC)/cstdint),)
            CXXFLAGS += -isystem $(SDK_CXX_INC)
        endif
    endif
endif

.PHONY: all foundation solver tools test check-data clean

all: foundation

# The unmodified baseline that students must run to obtain T_base.
foundation: dijkstra_foundation.cpp
	$(CXX) $(CXXFLAGS) -o foundation dijkstra_foundation.cpp

# Students rename/extend dijkstra_foundation.cpp into solver.cpp.
solver: solver.cpp
	$(CXX) $(CXXFLAGS) -o solver solver.cpp

# --- instructor / tooling targets ------------------------------------------

# The reference solver used to produce the answer keys. Multi-threaded on
# purpose: the competition's single-thread rule applies to submissions.
tools:
	$(MAKE) -C tools/ref

# Regression tests for grade.py (uncapped scoring, the wrong-answer accounting,
# the enforced limits, and baseline extrapolation vs a full run).
test: foundation
	python3 tools/test_grade.py

# Cross-check every committed instance's answer key against the unmodified
# foundation. This is the two-independent-implementations check: refsolve and
# dijkstra_foundation share no code.
check-data: foundation
	@set -e; for g in instances/*_dev.graph; do \
	  n=$${g%.graph}; \
	  ./foundation $$g $$n.queries /tmp/.chk.out; \
	  if cmp -s /tmp/.chk.out $$n.answers; then echo "ok    $$n"; \
	  else echo "WRONG $$n"; exit 1; fi; \
	done; rm -f /tmp/.chk.out; echo "all answer keys match the foundation"

clean:
	rm -f foundation solver
	$(MAKE) -C tools/ref clean

# asclepios — Abaqus-like lite FE solver
CXX      ?= clang++
OMP_DIR  ?= /opt/homebrew/opt/libomp
EIGEN    ?= /opt/homebrew/include/eigen3
CXXFLAGS ?= -O3 -march=native -std=c++17 -mmacosx-version-min=14.0 -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
CPPFLAGS += -Isrc -I$(EIGEN) -Xpreprocessor -fopenmp -I$(OMP_DIR)/include
LDFLAGS  += -L$(OMP_DIR)/lib -lomp

SRC := $(shell find src -name '*.cpp')
OBJ := $(SRC:src/%.cpp=build/%.o)
DEP := $(OBJ:.o=.d)

all: bin/asclepios

bin/asclepios: $(OBJ)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

elemtest: bin/elemtest

bin/elemtest: tools/elemtest.cpp $(filter-out build/main.o,$(OBJ))
	@mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

debug: CXXFLAGS = -O0 -g -std=c++17 -mmacosx-version-min=14.0 -Wall -fsanitize=address,undefined
debug: LDFLAGS += -fsanitize=address,undefined
debug: clean all

test: bin/asclepios
	python3 tests/run_tests.py

validate: bin/asclepios
	python3 validation/compare.py; status=$$?; python3 validation/coverage.py; exit $$status

parallel: bin/asclepios
	python3 validation/parallel.py

literature: bin/asclepios
	python3 validation/literature.py

clean:
	rm -rf build bin

-include $(DEP)
.PHONY: all clean test validate parallel literature elemtest debug

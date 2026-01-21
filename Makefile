
.PHONY: all
all: main

SOURCES = main.cpp modpoly.cpp Fp2k.cpp ec.cpp ecp.cpp isog.cpp utils.cpp conv_gmp_ntl.cpp id2iso.cpp choosetorsion.cpp costmodel.cpp interpolation.cpp crt.cpp klpt.cpp quaternions.cpp fast_quaternions.cpp fast_ff.cpp hashmap.cpp ordertojinvbigset.cpp getweber.cpp multivariates.cpp getresultants.cpp smallint.cpp mont.cpp
TESTSOURCES = tests.cpp modpoly.cpp Fp2k.cpp ec.cpp ecp.cpp isog.cpp utils.cpp conv_gmp_ntl.cpp id2iso.cpp choosetorsion.cpp costmodel.cpp interpolation.cpp crt.cpp klpt.cpp quaternions.cpp fast_quaternions.cpp fast_ff.hpp hashmap.cpp ordertojinvbigset.cpp getweber.cpp multivariates.cpp getresultants.cpp smallint.cpp mont.cpp
HEADERS = Fp2k.hpp ec.hpp ecp.hpp isog.hpp utils.hpp id2iso.hpp fast_quaternions.hpp fast_ff.hpp quaternions.hpp quatlatenum.hpp endring.hpp choosetorsion.hpp costmodel.hpp interpolation.hpp crt.hpp modpoly.hpp klpt.hpp hashmap.hpp ordertojinvbigset.hpp getweber.hpp multivariates.hpp getresultants.hpp smallint.hpp mont.hpp

CXXFLAGS = -std=c++20 -pedantic -Wall -Wextra -O3
# LDFLAGS = -lm -lntl -lgmp -lfplll -lmpfr
LDFLAGS = -lm -lntl -lgmp -lfplll -lmpfr -lprofiler

objs-main/%.o: CXXFLAGS += -march=native -DNDEBUG
objs-debug/%.o: CXXFLAGS += -g
objs-tests/%.o: CXXFLAGS += -g

objs-main/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

objs-debug/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

objs-tests/%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

main: $(patsubst %.cpp,objs-main/%.o,$(SOURCES))
	$(CXX) -DNDEBUG $^ -o $@ $(LDFLAGS)

debug: $(patsubst %.cpp,objs-debug/%.o,$(SOURCES))
	$(CXX) $^ -o $@ $(LDFLAGS) 

tests: $(patsubst %.cpp,objs-debug/%.o,$(TESTSOURCES))
	$(CXX) $^ -o $@ $(LDFLAGS) 

.PHONY: clean
clean:
	rm -f main objs-main/*
	rm -f debug objs-debug/*
	rm -f tests objs-debug/*


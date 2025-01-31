
.PHONY: all
all: main

SOURCES = main.cpp interpolation.cpp crt.cpp biglevel_hybrid.cpp
HEADERS = interpolation.hpp crt.hpp biglevel_hybrid.hpp

main: $(SOURCES) $(HEADERS)
	$(CXX) \
		-std=c++17 -pedantic \
		-Wall -Wextra \
		-O2 -march=native -DNDEBUG \
		-lntl \
		-lgmp -lm \
		$(SOURCES) \
		-o $@

debug: $(SOURCES) $(HEADERS)
	$(CXX) \
		-std=c++17 -pedantic \
		-Wall -Wextra \
		-O2 -g \
		-lntl \
		-lgmp -lm \
		$(SOURCES) \
		-o $@

.PHONY: clean
clean:
	rm -f main main-debug

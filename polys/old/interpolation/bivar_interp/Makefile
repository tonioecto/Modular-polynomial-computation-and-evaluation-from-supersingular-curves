.PHONY: all
all: interpolation 

%: %.cpp
	$(CXX) \
		-std=c++17 -pedantic \
		-Wall -Wextra \
		-O2 -march=native -DNDEBUG \
		-lntl \
		-lgmp -lm \
		$^ \
		-o $@

%-debug: %.cpp
	$(CXX) \
		-std=c++17 -pedantic \
		-Wall -Wextra \
		-O2 -g \
		-lntl \
		-lgmp -lm \
		$^ \
		-o $@

.PHONY: clean
clean:
	rm -f interpolation interpolation-debug
# Summary map
# Copyright (C) 2012 Daniel Beer <dlbeer@gmail.com>
#
# Permission to use, copy, modify, and/or distribute this software for
# any purpose with or without fee is hereby granted, provided that the
# above copyright notice and this permission notice appear in all
# copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
# AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
# DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
# PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
# TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.

CXX = g++
SM_CXXFLAGS = $(CXXFLAGS) -O2 -Wall -ggdb -std=c++0x
TESTS = \
    compile.test \
    rbt.test \
    sum.test \
    perf.test

all: $(TESTS)

test: all
	@for x in $(TESTS); do echo $$x; ./$$x > /dev/null; done

compile.test: test_compile.o
	$(CXX) -o $@ $^

rbt.test: test_rbt.o test_util.o
	$(CXX) -o $@ $^

sum.test: test_sum.o test_util.o
	$(CXX) -o $@ $^

perf.test: test_perf.o test_util.o
	$(CXX) -o $@ $^

%.o: %.cpp
	$(CXX) $(SM_CXXFLAGS) -o $*.o -c $*.cpp

clean:
	rm -f *.test
	rm -f *.o

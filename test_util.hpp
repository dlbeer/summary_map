// Summary map
// Copyright (C) 2012 Daniel Beer <dlbeer@gmail.com>
//
// Permission to use, copy, modify, and/or distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

#ifndef TEST_UTIL_HPP_
#define TEST_UTIL_HPP_

#include <vector>
#include <map>
#include <string>
#include <utility>

struct summary {
    int		sum;

    summary() : sum(0) { }

    summary(const std::pair<const std::string, int>& val) :
	sum(val.second) { }

    summary(const summary& a, const summary& b) :
	sum(a.sum + b.sum) { }
};

typedef summary_map<std::string, int, summary> test_map;

typedef std::vector<std::pair<std::string, int> > test_vector;

void load_words(const char *path, test_vector& vec, int limit = -1);
void shuffle(test_vector& vec, int seed);

#endif

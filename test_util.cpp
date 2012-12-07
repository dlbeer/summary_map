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

#include <cassert>
#include <fstream>
#include <random>

#include "summary_map.hpp"
#include "test_util.hpp"

static int x33a(const char *word)
{
    int sum = 0;

    while (*word)
	sum = sum * 33 + *(word++);

    return sum;
}

void load_words(const char *path, test_vector& vec, int limit)
{
    std::ifstream in(path);

    while (!in.eof()) {
	char buf[128];

	if (limit >= 0 && !(limit--))
	    break;

	assert(in.good());
	in.getline(buf, sizeof(buf));

	vec.push_back(std::make_pair(buf, x33a(buf)));
    }
}

void shuffle(test_vector& vec, int seed)
{
    std::mt19937 engine;

    engine.seed(seed);

    for (int i = vec.size() - 1; i > 0; i--)
	std::swap(vec[i], vec[engine() % i]);
}

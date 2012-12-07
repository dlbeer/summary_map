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

#include <iostream>
#include <algorithm>
#include <cassert>
#include <random>
#include <ctime>

#include "summary_map.hpp"
#include "test_util.hpp"

static void test_end(const char *label, clock_t start)
{
    clock_t end = clock();
    int ms = (long long)(end - start) * 1000LL / CLOCKS_PER_SEC;

    std::cout << "    " << label << ": " << ms << "ms\n";
}

template<class Map>
static void run_test(const char *label, test_vector& vec, Map& m)
{
    clock_t c;

    std::cout << label << '\n';

    shuffle(vec, 9);

    c = clock();
    m.insert(vec.begin(), vec.end());
    test_end("random insert", c);

    shuffle(vec, 10);

    c = clock();
    for (unsigned int i = 0; i < vec.size(); i++)
	m.find(vec[i].first);
    test_end("find", c);

    shuffle(vec, 11);

    c = clock();
    for (unsigned int i = 0; i < vec.size(); i++)
	m.erase(vec[i].first);
    test_end("random erase", c);

    std::sort(vec.begin(), vec.end());

    c = clock();
    m.insert(vec.begin(), vec.end());
    test_end("seq insert", c);

    c = clock();
    m.erase(m.begin(), m.end());
    test_end("seq erase", c);

    std::cout << '\n';
}

int main(int argc, char **argv)
{
    const char *path = "/usr/share/dict/words";

    if (argc > 1)
	path = argv[1];

    std::cout << "Loading words from " << path << '\n';

    test_vector vec;
    test_map m;
    std::map<std::string, int> sm;

    load_words(path, vec);
    assert(vec.size() >= 1000);
    std::cout << "Working with " << vec.size() << " items\n\n";

    run_test("summary_map", vec, m);
    run_test("std::map", vec, sm);

    return 0;
}

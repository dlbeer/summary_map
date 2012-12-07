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

#include "summary_map.hpp"
#include "test_util.hpp"

static void range_test(std::mt19937& rnd, test_map& m,
		       const test_vector& vec,
		       const std::vector<int>& series,
		       unsigned int limit)
{
    int a = rnd() % (limit + 1);
    int b = rnd() % (limit + 1);

    if (a > b)
	std::swap(a, b);

    test_map::const_iterator ia = m.find(vec[a].first);
    test_map::const_iterator ib = m.find(vec[b].first);

    const int sa = a ? series[a - 1] : 0;
    const int sb = b ? series[b - 1] : 0;

    const int vsum = sb - sa;
    const int msum = m.sum(ia, ib).sum;

    assert(vsum == msum);
}

int main(int argc, char **argv)
{
    const char *path = "/usr/share/dict/words";

    if (argc > 1)
	path = argv[1];

    std::cout << "Loading words from " << path << '\n';

    test_vector vec;
    std::vector<int> series;
    std::mt19937 rnd;
    test_map m;

    load_words(path, vec);
    assert(vec.size() >= 1000);
    std::cout << "Working with " << vec.size() << " items\n";

    std::sort(vec.begin(), vec.end());
    rnd.seed(5);

    int all_sum = 0;

    for (unsigned int i = 0; i < vec.size(); i++) {
	all_sum += vec[i].second;
	series.push_back(all_sum);
    }

    std::cout << "Inserting...\n";

    for (unsigned int i = 0; i < vec.size(); i++) {
	m.insert(vec[i]);
	range_test(rnd, m, vec, series, i + 1);
	assert(m.sum().sum == series[i]);
    }

    std::cout << "Sum: " << m.sum().sum << '\n';

    std::cout << "Removing...\n";

    for (int i = vec.size() - 1; i >= 0; i--) {
	assert(m.sum().sum == series[i]);
	m.erase(vec[i].first);
	range_test(rnd, m, vec, series, i);
    }

    assert(m.sum().sum == 0);
    return 0;
}

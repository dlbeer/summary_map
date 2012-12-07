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

#define SUMMARY_MAP_DEBUG
#include "summary_map.hpp"

#include "test_util.hpp"

static void check_present(const test_map& m, const test_vector& vec,
			  unsigned int start, unsigned int end)
{
    for (unsigned int i = start; i < end; i++) {
	test_map::const_iterator it = m.find(vec[i].first);

	assert(it != m.end());
	assert(it->first == vec[i].first);
	assert(it->second == vec[i].second);
    }
}

static void check_absent(const test_map& m, const test_vector& vec,
			 unsigned int start, unsigned int end)
{
    for (unsigned int i = start; i < end; i++) {
	test_map::const_iterator it = m.find(vec[i].first);

	assert(it == m.end());
    }
}

int main(int argc, char **argv)
{
    const char *path = "/usr/share/dict/words";

    if (argc > 1)
	path = argv[1];

    std::cout << "Loading words from " << path << '\n';

    test_vector vec;
    test_map m;

    load_words(path, vec, 1000);
    assert(vec.size() == 1000);
    check_absent(m, vec, 0, vec.size());

    std::cout << "Inserting...\n";
    shuffle(vec, 1);

    for (unsigned int i = 0; i < vec.size(); i++) {
	m.insert(vec[i]);
	assert(m.size() == i + 1);
	m.check();
    }

    check_present(m, vec, 0, vec.size());

    std::cout << "Removing some...\n";
    shuffle(vec, 2);

    for (int i = vec.size() - 1; i >= (signed)vec.size() / 2; i--) {
	m.erase(vec[i].first);
	assert(m.size() == (unsigned)i);
	m.check();
    }

    check_present(m, vec, 0, vec.size() / 2);
    check_absent(m, vec, vec.size() / 2, vec.size());

    std::cout << "Adding some...\n";

    for (unsigned int i = vec.size() / 2; i < vec.size(); i++) {
	m.insert(vec[i]);
	assert(m.size() == i + 1);
	m.check();
    }

    check_present(m, vec, 0, vec.size());

    std::cout << "Removing...\n";
    shuffle(vec, 3);

    for (int i = vec.size() - 1; i >= 0; i--) {
	m.erase(vec[i].first);
	assert(m.size() == (unsigned)i);
	m.check();
    }

    check_absent(m, vec, 0, vec.size());

    return 0;
}

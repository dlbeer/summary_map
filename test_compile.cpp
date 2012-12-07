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

#include "summary_map.hpp"

#include <map>
#include <string>

#include "test_util.hpp"

#ifdef SUMMARY_MAP_NO_CPP11
static void test_cpp11(test_map& m) { }
#else
static void test_cpp11(test_map& m)
{
    test_map other{{"foo", 5}, {"bar", 7}};
    test_map foo(std::move(other));

    m = std::move(foo);
}
#endif

static void test_mutable(test_map& m)
{
    std::map<std::string, int> source;
    test_map other;

    m.begin();
    m.clear();
    m.end();
    m.equal_range("");

    m = other;
    m[""] = 3;

    m.erase(m.begin());
    m.erase("");
    m.erase(m.begin(), m.end());

    m.find("");

    m.insert(std::make_pair("", 0));
    m.insert(m.end(), std::make_pair("", 0));
    m.insert(source.begin(), source.end());

    m.lower_bound("");

    m.rbegin();
    m.rend();
    m.swap(other);
    m.upper_bound("");

    m.sum();
    m.sum(m.begin(), m.end());

    swap(m, other);
}

static void test_const(test_map& m)
{
    m.begin();
    m.count("");
    m.empty();
    m.end();
    m.equal_range("");
    m.find("");
    m.get_allocator();
    m.key_comp();
    m.lower_bound("");
    m.max_size();
    m.rbegin();
    m.rend();
    m.size();
    m.upper_bound("");
    m.value_comp();
}

struct int_aggregator {
    int nothing() const
    {
	return 0;
    }

    int summarize(const std::pair<const std::string, int>& v)
    {
	return v.second;
    }

    int combine(int a, int b)
    {
	return a + b;
    }
};

int main()
{
    test_map m;
    summary_map<std::string, int, int,
		std::less<std::string>, int_aggregator> custom_map;

    test_mutable(m);
    test_const(m);
    test_cpp11(m);

    return 0;
}

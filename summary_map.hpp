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

#ifndef SUMMARY_MAP_HPP_
#define SUMMARY_MAP_HPP_

#include <functional>
#include <utility>
#include <limits>
#include <memory>

#ifdef SUMMARY_MAP_DEBUG
#include <cassert>
#endif

// Default aggregator. It assumes that the summary type contains three
// constructors:
//
//   - one which takes no value (summary over no values)
//   - one which takes a value (std::pair<k, d>, summary over one value)
//   - one which takes two other summaries (associative combinator)

template <class Value, class Summary>
struct aggregator {
    typedef Value	value_type;
    typedef Summary	summary_type;

    summary_type nothing() const
    {
	return summary_type();
    }

    summary_type summarize(const value_type& v) const
    {
	return summary_type(v);
    }

    summary_type combine(const summary_type& a, const summary_type& b) const
    {
	return summary_type(a, b);
    }
};

// Summary map.

template <class Key, class Data, class Summary,
	  class Compare = std::less<Key>,
	  class Aggregate = aggregator<std::pair<const Key, Data>, Summary>,
	  class Allocator = std::allocator<std::pair<const Key, Data> > >
class summary_map {
public:
    // Basic typedefs for all related types
    typedef Key				key_type;
    typedef Data			mapped_type;
    typedef Summary			summary_type;
    typedef Compare			key_compare;
    typedef Aggregate			aggregator_type;

    typedef std::pair<const Key, Data>	value_type;

    class value_compare :
	public std::binary_function<value_type, value_type, bool> {
    private:
	key_compare	key_cmp;

    public:
	value_compare() { }
	value_compare(const key_compare& k) : key_cmp(k) { }

	bool operator()(const value_type& a, const value_type& b) const
	{
	    return key_cmp(a.first, b.first);
	}
    };

    typedef Allocator				allocator_type;
    typedef typename Allocator::reference	reference;
    typedef typename Allocator::const_reference	const_reference;

    typedef size_t				size_type;
    typedef ptrdiff_t				difference_type;

    typedef typename Allocator::pointer		pointer;
    typedef typename Allocator::const_pointer	const_pointer;

private:
    struct node {
	static const int RED = 0x01;
	static const int SUM_VALID = 0x02;

	value_type	value;

	summary_type	sum;
	int		flags;

	node		*left;
	node		*right;
	node		*parent;

	node() :
	    flags(0), left(0), right(0), parent(0) { }
	node(const value_type& val) :
	    value(val), flags(0), left(0), right(0), parent(0) { }
    };

    typedef typename Allocator::template rebind<node>::other node_allocator;

    key_compare		comparer;
    aggregator_type	aggregator;
    node_allocator	allocator;
    node		header;
    size_type		item_count;

    // Construct/destroy a node using the given allocator
    static node *new_node(node_allocator& alloc, const value_type& val)
    {
	node *n = alloc.allocate(1);

	try {
	    n = new (n) node(val);
	}
	catch (...) {
	    alloc.deallocate(n, 1);
	    throw;
	}

	return n;
    }

    static void delete_node(node_allocator& alloc, node *n)
    {
	n->~node();
	alloc.deallocate(n, 1);
    }

    // Delete a tree and all of its children.
    static void delete_tree(node_allocator& alloc, node *n)
    {
	if (!n)
	    return;

	delete_tree(alloc, n->left);
	delete_tree(alloc, n->right);
	delete_node(alloc, n);
    }

    // Copy a subtree in an exception-safe fashion.
    static node *copy_subtree(node_allocator& alloc, node *new_parent,
			      node *src)
    {
	if (!src)
	    return 0;

	node *copy = new_node(alloc, src->value);

	try {
	    copy->flags = src->flags;
	    copy->sum = src->sum;
	    copy->parent = new_parent;
	    copy->left = copy_subtree(alloc, copy, src->left);
	    copy->right = copy_subtree(alloc, copy, src->right);
	} catch (...) {
	    delete_tree(alloc, copy);
	    throw;
	}

	return copy;
    }

    // Find the left-most and right-most children in a node.
    static node *leftmost_child(node *n)
    {
	if (!n)
	    return n;

	while (n->left)
	    n = n->left;

	return n;
    }

    static node *rightmost_child(node *n)
    {
	if (!n)
	    return n;

	while (n->right)
	    n = n->right;

	return n;
    }

    // Structural utilities. These perform no sanity checking, and
    // assume that the node is of sufficient depth.
    static bool is_header(node *n)
    {
	return !n->parent;
    }

    static bool is_root(node *n)
    {
	return is_header(n->parent);
    }

    static bool is_left_child(node *n)
    {
	return n == n->parent->left;
    }

    static bool is_right_child(node *n)
    {
	return n == n->parent->right;
    }

    static node *sibling(node *n)
    {
	node *p = n->parent;

	return (p->left == n) ? p->right : p->left;
    }

    static node *grandparent(node *n)
    {
	return n->parent->parent;
    }

    // Find the successor of this node. The successor of the header node
    // is itself.
    static node *next_in_order(node *n)
    {
	// If we have a right subtree, return the first item in it
	if (n->right)
	    return leftmost_child(n->right);

	// Otherwise, return the first ancestor from which we descended
	// left, or return the header as an end-of-sequence marker.
	while (!is_root(n) && is_right_child(n))
	    n = n->parent;

	// This will return the header if we ascended to the root in the
	// previous step.
	return n->parent;
    }

    // Find the predecessor of this node. The predecessor of the first
    // node is itself.
    static node *prev_in_order(node *n)
    {
	// If we have a left subtree, return the first item in it.
	// Happily, this is the correct thing to do even in the case of
	// the header node (which points left to the root of the real
	// tree).
	if (n->left)
	    return rightmost_child(n->right);

	// Otherwise, find the first ancestor from which we descended
	// right.
	while (!is_root(n) && is_left_child(n))
	    n = n->parent;

	return n->parent;
    }

    // Search for the node whose key equals the given one. If no such
    // node is found, the header is returned.
    static node *find_eq(const key_compare& cmp, const key_type& target,
			 node *header)
    {
	node *n = header->left;

	while (n) {
	    if (cmp(target, n->value.first))
		n = n->left;
	    else if (cmp(n->value.first, target))
		n = n->right;
	    else
		return n;
	}

	return header;
    }

    // Search for the first node whose key is greater than or equal to
    // the given one. If no such node is found, the header is returned.
    static node *find_ge(const key_compare& cmp, const key_type& target,
			 node *header)
    {
	node *n = header->left;
	node *best = header;

	while (n) {
	    if (cmp(n->value.first, target)) {
		n = n->right;
	    } else if (cmp(target, n->value.first)) {
		best = n;
		n = n->left;
	    } else {
		return n;
	    }
	}

	return best;
    }

    // Search for the last node whose key is less then or equal to the
    // given one. If no such node is found, the header is returned.
    static node *find_le(const key_compare& cmp, const key_type& target,
			 node *header)
    {
	node *n = header->left;
	node *best = header;

	while (n) {
	    if (cmp(target, n->value.first)) {
		n = n->left;
	    } else if (cmp(n->value.first, target)) {
		best = n;
		n = n->right;
	    } else {
		return n;
	    }
	}

	return best;
    }

    // Is this a black/red node? (include's placeholder leaves)
    static bool is_black(const node *n)
    {
	return !(n && (n->flags & node::RED));
    }

    static bool is_red(const node *n)
    {
	return n && (n->flags & node::RED);
    }

    static void make_black(node *n)
    {
	n->flags &= ~node::RED;
    }

    static void make_red(node *n)
    {
	n->flags |= node::RED;
    }

    static void invalidate(node *n)
    {
	n->flags &= ~node::SUM_VALID;
    }

    // Invalidate n and its anscestors
    static void invalidate_up(node *n)
    {
	for (node *p = n; !is_header(p); p = p->parent)
	    invalidate(p);
    }

    // Fix the downwards pointer from p to old by setting it to n.
    static void fix_downptr(node *p, node *old, node *n)
    {
	if (p->left == old)
	    p->left = n;
	else
	    p->right = n;
    }

    /* Rotating a node (n) left invalidates the summary for (i.e. changes
     * the set of descendants of):
     *
     *    - n itself
     *    - n's originally right child
     *
     *    B                 D
     *   / \       ==>     / \
     *  A   D             B   E
     *     / \           / \
     *    C   E         A   C
     */
    static void rotate_left(node *n)
    {
	node *p = n->parent; // p is above [B]
	node *r = n->right;  // r is [D]

	n->right = r->left;  // [B] now points right to [C]
	if (n->right)
	    n->right->parent = n;

	r->left = n;         // [D] now points left to [B]
	n->parent = r;

	r->parent = p;       // Make [D] the root of the subtree
	fix_downptr(p, n, r);
    }


    /* Rotating a node (n) right invalidates the summary for (i.e. changes
     * the set of descendants of):
     *
     *    - n itself
     *    - n's originally left child
     *
     *      D                 B
     *     / \       ==>     / \
     *    B   E             A   D
     *   / \                   / \
     *  A   C                 C   E
     */
    static void rotate_right(node *n)
    {
	node *p = n->parent; // p is above [D]
	node *l = n->left;   // l is [B]

	n->left = l->right;  // [D] now points left to [C]
	if (n->left)
	    n->left->parent = n;

	l->right = n;        // [B] now points right to [D]
	n->parent = l;

	l->parent = p;       // Make [B] the root of the subtree
	fix_downptr(p, n, l);
    }

    // Rebalance the tree after having added the given node. This
    // function modifies the tree, but only the flags and pointers. It
    // does not throw any exceptions.
    static void rebalance_after_insert(node *n)
    {
	make_red(n);

	// Because the current node is now red, it's possible that we
	// have violated the property that red nodes have only black
	// children (because n might have a red parent). Go up from the
	// newly inserted node and attempt to fix the property by
	// recolouring.
	//
	// No structural changes occur in this phase, so no summaries
	// are invalidated.
	for (;;) {
	    // If we've hit the root, we're done -- just make sure it's
	    // black.
	    if (is_root(n)) {
		make_black(n);
		return;
	    }

	    // If the parent is black, everything's fine.
	    if (is_black(n->parent))
		return;

	    // We definitely have grandparent, because our parent is
	    // red (and therefore is not the root). If our uncle is also
	    // red, we can recolour and propagate upwards. Otherwise (if
	    // the uncle is black), we need to continue to the next
	    // phase of fix-ups.
	    node *uncle = sibling(n->parent);
	    if (is_black(uncle))
		break;

	    // Recolour as described above and continue upwards.
	    make_black(n->parent);
	    make_black(uncle);
	    n = grandparent(n);
	    make_red(n);
	}

	// If we got to this stage, we have a grandparent, our parent is
	// red, and our uncle is black.
	//
	// First, if the node and its parent descend in different
	// directions (zig-zag ancestry), then rotate the parent in such
	// a way as to straighten the path.
	//
	// We take, as our new node, the parent which has now been
	// pushed down (it is red). Note that since both the initial
	// node and the parent were both red, we haven't affected the
	// black-count.
	//
	// Also note that the summaries which would be invalidated by
	// this rotation have already been invalidated by the insertion
	// operation.
	if (is_left_child(n) && is_right_child(n->parent)) {
	    rotate_right(n->parent);
	    n = n->right;
	} else if (is_right_child(n) && is_left_child(n->parent)) {
	    rotate_left(n->parent);
	    n = n->left;
	}

	// The node has straight ancestry to its parent. If we now
	// recolour and rotate the grandparent, we will have balanced
	// the tree.
	make_black(n->parent);
	make_red(grandparent(n));

	if (is_left_child(n))
	    rotate_right(grandparent(n));
	else
	    rotate_left(grandparent(n));
    }

    // Insert the given value into the tree. If necessary, a new node
    // will be created and the tree will be rebalanced. Returns a
    // pointer to the new or existing node, and a boolean telling us
    // whether a node was created.
    static std::pair<node *, bool>
    add_node(const key_compare& cmp, node_allocator& alloc,
	     node *header, const value_type& val)
    {
	node **nptr = &header->left;
	node *p = header;

	// Keep descending until we find a null down-pointer
	while (*nptr) {
	    node *c = *nptr;

	    invalidate(c);

	    if (cmp(val.first, c->value.first))
		nptr = &c->left;
	    else if (cmp(c->value.first, val.first))
		nptr = &c->right;
	    else
		return std::make_pair(c, false);

	    p = c;
	}

	// Attach a new node to the tree
	node *c = new_node(alloc, val);

	*nptr = c;
	c->parent = p;

	// Perform rebalancing
	rebalance_after_insert(c);

	return std::make_pair(c, true);
    }

    // Swap the position of the given node with the lexicographically
    // smallest child in the right subtree. Colours are also exchanged,
    // so as to preserve the red-black tree properties. This invalidates
    // summaries from n to its successor.
    //
    // This must only ever be called on a node with two children.
    static void swap_with_successor(node *n)
    {
	node *p = n->parent;
	node *s = leftmost_child(n->right);

	if (s->parent == n) {
	    // This is a tricky case: we can't just swap pointers,
	    // because we'd end up with a cycle.
	    node *left = n->left;
	    node *right = s->right;

	    s->left = left;
	    s->right = n;
	    s->parent = p;

	    n->right = right;
	    n->parent = s;

	    left->parent = s;
	} else {
	    s->left = n->left;

	    std::swap(n->right, s->right);
	    std::swap(n->parent, s->parent);

	    s->left->parent = s;
	    s->right->parent = s;

	    n->parent->left = n;
	}

	std::swap(n->flags, s->flags);
	n->left = 0;

	fix_downptr(p, n, s);

	if (n->right)
	    n->right->parent = n;
    }

    static void rebalance_after_remove(node *p)
    {
	node *n = 0;
	node *s;

	for (;;) {
	    // If we've hit the root, we're done
	    if (is_header(p))
		return;

	    // At each step here, we have that:
	    //
	    //   - n is black
	    //   - s is non-NULL
	    //   - the black count through n is one less than the count
	    //     through s
	    //
	    // Note that we can't use the sibling(n) to compute s here,
	    // because n might be null.
	    s = (p->left == n) ? p->right : p->left;

	    // If s is red, then p is black and s's children are also
	    // black. Swap colours of p and s, then rotate p towards us.
	    //
	    // This doesn't affect the black count through either path,
	    // but ensures that we have a black parent. It also
	    // invalidates the summary of s (rotation rules).
	    if (is_red(s)) {
		make_red(p);
		make_black(s);
		invalidate(s);

		if (is_right_child(s)) {
		    rotate_left(p);
		    s = p->right;
		} else {
		    rotate_right(p);
		    s = p->left;
		}
	    }

	    // Otherwise, if s is able to become red, colour it so. We
	    // have then balanced p's child subtrees. However, p will
	    // now have a count too low with respect to p's sibling, and
	    // we will have propagated the error upwards.
	    if (is_black(p) && is_black(s->left) && is_black(s->right)) {
		make_red(s);
		n = p;
		p = n->parent;
		continue;
	    }

	    // Otherwise, we will have to fix this by rotations.
	    break;
	}

	// If s and its children are black, but the parent is red, we
	// can just swap colours. This doesn't affect the count through
	// s, but it increases the count through n, this fixing the
	// tree.
	if (is_red(p) && is_black(s->left) && is_black(s->right)) {
	    make_red(s);
	    make_black(p);
	    return;
	}

	// We now know that s has at least one red child. Perform a
	// rotation to ensure that s's red child (or one of them)
	// descends in a different direction to n.
	//
	// We do this, if necessary, by a recolouring and rotation which
	// does not affect black counts.
	if (is_right_child(s) && is_black(s->right)) {
	    make_red(s);
	    make_black(s->left);
	    invalidate(s);
	    invalidate(s->left);
	    rotate_right(s);
	    s = s->parent;
	} else if (is_left_child(s) && is_black(s->left)) {
	    make_red(s);
	    make_black(s->right);
	    invalidate(s);
	    invalidate(s->right);
	    rotate_left(s);
	    s = s->parent;
	}

	// n now has a sibling with a red child whose descent is
	// opposite to that of n's. We don't know what colour p is.
	//
	// Make p black, and then rotate it towards n. Make s take over
	// the colour of p. This increases the black count through n by
	// 1.
	//
	// The far child has had its black parent cut from the path, so
	// recolour it black to balance the count.
	//
	// s is invalidated by this step (p is already invalidated, so
	// we can just copy its flags).
	s->flags = p->flags;
	make_black(p);

	if (is_right_child(s)) {
	    make_black(s->right);
	    rotate_left(p);
	} else {
	    make_black(s->left);
	    rotate_right(p);
	}
    }

    // Unlink the given node from the tree and rebalance. The node is
    // not deallocated -- you need to do this yourself.
    static void unlink_node(node *n)
    {
	// If we have two children, swap and reduce the problem to that
	// of removing a node with a single child. After removing n,
	// ordering will be preserved.
	if (n->left && n->right)
	    swap_with_successor(n);

	invalidate_up(n);

	// Replace the node with its child (if any).
	node *s = n->left ? n->left : n->right;
	fix_downptr(n->parent, n, s);

	// If n has a child, then n is black and s is red (otherwise it
	// wouldn't be possible to be black-balanced with one non-null
	// child).
	//
	// In this case, we can just recolour it and return immediately.
	if (s) {
	    s->parent = n->parent;
	    make_black(s);
	    return;
	}

	// If we removed a red node (which had no children), or the
	// root, we haven't altered the black balance.
	if (is_red(n))
	    return;

	// Otherwise, we have things to fix
	rebalance_after_remove(n->parent);
    }

    // Return the summary of the given subtree. Returns either the
    // cached version, or computes it recursively (and caches the
    // result).
    static summary_type sum_tree(const aggregator_type& ag, node *n)
    {
	if (!n)
	    return ag.nothing();

	if (!(n->flags & node::SUM_VALID)) {
	    n->sum =
		ag.combine(
		    ag.combine(sum_tree(ag, n->left), ag.summarize(n->value)),
		    sum_tree(ag, n->right));
	    n->flags |= node::SUM_VALID;
	}

	return n->sum;
    }

    // Determine the depth of a node. The root node has a depth of 1.
    static int depth_of(node *n)
    {
	int r = 0;

	while (!is_header(n)) {
	    r++;
	    n = n->parent;
	}

	return r;
    }

    // Determine the lowest common ancestor of two nodes.
    static node *concestor(node *a, node *b)
    {
	int da = depth_of(a);
	int db = depth_of(b);

	// Make sure b is the deeper of the two
	if (da > db) {
	    std::swap(da, db);
	    std::swap(a, b);
	}

	// Climb up b until we find an ancestor of the same depth as a
	while (da < db) {
	    b = b->parent;
	    db--;
	}

	// Climb both simultaneously until we meet
	while (a != b) {
	    a = a->parent;
	    b = b->parent;
	}

	return a;
    }

    // Return the summary over all nodes strictly between start and
    // root. Root comes lexicographically after start, and is also its
    // ancestor.
    static summary_type sum_from(const aggregator_type& ag,
				 node *start, node *root)
    {
	summary_type sum = ag.nothing();
	node *n = start;

	while (n != root) {
	    if (n != start)
		sum = ag.combine(sum, n->value);

	    sum = ag.combine(sum, sum_tree(ag, n->right));

	    // Ascend to an ancestor which is lexicographically after n.
	    while (is_right_child(n) && n->parent != root)
		n = n->parent;

	    n = n->parent;
	}

	return sum;
    }

    // Return the sum over all nodes strictly between root and end. Root
    // comes lexicographically before end, and is also its ancestor.
    static summary_type sum_to(const aggregator_type& ag,
			       node *end, node *root)
    {
	summary_type sum = ag.nothing();
	node *n = end;

	while (n != root) {
	    if (n != end)
		sum = ag.combine(n->value, sum);

	    sum = ag.combine(sum_tree(ag, n->left), sum);

	    // Ascend to an ancestor which is lexicographically before
	    // n.
	    while (is_left_child(n) && n->parent != root)
		n = n->parent;

	    n = n->parent;
	}

	return sum;
    }

    // Return the summary between the given nodes (including the first,
    // but not including the last). We assume that first <= last.
    static summary_type sum_between(const aggregator_type& ag,
				    node *first, node *last)
    {
	if (is_header(first))
	    return ag.nothing();

	// This special case works, because last is the parent of the
	// root. The sum will include the root's right subtree iff the
	// path to first is left of the root
	if (is_header(last))
	    return ag.combine(ag.summarize(first->value),
			      sum_from(ag, first, last));

	if (first == last)
	    return ag.nothing();

	node *c = concestor(first, last);

	summary_type sum = ag.summarize(first->value);

	sum = ag.combine(sum, sum_from(ag, first, c));

	if ((c != last) && (c != first))
	    sum = ag.combine(sum, c->value);

	sum = ag.combine(sum, sum_to(ag, last, c));

	return sum;
    }

#ifdef SUMMARY_MAP_DEBUG
    // Recursive invariant check. The function asserts on each invariant
    // that the red-black tree is supposed to satisfy, and returns a
    // pair of counts: the total number of nodes in the subtree, and the
    // black-height of the subtree.
    static std::pair<size_type, size_type>
    check_subtree(const key_compare& cmp,
		  const node *parent, const node *n,
		  const key_type *lower_bound,
		  const key_type *upper_bound)
    {
	if (!n)
	    return std::make_pair(0, 0);

	// Check that the parent pointer points to the node we came
	// from.
	assert(n->parent == parent);

	// Check that the key is within the expected range
	assert(!upper_bound || cmp(n->value.first, *upper_bound));
	assert(!lower_bound || cmp(*lower_bound, n->value.first));

	// Red nodes always have black parents
	assert(is_black(n) || is_black(parent));

	// Recursively check subtrees
	std::pair<size_type, size_type> left =
	    check_subtree(cmp, n, n->left, lower_bound, &n->value.first);
	std::pair<size_type, size_type> right =
	    check_subtree(cmp, n, n->right, &n->value.first, upper_bound);

	// Verify that this subtree is black-balanced
	assert(left.second == right.second);

	// Compile totals for this subtree
	const size_type black_self = is_black(n) ? 1 : 0;

	return std::make_pair(left.first + right.first + 1,
			      left.second + black_self);
    }
#endif

public:
    // Basic constructors
    summary_map() : item_count(0) { }

    summary_map(const key_compare& cmp,
		const aggregator_type& ag,
		const node_allocator& alloc) :
	comparer(cmp), aggregator(ag), allocator(alloc), item_count(0) { }

    ~summary_map()
    {
	delete_tree(allocator, header.left);
    }

#ifndef SUMMARY_MAP_NO_CPP11
    // Move construction and assignment
    summary_map(summary_map&& r) : item_count(0)
    {
	swap(r);
    }

    summary_map& operator=(summary_map&& r)
    {
	summary_map m(std::move(r));

	swap(m);
	return *this;
    }

    summary_map(std::initializer_list<value_type> init) :
	item_count(0)
    {
	try {
	    insert(init.begin(), init.end());
	}
	catch (...) {
	    delete_tree(allocator, header.left);
	    throw;
	}
    }
#endif

    // Copy construction and assignment
    summary_map(const summary_map& r) :
	comparer(r.comparer), aggregator(r.aggregator), allocator(r.allocator)
    {
	header.left = copy_subtree(allocator, &header, r.header.left);
	item_count = r.item_count;
    }

    summary_map& operator=(const summary_map& r)
    {
	summary_map copy(r);

	swap(copy);
	return *this;
    }

    void swap(summary_map& r)
    {
	using std::swap;

	swap(comparer, r.comparer);
	swap(aggregator, r.aggregator);
	swap(allocator, r.allocator);

	std::swap(item_count, r.item_count);
	std::swap(header.left, r.header.left);

	if (header.left)
	    header.left->parent = &header;
	if (r.header.left)
	    r.header.left->parent = &r.header;
    }

    // Access to function objects
    key_compare key_comp() const
    {
	return comparer;
    }

    value_compare value_comp() const
    {
	return value_compare(key_comp());
    }

    allocator_type get_allocator() const
    {
	return allocator_type(allocator);
    }

    aggregator_type get_aggregator() const
    {
	return aggregator;
    }

    // Size, clear/empty
    void clear()
    {
	delete_tree(allocator, header.left);
	header.left = 0;
	item_count = 0;
    }

    bool empty() const
    {
	return item_count > 0;
    }

    size_type size() const
    {
	return item_count;
    }

    size_type max_size() const
    {
	return std::numeric_limits<size_type>::max();
    }

    // Iterator types
    struct iterator {
	typedef std::pair<const Key, Data>	value_type;
	typedef value_type&			reference;
	typedef value_type			*pointer;

	typedef ptrdiff_t		difference_type;

	typedef std::bidirectional_iterator_tag iterator_category;

	node		*m_cur;

	iterator() : m_cur(0) { }
	iterator(node *n) : m_cur(n) { }

	reference operator*() const
	{
	    invalidate_up(m_cur);
	    return m_cur->value;
	}

	pointer operator->() const
	{
	    invalidate_up(m_cur);
	    return &m_cur->value;
	}

	bool operator==(const iterator& r) const
	{
	    return m_cur == r.m_cur;
	}

	bool operator!=(const iterator& r) const
	{
	    return m_cur != r.m_cur;
	}

	iterator& operator--()
	{
	    m_cur = prev_in_order(m_cur);
	    return *this;
	}

	iterator& operator++()
	{
	    m_cur = next_in_order(m_cur);
	    return *this;
	}

	iterator operator--(int)
	{
	    iterator tmp = *this;

	    m_cur = prev_in_order(m_cur);
	    return tmp;
	}

	iterator operator++(int)
	{
	    iterator tmp = *this;

	    m_cur = next_in_order(m_cur);
	    return tmp;
	}
    };

    // Constant iterator type
    struct const_iterator {
	typedef std::pair<const Key, Data>	value_type;
	typedef const value_type&		reference;
	typedef const value_type		*pointer;

	typedef ptrdiff_t		difference_type;

	typedef std::bidirectional_iterator_tag iterator_category;

	const node		*m_cur;

	const_iterator() : m_cur(0) { }
	const_iterator(iterator r) : m_cur(r.m_cur) { }
	const_iterator(const node *n) : m_cur(n) { }

	reference operator*() const
	{
	    return m_cur->value;
	}

	pointer operator->() const
	{
	    return &m_cur->value;
	}

	bool operator==(const const_iterator& r) const
	{
	    return m_cur == r.m_cur;
	}

	bool operator!=(const const_iterator& r) const
	{
	    return m_cur != r.m_cur;
	}

	const_iterator& operator--()
	{
	    m_cur = prev_in_order((node *)m_cur);
	    return *this;
	}

	const_iterator& operator++()
	{
	    m_cur = next_in_order((node *)m_cur);
	    return *this;
	}

	const_iterator operator--(int)
	{
	    const_iterator tmp = *this;

	    m_cur = prev_in_order(m_cur);
	    return tmp;
	}

	const_iterator operator++(int)
	{
	    const_iterator tmp = *this;

	    m_cur = next_in_order(m_cur);
	    return tmp;
	}
    };

    // Iterators
    iterator begin()
    {
	return iterator(leftmost_child(&header));
    }

    iterator end()
    {
	return iterator(&header);
    }

    const_iterator begin() const
    {
	return const_iterator(leftmost_child((node *)&header));
    }

    const_iterator end() const
    {
	return const_iterator(&header);
    }

    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    reverse_iterator rbegin()
    {
	return reverse_iterator(end());
    }

    reverse_iterator rend()
    {
	return reverse_iterator(begin());
    }

    const_reverse_iterator rbegin() const
    {
	return const_reverse_iterator(end());
    }

    const_reverse_iterator rend() const
    {
	return const_reverse_iterator(begin());
    }

    // Search functions
    iterator find(const key_type& target)
    {
	return find_eq(comparer, target, &header);
    }

    const_iterator find(const key_type& target) const
    {
	return const_iterator(find_eq(comparer, target, (node *)&header));
    }

    iterator lower_bound(const key_type& target)
    {
	return find_ge(comparer, target, &header);
    }

    const_iterator lower_bound(const key_type& target) const
    {
	return const_iterator(find_ge(comparer, target, (node *)&header));
    }

    iterator upper_bound(const key_type& target)
    {
	return find_le(comparer, target, &header);
    }

    const_iterator upper_bound(const key_type& target) const
    {
	return const_iterator(find_le(comparer, target, (node *)&header));
    }

    // Miscellaneous search wrappers
    size_type count(const key_type& target) const
    {
	return (find(target) == end()) ? 0 : 1;
    }

    std::pair<iterator, iterator> equal_range(const key_type& target)
    {
	iterator b = find(target);
	iterator e = b;

	if (e != end())
	    ++e;

	return std::make_pair(b, e);
    }

    std::pair<const_iterator, const_iterator>
	equal_range(const key_type& target) const
    {
	const_iterator b = find(target);
	const_iterator e = b;

	if (e != end())
	    ++e;

	return std::make_pair(b, e);
    }

    std::pair<iterator, bool> insert(const value_type& val)
    {
	std::pair<node *, bool> r =
	    add_node(comparer, allocator, &header, val);

	if (r.second)
	    item_count++;

	return std::make_pair(r.first, r.second);
    }

    iterator insert(iterator pos, const value_type& val)
    {
	// Special-case optimization: insertion of sequential data (with
	// end() supplied as a hint).
	if ((pos == end()) && header.left) {
	    node *last = rightmost_child(header.left);

	    // Does it come after the last item?
	    if (comparer(last->value.first, val.first)) {
		node *n = new_node(allocator, val);

		last->right = n;
		n->parent = last;

		rebalance_after_insert(n);
		return n;
	    }
	}

	return insert(val).first;
    }

    template<class InputIterator>
    void insert(InputIterator first, InputIterator last)
    {
	while (first != last)
	    insert(end(), *(first++));
    }

    mapped_type& operator[](const key_type& key)
    {
	iterator i = find(key);

	if (i == end())
	    i = insert(std::make_pair(key, mapped_type())).first;

	return i->second;
    }

    void erase(iterator pos)
    {
	node *n = pos.m_cur;

	unlink_node(n);
	delete_node(allocator, n);
	item_count--;
    }

    size_type erase(const key_type& key)
    {
	iterator i = find(key);

	if (i == end())
	    return 0;

	erase(i);
	return 1;
    }

    void erase(iterator first, iterator last)
    {
	while (first != last)
	    erase(first++);
    }

    summary_type sum()
    {
	return sum_tree(aggregator, header.left);
    }

    summary_type sum(const_iterator first, const_iterator last)
    {
	return sum_between(aggregator,
		(node *)first.m_cur, (node *)last.m_cur);
    }

#ifdef SUMMARY_MAP_DEBUG
    void check()
    {
	std::pair<size_type, size_type> r =
	    check_subtree(comparer, &header, header.left, 0, 0);

	assert(r.first == item_count);
	assert(is_black(header.left));
    }
#endif
};

template <class K, class D, class S, class C, class A, class Alloc>
inline void swap(summary_map<K, D, S, C, A, Alloc>& a,
		 summary_map<K, D, S, C, A, Alloc>& b)
{
    a.swap(b);
}

#endif

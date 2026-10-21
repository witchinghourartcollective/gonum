// Copyright ©2015 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package topo

import (
	"slices"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/internal/set"
)

// DegeneracyOrdering returns the degeneracy ordering and the k-cores of
// the undirected graph g.
func DegeneracyOrdering(g graph.Undirected) (order []graph.Node, cores [][]graph.Node) {
	order, offsets := degeneracyOrdering(g)

	slices.Reverse(order)
	cores = make([][]graph.Node, len(offsets))
	offset := len(order)
	for i, n := range offsets {
		cores[i] = order[offset-n : offset]
		offset -= n
	}
	return order, cores
}

// KCore returns the k-core of the undirected graph g with nodes in an
// optimal ordering for the coloring number.
func KCore(k int, g graph.Undirected) []graph.Node {
	order, offsets := degeneracyOrdering(g)

	var offset int
	for _, n := range offsets[:k] {
		offset += n
	}
	core := make([]graph.Node, len(order)-offset)
	copy(core, order[offset:])
	return core
}

// degeneracyOrdering is the common code for DegeneracyOrdering and KCore. It
// returns l, the nodes of g in optimal ordering for coloring number and
// s, a set of relative offsets into l for each k-core, where k is an index
// into s.
func degeneracyOrdering(g graph.Undirected) (l []graph.Node, s []int) {
	nodes := graph.NodesOf(g.Nodes())

	// The algorithm used here is essentially as described at
	// http://en.wikipedia.org/w/index.php?title=Degeneracy_%28graph_theory%29&oldid=640308710

	// Initialize an output list L in return parameters.

	// Compute a number d_v for each vertex v in G,
	// the number of neighbors of v that are not already in L.
	// Initially, these numbers are just the degrees of the vertices.
	dv := make(map[int64]int, len(nodes))
	var (
		maxDegree  int
		neighbours = make(map[int64][]graph.Node)
	)
	for _, n := range nodes {
		id := n.ID()
		adj := graph.NodesOf(g.From(id))
		neighbours[id] = adj
		dv[id] = len(adj)
		if len(adj) > maxDegree {
			maxDegree = len(adj)
		}
	}

	// Initialize an array D such that D[i] contains a list of the
	// vertices v that are not already in L for which d_v = i.
	d := make([][]graph.Node, maxDegree+1)
	for _, n := range nodes {
		deg := dv[n.ID()]
		d[deg] = append(d[deg], n)
	}

	// Initialize k to 0.
	k := 0
	// Repeat n times:
	s = []int{0}
	for range nodes {
		// Scan the array cells D[0], D[1], ... until
		// finding an i for which D[i] is nonempty.
		var (
			i  int
			di []graph.Node
		)
		for i, di = range d {
			if len(di) != 0 {
				break
			}
		}

		// Set k to max(k,i).
		if i > k {
			k = i
			s = append(s, make([]int, k-len(s)+1)...)
		}

		// Select a vertex v from D[i]. Add v to the
		// beginning of L and remove it from D[i].
		var v graph.Node
		v, d[i] = di[len(di)-1], di[:len(di)-1]
		l = append(l, v)
		s[k]++
		delete(dv, v.ID())

		// For each neighbor w of v not already in L,
		// subtract one from d_w and move w to the
		// cell of D corresponding to the new value of d_w.
		for _, w := range neighbours[v.ID()] {
			dw, ok := dv[w.ID()]
			if !ok {
				continue
			}
			for i, n := range d[dw] {
				if n.ID() == w.ID() {
					d[dw][i], d[dw] = d[dw][len(d[dw])-1], d[dw][:len(d[dw])-1]
					dw--
					d[dw] = append(d[dw], w)
					break
				}
			}
			dv[w.ID()] = dw
		}
	}

	return l, s
}

// BronKerbosch returns the set of maximal cliques of the undirected graph g.
//
// The algorithm used is Bron–Kerbosch with an outer degeneracy ordering
// (Eppstein, Löffler and Strash) and the pivot rule of Tomita, Tanaka and
// Takahashi, which selects the pivot u ∈ P ∪ X maximising |P ∩ N(u)|. This
// gives a worst-case running time of O(3^(n/3)), which is optimal for listing
// all maximal cliques.
//
// See https://doi.org/10.1016/j.tcs.2006.06.015 and
// https://doi.org/10.1007/978-3-642-17517-6_36.
func BronKerbosch(g graph.Undirected) [][]graph.Node {
	nodes := graph.NodesOf(g.Nodes())

	// Build the adjacency sets once; the recursion would otherwise
	// repeatedly materialise neighbour lists from g.
	adj := make(map[int64]set.Nodes, len(nodes))
	for _, n := range nodes {
		nid := n.ID()
		to := g.From(nid)
		nb := set.NewNodesSize(to.Len())
		for to.Next() {
			u := to.Node()
			if u.ID() == nid {
				// Self loops do not contribute to cliques.
				continue
			}
			nb.Add(u)
		}
		adj[nid] = nb
	}

	p := set.NewNodesSize(len(nodes))
	for _, n := range nodes {
		p.Add(n)
	}
	x := set.NewNodes()
	bk := bronKerbosch{adj: adj}
	order, _ := degeneracyOrdering(g)
	slices.Reverse(order)
	for _, v := range order {
		nv := adj[v.ID()]
		bk.maximalCliquePivot([]graph.Node{v}, intersect(p, nv), intersect(x, nv))
		p.Remove(v)
		x.Add(v)
	}
	return bk.cliques
}

type bronKerbosch struct {
	adj     map[int64]set.Nodes
	cliques [][]graph.Node
}

func (bk *bronKerbosch) maximalCliquePivot(r []graph.Node, p, x set.Nodes) {
	if len(p) == 0 && len(x) == 0 {
		bk.cliques = append(bk.cliques, r)
		return
	}

	nu := bk.adj[bk.choosePivotFrom(p, x).ID()]

	// Collect P \ N(u) before mutating p.
	candidates := make([]graph.Node, 0, len(p))
	for id, v := range p {
		if _, ok := nu[id]; !ok {
			candidates = append(candidates, v)
		}
	}
	for _, v := range candidates {
		nv := bk.adj[v.ID()]
		sr := append(r[:len(r):len(r)], v)
		bk.maximalCliquePivot(sr, intersect(p, nv), intersect(x, nv))
		p.Remove(v)
		x.Add(v)
	}
}

// choosePivotFrom returns the node u in p ∪ x that maximises |p ∩ N(u)|,
// following Tomita, Tanaka and Takahashi.
func (bk *bronKerbosch) choosePivotFrom(p, x set.Nodes) graph.Node {
	var (
		best  = -1
		pivot graph.Node
	)
	consider := func(s set.Nodes) {
		for _, u := range s {
			nb := bk.adj[u.ID()]
			if len(nb) <= best {
				// Cannot beat the current best.
				continue
			}
			c := intersectionCount(p, nb, best)
			if c > best {
				best = c
				pivot = u
				if best == len(p) {
					return
				}
			}
		}
	}
	consider(p)
	if best < len(p) {
		consider(x)
	}
	if pivot == nil {
		panic("bronKerbosch: empty set")
	}
	return pivot
}

// intersectionCount returns |a ∩ b|. It may return early with a value less
// than or equal to floor once the count can no longer exceed floor.
func intersectionCount(a, b set.Nodes, floor int) int {
	if len(a) > len(b) {
		a, b = b, a
	}
	var c int
	remaining := len(a)
	for id := range a {
		remaining--
		if _, ok := b[id]; ok {
			c++
		}
		if c+remaining <= floor {
			return c
		}
	}
	return c
}

// intersect returns a new set holding a ∩ b, iterating over the smaller set.
func intersect(a, b set.Nodes) set.Nodes {
	if len(a) > len(b) {
		a, b = b, a
	}
	dst := make(set.Nodes, len(a))
	for id, n := range a {
		if _, ok := b[id]; ok {
			dst[id] = n
		}
	}
	return dst
}

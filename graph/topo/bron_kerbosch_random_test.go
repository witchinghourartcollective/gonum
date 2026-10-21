// Copyright ©2026 The Gonum Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package topo

import (
	"fmt"
	"slices"
	"testing"

	"math/rand/v2"

	"gonum.org/v1/gonum/graph"
	"gonum.org/v1/gonum/graph/simple"
)

func randomUndirected(n int, p float64, seed uint64) *simple.UndirectedGraph {
	rnd := rand.New(rand.NewPCG(seed, seed))
	g := simple.NewUndirectedGraph()
	for i := 0; i < n; i++ {
		g.AddNode(simple.Node(i))
	}
	for i := 0; i < n; i++ {
		for j := i + 1; j < n; j++ {
			if rnd.Float64() < p {
				g.SetEdge(simple.Edge{F: simple.Node(i), T: simple.Node(j)})
			}
		}
	}
	return g
}

// bruteMaximalCliques enumerates maximal cliques by checking every vertex
// subset. It is only usable for small n.
func bruteMaximalCliques(g graph.Undirected, n int) map[string]bool {
	isClique := func(mask uint32) bool {
		for i := 0; i < n; i++ {
			if mask&(1<<i) == 0 {
				continue
			}
			for j := i + 1; j < n; j++ {
				if mask&(1<<j) != 0 && !g.HasEdgeBetween(int64(i), int64(j)) {
					return false
				}
			}
		}
		return true
	}
	out := make(map[string]bool)
	for mask := uint32(1); mask < 1<<n; mask++ {
		if !isClique(mask) {
			continue
		}
		maximal := true
		for k := 0; k < n; k++ {
			if mask&(1<<k) == 0 && isClique(mask|1<<k) {
				maximal = false
				break
			}
		}
		if maximal {
			var ids []int64
			for i := 0; i < n; i++ {
				if mask&(1<<i) != 0 {
					ids = append(ids, int64(i))
				}
			}
			out[fmt.Sprint(ids)] = true
		}
	}
	return out
}

func TestBronKerboschRandomAgainstBruteForce(t *testing.T) {
	t.Parallel()
	for seed := uint64(1); seed <= 60; seed++ {
		n := 4 + int(seed%11)
		p := []float64{0.1, 0.3, 0.5, 0.7, 0.9}[seed%5]
		g := randomUndirected(n, p, seed)

		want := bruteMaximalCliques(g, n)
		got := make(map[string]bool)
		for _, c := range BronKerbosch(g) {
			ids := make([]int64, len(c))
			for i, v := range c {
				ids[i] = v.ID()
			}
			slices.Sort(ids)
			k := fmt.Sprint(ids)
			if got[k] {
				t.Errorf("seed %d: duplicate clique %s", seed, k)
			}
			got[k] = true
		}
		if len(got) != len(want) {
			t.Errorf("seed %d (n=%d p=%v): got %d cliques, want %d", seed, n, p, len(got), len(want))
		}
		for k := range want {
			if !got[k] {
				t.Errorf("seed %d: missing clique %s", seed, k)
			}
		}
	}
}

func BenchmarkBronKerboschRandom(b *testing.B) {
	for _, c := range []struct {
		n int
		p float64
	}{{60, 0.3}, {100, 0.3}, {150, 0.3}, {200, 0.3}, {1000, 0.02}, {5000, 0.002}} {
		g := randomUndirected(c.n, c.p, 1)
		b.Run(fmt.Sprintf("n=%d_p=%v", c.n, c.p), func(b *testing.B) {
			var k int
			for b.Loop() {
				k = len(BronKerbosch(g))
			}
			b.ReportMetric(float64(k), "cliques")
		})
	}
}

package ancs

import (
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fasta"
	"sort"
)

type Seg struct {
	s int
	l int
}

func (seg *Seg) End() int {
	return seg.s + seg.l - 1
}
func NewSeg() Seg {
	return Seg{s: 0, l: 0}
}
func SortByStart(s []Seg) []Seg {
	sort.Slice(s, func(i, j int) bool {
		return s[i].s < s[j].s
	})
	return s
}
func MinAncLen(l int, g float64, t float64) int {
	x := 1
	cq := 0.0
	for cq < t {
		x++
		cq = cq + sus.Prob(l, g, x)
	}
	return x
}
func FindHomologies(
	query *fasta.Sequence,
	e *esa.Esa,
	subjectLen int,
	a int) []Seg {
	var qc, qp int
	var mc *esa.Minterval
	var prevStartS, prevEndS int
	// seg := NewSeg()
	n := make(map[int]bool)
	var seg Seg
	var h []Seg
	areEquidist, rightAnchorFound := false, false
	queryLen := query.Length()
	for qc < queryLen {
		queryPrefix := query.Data()[qc:queryLen]
		mc = e.MatchPref(queryPrefix)
		if mc.L > a && isUnique(mc) {
			areEquidist = qc-qp == e.Sa[mc.I]-prevStartS
			if qc > qp && areEquidist {
				seg.l = seg.l + qc - prevEndS + mc.L
				for pos := prevEndS + 1; pos < e.Sa[mc.I]; pos++ {
					n[pos] = true
				}
				rightAnchorFound = true
			} else {
				if rightAnchorFound {
					if seg.s > subjectLen {
						seg.s = seg.s - subjectLen - 1
					}
					h = append(h, seg)
				}
				seg.s = e.Sa[mc.I]
				seg.l = mc.L
				rightAnchorFound = false
			}
			qp = qc
		}
		qc = qc + mc.L + 1
	}
	SortByStart(h)
	predecessor := make(map[int]int)
	score := make(map[int]int)
	visited := make(map[int]bool)
	for i := -1; i < len(h); i++ {
		predecessor[i] = -1
		score[i] = 0
		visited[i] = false
	}
	score[0] = h[0].l
	maxScore, maxIndex := 0, -1
	for i := 1; i < len(h); i++ {
		maxScore = 0
		maxIndex = -1
		for k := 0; k < i-1; k++ {
			if h[k].End() < h[k].l {
				if score[k] > maxScore {
					maxScore = score[k]
					maxIndex = k
				}
			}
		}
		predecessor[i] = maxIndex
		score[i] = score[maxIndex] + h[i].l
	}
	s := argmaxMapInt(score)
	for s > -1 {
		visited[s] = true
		s = predecessor[s]
	}
	var hred []Seg
	for i := 0; i < len(h); i++ {
		if visited[i] {
			hred = append(hred, h[i])
		}
	}
	return hred
}
func isUnique(m *esa.Minterval) bool {
	if m.J == m.I {
		return true
	} else {
		return false
	}
}
func argmaxMapInt(m map[int]int) int {
	var maxKey, maxValue int
	for key, value := range m {
		if value > maxValue {
			maxValue = value
			maxKey = key
		}
	}
	return maxKey
}

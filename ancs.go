package ancs

import (
	"fmt"
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fasta"
	"os"
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
	a int) ([]Seg, map[int]bool) {
	var qc, qp int
	var currLen, currStartS int
	var prevLen, prevStartS, prevEndS int
	var seg Seg
	var h []Seg
	n := make(map[int]bool)
	rightAnchorFound := false
	subjectStrandLen := subjectLen / 2
	queryLen := query.Length()
	for qc < queryLen {
		queryPrefix := query.Data()[qc:queryLen]
		if anchorLongMatch(&currStartS, &currLen,
			n, a, queryPrefix, e) {
			prevEndQ := qp + prevLen
			prevEndS = prevStartS + prevLen
			afterPrev := currStartS > prevEndS
			areEquidist := qc-prevEndQ == currStartS-prevEndS
			onSameStrand := (currStartS < subjectStrandLen) ==
				(prevStartS < subjectStrandLen)
			segCanBeExtended := afterPrev && areEquidist && onSameStrand
			if segCanBeExtended {
				seg.l = seg.l + qc - prevEndQ + currLen
				rightAnchorFound = true
			} else {
				if rightAnchorFound || prevLen/2 >= a {
					if seg.s > subjectStrandLen {
						seg.s = subjectLen + 1 - seg.s - seg.l
					}
					h = append(h, seg)
				}
				seg.s = currStartS
				seg.l = currLen
				rightAnchorFound = false
			}
			qp = qc
			prevLen = currLen
			prevStartS = currStartS
		}
		qc = qc + currLen + 1
		fmt.Println()
	}
	//Close the last segment if open:
	if rightAnchorFound || prevLen/2 >= a {
		if seg.s > subjectStrandLen {
			seg.s = subjectLen + 1 - seg.s - seg.l
		}
		h = append(h, seg)
	}

	if len(h) == 0 {
		fmt.Fprintln(os.Stderr, "No homologous regions found\n")
		os.Exit(0)
	}
	return h, n
}
func ReduceOverlaps(h []Seg) []Seg {
	hlen := len(h)
	if hlen < 2 {
		return h
	}
	predecessor := make([]int, hlen)
	score := make([]int, hlen)
	visited := make([]bool, hlen)
	score[0] = h[0].l
	predecessor[0] = -1
	for i := 1; i < hlen; i++ {
		maxScore := 0
		maxIndex := -1
		for k := 0; k < i; k++ {
			if h[k].End() < h[i].s {
				if score[k] > maxScore {
					maxScore = score[k]
					maxIndex = k
				}
			}
		}
		predecessor[i] = maxIndex
		if maxIndex != -1 {
			score[i] = h[i].l + score[maxIndex]
		} else {
			score[i] = h[i].l
		}
	}

	// Debug messages. Will be removed in the future
	//fmt.Println("***Homologies:")
	//for _, el := range(h) {
	//      fmt.Printf("***(%d, %d)\n", el.s, el.End())
	//}
	//fmt.Println("***Scores:", score)
	s := argmax(score)
	for s != -1 {
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
func TotalSegLen(segments []Seg) int {
	sumlen := 0
	for _, s := range segments {
		sumlen += s.l
	}
	return sumlen
}
func PrintSegsiteRanges(m map[int]bool, file *os.File) {
	if len(m) == 0 {
		fmt.Fprintf(file, "No segregation sites found\n")
	} else {
		k := append([]int{-1}, getSortedIntKeys(m)...)
		k = append(k, -1)
		for i := 1; i < len(k)-1; i++ {
			prev := k[i] == k[i-1]+1
			next := k[i] == k[i+1]-1
			if prev && next {
				continue
			}
			if next {
				fmt.Fprintf(file, "[%d", k[i]+1)
			} else if prev {
				fmt.Fprintf(file, ":%d] ", k[i]+1)
			} else {
				fmt.Fprintf(file, "%d ", k[i]+1)
			}
		}
		fmt.Fprintf(file, "\n")
	}
}
func SegToFasta(segments []Seg,
	e *esa.Esa,
	n map[int]bool,
	printNs bool) []*fasta.Sequence {
	var segfasta []*fasta.Sequence
	for i, s := range segments {
		start := s.s
		end := s.End()
		var data []byte
		for j := start; j < end+1; j++ {
			if printNs && n[j] {
				data = append(data, 'N')
			} else {
				data = append(data, e.T[j])
			}
		}
		segname := fmt.Sprintf("Segment_%d (%d..%d)", i+1, start+1, end+1)
		converted := fasta.NewSequence(segname, data)
		segfasta = append(segfasta, converted)
	}
	return segfasta
}
func anchorLongMatch(
	currStartS, currLen *int,
	n map[int]bool,
	a int,
	queryPrefix []byte,
	e *esa.Esa) bool {
	mc := e.MatchPref(queryPrefix)
	newStartS := e.Sa[mc.I]
	newCurrLen := mc.L
	*currStartS = newStartS
	*currLen = newCurrLen
	//Add the guaranteed mismatch to the segsite map
	n[newStartS+newCurrLen] = true
	lu := (mc.J == mc.I) && (newCurrLen >= a)
	return lu
}
func argmax(x []int) int {
	maxIdx := 0
	for i := 1; i < len(x); i++ {
		if x[i] > x[maxIdx] {
			maxIdx = i
		}
	}
	return maxIdx
}
func getSortedIntKeys(m map[int]bool) []int {
	keys := make([]int, 0, len(m))
	for k, _ := range m {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	return keys
}

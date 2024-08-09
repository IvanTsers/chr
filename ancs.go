package ancs

import (
	"fmt"
	"github.com/evolbioinf/esa"
	"github.com/ivantsers/fasta"
	"os"
	"sort"
)

// Data type Seg contains a zero-based start and length of a segment.
type Seg struct {
	s int
	l int
}

// Method End() returns an inclusive coordinate of the end of a segment.
func (seg *Seg) End() int {
	return seg.s + seg.l
}

// Constructor method NewSeg() returns a segment of specified start and length.
func NewSeg(x, y int) Seg {
	return Seg{s: x, l: y}
}

// SortByStart accepts a slice of segments s and sorts the segments by their start positions in ascending order.
func SortByStart(s []Seg) []Seg {
	sort.Slice(s, func(i, j int) bool {
		return s[i].s < s[j].s
	})
	return s
}

// The function FindHomologies accepts query and subject sequences, an enhanced suffix array of the subject, and the minimum anchor length a. The function returns a slice of segments (homologous regions, or homologies) and a bool map of segregation sites found within the homologies. If no homologies have been found, an empty slice of segments and empty map of segregation sites are returned.
func FindHomologies(
	query *fasta.Sequence,
	subject *fasta.Sequence,
	e *esa.Esa,
	a int) ([]Seg, map[int]bool) {
	var qc, qp int
	var currLen, currStartS int
	var prevLen, prevStartS, prevEndS int
	var seg Seg
	var h []Seg
	n := make(map[int]bool)
	rightAnchor := false
	subjectLen := subject.Length()
	subjectStrandLen := subjectLen / 2
	queryLen := query.Length()
	for qc < queryLen {
		queryPrefix := query.Data()[qc:queryLen]
		if lcpAnchor(&currStartS, &currLen,
			prevStartS, prevLen,
			subjectLen, queryLen, qc, qp,
			a, queryPrefix, subject) || esaAnchor(&currStartS, &currLen, a, queryPrefix, e) {
			prevEndQ := qp + prevLen
			prevEndS = prevStartS + prevLen
			afterPrev := currStartS > prevEndS

			areEquidist := qc-prevEndQ == currStartS-prevEndS

			onSameStrand := (currStartS < subjectStrandLen) ==
				(prevStartS < subjectStrandLen)

			segCanBeExtended := afterPrev &&
				areEquidist &&
				onSameStrand
			if segCanBeExtended {
				prevSegEnd := seg.End()
				gap := qc - prevEndQ
				seg.l = seg.l + gap + currLen
				a := subject.Data()[prevSegEnd : prevSegEnd+gap]
				b := query.Data()[prevEndQ : prevEndQ+gap]

				for i := 0; i < gap; i++ {
					if a[i] != b[i] {
						n[prevSegEnd+i] = true
					}
				}
				rightAnchor = true
			} else {
				if rightAnchor || prevLen/2 >= a {
					if seg.s > subjectStrandLen {
						seg.s = subjectLen + 1 - seg.s - seg.l
					}
					h = append(h, seg)
				}
				seg.s = currStartS
				seg.l = currLen
				rightAnchor = false
			}
			qp = qc
			prevLen = currLen
			prevStartS = currStartS
		}
		qc = qc + currLen + 1
	}
	//Close the last segment if open:
	if rightAnchor || prevLen/2 >= a {
		if seg.s > subjectStrandLen {
			seg.s = subjectLen + 1 - seg.s - seg.l
		}
		h = append(h, seg)
	}
	return h, n
}

// ReduceOverlaps() accepts a sorted slice of segments
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

// \ty(TotalSegLen()) accepts a slice of segments and returns their total length.
func TotalSegLen(segments []Seg) int {
	sumlen := 0
	for _, s := range segments {
		sumlen += s.l
	}
	return sumlen
}

// PrintSegSiteRanges() accepts a bool map of Ns (segregation sites), a slice of segments, and a pointer to an output file, and prints segregation site coordinate ranges.
func PrintSegsiteRanges(n map[int]bool,
	h []Seg, file *os.File) {
	if len(n) == 0 {
		fmt.Fprintf(file, "No segregation sites found\n")
	} else {
		for _, seg := range h {
			k := []int{-1}
			for i := seg.s; i < seg.End(); i++ {
				if n[i] {
					k = append(k, i-seg.s+1)
				}
			}
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
}

// SegToFasta() converts a slice of segments into actual fasta sequences. It accepts a slice of segments, a pointer to the corresponding ESA, a map of Ns, and a bool toggle for printing Ns. It returns a slice of pointers to fasta entries (type fasta.Sequence).
func SegToFasta(segments []Seg,
	e *esa.Esa,
	n map[int]bool,
	printNs bool) []*fasta.Sequence {
	var segfasta []*fasta.Sequence
	for i, s := range segments {
		start := s.s
		end := s.End()
		var data []byte
		for j := start; j < end; j++ {
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

// The function lcpAnchor accepts the following inputs: 1) pointer to the current match length; 2) start of the current match in the subject; 3) a map of segregation sites; 4) the minimum anchor length; 5) the current query prefix; 6) a pointer to the subject. The function returns a boolean. Regardless of the significance of the match, the function updates the start and the length of the current match.
func lcpAnchor(currStartS, currLen *int,
	prevStartS, prevLen int,
	subjectLen, queryLen int,
	qc, qp int,
	a int,
	queryPrefix []byte,
	subject *fasta.Sequence) bool {
	advance := qc - qp
	gap := advance - prevLen
	tryS := prevStartS + advance
	if tryS >= subjectLen || gap > a {
		return false
	}
	*currStartS = tryS
	newCurrLen := lcp(queryLen,
		queryPrefix, subject.Data()[tryS:])
	*currLen = newCurrLen
	return newCurrLen >= a
}
func lcp(max int, a, b []byte) int {
	count := 0
	for i := 0; i < max; i++ {
		if i >= len(a) || i >= len(b) || a[i] != b[i] {
			break
		}
		count++
	}
	return count
}

// The function anchorEsa() accepts: 1) a pointer to the current match length; 2) start of the match the subject; 3) a map of segregation sites; 4) the minimum anchor length; 5) the current query prefix; 6) a pointer to the subject ESA. The function returns a boolean. Regardless of the significance of the match, the function updates the start and the length of the current match.
func esaAnchor(
	currStartS, currLen *int,
	a int,
	queryPrefix []byte,
	e *esa.Esa) bool {
	mc := e.MatchPref(queryPrefix)
	newStartS := e.Sa[mc.I]
	newCurrLen := mc.L
	*currStartS = newStartS
	*currLen = newCurrLen
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

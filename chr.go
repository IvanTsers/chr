package chr

import (
	"fmt"
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fasta"
	"math"
	"os"
	"sort"
)

// Data structure Homologs describes homologous regions of a subject found in queries. This data type contains a slice of segments of the subject S and a map of segregation sites N.
type Homologs struct {
	S []seg
	N map[int]bool
}
type seg struct {
	s int
	l int
}
type subject struct {
	esa            *esa.Esa
	totalL         int
	strandL        int
	a              int
	contigHeaders  []string
	contigSegments []seg
}
type query struct {
	seq    []byte
	l      int
	suffix []byte
}
type match struct {
	l      int
	startS int
	startQ int
	endS   int
	endQ   int
}

// Fields of this data structure contain parameters used to call Intersect(). The parameters include: 1) a reference; 2) path to the directory of target genomes minus the reference; 3) threshold, the minimum fraction of intersecting genomes; 4) p-value of the shustring length (needed for sus.Quantile); 5) a switch to clean* subject's sequence; 6) a switch to clean* query's sequences; 7) a switch to print positions of segregation sites in output's headers; 8) a switch to print N at the positions of mismatches; 8) a switch to print one-based coordinates. *To clean a sequence is to remove non-ATGC nucleotides.
type Parameters struct {
	Reference       []*fasta.Sequence
	TargetDir       string
	Threshold       float64
	ShustrPval      float64
	CleanSubject    bool
	CleanQuery      bool
	PrintSegSitePos bool
	PrintN          bool
	PrintOneBased   bool
}

func (h *Homologs) reduceOverlaps() {
	slen := len(h.S)
	if slen < 2 {
		return
	}
	h.sort()
	segs := h.S
	predecessor := make([]int, slen)
	score := make([]int, slen)
	visited := make([]bool, slen)
	score[0] = segs[0].l
	predecessor[0] = -1
	for i := 1; i < slen; i++ {
		maxScore := 0
		maxIndex := -1
		for k := 0; k < i; k++ {
			if segs[k].end() <= segs[i].s {
				if score[k] > maxScore {
					maxScore = score[k]
					maxIndex = k
				}
			}
		}
		predecessor[i] = maxIndex
		if maxIndex != -1 {
			score[i] = segs[i].l + score[maxIndex]
		} else {
			score[i] = segs[i].l
		}
	}
	s := argmax(score)
	for s != -1 {
		visited[s] = true
		s = predecessor[s]
	}
	var segred []seg
	for i := 0; i < slen; i++ {
		if visited[i] {
			segred = append(segred, segs[i])
		}
	}
	h.S = segred
}
func (seg *seg) end() int {
	return seg.s + seg.l
}
func newSeg(x, y int) seg {
	return seg{s: x, l: y}
}
func (q *query) updSuffix(x int) {
	q.suffix = q.seq[x:]
}
func (h *Homologs) sort() *Homologs {
	sort.Slice(h.S, func(i, j int) bool {
		return h.S[i].s < h.S[j].s
	})
	return h
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

// The function Intersect accepts a struct of Parameters and returns sequences of homologous regions, common for subject (reference) and query sequences.
func Intersect(parameters Parameters) []*fasta.Sequence {
	r := parameters.Reference
	d := parameters.TargetDir
	if parameters.ShustrPval == 0.0 {
		parameters.ShustrPval = 0.95
	}
	numFiles := 0
	dirEntries, err := os.ReadDir(d)
	if err != nil {
		fmt.Fprintf(os.Stderr,
			"chr.Intersect: error reading %v: %v", d, err)
		os.Exit(1)
	}
	numFiles = len(dirEntries)
	if numFiles == 0 {
		return r
	} else {
		var subject subject
		for i, _ := range r {
			if parameters.CleanSubject {
				r[i].Clean()
			}
			r[i].DataToUpper()
		}
		subjectHeader := r[0].Header()
		subjectData := r[0].Data()
		contigHeaders := []string{subjectHeader}
		contigSegs := []seg{newSeg(0, len(subjectData))}
		cL := len(subjectData)
		if len(r) > 1 {
			for i := 1; i < len(r); i++ {
				seq := r[i]
				seqH := seq.Header()
				seqD := seq.Data()
				seqL := len(seqD)
				contigHeaders = append(contigHeaders, seqH)
				subjectData = append(subjectData, '!')
				cL += 1
				cseg := newSeg(cL, seqL)
				contigSegs = append(contigSegs, cseg)
				subjectData = append(subjectData, seqD...)
				cL += seqL
			}
		}
		atgc := 0.0
		gc := 0.0
		for _, c := range subjectData {
			if c == 'A' || c == 'C' || c == 'G' || c == 'T' {
				atgc++
				if c == 'C' || c == 'G' {
					gc++
				}
			}
		}
		gcContent := gc / atgc
		pval := parameters.ShustrPval
		minAncLen := sus.Quantile(cL, gcContent, pval)
		rev := fasta.NewSequence("reverse", subjectData)
		rev.ReverseComplement()
		subjectData = append(subjectData, rev.Data()...)
		sa := esa.MakeEsa(subjectData)
		subject.esa = sa
		subject.totalL = len(subjectData)
		subject.strandL = len(subjectData) / 2
		subject.a = minAncLen
		subject.contigHeaders = contigHeaders
		subject.contigSegments = contigSegs
		homologs := Homologs{S: []seg{}, N: make(map[int]bool)}
		for _, entry := range dirEntries {
			var query query
			filePath := d + "/" + entry.Name()
			f, _ := os.Open(filePath)
			queryData := fasta.ReadAll(f)
			f.Close()
			qSeq := fasta.Concatenate(queryData, 0)
			if parameters.CleanQuery {
				qSeq.Clean()
			}
			qSeq.DataToUpper()
			query.seq = qSeq.Data()
			query.l = qSeq.Length()
			h := findHomologs(query, subject)
			homologs.S = append(homologs.S, h.S...)
			for pos, _ := range h.N {
				homologs.N[pos] = true
			}
		}
		f := parameters.Threshold
		g := numFiles
		t := int(math.Floor(f * float64(g)))
		if t == 0 {
			t = 1
		}
		p := pileHeights(homologs, subject.strandL)
		isAdj := makeMapAdj(homologs)
		intersection := pileToSeg(p, t, isAdj)
		homologs.S = intersection
		printN := parameters.PrintN
		printOneBased := parameters.PrintOneBased
		printSegSitePos := parameters.PrintSegSitePos
		result := homologsToFasta(homologs, subject, printN,
			printOneBased, printSegSitePos)
		return result
	}
}
func findHomologs(query query, subject subject) Homologs {
	h := Homologs{S: []seg{}, N: make(map[int]bool)}
	var qc, qp int
	var c, p match
	var seg seg
	rightAnchor := false
	for qc < query.l {
		query.updSuffix(qc)
		if lcpAnchor(&c, &p, query, subject, qc, qp) || esaAnchor(&c, query, subject) {
			p.endQ = qp + p.l
			p.endS = p.startS + p.l
			afterPrev := c.startS > p.endS

			areEquidist := qc-p.endQ == c.startS-p.endS

			onSameStrand := (c.startS < subject.strandL) ==
				(p.startS < subject.strandL)

			segCanBeExtended := afterPrev &&
				areEquidist &&
				onSameStrand
			if segCanBeExtended {
				prevSegEnd := seg.end()
				gapLen := qc - p.endQ
				seg.l = seg.l + gapLen + c.l
				gapSeqSubject := subject.esa.T[prevSegEnd : prevSegEnd+gapLen]
				gapSeqQuery := query.seq[p.endQ : p.endQ+gapLen]
				for i := 0; i < gapLen; i++ {
					if gapSeqSubject[i] != gapSeqQuery[i] {
						h.N[prevSegEnd+i] = true
					}
				}
				rightAnchor = true
			} else {
				if rightAnchor || p.l/2 >= subject.a {
					if seg.s > subject.strandL {
						seg.s = subject.totalL + 1 - seg.s - seg.l

					}
					h.S = append(h.S, seg)
				}
				seg.s = c.startS
				seg.l = c.l
				rightAnchor = false
			}
			qp = qc
			p.l = c.l
			p.startS = c.startS
		}
		qc = qc + c.l + 1
	}
	//Close the last segment if open:
	if rightAnchor || p.l/2 >= subject.a {
		if seg.s > subject.strandL {
			seg.s = subject.totalL + 1 - seg.s - seg.l

		}
		h.S = append(h.S, seg)
	}
	h.reduceOverlaps()
	return h
}
func lcpAnchor(c *match, p *match,
	query query, subject subject,
	qc, qp int) bool {
	advance := qc - qp
	gap := advance - p.l
	tryS := p.startS + advance
	if tryS >= subject.totalL || gap > subject.a {
		return false
	}
	c.startS = tryS
	newL := lcp(query.l, query.suffix, subject.esa.T[tryS:])
	c.l = newL
	return newL >= subject.a
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
func esaAnchor(c *match, query query, subject subject) bool {
	mc := subject.esa.MatchPref(query.suffix)
	newStartS := subject.esa.Sa[mc.I]
	newL := mc.L
	c.startS = newStartS
	c.l = newL
	lu := (mc.J == mc.I) && (newL >= subject.a)
	return lu
}
func pileHeights(h Homologs, strandL int) []int {
	pile := make([]int, strandL)
	for i := 0; i < len(h.S); i++ {
		seg := h.S[i]
		for j := seg.s; j < seg.end(); j++ {
			pile[j] += 1
		}
	}
	return pile
}
func makeMapAdj(h Homologs) map[int]bool {
	starts := make(map[int]bool)
	ends := []int{}
	for i := 0; i < len(h.S); i++ {
		seg := h.S[i]
		starts[seg.s] = true
		ends = append(ends, seg.end())
	}
	isAdj := make(map[int]bool)
	for _, e := range ends {
		if starts[e] {
			isAdj[e] = true
		}
	}
	return isAdj
}
func pileToSeg(p []int, t int, isAdj map[int]bool) []seg {
	var segs []seg
	var seg seg
	segIsOpen := false
	for k, v := range p {
		if segIsOpen {
			if v < t {
				segs = append(segs, seg)
				segIsOpen = false
			} else {
				seg.l += 1
				if isAdj[k+1] {
					segs = append(segs, seg)
					segIsOpen = false
				}
			}
		} else {
			if v >= t {
				seg.s = k
				seg.l = 1
				segIsOpen = true
			}
		}
	}
	if segIsOpen {
		segs = append(segs, seg)
	}
	return segs
}
func homologsToFasta(h Homologs, subject subject,
	printN bool, printOneBased bool,
	printSegSitePos bool) []*fasta.Sequence {
	var sequences []*fasta.Sequence
	segs := h.S
	ns := h.N
	for _, seg := range segs {
		start := seg.s
		end := seg.end()
		data := make([]byte, 0, seg.l)
		for j := start; j < end; j++ {
			if printN && ns[j] {
				data = append(data, 'N')
			} else {
				data = append(data, subject.esa.T[j])
			}
		}
		ch, cs, ce := findSegment(seg, subject)
		if printOneBased {
			cs += 1
			ce += 1
		}
		header := fmt.Sprintf("%s (%d..%d)", ch, cs, ce)
		if printSegSitePos {
			segsites := buildSegSiteStr(seg, ns, printOneBased)
			header += " " + segsites
		}
		seq := fasta.NewSequence(header, data)
		sequences = append(sequences, seq)
	}
	return sequences
}
func findSegment(seg seg, subject subject) (string, int, int) {
	var ch string
	var cs, ce int
	contigHeaders := subject.contigHeaders
	contigSegments := subject.contigSegments
	for i, contigSeg := range contigSegments {
		if isWithin(seg, contigSeg) {
			ch = contigHeaders[i]
			cs = seg.s - contigSeg.s + 1
			ce = cs + seg.l
			break
		}
	}
	return ch, cs, ce
}
func isWithin(in seg, out seg) bool {
	return in.s >= out.s && in.end() <= out.end()
}
func buildSegSiteStr(seg seg, ns map[int]bool,
	printOneBased bool) string {
	var segSiteStr string
	if len(ns) == 0 {
		segSiteStr = "0;"
	} else {
		segSiteStr += fmt.Sprintf("%d; ", len(ns))
		k := []int{-1}
		for i := seg.s; i < seg.end(); i++ {
			if ns[i] {
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
			coord := k[i]
			if printOneBased {
				coord = k[i] + 1
			}
			if next {
				segSiteStr += fmt.Sprintf("[%d", coord)
			} else if prev {
				segSiteStr += fmt.Sprintf(":%d] ", coord)
			} else {
				segSiteStr += fmt.Sprintf("%d ", coord)
			}
		}
	}
	return segSiteStr
}

package chr

import (
	"fmt"
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/fasta"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fastautils"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
)

// Data structure Homologs describes homologous regions of a subject found in queries. This data type contains a slice of segments of the subject S and a map of segregating sites N.
type Homologs struct {
	S []seg
	N map[int]bool
}
type seg struct {
	s int
	l int
	n map[int]bool
}
type subject struct {
	esa            *esa.Esa
	totalL         int
	strandL        int
	a              int
	contigHeaders  []string
	contigSegments []seg
	contigShifts   map[string]int
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

// Fields of this data structure contain parameters used to call Intersect(). The parameters include:  1) a reference; 2) a switch to shift reference contig start coordinates to the rigth; 3) a path to the directory of target genomes minus the reference; 4) a threshold, the minimum fraction of intersecting genomes; 5) p-value of the shustring length (needed for sus.Quantile); 6) a switch to clean* subject sequence; 7) a switch to clean* query sequences; 8) a switch to print positions of segregating sites in output headers; 9) a switch to print N at the positions of mismatches; 10) a switch to print one-based coordinates; 11) a number of threads.  *To clean a sequence is to remove non-ATGC nucleotides.
type Parameters struct {
	Reference       []*fasta.Sequence
	ShiftRefRight   bool
	TargetDir       string
	Threshold       float64
	ShustrPval      float64
	CleanSubject    bool
	CleanQuery      bool
	PrintSegSitePos bool
	PrintN          bool
	PrintOneBased   bool
	NumThreads      int
}

func (h *Homologs) filterOverlaps() {
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

// The function Intersect accepts a struct of Parameters and returns sequences of homologous regions, common to subject (reference) and query sequences.
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
		subject.contigShifts = make(map[string]int)
		for i, _ := range r {
			if parameters.CleanSubject {
				fastautils.Clean(r[i])
			}
			fastautils.DataToUpper(r[i])
		}
		shiftRefRight := parameters.ShiftRefRight
		subjectHeader := r[0].Header()
		subjectData := r[0].Data()
		if shiftRefRight {
			err := extractShiftField(subjectHeader, &subject)
			if err != nil {
				fmt.Fprint(os.Stderr, err)
				os.Exit(1)
			}
		}
		contigHeaders := []string{subjectHeader}
		contigSegs := []seg{newSeg(0, len(subjectData))}
		cL := len(subjectData)
		if len(r) > 1 {
			for i := 1; i < len(r); i++ {
				seq := r[i]
				seqH := seq.Header()
				seqD := seq.Data()
				seqL := len(seqD)
				var cseg seg // initialize a segment
				if shiftRefRight {
					err := extractShiftField(seqH, &subject)
					if err != nil {
						fmt.Fprint(os.Stderr, err)
						os.Exit(1)
					}
				}
				contigHeaders = append(contigHeaders, seqH)
				subjectData = append(subjectData, '!')
				cL += 1
				cseg = newSeg(cL+1, seqL-1)
				contigSegs = append(contigSegs, cseg)
				subjectData = append(subjectData, seqD...)
				cL += seqL
			}
		}
		numBases := float64(len(subjectData))
		numGC := 0.0
		for _, c := range subjectData {
			if c == 'C' || c == 'G' {
				numGC++
			}
		}
		gc := numGC / numBases
		pval := parameters.ShustrPval
		minAncLen := sus.Quantile(cL, gc, pval)
		rev := fasta.NewSequence("reverse", subjectData)
		rev.ReverseComplement()
		subjectData = append(subjectData, '#')
		subjectData = append(subjectData, rev.Data()...)
		sa := esa.MakeEsa(subjectData)
		subject.esa = sa
		subject.totalL = len(subjectData)
		subject.strandL = subject.totalL / 2
		subject.a = minAncLen
		subject.contigHeaders = contigHeaders
		subject.contigSegments = contigSegs
		// Prepare to process queries concurrently>>
		var wg sync.WaitGroup
		var mu sync.Mutex
		fileChan := make(chan string)
		homologs := Homologs{S: []seg{}, N: make(map[int]bool)}
		worker := func() {
			defer wg.Done()
			for file := range fileChan {
				var query query
				f, _ := os.Open(file)
				queryData := fastautils.ReadAll(f)
				f.Close()
				qSeq, err := fastautils.Concatenate(queryData, '?')
				if err != nil {
					fmt.Fprint(os.Stderr, err)
					os.Exit(1)
				}
				if parameters.CleanQuery {
					fastautils.Clean(qSeq)
				}
				fastautils.DataToUpper(qSeq)
				query.seq = qSeq.Data()
				query.l = len(qSeq.Data())
				h := findHomologs(query, subject)
				mu.Lock()
				homologs.S = append(homologs.S, h.S...)
				homologs.N = appendKeys(homologs.N, h.N)
				mu.Unlock()
			}
		}
		numThreads := parameters.NumThreads
		if numThreads <= 0 {
			numThreads = 1
		}
		// Start workers
		for i := 0; i < numThreads; i++ {
			wg.Add(1)
			go worker()
		}
		go func() {
			for _, entry := range dirEntries {
				filePath := d + "/" + entry.Name()
				fileChan <- filePath
			}
			close(fileChan)
		}()
		// Wait for all workers to finish
		wg.Wait()
		f := parameters.Threshold
		g := numFiles
		t := int(math.Floor(f * float64(g)))
		if t == 0 {
			t = 1
		}
		p := pileHeights(homologs, subject.strandL)
		isAdj := make(map[int]bool)
		if parameters.CleanQuery {
			isAdj = makeMapAdj(homologs)
		}
		contigBounds := make(map[int]bool)
		for _, contig := range subject.contigSegments {
			contigBounds[contig.end()+1] = true
		}
		intersection := pileToSeg(p, t, isAdj, contigBounds)
		homologs.S = intersection
		printN := parameters.PrintN
		printOneBased := parameters.PrintOneBased
		printSegSitePos := parameters.PrintSegSitePos
		result := homologsToFasta(homologs, subject, printN,
			printOneBased, printSegSitePos, shiftRefRight)
		return result
	}
}
func extractShiftField(header string, subject *subject) error {
	var err error
	if header == "" {
		err = fmt.Errorf("chr.Intersect: failed " +
			"extractShiftField: header is empty\n")
		return err
	}
	if subject == nil {
		err = fmt.Errorf("chr.Intersect: failed " +
			"extractShiftField: subject is nil\n")
		return err
	}
	headerFields := strings.Split(header, "$")

	if len(headerFields) < 3 || headerFields[1] == "" || headerFields[2] == "" {
		err = fmt.Errorf(
			"chr.Intersect: error reading "+
				"shift fields in the reference "+
				"header %v\n", header)
	}

	shift, errConv := strconv.Atoi(headerFields[2])
	if errConv != nil {
		err = fmt.Errorf(
			"chr.Intersect: error converting "+
				"a shift field into an integer:"+
				"\n\t%v\n", errConv)
	}

	subject.contigShifts[header] = shift

	return err
}

func appendKeys(a map[int]bool, b map[int]bool) map[int]bool {
	for key, _ := range b {
		a[key] = true
	}
	return a
}
func findHomologs(query query, subject subject) Homologs {
	h := Homologs{S: []seg{}, N: make(map[int]bool)}
	var qc, qp int
	var c, p match
	seg := seg{s: 0, l: 0, n: make(map[int]bool)}
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
					isSegsite := compare(gapSeqSubject[i], gapSeqQuery[i])
					ssPos := -1
					if isSegsite {
						if prevSegEnd > subject.strandL {
							ssPos = subject.totalL - prevSegEnd - i - 1
							seg.n[ssPos] = true
						} else {
							ssPos = prevSegEnd + i
							seg.n[ssPos] = true
						}
					}
				}
				rightAnchor = true
			} else {
				if rightAnchor || p.l/2 >= subject.a {
					if seg.s > subject.strandL {
						seg.s = subject.totalL - seg.s - seg.l
					}
					h.S = append(h.S, seg)
				}
				seg.s = c.startS
				seg.l = c.l
				seg.n = make(map[int]bool)
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
			seg.s = subject.totalL - seg.s - seg.l
		}
		h.S = append(h.S, seg)
	}
	h.filterOverlaps()
	for _, seg := range h.S {
		h.N = appendKeys(h.N, seg.n)
	}
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
	newL := lcpLen(query.l, query.suffix, subject.esa.T[tryS:])
	c.l = newL
	return newL >= subject.a
}
func lcpLen(max int, a, b []byte) int {
	limit := min(len(a), len(b), max)
	count := 0
	for i := 0; i < limit; i++ {
		if a[i] != b[i] {
			break
		}
		count++
	}
	return count
}

func min(vals ...int) int {
	m := vals[0]
	for _, v := range vals {
		if v < m {
			m = v
		}
	}
	return m
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
func compare(s byte, q byte) bool {
	notEqual := true
	if s == q {
		notEqual = false
	}
	return notEqual
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
func pileToSeg(p []int, t int,
	isAdj map[int]bool,
	contigBounds map[int]bool) []seg {
	var segs []seg
	var seg seg
	segIsOpen := false
	for k, v := range p {
		if contigBounds[k] {
			continue
		}
		if segIsOpen {
			if v < t || contigBounds[k+1] {
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
	printN bool,
	printOneBased bool,
	printSegSitePos bool,
	shiftRefRight bool) []*fasta.Sequence {

	var sequences []*fasta.Sequence
	segs := h.S
	ns := h.N
	for _, seg := range segs {
		start := seg.s
		end := seg.end()
		data := make([]byte, seg.l)
		copy(data, subject.esa.T[start:end])
		if printN {
			for j := 0; j < seg.l; j++ {
				if ns[start+j] {
					data[j] = 'N'
				}
			}
		}
		ch, cs, ce := findSegment(seg, subject, shiftRefRight)
		if printOneBased {
			cs += 1
		}
		header := fmt.Sprintf("%s_(%d..%d)", ch, cs, ce)
		if printSegSitePos {
			segsites := buildSegSiteStr(seg, ns, printOneBased)
			header += " " + segsites
		}
		seq := fasta.NewSequence(header, data)
		sequences = append(sequences, seq)
	}
	return sequences
}
func findSegment(seg seg, subject subject,
	shiftRefRight bool) (string, int, int) {

	var ch string
	var cs, ce, shift int
	contigHeaders := subject.contigHeaders
	contigSegments := subject.contigSegments
	contigShifts := subject.contigShifts

	for i, contigSeg := range contigSegments {
		if isWithin(seg, contigSeg) {
			ch = contigHeaders[i]
			if shiftRefRight {
				shift = contigShifts[ch]
				chFields := strings.Split(ch, "$")
				ch = chFields[0]
			}
			cs = seg.s - contigSeg.s + shift
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
	var positions []int
	for i := seg.s; i < seg.end(); i++ {
		if ns[i] {
			pos := i - seg.s
			if printOneBased {
				pos++
			}
			positions = append(positions, pos)
		}
	}
	var builder strings.Builder
	numSegSites := len(positions)
	builder.WriteString(strconv.Itoa(numSegSites))

	for _, pos := range positions {
		builder.WriteString(" ")
		builder.WriteString(strconv.Itoa(pos))
	}

	return builder.String()
}

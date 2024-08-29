package chr

import (
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fasta"
	"math"
	"math/rand"
	"os"
	"reflect"
	"sort"
	"testing"
	"time"
)

func readFasta(path string) []*fasta.Sequence {
	f, _ := os.Open(path)
	contigs := fasta.ReadAll(f)
	f.Close()
	return contigs
}
func prepareSubject(path string) subject {
	var subject subject
	r := readFasta(path)
	for i, _ := range r {
		r[i].Clean()
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
	minAncLen := sus.Quantile(cL, gcContent, 0.95)
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
	return subject
}
func prepareQuery(path string) query {
	var query query
	q := readFasta(path)
	qSeq := fasta.Concatenate(q, 0)
	qSeq.Clean()
	query.seq = qSeq.Data()
	query.l = qSeq.Length()
	return query
}
func containsData(slice []*fasta.Sequence,
	el *fasta.Sequence) bool {
	for _, item := range slice {
		if reflect.DeepEqual(item.Data(), el.Data()) {
			return true
		}
	}
	return false
}
func TestSort(t *testing.T) {
	var starts []int
	var segments []seg
	for i := 0; i < 12; i++ {
		rand.Seed(time.Now().UnixNano())
		rs := rand.Intn(1000000000)
		rl := 1 + rand.Intn(999999999)
		starts = append(starts, rs)
		segments = append(segments, newSeg(rs, rl))
	}
	h := Homologs{S: segments, N: make(map[int]bool)}
	hSorted := h.sort()
	sort.Ints(starts)
	for _, i := range [3]int{0, 5, 11} {
		want := starts[i]
		get := hSorted.S[i].s
		if want != get {
			t.Errorf("want:\n%d\nget:\n%d\n",
				want, get)
		}
	}
	x := newSeg(12, 221)
	testCases := []struct {
		name  string
		input []seg
		want  []seg
	}{
		{"Empty", []seg{}, []seg{}},
		{"Single element", []seg{x}, []seg{x}},
		{"Identical", []seg{x, x, x, x}, []seg{x, x, x, x}},
		{"Already sorted", hSorted.S, hSorted.S},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			h.S = tc.input
			hSort := h.sort()
			get := hSort.S
			if !reflect.DeepEqual(get, tc.want) {
				t.Errorf("want:\n%v\nget:\n%v\n",
					tc.want, get)
			}
		})
	}
}
func TestReduceOverlaps(t *testing.T) {
	x := newSeg(0, 10)
	y := newSeg(5, 10)
	z := newSeg(5, 20)
	a := newSeg(21, 2)
	b := newSeg(25, 5)
	c := newSeg(27, 1)
	d := newSeg(30, 5)
	e := newSeg(35, 10)
	f := newSeg(50, 2)
	g := newSeg(0, 5)
	testCases := []struct {
		name  string
		input []seg
		want  []seg
	}{
		{"Empty", []seg{}, []seg{}},
		{"Sinlge", []seg{x}, []seg{x}},
		{"No overlap", []seg{x, f}, []seg{x, f}},
		{"Perfect overlap", []seg{x, x}, []seg{x}},
		{"Partial overlap same len", []seg{x, y}, []seg{x}},
		{"Partial overlap diff len", []seg{x, z}, []seg{z}},
		{"Various", []seg{g, x, y, z, a, b, c, d, e, f},
			[]seg{g, z, b, d, e, f}},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			h := Homologs{S: tc.input, N: make(map[int]bool)}
			h.reduceOverlaps()
			get := h.S
			if !reflect.DeepEqual(get, tc.want) {
				t.Errorf("\nwant:\n%v\nget:\n%v\n", tc.want, get)
			}
		})
	}
}
func TestFindHomologs(t *testing.T) {
	identical := readFasta("data/o/identical.fasta")
	tooDistant := readFasta("data/o/tooDistant.fasta")
	s1q1 := readFasta("data/o/s1q1.fasta")
	s1q2a := readFasta("data/o/s1q2a.fasta")
	s1q2b := readFasta("data/o/s1q2b.fasta")
	s2q1a := readFasta("data/o/s2q1a.fasta")
	s2q1b := readFasta("data/o/s2q1b.fasta")
	s2q2 := readFasta("data/o/s2q2.fasta")
	onEnds := readFasta("data/o/onEnds.fasta")
	mutations := readFasta("data/o/mutations.fasta")
	testCases := []struct {
		name  string
		input string
		want  []*fasta.Sequence
	}{
		{"identical", "identical.fasta", identical},
		{"tooDistant", "tooDistant.fasta", tooDistant},
		{"s1q1", "s1q1.fasta", s1q1},
		{"s1q2a", "s1q2a.fasta", s1q2a},
		{"s1q2b", "s1q2b.fasta", s1q2b},
		{"s2q1a", "s2q1a.fasta", s2q1a},
		{"s2q1b", "s2q1b.fasta", s2q1b},
		{"s2q2", "s2q2.fasta", s2q2},
		{"onEnds", "onEnds.fasta", onEnds},
		{"mutations", "mutations.fasta", mutations},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			subject := prepareSubject("data/s/" + tc.input)
			query := prepareQuery("data/q/" + tc.input)
			h := findHomologs(query, subject)
			h.reduceOverlaps()
			get := homologsToFasta(h, subject, false, true)
			PrintSegsiteRanges(h, os.Stdout, false)
			for i := range get {
				if !reflect.DeepEqual(get[i].Data(),
					tc.want[i].Data()) {
					t.Errorf("\nwant:\n%v\nget:\n%v\n",
						string(tc.want[i].Data()),
						string(get[i].Data()))
				}
			}
		})
	}
}
func TestPileToSeg(t *testing.T) {
	slen := 10
	a := newSeg(0, 2)
	b := newSeg(0, 3)
	c := newSeg(1, 2)
	d := newSeg(3, 5)
	e := newSeg(5, 2)
	f := newSeg(4, 4)
	g := newSeg(8, 1)
	segments := []seg{a, b, a, a, c,
		d, f, e, e, f, g,
	}
	h := Homologs{S: segments, N: make(map[int]bool)}
	t.Run("pileHeights", func(t *testing.T) {
		want := []int{4, 5, 2, 1, 3, 5, 5, 3, 1, 0}
		get := pileHeights(h, slen)
		if !reflect.DeepEqual(want, get) {
			t.Errorf("\nwant:\n%v\nget:\n%v\n", want, get)
		}
	})
	zero := []seg{newSeg(0, 3), newSeg(3, 5), newSeg(8, 1)}
	forty := []seg{newSeg(0, 3), newSeg(4, 4)}
	sixty := []seg{newSeg(0, 2), newSeg(4, 4)}
	eighty := []seg{newSeg(0, 2), newSeg(5, 2)}
	hundred := []seg{newSeg(1, 1), newSeg(5, 2)}
	testCases := []struct {
		name      string
		threshold float64
		want      []seg
	}{
		{"0%", 0.0, zero},
		{"10%", 0.1, zero},
		{"20%", 0.2, zero},
		{"40%", 0.4, forty},
		{"60%", 0.6, sixty},
		{"80%", 0.8, eighty},
		{"90%", 0.9, eighty},
		{"100%", 1.0, hundred},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			want := tc.want
			p := pileHeights(h, slen)
			isAdj := makeMapAdj(h)
			threshold := int(math.Floor(tc.threshold * 5.0))
			if threshold == 0 {
				threshold = 1
			}
			get := pileToSeg(p, threshold, isAdj)
			if !reflect.DeepEqual(want, get) {
				t.Errorf("\nwant:\n%v\nget:\n%v\n", want, get)
			}
		})
	}
}
func TestIntersect(t *testing.T) {
	stan := readFasta("data/o/stan.fasta")
	sars := readFasta("data/o/sars.fasta")

	testCases := []struct {
		name  string
		input string
		want  []*fasta.Sequence
	}{
		{"stan", "data/i/stan/*", stan},
		{"sars", "data/i/sars/*", sars},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			want := tc.want
			r := readFasta("data/i/" + tc.name + "/t1.fasta")
			parameters := Parameters{
				Reference:     r,
				TargetDir:     "data/i/" + tc.name + "/t",
				Threshold:     1.0,
				PrintN:        false,
				PrintOneBased: true,
			}
			get := Intersect(parameters)
			wL := len(want)
			gL := len(get)
			minLen := gL
			if gL != wL {
				t.Errorf("\nthe result has %d sequences, expected %d",
					gL, wL)
				if wL < gL {
					minLen = wL
				}

			}
			wConcat := fasta.Concatenate(want, 0)
			gConcat := fasta.Concatenate(get, 0)
			if len(wConcat.Data()) != len(gConcat.Data()) {
				t.Errorf("\nthe result has %d nucleotides, expected %d",
					len(gConcat.Data()), len(wConcat.Data()))
			}
			for i := 0; i < minLen; i++ {
				if !containsData(get, want[i]) {
					t.Errorf("\n'%v' is absent from results:\n%v\n",
						want[i].Header(), string(want[i].Data()))
				}
			}
		})
	}
}

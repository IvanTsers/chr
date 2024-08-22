package ancs

import (
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fasta"
	"math/rand"
	"os"
	"path/filepath"
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
func prepareSeq(path string) *fasta.Sequence {
	contigs := readFasta(path)
	for i, _ := range contigs {
		contigs[i].Clean()
	}
	seq := fasta.Concatenate(contigs, '!')
	return seq
}
func prepareSubject(path string) (*esa.Esa, int) {
	seq := prepareSeq(path)
	rev := fasta.NewSequence("reverse", seq.Data())
	rev.ReverseComplement()
	seq.SetData(append(seq.Data(), rev.Data()...))
	sa := esa.MakeEsa(seq.Data())
	gc := seq.GC()
	ma := sus.Quantile(seq.Length()/2, gc, 0.95)
	return sa, ma
}
func approxAlignment(ref string, inFiles string) []*fasta.Sequence {
	allFiles, _ := getFiles(inFiles)
	queryNames := []string{}
	for _, fileName := range allFiles {
		if fileName != ref {
			queryNames = append(queryNames, fileName)
		}
	}
	numQueries := len(queryNames)
	sa, ma := prepareSubject(ref)

	allHomologies := []Seg{}
	allSegsites := make(map[int]bool)

	for _, q := range queryNames {
		query := prepareSeq(q)
		h, n := FindHomologies(query, sa, ma)
		h = SortByStart(h)
		h = ReduceOverlaps(h)
		allHomologies = append(allHomologies, h...)
		for pos, _ := range n {
			allSegsites[pos] = true
		}
	}
	intersection := Intersect(allHomologies,
		numQueries, 1.0, len(sa.T)/2)

	result := SegToFasta(intersection, sa, allSegsites, false)
	return result
}
func getFiles(pattern string) ([]string, error) {
	files, err := filepath.Glob(pattern)
	if err != nil {
		return nil, err
	}
	return files, nil
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
func TestSeg(t *testing.T) {
	seg := NewSeg(1, 10)
	want := 11
	get := seg.End()
	if get != want {
		t.Errorf("want:\n%d\nget:\n%d\n",
			want, get)
	}
	seg.l += 10
	want += 10
	get = seg.End()
	if get != want {
		t.Errorf("want:\n%d\nget:\n%d\n",
			want, get)
	}
}
func TestSortByStart(t *testing.T) {
	var starts []int
	var segments []Seg
	for i := 0; i < 12; i++ {
		rand.Seed(time.Now().UnixNano())
		rs := rand.Intn(1000000000)
		rl := 1 + rand.Intn(999999999)
		starts = append(starts, rs)
		segments = append(segments, NewSeg(rs, rl))
	}
	sorted := SortByStart(segments)
	sort.Ints(starts)
	for _, i := range [3]int{0, 5, 11} {
		want := starts[i]
		get := sorted[i].s
		if want != get {
			t.Errorf("want:\n%d\nget:\n%d\n",
				want, get)
		}
	}
	x := NewSeg(12, 221)
	testCases := []struct {
		name  string
		input []Seg
		want  []Seg
	}{
		{"Empty", []Seg{}, []Seg{}},
		{"Single element", []Seg{x}, []Seg{x}},
		{"Identical", []Seg{x, x, x, x}, []Seg{x, x, x, x}},
		{"Already sorted", sorted, sorted},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			inputCopy := make([]Seg, len(tc.input))
			copy(inputCopy, tc.input)
			get := SortByStart(inputCopy)
			if !reflect.DeepEqual(get, tc.want) {
				t.Errorf("want:\n%v\nget:\n%v\n",
					tc.want, get)
			}
		})
	}
}
func TestReduceOverlaps(t *testing.T) {
	x := NewSeg(0, 10)
	y := NewSeg(5, 10)
	z := NewSeg(5, 20)
	a := NewSeg(21, 2)
	b := NewSeg(25, 5)
	c := NewSeg(27, 1)
	d := NewSeg(30, 5)
	e := NewSeg(35, 10)
	f := NewSeg(50, 2)
	g := NewSeg(0, 5)
	testCases := []struct {
		name  string
		input []Seg
		want  []Seg
	}{
		{"Empty", []Seg{}, []Seg{}},
		{"Sinlge", []Seg{x}, []Seg{x}},
		{"Perfect overlap", []Seg{x, x}, []Seg{x}},
		{"Partial overlap same len", []Seg{x, y}, []Seg{x}},
		{"Partial overlap diff len", []Seg{x, z}, []Seg{z}},
		{"Various", []Seg{g, x, y, z, a, b, c, d, e, f},
			[]Seg{g, z, b, d, e, f}},
	}
	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			inputCopy := make([]Seg, len(tc.input))
			copy(inputCopy, tc.input)
			get := ReduceOverlaps(inputCopy)
			if !reflect.DeepEqual(get, tc.want) {
				t.Errorf("\nwant:\n%v\nget:\n%v\n", tc.want, get)
			}
		})
	}
}
func TestFindHomologies(t *testing.T) {
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
			sa, ma := prepareSubject("data/s/" + tc.input)
			q := prepareSeq("data/q/" + tc.input)
			h, n := FindHomologies(q, sa, ma)
			h = SortByStart(h)
			h = ReduceOverlaps(h)
			get := SegToFasta(h, sa, n, false)
			PrintSegsiteRanges(n, h, os.Stdout)
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
func TestIntersect(t *testing.T) {
	slen := 10
	a := NewSeg(0, 2)
	b := NewSeg(0, 3)
	c := NewSeg(1, 2)
	d := NewSeg(3, 5)
	e := NewSeg(5, 2)
	f := NewSeg(4, 4)
	g := NewSeg(8, 1)
	h := []Seg{a, b, a, a, c,
		d, f, e, e, f, g,
	}
	t.Run("pileHeights", func(t *testing.T) {
		want := []int{4, 5, 2, 1, 3, 5, 5, 3, 1, 0}
		get := pileHeights(h, slen)
		if !reflect.DeepEqual(want, get) {
			t.Errorf("\nwant:\n%v\nget:\n%v\n", want, get)
		}
	})
	zero := []Seg{NewSeg(0, 3), NewSeg(3, 5), NewSeg(8, 1)}
	forty := []Seg{NewSeg(0, 3), NewSeg(4, 4)}
	sixty := []Seg{NewSeg(0, 2), NewSeg(4, 4)}
	eighty := []Seg{NewSeg(0, 2), NewSeg(5, 2)}
	hundred := []Seg{NewSeg(1, 1), NewSeg(5, 2)}
	testCases := []struct {
		name      string
		threshold float64
		want      []Seg
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
			get := Intersect(h, 5, tc.threshold, slen)
			if !reflect.DeepEqual(want, get) {
				t.Errorf("\nwant:\n%v\nget:\n%v\n", want, get)
			}
		})
	}
}
func TestAlignment(t *testing.T) {
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
			get := approxAlignment("data/i/"+tc.name+"/t1.fasta", tc.input)
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

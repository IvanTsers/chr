package ancs

import (
	"github.com/evolbioinf/esa"
	"github.com/evolbioinf/sus"
	"github.com/ivantsers/fasta"
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
func prepareSeq(path string) *fasta.Sequence {
	contigs := readFasta(path)
	for i, _ := range contigs {
		contigs[i].Clean()
	}
	seq := fasta.Concatenate(contigs, '!')
	return seq
}
func prepareSubject(path string) (*fasta.Sequence, *esa.Esa, int) {
	seq := prepareSeq(path)
	rev := fasta.NewSequence("reverse", seq.Data())
	rev.ReverseComplement()
	seq.SetData(append(seq.Data(), rev.Data()...))
	sa := esa.MakeEsa(seq.Data())
	gc := seq.GC()
	ma := sus.Quantile(seq.Length()/2, gc, 0.95)
	return seq, sa, ma
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
			[]Seg{z, c, e, f}},
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
			s, sa, ma := prepareSubject("data/s/" + tc.input)
			q := prepareSeq("data/q/" + tc.input)
			h, n := FindHomologies(q, s, sa, ma)
			h = SortByStart(h)
			h = ReduceOverlaps(h)
			get := SegToFasta(h, sa, n, false)
			PrintSegsiteRanges(n, h, os.Stdout)
			for i := range get {
				if !reflect.DeepEqual(get[i].Data(), tc.want[i].Data()) {
					t.Errorf("\nwant:\n%v\nget:\n%v\n",
						string(tc.want[i].Data()),
						string(get[i].Data()))
				}
			}
		})
	}
}

# [`chr`](https://owncloud.gwdg.de/index.php/s/rUiY9yfvhxMp0ja): common homologous regions
This package contains functions for fast search of homologous regions shared by multiple closely related genomes. The interface of `chr` consists of the `Intersect` function and the `Parameters` struct:
```go
chr.Intersect(chr.Parameters)
```
The function returns a slice of pointers to fasta sequences `[]*fasta.Sequence`. The type `fasta` is defined in the [`fasta`](https://github.com/evolbioinf/fasta) package.

### `Parameters`
```go
type Parameters struct {
        Reference       []*fasta.Sequence
        QueryPaths      []string
        Threshold       float64
        ShustrPval      float64
        PrintN          bool
        NumThreads   	int
}
```
 Fields of this data structure contain parameters used to call Intersect().
    The parameters include: 

 1. reference slice of `fasta` sequences
 3. paths to queries
 4. a threshold, the minimum fraction of intersecting genomes
 5. a p-value of the shustring length (needed for [`sus`](https://github.com/evolbioinf/sus).Quantile)
 6. number of threads (maximum threads used for parallel query processing)

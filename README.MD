# [`chr`](https://owncloud.gwdg.de/index.php/s/rUiY9yfvhxMp0ja): common homologous regions
This package contains functions for fast search of homologous regions shared by multiple closely related genomes. The interface of `chr` consists of the `Intersect` function and the `Parameters` struct:
```go
chr.Intersect(chr.Parameters)
```
The function returns a slice of pointers to fasta sequences `[]*fasta.Sequence`. The type `fasta` is defined in the [`fasta`](https://github.com/IvanTsers/fasta) package.

### `Parameters`
```go
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
```
 Fields of this data structure contain parameters used to call Intersect().
    The parameters include: 

 1. reference slice of [`fasta`](https://github.com/IvanTsers/fasta) sequences
 2. path to the directory of target genomes minus the reference
 3. a threshold, the minimum fraction of intersecting genomes
 4. a p-value of the shustring length (needed for [`sus`](https://github.com/evolbioinf/sus).Quantile)
 5. a switch to clean subject's sequence. To clean a sequence is to remove non-ATGC nucleotides.
 6. a switch to clean query sequences
 7. a switch to print positions of segregation sites in output's headers
 8. a switch to print N at the positions of mismatches
 9. a switch to print one-based coordinates.
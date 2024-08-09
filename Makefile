EXE = ancs

all: $(EXE)

$(EXE): $(EXE).go
	go build $(EXE).go
$(EXE).go: $(EXE).org
	awk -f scripts/preTangle.awk $(EXE).org | bash scripts/org2nw | notangle -R$(EXE).go | gofmt > $(EXE).go
test: $(EXE)_test.go $(EXE)
	go test -v
$(EXE)_test.go: $(EXE)_test.org
	awk -f scripts/preTangle.awk $(EXE)_test.org | bash scripts/org2nw | notangle -R$(EXE)_test.go | gofmt > $(EXE)_test.go

.PHONY: doc clean

doc:
	make -C doc

clean:
	rm -f $(EXE) *.go
	make clean -C doc

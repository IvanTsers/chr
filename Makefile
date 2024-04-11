EXE = ancs

all: $(EXE)

$(EXE): $(EXE).go
	go build $(EXE).go
$(EXE).go: $(EXE).org
	bash scripts/org2nw $(EXE).org | notangle -R$(EXE).go | gofmt > $(EXE).go
#test: $(EXE) $(EXE)_test.go
#	go test -v
#$(EXE)_test.go: $(EXE).org
#	bash scripts/org2nw $(EXE).org | notangle -R$(EXE)_test.go | gofmt > $(#EXE)_test.go

.PHONY: doc clean

doc:
	make -C doc

clean:
	rm -f $(EXE) *.go
	make clean -C doc

package main

import (
	"os"
	"reflect"
	"testing"
)

func TestWritePBS(t *testing.T) {
	p := Job{MakeName(Input[Geometry]), "opt.inp", 35, "", "", numJobs}
	write := "testfiles/write/mp.pbs"
	right := "testfiles/right/mp.pbs"
	WritePBS(write, &p, pbsSequoia)
	if !compareFile(write, right) {
		t.Errorf("mismatch between %s and %s\n", right, write)
	}
}

func TestSubmit(t *testing.T) {
	got := Submit("opt/mp.pbs")
	want := "775241"
	if got != want {
		t.Errorf("got %v, wanted %v\n", got, want)
	}
}

func TestReadPBSNodes(t *testing.T) {
	// cn074 has 6 jobs
	f, _ := os.Open("testfiles/read/pbsnodes")
	got := readPBSnodes(f)
	want := []string{"workq:cn064", "workq:cn065", "workq:cn066", "workq:cn067"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("got %q, wanted %q\n", got, want)
	}
}

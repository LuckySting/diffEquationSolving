package main

import (
	"fmt"
	"github.com/cpmech/gosl/la"
	"github.com/james-bowman/sparse"
	"github.com/james-bowman/sparse/blas"
	"gonum.org/v1/gonum/mat"
	"math"
)

type CSRProxy struct {
	sparse.CSR
}

func (c *CSRProxy) MulVecTo(dst *mat.VecDense, trans bool, x mat.Vector) {
	ar, ac := c.Dims()
	if trans {
		ar, ac = ac, ar
	}
	if ac != x.Len() || ar != dst.Len() {
		panic(mat.ErrShape)
	}
	blas.Dusmv(trans, 1, c.RawMatrix(), x.(*mat.VecDense).RawVector().Data, 1, dst.RawVector().Data, 1)
}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

func boundary(x float64, z float64) float64 {
	return 300
}

func function(x float64, z float64) float64 {
	return 1000 * math.Exp(-math.Pow(x-0, 2)*math.Pow(z-0, 2))
}

func coefficient(x float64, z float64) float64 {
	return 0.536
}

func getPIndex(i int, j int, nX int) int {
	return j*(nX+1) + i
}

func solve(xBound float64, zBound float64, hStep float64) {
	nX := int(math.Round(xBound / hStep))
	nZ := int(math.Round(zBound / hStep))
	n := (nX + 1) * (nZ + 1)
	bMatrix := la.NewVector(n)
	eqMatrixTriplet := la.NewTriplet(n, n, 0)
	for i := 0; i < nX+1; i++ {
		for j := 0; j < nZ+1; j++ {
			pIndex := getPIndex(i, j, nX)
			if i == 0 || j == 0 || i == nX || j == nZ {
				eqMatrixTriplet.Put(pIndex, pIndex, 1)
				bMatrix[pIndex] = boundary(hStep*float64(i), hStep*float64(j))
			} else {
				coefA := coefficient((float64(i)-0.5)*hStep, float64(j)*hStep)
				coefB := coefficient((float64(i)+0.5)*hStep, float64(j)*hStep)
				coefC := coefficient(float64(i)*hStep, (float64(j)-0.5)*hStep)
				coefD := coefficient(float64(i)*hStep, (float64(j)+0.5)*hStep)
				coefF := -coefA - coefB - coefC - coefD

				eqMatrixTriplet.Put(pIndex, getPIndex(i-1, j, nX), coefA)

				eqMatrixTriplet.Put(pIndex, getPIndex(i+1, j, nX), coefB)

				eqMatrixTriplet.Put(pIndex, getPIndex(i, j-1, nX), coefC)

				eqMatrixTriplet.Put(pIndex, getPIndex(i, j+1, nX), coefD)

				eqMatrixTriplet.Put(pIndex, pIndex, coefF)

				bMatrix[pIndex] = -function(hStep*float64(i), hStep*float64(j)) * hStep * hStep
			}
		}
	}
	sparseSolver := la.NewSparseSolver("Umfpack")
	defer sparseSolver.Free()
	sparseSolver.Init(eqMatrixTriplet, nil)
	x := la.NewVector(len(bMatrix))
	sparseSolver.Solve(x, bMatrix)
	matPrint(mat.NewDense(nX, nZ, x))
}

func main() {
	solve(1, 1, 0.25)
}

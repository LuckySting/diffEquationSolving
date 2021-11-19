package helpers

import (
	"fmt"
	"gonum.org/v1/gonum/mat"
)

func FillSlice(from int, to int, step float64, slice []float64, f func(float64) float64) {
	for i := from; i < to; i++ {
		slice[i] = f(step * float64(i))
	}
}

func ThomasAlgorithm(botDiag []float64, midDiag []float64, topDiag []float64, bVector []float64) []float64 {
	n := len(bVector)
	vVector := make([]float64, n)
	uVector := make([]float64, n)
	vVector[0] = topDiag[0] / -midDiag[0]
	uVector[0] = -bVector[0] / -midDiag[0]
	for i := 1; i < n-1; i++ {
		alpha := -midDiag[i] - botDiag[i]*vVector[i-1]
		vVector[i] = topDiag[i] / alpha
		uVector[i] = (botDiag[i]*uVector[i-1] - bVector[i]) / alpha
	}
	vVector[n-1] = 0
	uVector[n-1] = (botDiag[n-1]*uVector[n-1] - bVector[n-1]) / (-midDiag[n-1] - botDiag[n-1]*vVector[n-1])
	xVector := make([]float64, n)
	xVector[n-1] = uVector[n-1]
	for i := n - 2; i > -1; i-- {
		xVector[i] = vVector[i]*xVector[i+1] + uVector[i]
	}
	return xVector
}

// *  1  0  0  0  0		1
// *  1  2  1  0  0		4
// *  0  1  2  1  0		4
// *  0  0  1  2  1		4
// *  0  0  0  0  1		1

// v0 =    0 ; u0 = 1
// alpha = -2 v1 = -0.5 ; u1 = 1.5

func Test() {
	topDiag := []float64{0, 1, 1, 1, 0}
	midDiag := []float64{1, 2, 2, 2, 1}
	botDiag := []float64{0, 1, 1, 1, 0}
	bVector := []float64{1, 4, 4, 4, 1}
	rawData := make([]float64, 5*3)
	for i := 0; i < 5; i++ {
		rawData[3*i] = botDiag[i]
		rawData[3*i+1] = midDiag[i]
		rawData[3*i+2] = topDiag[i]
	}
	A := mat.NewBandDense(5, 5, 1, 1, rawData)

	x := mat.NewVecDense(5, nil)
	x.SolveVec(A, mat.NewVecDense(5, bVector))
	fx := mat.Formatted(x, mat.FormatPython())
	fmt.Printf("Original = %v\n", fx)
	x2 := ThomasAlgorithm(botDiag, midDiag, topDiag, bVector)
	fmt.Println(x2)
}

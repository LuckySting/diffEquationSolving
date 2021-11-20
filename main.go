package main

import (
	"C"
	"diffEquationSolving/helpers"
	"fmt"
	"gonum.org/v1/gonum/mat"
	"math"
	"time"
)

func boundary(x float64, z float64) float64 {
	return 300
}

func mainFunc(x float64, z float64) float64 {
	return 1000 * math.Exp(-math.Pow(x-2, 2)*math.Pow(z-2, 2))
}

func kFunc(x float64, z float64) float64 {
	return 0.536
}

type solvedLineStruct struct {
	index  int
	result []float64
}

func solveLine(prevULine []float64, funcLine []float64, kVal float64, startVal float64, endVal float64, hStep float64, tStep float64) []float64 {
	n := len(prevULine)
	topDiag := make([]float64, n)
	midDiag := make([]float64, n)
	botDiag := make([]float64, n)
	bVector := make([]float64, n)
	midDiag[0] = 1
	midDiag[len(midDiag)-1] = 1
	bVector[0] = startVal
	bVector[len(bVector)-1] = endVal
	for idx := 1; idx < n-1; idx++ {
		topDiag[idx] = -1
		midDiag[idx] = 2 + (hStep*hStep)/tStep/kVal
		botDiag[idx] = -1
		bVector[idx] = prevULine[idx-1]
		bVector[idx] += -prevULine[idx] * 2
		bVector[idx] += prevULine[idx+1]
		bVector[idx] += prevULine[idx] * (hStep * hStep) / tStep / kVal
		bVector[idx] += funcLine[idx] / kVal * (hStep * hStep)
	}
	res := helpers.ThomasAlgorithm(botDiag, midDiag, topDiag, bVector)
	return res
}

func parallelSolve(xBound float64, zBound float64, hStep float64, tStep float64, iterations int) [][]float64 {
	nX := int(xBound / hStep)
	nZ := int(zBound / hStep)

	originU := make([][]float64, nX+1)
	for idx, _ := range originU {
		originU[idx] = make([]float64, nZ+1)
	}
	prevU := make([][]float64, nX+1)
	for idx, _ := range prevU {
		prevU[idx] = make([]float64, nZ+1)
	}
	currentU := make([][]float64, nX+1)
	for idx, _ := range currentU {
		currentU[idx] = make([]float64, nZ+1)
	}
	for iteration := 0; iteration < iterations; iteration++ {
		for i := 0; i < nZ; i++ { //set boundaries
			currentU[i][0] = boundary(0, 0)
			currentU[i][len(currentU[i])-1] = boundary(0, 0)
		}
		for j := 0; j < nX; j++ { //set boundaries
			currentU[0][j] = boundary(0, 0)
			currentU[len(currentU)-1][j] = boundary(0, 0)
		}

		resChan := make(chan solvedLineStruct, nX)

		for j := 1; j < nX; j++ { // by cols
			prevUCol := make([]float64, nX+1)
			for i, v := range prevU {
				prevUCol[i] = v[j]
			}
			funcValCol := make([]float64, nZ)
			helpers.FillSlice(0, nZ, hStep, funcValCol, func(z float64) float64 {
				return mainFunc(float64(j)*hStep, z)
			})
			kVal := kFunc(0, 0)
			startVal := boundary(0, 0)
			endVal := boundary(0, 0)
			go func(index int) {
				res := solveLine(prevUCol, funcValCol, kVal, startVal, endVal, hStep, tStep)
				resChan <- solvedLineStruct{
					index:  index,
					result: res,
				}
			}(j)
		}
		for p := 0; p < nX-1; p++ {
			solved := <-resChan
			for idx, val := range solved.result {
				currentU[idx][solved.index] = val
			}
		}

		for i, _ := range currentU { // copy currentU to prevU
			for j, _ := range currentU[i] {
				prevU[i][j] = currentU[i][j]
			}
		}

		resChan = make(chan solvedLineStruct, nZ)

		for i := 1; i < nZ; i++ { // by rows
			prevURow := make([]float64, nX+1)
			for j, v := range prevU[i] {
				prevURow[j] = v
			}
			funcValRow := make([]float64, nX+1)
			helpers.FillSlice(0, nX+1, hStep, funcValRow, func(x float64) float64 {
				return mainFunc(x, float64(i)*hStep)
			})
			kVal := kFunc(0, 0)
			startVal := boundary(0, 0)
			endVal := boundary(0, 0)
			go func(index int) {
				res := solveLine(prevURow, funcValRow, kVal, startVal, endVal, hStep, tStep)
				resChan <- solvedLineStruct{
					index:  index,
					result: res,
				}
			}(i)
		}
		for p := 0; p < nZ-1; p++ {
			solved := <-resChan
			for idx, val := range solved.result {
				currentU[solved.index][idx] = val
			}
		}

		var maxErr float64 // calc max error
		for i, _ := range currentU {
			for j, _ := range currentU[i] {
				err := math.Abs(originU[i][j] - currentU[i][j])
				if err > maxErr {
					maxErr = err
				}
			}
		}

		fmt.Printf("\rIteration: %d/%d; U(s+1) - U(s): %0.5f", iteration, iterations, maxErr)

		for i, _ := range currentU { // copy currentU to prevU
			for j, _ := range currentU[i] {
				prevU[i][j] = currentU[i][j]
			}
		}
		for i, _ := range currentU { // copy currentU to originU
			for j, _ := range currentU[i] {
				originU[i][j] = currentU[i][j]
			}
		}
	}
	return currentU
}

//export solver
func solver() *C.char {
	hStep := 0.05
	start := time.Now()
	res := parallelSolve(10, 10, hStep, 0.01, 2000)
	elapsed := time.Since(start)
	fmt.Printf("\nSolving took %s", elapsed)
	n := len(res)
	m := len(res[0])
	resMat := mat.NewDense(n, m, nil)
	for i := 0; i < n; i++ {
		for j := 0; j < m; j++ {
			resMat.Set(i, j, res[i][j])
		}
	}

	fa := mat.Formatted(resMat, mat.FormatPython())
	return C.CString(fmt.Sprintf("%#v\n", fa))
}

func main() {
	solver()
}

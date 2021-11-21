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

func solveLine(prevU3Lines [][3]float64, funcLine []float64, kVal float64, startVal float64, endVal float64, hStep float64, tStep float64) []float64 {
	n := len(prevU3Lines)
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
		bVector[idx] = prevU3Lines[idx][0]
		bVector[idx] += -prevU3Lines[idx][1] * 2
		bVector[idx] += prevU3Lines[idx][2]
		bVector[idx] += prevU3Lines[idx][1] * (hStep * hStep) / tStep / kVal
		bVector[idx] += funcLine[idx] / kVal * (hStep * hStep)
	}
	res := helpers.ThomasAlgorithm(botDiag, midDiag, topDiag, bVector)
	return res
}

func parallelSolve(xBound float64, zBound float64, hStep float64, tStep float64, iterations int, processes int) [][]float64 {
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
		for i := 0; i < nZ+1; i++ { //set boundaries
			currentU[i][0] = boundary(0, 0)
			currentU[i][len(currentU[i])-1] = boundary(0, 0)
		}
		for j := 0; j < nX+1; j++ { //set boundaries
			currentU[0][j] = boundary(0, 0)
			currentU[len(currentU)-1][j] = boundary(0, 0)
		}

		resChan := make(chan solvedLineStruct, nX)
		goroutinesCount := processes
		batchSize := nX/goroutinesCount + 1
		for g := 0; g < goroutinesCount; g++ {
			startIndex := g * batchSize
			endIndex := startIndex + batchSize
			if startIndex == 0 {
				startIndex = 1
			}
			if endIndex > nX {
				endIndex = nX
			}
			if endIndex < startIndex {
				continue
			}
			go func() {
				buffer := make([]solvedLineStruct, endIndex-startIndex)
				bufIdx := 0
				for j := startIndex; j < endIndex; j++ { // by cols
					prevU3Cols := make([][3]float64, nX+1)
					for i := 0; i < len(prevU3Cols); i++ {
						prevU3Cols[i][0] = prevU[i][j-1]
						prevU3Cols[i][1] = prevU[i][j]
						prevU3Cols[i][2] = prevU[i][j+1]
					}
					funcValCol := make([]float64, nZ)
					helpers.FillSlice(0, nZ, hStep, funcValCol, func(z float64) float64 {
						return mainFunc(float64(j)*hStep, z)
					})
					kVal := kFunc(0, 0)
					startVal := boundary(0, 0)
					endVal := boundary(0, 0)
					res := solveLine(prevU3Cols, funcValCol, kVal, startVal, endVal, hStep, tStep)
					buffer[bufIdx] = solvedLineStruct{
						index:  j,
						result: res,
					}
					bufIdx++
				}
				for _, v := range buffer {
					resChan <- v
				}
			}()
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
		batchSize = nZ/goroutinesCount + 1
		for g := 0; g < goroutinesCount; g++ {
			startIndex := g * batchSize
			endIndex := startIndex + batchSize
			if startIndex == 0 {
				startIndex = 1
			}
			if endIndex > nZ {
				endIndex = nZ
			}
			if endIndex < startIndex {
				continue
			}
			go func() {
				buffer := make([]solvedLineStruct, endIndex-startIndex)
				bufIdx := 0
				for i := startIndex; i < endIndex; i++ { // by rows
					prevU3Rows := make([][3]float64, nX+1)
					for j := 0; j < len(prevU3Rows); j++ {
						prevU3Rows[j][0] = prevU[i-1][j]
						prevU3Rows[j][1] = prevU[i][j]
						prevU3Rows[j][2] = prevU[i+1][j]
					}
					funcValRow := make([]float64, nX+1)
					helpers.FillSlice(0, nX+1, hStep, funcValRow, func(x float64) float64 {
						return mainFunc(x, float64(i)*hStep)
					})
					kVal := kFunc(0, 0)
					startVal := boundary(0, 0)
					endVal := boundary(0, 0)
					res := solveLine(prevU3Rows, funcValRow, kVal, startVal, endVal, hStep, tStep)
					buffer[bufIdx] = solvedLineStruct{
						index:  i,
						result: res,
					}
					bufIdx++
				}
				for _, v := range buffer {
					resChan <- v
				}
			}()
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
func solver(xBound float64, zBound float64, hStep float64, tStep float64, iterations int, processes int) *C.char {
	start := time.Now()
	res := parallelSolve(xBound, zBound, hStep, tStep, iterations, processes)
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
	solver(10, 10, 0.01, 0.01, 300, 8)
}

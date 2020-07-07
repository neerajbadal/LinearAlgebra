'''
Created on 04-Sep-2019

@author: Neeraj Badal
'''
def readFile(filePath):
    mat_A = list()
    with open(filePath) as f: 
        for line in f:
             
            rowA = str(line).splitlines()[0].split(sep=" ")
            rowA = [ int(x) for x in rowA ]
            mat_A.append(rowA)

    
    return mat_A
def getDimensionOfMatrix(mat_A):
    noOfRows = len(mat_A)
    if noOfRows > 0:
        noOfColumns = len(mat_A[0])
    else:
        noOfColumns = 0
    return [noOfRows,noOfColumns]

def generateIdentityMatrix(rows,columns):
    mat_i = list()
    for i in range(0,rows):
        row_i = list()
        for j in range(0,columns):
            if(i==j):
                row_i.append(1)
            else:
                row_i.append(0)
        mat_i.append(row_i)
    return mat_i
def getElementaryMatrix(mat_A,rowIndex,colIndex,pivotRow,pivotColumn):
    dim_A = getDimensionOfMatrix(mat_A)
    mat_e = generateIdentityMatrix(dim_A[0],dim_A[0])
    changeVal = mat_A[rowIndex][colIndex] / mat_A[pivotRow][pivotColumn]
    mat_e[rowIndex][pivotRow] = changeVal
    mat_e[rowIndex][rowIndex] = -1
    
    return mat_e

def getElementaryMatrixBackward(mat_A,rowIndex,colIndex,pivotRow,pivotColumn):
    dim_A = getDimensionOfMatrix(mat_A)
    mat_e = generateIdentityMatrix(dim_A[0],dim_A[0])
    changeVal = mat_A[rowIndex][colIndex]
#     if changeVal < 0:
#         mat_e[rowIndex][pivotRow] = -changeVal
#     else:
    mat_e[rowIndex][pivotRow] = -changeVal
    mat_e[rowIndex][rowIndex] = 1
    
    return mat_e


def getElementaryMatrixForSelfChange(mat_A,rowIndex,colIndex):
    dim_A = getDimensionOfMatrix(mat_A)
    mat_e = generateIdentityMatrix(dim_A[0],dim_A[0])
    changeVal = 1 / mat_A[rowIndex][colIndex]
    mat_e[rowIndex][rowIndex] = changeVal
   
    return mat_e


def matrixRowOperation(mat_A,rowIndex,rowChangeVector):
    dim_A = getDimensionOfMatrix(mat_A)
#     print("row change vector ",rowChangeVector)
    mat_A[rowIndex] = [i * rowChangeVector[rowIndex] for i in mat_A[rowIndex]]
#     print(mat_A[rowIndex])
    for row in range(dim_A[0]):
            if row != rowIndex:
                temp_row = [i * rowChangeVector[row] for i in mat_A[row]]
                mat_A[rowIndex] = [round(sum(x),5) for x in zip(temp_row, mat_A[rowIndex])]
            
#     print(temp_row)
    return mat_A

# from scipy.linalg import lu
# import numpy as np
# # M = np.array([[5,3,1], [6,15,13], [0,25,18]])
# M = np.array([[1,3,1], [5,10,2], [1,8,9]])
# p,l,u = lu(M)
# print(u)
#    
# exit(0)
if __name__ == "__main__":
        mat_A = readFile("task1_input.dat")
#         mat_A = [[5,3,1], [6,15,13], [0,25,18]]
        dim_A = getDimensionOfMatrix(mat_A)
        print(mat_A)
        print("dimension of matrix ",dim_A)
        mat_i = generateIdentityMatrix(dim_A[0],dim_A[1])
        
        print(mat_i)
        prevPivotRow = 0
        prevPivotCol = 0
        pivotRow = 0
        pivotCol = 0
        for col in range(0,dim_A[1]):
            if col < dim_A[0]:
                nonZeroCount = 0
                print(" col pos ",col)
                for row in range(prevPivotRow,dim_A[0]):
                    if(mat_A[row][col] != 0):
                        pivotRow = row
                        nonZeroCount = nonZeroCount + 1
                        break;
                if nonZeroCount == 0:
                    continue
                
                if pivotRow != prevPivotCol:
                    tempRow = mat_A[prevPivotCol]
                    mat_A[prevPivotCol] = mat_A[pivotRow]
                    mat_A[pivotRow] = tempRow
                pivotRow = prevPivotRow
               
                for row in range(pivotRow+1,dim_A[0]):
                    if mat_A[row][col] !=0:
                        elem_mat = getElementaryMatrix(mat_A,row,col,pivotRow,col)
                        print(row,col," ",elem_mat)
                        print("before ",mat_A)
                        mat_A = matrixRowOperation(mat_A, row,elem_mat[row])
                        print("after",mat_A)
                prevPivotRow = prevPivotRow + 1
                prevPivotCol = prevPivotCol + 1
            print(pivotRow,"  ",pivotCol)
            for row in range(0,dim_A[0]):
                for col in range(0,dim_A[1]):
                    print(mat_A[row][col],end='\t')
                print()
            print("----------------------------------------")
#             break
#             print(mat_A)     
#         print(mat_A)
        
        pivotColumnList = list()
        for row in range(0,dim_A[0]):
            for col in range(0,dim_A[1]):
                if mat_A[row][col] !=0:
                    elem_mat = getElementaryMatrixForSelfChange(mat_A, row, col)
                    mat_A = matrixRowOperation(mat_A, row,elem_mat[row])
                    pivotColumnList.append(col)
                    break
        
        print("Pivot Column List ",pivotColumnList)
        for row in range(0,dim_A[0]):
                for col in range(0,dim_A[1]):
                    print(mat_A[row][col],end='\t')
                print()
        pivotCount = 0        
        
        for col in pivotColumnList:
            for row in range(0,pivotCount):
                if mat_A[row][col] !=0:
                        elem_mat = getElementaryMatrixBackward(mat_A,row,col,pivotCount,col)
                        print(row,col," ",elem_mat)
#                         print("before ",mat_A)
                        mat_A = matrixRowOperation(mat_A, row,elem_mat[row])
#                         print("after",mat_A)
                        print("------------------------------")
                        for rowi in range(0,dim_A[0]):
                                for coli in range(0,dim_A[1]):
                                    print(mat_A[rowi][coli],end='\t')
                                print()
            
            pivotCount = pivotCount + 1
        print("------------------------------")
        for row in range(0,dim_A[0]):
                for col in range(0,dim_A[1]):
                    print(round(mat_A[row][col],2),end='\t')
                print()
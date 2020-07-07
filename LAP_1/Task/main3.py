'''
Created on 04-Sep-2019

@author: Neeraj Badal
'''
'''
Created on 20-Aug-2019

@author: neerajbadal
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
def getElementaryMatrix(mat_A,rowIndex,colIndex):
    dim_A = getDimensionOfMatrix(mat_A)
    mat_e = generateIdentityMatrix(dim_A[0],dim_A[1])
    changeVal = mat_A[rowIndex][colIndex] / mat_A[colIndex][colIndex]
    mat_e[rowIndex][colIndex] = changeVal
    mat_e[rowIndex][rowIndex] = -1
    
    return mat_e

def matrixRowOperation(mat_A,rowIndex,rowChangeVector):
    dim_A = getDimensionOfMatrix(mat_A)
    mat_A[rowIndex] = [i * rowChangeVector[rowIndex] for i in mat_A[rowIndex]]
    print(mat_A[rowIndex])
    for col in range(dim_A[1]):
            if col != rowIndex:
                temp_row = [i * rowChangeVector[col] for i in mat_A[col]]
                mat_A[rowIndex] = [sum(x) for x in zip(temp_row, mat_A[rowIndex])]
            
#     print(temp_row)
    return mat_A

# from scipy.linalg import lu
# import numpy as np
# M = np.array([[5,3,1], [6,15,13], [0,25,18]])
# p,l,u = lu(M)
# print(u)
#  
# exit(0)
if __name__ == "__main__":
#         mat_A = readFile("task1_input.dat")
        mat_A = [[5,3,1], [6,15,13], [0,25,18]]
        dim_A = getDimensionOfMatrix(mat_A)
        print(mat_A)
        print("dimension of matrix ",dim_A)
        mat_i = generateIdentityMatrix(dim_A[0],dim_A[1])
        
        print(mat_i)
        
        for j in range(dim_A[1]):
            for i in range(dim_A[0]):
                if (i>j ):
                    
                    elem_mat = getElementaryMatrix(mat_A,i,j)
                    print(i,j," ",elem_mat)
                    mat_A = matrixRowOperation(mat_A, i,elem_mat[i])
                    print(mat_A)
#                     print(matrixRowOperation(mat_A, i, getElementaryMatrix(mat_i,mat_A,i,j)[i]))
#                     exit(0)
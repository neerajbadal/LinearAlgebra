'''
Created on 04-Sep-2019

@author: Neeraj Badal
'''
import copy
import sys
def readFile(filePath):
    '''
    reads input file from filePath and returns matrix correspondingly
    filePath : absolute file path containing matrix/system of linear equations\
    returns:
    mat_A : 2-D corresponding matrix
    '''
    mat_A = list()
    with open(filePath) as f: 
        for line in f:
             
            rowA = str(line).splitlines()[0].split(sep=" ")
            rowA = [ float(x) for x in rowA ]
            mat_A.append(rowA)

    
    return mat_A

def getLCM(x_,y_):
    '''
    computes LCM of two numbers
    '''
    a = x_
    b = y_
    while(b!=0):
        t = b
        b = a%b
        a = t
    gcd = a
    return (x_*y_)/gcd
        


def printMatrix(mat_A_):
    '''
    provides custom print formatting to print the matrix
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    print("----------------------------------------")
    for row in range(0,dim_A_[0]):
#         mat_A[row] = [round(x,3) for x in mat_A[row]]
#         print(mat_A[row],"|")        
                for col in range(0,dim_A_[1]):
                    print(truncate(mat_A_[row][col],3),end='\t')
                print(end='|')    
                print()
    print("----------------------------------------")

def getDimensionOfMatrix(mat_A_):
    '''
    returns dimensions of the input matrix,
    as a list having two elemnts [noOfRows,noOfColumns]
    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    noOfRows = len(mat_A_)
    if noOfRows > 0:
        noOfColumns = len(mat_A_[0])
    else:
        noOfColumns = 0
    return [noOfRows,noOfColumns]

def truncate(val,decimals=0):
    multiplier = 10 ** decimals
    return int(val * multiplier) / multiplier

def generateIdentityMatrix(rows,columns):
    '''
    Generates Identity Matrix
    returns:
    Identity Matrix
    '''
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
def getElementaryMatrix(mat_A_,rowIndex,colIndex,pivotRow,pivotColumn):
    '''
    computes elementary matrix for a single elementary row operation
    mat_A_ : input matrix
    rowIndex : index of the row being modified
    colIndex : used for matching element with pivotColumn Index
    pivotRow : Row index of the pivot element
    pivotColumn : Column index of the pivot element 
    
    returns:
    permutation matrix signifying elementary row operation    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    mat_e = generateIdentityMatrix(dim_A_[0],dim_A_[0])
#     lcm = getLCM(mat_A_[rowIndex][colIndex], mat_A_[pivotRow][pivotColumn])
    changeVal = mat_A_[rowIndex][colIndex] / mat_A_[pivotRow][pivotColumn]
#     changeVal = round(changeVal,3)
#     mat_e[rowIndex][pivotRow] = changeVal
#     mat_e[rowIndex][rowIndex] = -1
    
    if mat_A_[pivotRow][pivotColumn] < mat_A_[rowIndex][colIndex]:
#         mat_e[rowIndex][pivotRow] = lcm / mat_A_[pivotRow][pivotColumn]
#         mat_e[rowIndex][rowIndex] = -1*(lcm/mat_A_[rowIndex][colIndex])
        mat_e[rowIndex][pivotRow] = changeVal
        mat_e[rowIndex][rowIndex] = -1  
    else:
#         mat_e[rowIndex][pivotRow] = -1 * lcm / mat_A_[pivotRow][pivotColumn]
#         mat_e[rowIndex][rowIndex] = (lcm/mat_A_[rowIndex][colIndex])
        mat_e[rowIndex][pivotRow] = -1 * changeVal
        mat_e[rowIndex][rowIndex] = 1
    return mat_e

def getElementaryMatrixBackward(mat_A_,rowIndex,colIndex,pivotRow,pivotColumn):
    '''
    Convert RowEchleon form matrix to row-reduced Echleon Matrix 
    using BackPropagation
    mat_A_ : input matrix
    rowIndex : index of the row being modified
    colIndex : used for matching element with pivotColumn Index
    pivotRow : Row index of the pivot element
    pivotColumn : Column index of the pivot element 
    
    returns:
    permutation matrix signifying elementary row operation
    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    mat_e = generateIdentityMatrix(dim_A_[0],dim_A_[0])
    changeVal = mat_A_[rowIndex][colIndex]
#     if changeVal < 0:
#         mat_e[rowIndex][pivotRow] = -changeVal
#     else:
    mat_e[rowIndex][pivotRow] = -changeVal
    mat_e[rowIndex][rowIndex] = 1
    
    return mat_e


def getElementaryMatrixForSelfChange(mat_A_,rowIndex,colIndex):
    '''
    convert all the pivot element to be equal to 1.
    mat_A_ : input matrix
    rowIndex : index of the row being modified
    colIndex : pivotColumn Index
    
    returns:
    permutation matrix signifying elementary row operation
    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    mat_e = generateIdentityMatrix(dim_A_[0],dim_A_[0])
    changeVal = 1 / mat_A_[rowIndex][colIndex]
#     changeVal = truncate(changeVal,3)
    mat_e[rowIndex][rowIndex] = changeVal
   
    return mat_e


def matrixRowOperation(mat_A_,rowIndex,rowChangeVector):
    '''
    computes elementary row operation on the specified row 
    mat_A_ : input matrix
    rowIndex : index of the row being modified
    rowChangeVector: implied in the name
    returns:
        modified mat_A_ based on the row operation
    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    mat_A_[rowIndex] = [i * rowChangeVector[rowIndex] for i in mat_A_[rowIndex]]
    rowToBeModified = []
    for col in range(0,len(rowChangeVector)):
        if rowChangeVector[col] != 0:
            rowToBeModified.append(col)
    for row in rowToBeModified:
            if row != rowIndex:
                temp_row = [i * rowChangeVector[row] for i in mat_A_[row]]
#                 mat_A_[rowIndex] = [round(sum(x),5) for x in zip(temp_row, mat_A_[rowIndex])]
                mat_A_[rowIndex] = [truncate(sum(x),3) for x in zip(temp_row, mat_A_[rowIndex])]
    return mat_A_


def getAbsolute(val):
    if val >= 0:
        return val
    else:
        return (val * -1)

def getRowEchleonForm(mat_A_,elementaryMatList):
    '''
    Perform Row Echleon Form on the input Matrix
    mat_A_ : input matrix
    elemetaryMatList : empty list in which series of elemtary row operations
                      performed will be appended.
    
    returns:
        2-element tuple where first element is the modified matrix and,
        second element is the filled elemtaryMatList
    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    prevPivotRow = 0
    prevPivotCol = 0
    pivotRow = 0
    pivotCol = 0
    
    for col in range(0,dim_A_[1]):
        nonZeroCount = 0
        maxPivotValue = -1
        for row in range(prevPivotRow,dim_A_[0]):
            if(mat_A_[row][col] != 0 and getAbsolute(mat_A_[row][col]) > maxPivotValue):
                pivotRow = row
                maxPivotValue = getAbsolute(mat_A_[row][col])
                nonZeroCount = nonZeroCount + 1
#                 break;
        if nonZeroCount == 0:
            '''
            current being zero. Algorithm will
            shift to next column.
            '''
            continue
         
        if pivotRow != prevPivotCol:
            '''
            pivot element of this column is not the
            leading term, so swap it to make it the 
            leading term
            '''
            tempRow = mat_A_[prevPivotRow]
            mat_A_[prevPivotRow] = mat_A_[pivotRow]
            mat_A_[pivotRow] = tempRow
            permMatrix = generateIdentityMatrix(dim_A_[0],dim_A_[0])
            permMatrix[prevPivotRow][prevPivotRow] = 0
            permMatrix[prevPivotRow][pivotRow] = 1
            
            permMatrix[pivotRow][prevPivotRow] = 1
            permMatrix[pivotRow][pivotRow] = 0
              
            elementaryMatList.append(permMatrix)
        '''
        Reset row to start from pivot row + 1
        '''
        pivotRow = prevPivotRow
        for row in range(pivotRow+1,dim_A_[0]):
            if mat_A_[row][col] !=0:
                elem_mat = getElementaryMatrix(mat_A_,row,col,pivotRow,col)
                elementaryMatList.append(elem_mat)
#                 print(mat_A_)
                mat_A_ = matrixRowOperation(mat_A_, row,elem_mat[row])
        prevPivotRow = prevPivotRow + 1
        prevPivotCol = prevPivotCol + 1
    return (mat_A_,elementaryMatList)
         
def convertLeadingTermsToUnity(mat_A_):
    '''
    Driver to make pivot element to equal to unity of
    the input row echleon matrix mat_A_
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    for row in range(0,dim_A_[0]):
        for col in range(0,dim_A_[1]):
            if mat_A_[row][col] !=0:
                elem_mat = getElementaryMatrixForSelfChange(mat_A_, row, col)
                mat_A_ = matrixRowOperation(mat_A_, row,elem_mat[row])
                break
    return mat_A_
def getPivotColumnList(mat_A_):
    '''
    returns a list having column index of the pivot element in each row.
    for a row not having pivot -1 is assigned.
    
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    pivotColumnList = list()
    dim_A_ = getDimensionOfMatrix(mat_A_)
    for row in range(0,dim_A_[0]):
        nonZeroCounter = 0
        for col in range(0,dim_A_[1]):
            if mat_A_[row][col] !=0:
#                 elem_mat = getElementaryMatrixForSelfChange(mat_A_, row, col)
#                 mat_A_ = matrixRowOperation(mat_A_, row,elem_mat[row])
                nonZeroCounter = nonZeroCounter + 1
                pivotColumnList.append(col)
                break
        if nonZeroCounter == 0:
            pivotColumnList.append(-1)
    return pivotColumnList
  
def getMatrixRank(mat_A_):
    '''
    returns rank of the matrix
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    pivotColumnList = getPivotColumnList(mat_A_)
    x = pivotColumnList.count(-1)
    return (len(pivotColumnList) - x)

def getSolution(mat_A_):
    '''
    Do back-substitution on the row-reduced echleon form
    matrix to get solution for the system.
    returns:
    solution vector in both the cases for Unique and Infinite 
    solution 
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    solutionVector = list()
    dim_A_ = getDimensionOfMatrix(mat_A_)
#     pivotColumnList = getPivotColumnList(mat_A_)
    for index in range(0,dim_A_[1]-1):
        solutionVector.append(0)
    
    freeColumns = [x for x in range(0,dim_A_[1]) if x not in pivotColumnList]
    freeColumnUsed = freeColumns[0]
    for row in range(dim_A_[0]-1,-1,-1):
        value = mat_A_[row][dim_A_[1]-1]
        if pivotColumnList[row] != -1:
            for col in range(pivotColumnList[row]+1,dim_A_[1]-1):
                
                if col in pivotColumnList:
                    value = value - (mat_A_[row][col] * solutionVector[col])
                elif col == freeColumnUsed:
                    value = value - mat_A_[row][col]
                    solutionVector[col] = 1
                    break;
        
            solutionVector[pivotColumnList[row]] = truncate(value,3)    
    return solutionVector

if __name__ == "__main__":
        if len (sys.argv) != 3 :
            print("Usage: python main.py problem# <absoluteInputDataPath>")
            sys.exit(1)
        problemIndex = sys.argv[1]
        inputFile = sys.argv[2]
        mat_A = readFile(inputFile)
#         mat_A = [[5,3,1], [6,15,13], [0,25,18]]
        dim_A = getDimensionOfMatrix(mat_A)
        
        if problemIndex == "problem1":
            elemetaryList = list()
            modified_mat_A,elemetaryList = getRowEchleonForm(mat_A,elemetaryList)
            print("RANK OF THE MATRIX: ",getMatrixRank(modified_mat_A))
            print("ROW ECHELON FORM:")
            printMatrix(modified_mat_A)
            print("SEQUENCE OF ELEMENTARY MATRICES USED:")
#             for row in range(0,dim_A[0]):
#                 for listIndex in range(len(elemetaryList)-1,-1,-1):
#                     print(end='|')
#                     for col in range(0,dim_A[0]):    
#                         print("%.3f" % round(elemetaryList[listIndex][row][col],3),end="\t")
#                     print(end="| ")
#                 print()
            for listIndex in range(len(elemetaryList)-1,-1,-1):
                print(elemetaryList[listIndex],end="\t")
            
        elif problemIndex == "problem2":
            elemetaryList = list();
            modified_mat_A,elemetaryList = getRowEchleonForm(mat_A,elemetaryList)
            modified_mat_A = convertLeadingTermsToUnity(modified_mat_A)
            pivotColumnList = getPivotColumnList(modified_mat_A)
            mat_A = modified_mat_A
            pivotCount = 0
            
            for col in pivotColumnList:
                if col != -1:
                    pivotCount = pivotCount + 1
                    for row in range(0,pivotCount):
                        if mat_A[row][col] !=0:
                                elem_mat = getElementaryMatrixBackward(mat_A,row,col,pivotCount-1,col)
                                mat_A = matrixRowOperation(mat_A, row,elem_mat[row])
            stripedMat_A = copy.deepcopy(mat_A)               
            for row in range(0,dim_A[0]):
                del stripedMat_A[row][dim_A[1]-1]
            
            origMat_A_Rank = getMatrixRank(stripedMat_A)
            augmentedMatrixRank = getMatrixRank(mat_A)
            origMat_ADim = getDimensionOfMatrix(stripedMat_A)
            if  augmentedMatrixRank > origMat_A_Rank:
                print("NO SOLUTION EXISTS !")
            elif augmentedMatrixRank == origMat_A_Rank and augmentedMatrixRank == origMat_ADim[1]:
                print("UNIQUE SOLUTION EXISTS !")
                print(getSolution(mat_A))
            elif augmentedMatrixRank == origMat_A_Rank and augmentedMatrixRank < origMat_ADim[1]:
                print("MANY SOLUTIONS EXISTS !")
                print(getSolution(mat_A))
                
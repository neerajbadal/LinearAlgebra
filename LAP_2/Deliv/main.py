'''
Created on 10-Nov-2019

@author: Neeraj Badal
'''
import numpy as np
import copy
import sys
import networkx as nx
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import collections
from sklearn.manifold import TSNE

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

def getDotProduct(vec_1_,vec_2_):
    '''
    returns dot product of 2 given vectors
    '''
    dot_sum = 0.0
    vec_1 = copy.deepcopy(vec_1_)
    vec_2 = copy.deepcopy(vec_2_)
    
    dot_sum = sum((float(vec_1[index_i]) * float(vec_2[index_i])) for index_i in range(0,len(vec_1)))
    return dot_sum  

def getNorm2(vec_1):
    '''
    returns L2 norm of the given vector.
    '''
    norm_square = getDotProduct(vec_1, vec_1)
    return (norm_square**(0.5))

def getMatrixVectorProduct(mat_A_,vec_1,n_rows=3,n_cols=3):
    '''
    returns vector obtained from the matrix vector product
    '''
    return [getDotProduct(row_, vec_1) for row_ in mat_A_]
            
    

def getDimensionOfMatrix(mat_A_):
    '''
    returns dimensions of the input matrix,
    as a list having two elemnts [noOfRows,noOfColumns]
    
    '''
#     mat_A_ = copy.deepcopy(mat_A_)
    noOfRows = len(mat_A_)
    if noOfRows > 0:
        noOfColumns = len(mat_A_[0])
    else:
        noOfColumns = 0
    return [noOfRows,noOfColumns]

def getMaxAbsoluteNorm(vec_1):
    maxVal = -9999999.9
    for vec_in in range(0,len(vec_1)):
        if maxVal < abs(vec_1[vec_in]):
            maxVal = abs(vec_1[vec_in])
    return maxVal
    
def getLargestEigenPair(mat_A_):
    '''
    returns largest eigen pair computed using power method
    output : [eigenvalue,eigenvector]
    '''
    initial_vec = list()
    dim_A_ = getDimensionOfMatrix(mat_A_)
    n_cols = dim_A_[1]
    for i in range(n_cols):
        initial_vec.append(1.0)#n_cols - i)

    eigVal = 1
    iteration_ = 1000
    tolerance = 0.0000000000001
    new_vec = initial_vec
    signMultiplier = 1
    for iter_i in range(0,iteration_):
        old_vec = copy.deepcopy(new_vec)
        old_eigVal = eigVal
        matDim = getDimensionOfMatrix(mat_A_)
        
        new_vec = getMatrixVectorProduct(mat_A_, new_vec, matDim[0], matDim[1])#np.dot(mat_A_,new_vec)
        eigVal =  getMaxAbsoluteNorm(new_vec)#max(abs(new_vec))
        
        for index_ in range(0,len(new_vec)):
            if (new_vec[index_] * old_vec[index_] >= 0):
                signMultiplier = 1
            elif(new_vec[index_] * old_vec[index_] < 0):
                signMultiplier = -1
        
        for index_ in range(0,len(new_vec)):
            new_vec[index_] = new_vec[index_] / float(eigVal)
        
        eigValChange = abs(eigVal-old_eigVal)/float(eigVal)
        if iter_i > 1:
            if eigValChange <= tolerance:
                break
        
        
    new_vec_norm = getNorm2(new_vec)#np.sqrt(np.dot(new_vec.T,new_vec))
    
    new_vec = [(elem/new_vec_norm) for elem in new_vec]
    
    return [eigVal*signMultiplier,new_vec]

def performRank1Mul(vec_1,vec_2):
    '''
    returns matrix generated from 2 rank-1 matrix multiplication
    '''
    vec_1_len = len(vec_1)
    vec_2_len = len(vec_2)
    mod_mat_B = list()
    for row_i in range(0,vec_1_len):
        row_list = list()
        for col_i in range(0,vec_2_len):
            
            temVal = vec_1[row_i] * vec_2[col_i]
            row_list.append(temVal)
        mod_mat_B.append(row_list)
    
    return mod_mat_B

def truncate(val,decimals=0):
    multiplier = 10 ** decimals
    return int(val * multiplier) / multiplier



def all_pair_shortest_floyd(mat_c):
    '''
    returns the shortest path distance between each pair of vertices in the graph
    '''
    mat_c = copy.deepcopy(mat_c)
    dim_c = getDimensionOfMatrix(mat_c)
    dist = mat_c#map(lambda i : map(lambda j : j , i) , mat_c)
    INF = 99999
    for l_rows in range(0,dim_c[0]):
        for l_cols in range(0,dim_c[1]):
            if dist[l_rows][l_cols] == 0 and l_rows != l_cols :
                dist[l_rows][l_cols] = INF 
    
    
                
    V = dim_c[0]
    for k in range(0,V): 
        for i in range(0,V): 
            for j in range(0,V): 
                dist[i][j] = min(dist[i][j],dist[i][k]+ dist[k][j]) 

    return dist


def generatePath(pathDict,key):
    '''
    a recursive function to trace the path for all shortest path between a source to
    all vertices in the graph.
    returns the list of vertices that come in the shortest-path between the source and
    all other vertex
    
    '''
    
    newVertexList = []
    
    vals = pathDict.get(key)
    if vals is not None:
        for iter_ii in range(0,len(vals)):
            newVertexList.extend([key[1]])
        for new_val in vals:
            sub_val = pathDict.get((key[0],new_val))
            if sub_val is not None:
                for iter_iii in range(0,len(sub_val)-1):
                    newVertexList.extend([key[1]])
            newVertexList.extend(generatePath(pathDict, (key[0],new_val)))
    
    return newVertexList
    
def modBFS(mat_c):
    '''
    Modified version of BFS to calculate betweeness centrality measure of 
    vertices in the given graph
    
    returns the between centrality measure for each vertex in the graph
    
    '''
    dim_c = getDimensionOfMatrix(mat_c)
    betweenNessCentrality = list()
    
    vertexInPathList = list()
    shortestPathCounts = list()
    for src in range(0,dim_c[0]):
        visited = list()
        INF = 99999
        dist = list()
        paths = list()
        im_vertex = list()
        d = collections.defaultdict(list)
        for i_c in range(0,dim_c[0]):
            visited.append(0) 
            dist.append(INF)
            paths.append(0)
            im_vertex.append(0)
        dist[src] = 0 
        paths[src] = 1
      
        q = list() 
        q.append(src); 
        visited[src] = 1
        while len(q)!=0: 
        
            curr = q[0]
            q.pop(0) 
      
            for curr_i in range(0,dim_c[0]):
                if mat_c[curr][curr_i] !=0:
                 
                    if (visited[curr_i] == 0): 
                        q.append(curr_i); 
                        visited[curr_i] = 1
                    
                    if (dist[curr_i] > dist[curr] + 1): 
                        dist[curr_i] = dist[curr] + 1; 
                        paths[curr_i] = paths[curr]
                        d[src,curr_i].append(curr) 
                    
                    elif (dist[curr_i] == dist[curr] + 1): 
                        paths[curr_i] += paths[curr]
                        d[src,curr_i].append(curr)
                    
                   
        pathVertexDict = copy.deepcopy(d)
        for d_keys in pathVertexDict.keys():
            pathVertexDict[d_keys].clear()
       
        for d_keys in d.keys():
            pathVertexDict[d_keys].extend(generatePath(d, d_keys))

        tempVertexList = list()
        for index_ in range(0,dim_c[0]):
            intermVertex = pathVertexDict.get((src,index_))
            vertexCount = dict()
            vertexCount[index_] = 1
            if intermVertex is not None:
                for vert_ in intermVertex:
                    if vertexCount.get(vert_) is not None:
                        vertexCount[vert_] += 1
                    else:
                        vertexCount[vert_] = 1
            
            if vertexCount.get(src) is None:
                vertexCount[src] = 1
            tempVertexList.append(vertexCount)

        shortestPathCounts.append(paths)
        vertexInPathList.append(tempVertexList)

    '''
    computing betweeness centrality from the vertex count and list obtained from 
    modified BFS
    '''
        
    for v_i in range(0,dim_c[0]):
        tempB = 0.0
        for v_j in range(0,dim_c[0]):
            for v_k in range(0,dim_c[0]):
                if v_i != v_j and v_i != v_k:
                    tempB += (vertexInPathList[v_j][v_k].get(v_i,0))/float(shortestPathCounts[v_j][v_k])
        betweenNessCentrality.append(tempB/2.0)
    return betweenNessCentrality

def generateEmptyMat(row_dim,col_dim):
    '''
    returns a zero matrix of specified rows and columns
    '''
    
    mat_empty = [[0.0 for cols_l in range(0,col_dim)] for rows_l in range(0,row_dim)]
    return mat_empty
def getTranspose(mat_A_):
    '''
    returns transposed list of the given matrix
    '''
    dim_A_ = getDimensionOfMatrix(mat_A_)
    mat_B_ = [[row[i] for row in mat_A_] for i in range(0,dim_A_[1])]
                
    return mat_B_
def getColumnsOfMat(mat_A):
    '''
    returns columns of a matrix as list
    '''
    mat_B_ = getTranspose(mat_A)
    dim_B_ = getDimensionOfMatrix(mat_B_)
    colList = list()
    for rows_i in range(0,dim_B_[0]):
        colList.append(mat_B_[rows_i])
    return colList
    

def scalarMul(scalar,vec):
    '''
    perform scalar multiplication with a vector and returns
    the scaled vector
    '''
    vec = [elem * scalar for elem in vec]
    return vec

def vecAdd(vec1,vec2):
    '''
    performs addition between 2 vectors
    '''
    vec3 = [ vec1[add_i]+vec2[add_i] for add_i in range(0,len(vec1))]
    return vec3

def vecSub(vec1,vec2):
    '''
    performs subtraction between 2 vectors
    '''
    vec3 = [ float(vec1[add_i])-float(vec2[add_i]) for add_i in range(0,len(vec1))]
    return vec3

    
     



def orthogonalizeMat(mat_A_r):

    '''
    orthogonalize the columns of the input matrix.
    process used : Gram-Schmidt
    '''
      
    mat_A_ = mat_A_r#copy.deepcopy(mat_A_r)
    cols_A = getColumnsOfMat(mat_A_)
    
    col_norms = [ getNorm2(cols_A[ind_]) for ind_ in range(0,len(cols_A)) ]
    cols_A = [ [elem/col_norms[ind_] for elem in cols_A[ind_]] for ind_ in range(0,len(cols_A)) ]
    
    orthogonal_col = list()
    orthogonal_col.append(cols_A[0])
    for col_i in range(1,len(cols_A)):
        tempCol = cols_A[col_i]
        for col_j in range(0,col_i):
#             numerator = getDotProduct(cols_A[col_i], orthogonal_col[col_j])
#             denom = getDotProduct(orthogonal_col[col_j], orthogonal_col[col_j])
            
            numerator = np.dot(cols_A[col_i], orthogonal_col[col_j])
            denom = np.dot(orthogonal_col[col_j], orthogonal_col[col_j])
            
            scalarFactor = float(numerator)/denom 
#             tempSum = scalarMul(scalarFactor, orthogonal_col[col_j])
#             tempCol = vecSub(tempCol,tempSum)
            tempSum = np.array(orthogonal_col[col_j])*scalarFactor
            tempCol = tempCol-tempSum
            
            
        orthogonal_col.append(tempCol)
        print("orthogonalizing column no. : ",col_i)
    
    print("All columns orthogonalized")
    
#     for col_i in range(0,len(orthogonal_col)):
#         tempNorm = getNorm2(orthogonal_col[col_i])
#         rec_norm = 1.0 / tempNorm
#         orthogonal_col[col_i] = scalarMul(rec_norm, orthogonal_col[col_i])
    col_norms = [ getNorm2(orthogonal_col[col_i]) for col_i in range(0,len(orthogonal_col)) ]
    orthogonal_col = [ [elem/col_norms[ind_] for elem in orthogonal_col[ind_]] for ind_ in range(0,len(orthogonal_col)) ]  
        
        
        
#     for col_i in range(0,len(orthogonal_col)):
#         tempNorm = getNorm2(orthogonal_col[col_i])
#         if tempNorm > 1.2:
#             print(tempNorm,orthogonal_col[col_i],col_i)
    
    orthogonalizedMat = getTranspose(orthogonal_col)
    
    return orthogonalizedMat



def getMatMul(mat_A_,mat_B_):
    '''
    returns matrix multiplication product of two matrices
    '''
    
    result = [[sum(a*b for a,b in zip(X_row,Y_col)) for Y_col in zip(*mat_B_)] for X_row in mat_A_]

    return result


def checkTolerance(mat_old,mat_new,tolerance):
    '''
    returns True if error is outside the tolerance limit specified
    '''
    dim_A_ = getDimensionOfMatrix(mat_old)
    breakFlag = False
    for r_i in range(0,dim_A_[0]):
        for col_i in range(0,dim_A_[1]):
            if r_i == col_i:
                current_ratio = float(abs(mat_old[r_i][col_i] - mat_new[r_i][col_i]))
                current_ratio /= float(mat_new[r_i][col_i])
                if current_ratio <= tolerance:
                    breakFlag = False
                else:
                    breakFlag = True
                    break
        if breakFlag == True:
            break
    return breakFlag

def clipDiagonalEntries(mat_A_):
    '''
    returns the diagonal elements of a matrix
    '''
    dim_A_ = getDimensionOfMatrix(mat_A_)
    
    diagVal = list()
    
    for r_i in range(0,dim_A_[0]):
        for col_i in range(0,dim_A_[1]):
            if r_i == col_i:
                diagVal.append(truncate(mat_A_[r_i][col_i],3))
    return diagVal
def performQRItertation(mat_A):
    '''
    QR implementation for computation of eigen values and
    eigen vectors
    returns :
    list : [eigen_value_list,eigen_vector_list]
    The eigen-vector corresponds to the order of eigen values in eigen_value_list
    
    Orthogonalization is done through GivensRotation
    
    '''
    nIter = 1000
    tolerance = 10e-6
    mat_A_ = copy.deepcopy(mat_A)
    prodQ = None
    for iter_ in range(0,nIter):

        qr_pair = qrFactorGivens(mat_A_)
        matQ = qr_pair[0]
#         matQ_T = getTranspose(matQ)
        matR = qr_pair[1]

        mat_A_old = copy.deepcopy(mat_A_) 
        mat_A_ = getMatMul(matR, matQ)
        
        if iter_ == 0:
            prodQ = matQ
        if iter_ > 0:
            prodQ = getMatMul(prodQ,matQ)
        if not checkTolerance(mat_A_old, mat_A_, tolerance) and (iter_ > 100):
            break
        
    cols_Q = getColumnsOfMat(prodQ) 
    
#     eigVectors = list()
#     for col_i in range(0,len(cols_Q)):
#         tempNorm = getNorm2(cols_Q[col_i])
#         recNorm = 1.0 / tempNorm
#         eigVectors.append(scalarMul(recNorm,cols_Q[col_i]))
        
    col_norms = [ getNorm2(cols_Q[col_i]) for col_i in range(0,len(cols_Q)) ]
    eigVectors = [ [elem/col_norms[ind_] for elem in cols_Q[ind_]] for ind_ in range(0,len(cols_Q)) ]      
    
    
    eigValues = clipDiagonalEntries(mat_A_)
    return [eigValues,eigVectors]

def combineListToDict(key_list,val_list):
    '''
    combines list to key names and returns a dict
    '''
    nameDict = dict()
    for list_i in range(0,len(key_list)):
        nameDict[key_list[list_i]] = val_list[list_i]

    return nameDict

def getLaplacianMatrix(diagonalMat,mat_A_):
    '''
    returns laplacian of the input adjacency matrix
    '''
    dim_A_ = getDimensionOfMatrix(mat_A_)
    L_mat = list()
    
    for r_i in range(0,dim_A_[0]):
        r_list = list()
        for c_i in range(0,dim_A_[1]):
            r_list.append(float(diagonalMat[r_i][c_i]) - float(mat_A_[r_i][c_i]))
        L_mat.append(r_list)
    
    return L_mat


def generatePseudoIdentity(rows,columns):
    '''
    Generates Pseudo Identity Matrix
    returns:
    Identity Matrix
    '''
    mat_i = list()
    for i in range(0,rows):
        row_i = list()
        for j in range(0,columns):
            if(i==j):
                row_i.append(1.0)
            else:
                row_i.append(0.0)
        mat_i.append(row_i)
    return mat_i

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

def signOfScalar(scalar_):
    '''
    returns sign of a scalar value
    '''
    sign_ = scalar_ / abs(scalar_)
    return sign_

def getGivenRotMat(row_i,row_j,a,b,dim_mat):
    '''
    returns given rotation matrix needed to convert the
    upper triangular formation of the matrix
    '''
    
    if b == 0.0:
        c_l = 1.0
        s_l = 0.0
        r_l = 1.0 #a
        rho_ = 0.0
    elif a == 0.0:
        c_l = 0.0
        s_l = 1.0
        r_l = 1.0 #b
        rho_ = 1.0 
    else:
        if abs(b) > abs(a):
            r_l = -a/(b)
            s_l = 1 + (r_l**2)
            s_l = s_l**(0.5)
            s_l = 1/s_l
            c_l = s_l * r_l
            rho_ = 2.0/c_l
        else:
            r_l = -b/(a)
            c_l = 1 + (r_l**2)
            c_l = c_l**(0.5)
            c_l = 1/c_l
            s_l = c_l * r_l
            rho_ = s_l / 2.0
    
    givens_mat_T = [[c_l,-s_l],[s_l,c_l]]
    givens_mat = [[c_l,s_l],[-s_l,c_l]]
    return [givens_mat,givens_mat_T]#identity_mat


def subMatrixMul(mat_A_,mat_B_,start_r,end_r,n_rows,n_cols):
    '''
    A utility function to perform matrix multiplication on specified
    rows and columns
    '''
    tempMat = list()
    for r_i in range(start_r,end_r+1):
        tempMat.append(mat_A_[r_i])
# # 
    slicedMat = getMatMul(mat_B_,tempMat)
#     
    for r_i in range(start_r,end_r+1):
        mat_A_[r_i] = slicedMat[r_i-start_r]
    
    return mat_A_

def qrFactorGivens(mat_A_):
    '''
    returns QR Factorization of the
    input matrix using Givens Rotation
    
    output : list : [Q_Matrix,R_Matrix] 
    '''
    mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    
    max_r = max(dim_A_[0],dim_A_[1])
#     identity_mat = generatePseudoIdentity(dim_A_[0],dim_A_[1])
    identity_mat = generateIdentityMatrix(max_r,max_r)
    for m_c_i in range(0,dim_A_[1]):
        for m_r_i in range(dim_A_[0]-1,m_c_i,-1):
            if m_r_i != m_c_i:
                givens_pair = getGivenRotMat(m_r_i, m_r_i-1, mat_A_[m_r_i-1][m_c_i],
                                           mat_A_[m_r_i][m_c_i] ,
                                            dim_A_)
                identity_mat = subMatrixMul(identity_mat, givens_pair[1],
                                             m_r_i-1, m_r_i,max_r,max_r)
                
                
                mat_A_ = subMatrixMul(mat_A_, givens_pair[1], m_r_i-1, m_r_i, dim_A_[0],dim_A_[1])
                
    
    q_mat = getTranspose(identity_mat)
    r = mat_A_
    return [q_mat,r]

def getNonZeroValCount(vec_):
    '''
    returns count of non-zero values in the input vector
    '''
    non_zero = sum(1  for val in vec_ if abs(val) > 10e-6)
    return non_zero 

def extendBasis(basis_dim_required,basis):
    '''
    performs basis extension, by calling orthogonalize 
    through GramShmidt process
    '''
    basis = copy.deepcopy(basis)
    noOfVectorsReq = abs(len(basis) - basis_dim_required)
    basis_vec_len = len(basis[0])
    extra_vecs = [[1.0 if  r_i==c_i else 0.0 for r_i in range(0,basis_vec_len)] for c_i in range(len(basis),len(basis)+noOfVectorsReq)]
    
    basis.extend(extra_vecs)
    basis = getTranspose(basis)
    basis = orthogonalizeMat(basis)
    return basis
    

def getSVDFactors(mat_A_,n_components = None):
    '''
    returns svd factors of the given input matrix,
    depending number of significant components specified
    in n_components parameter
    
    '''
    print("----------------------------------")
    
#     mat_A_ = copy.deepcopy(mat_A_)
    dim_A_ = getDimensionOfMatrix(mat_A_)
    if n_components is None:
        n_components = dim_A_[1]
        print("performing SVD factorization")
    else:
        print("performing truncated for SVD for ",n_components," components")
    
    
#     print("no. of components ",n_components)
    mat_A_T = np.transpose(mat_A_)
    mat_A_T_A = np.matmul(mat_A_T,mat_A_)
#     print("A_T_A generated",mat_A_T_A.shape)
    e_val,e_vector = np.linalg.eig(mat_A_T_A) 
    e_val = e_val.astype(np.float32)
    e_vector = e_vector.astype(np.float32)
    
    e_vector = e_vector.tolist()
    
    e_vector = getTranspose(e_vector)
    
    
    abs_e_val = [abs(x) for x in e_val]
    
    
    sortIndex = [i[0] for i in sorted(enumerate(abs_e_val), key=lambda x:x[1])]
    
    
    
    sigma_val = [ abs(e_val[sortIndex[dim_A_[1]-1-c_i]])**0.5 for c_i in range(0,dim_A_[1]) ]
    
#     sigma_val = [ 0.0 if ind_ >= n_components else sigma_val[ind_] for ind_ in range(0,dim_A_[1])]
    
    sigma_val = [sigma_val[ind_] for ind_ in range(0,n_components)]
    
    e_vector = [e_vector[sortIndex[dim_A_[1]-1-c_i]] for c_i in range(0,dim_A_[1])] 
    
    e_vector = [e_vector[ind_] for ind_ in range(0,n_components)]
    
    
    mat_sigma = [[ sigma_val[c_i]  if  r_i==c_i else 0.0 for c_i in range(0,n_components)] for r_i in range(0,n_components)]
    
    print(" Sigma Matrix generated ")
    
#     e_vect_norms = [getNorm2(vec) for vec in e_vector]
    
    e_vect_norms = [np.linalg.norm(vec) for vec in e_vector]
    
    
    eig_vector_n = [ [(elem/e_vect_norms[index_]) for elem in e_vector[index_]] for index_ in range(0,len(e_vect_norms))]
    
    u_vector_count = getNonZeroValCount(sigma_val)

    u_cols_size = min(len(eig_vector_n),u_vector_count)
    
#     A_vi = [ getMatrixVectorProduct(mat_A_,eig_vector_n[ind_],3,3) for ind_ in range(0,u_cols_size) ]
    
#     print("computing A_vi ")
    
    A_vi = np.matmul(mat_A_,np.transpose(eig_vector_n))
    
#     A_vi = [ np.dot(mat_A_,eig_vector_n[ind_]) for ind_ in range(0,u_cols_size) ]
    
    A_vi = np.transpose(A_vi)
    
    mat_u = [ [elem/sigma_val[ind_] for elem in A_vi[ind_]] for ind_ in range(0,len(A_vi)) ]
    
#     mat_u = [ [elem/sigma_val[ind_] for elem in A_vi[ind_]] for ind_ in range(0,len(A_vi)) ] 
    
    
    
    mat_V_T = eig_vector_n#getTranspose(eig_vector_n)
    
    
    print(" generating U matrix") 
    
    if len(mat_u) < n_components:
        mat_u = extendBasis(n_components,mat_u)
    else:
        mat_u = getTranspose(mat_u)
        mat_u = orthogonalizeMat(mat_u)
        
        
    print(" matrix U generated ")
    
    rec_a = np.matmul(mat_sigma, mat_V_T)
    rec_a = np.matmul(mat_u, rec_a)
    print("SVD Factorization Completed for ",n_components," components")
    print("----------------------------------")
    
    return [mat_u,mat_sigma,mat_V_T]
    

if __name__ == '__main__':

    if len (sys.argv) != 3 :
            print("Usage: python main.py problem# <absoluteInputDataPath>")
            sys.exit(1)
    problemIndex = sys.argv[1]
        
    if problemIndex == "problem1":
#             mat_A = [[-2.0,-4.0,2],[-2,1,2],[4,2,5]]
#             mat_A = [[2,-12],[1,-5]]
#             mat_A = [[12,-51,4],[6,167,-68],[-4,24,-41]]
            inputFile = sys.argv[2]
            mat_A = readFile(inputFile)
            dim_A = getDimensionOfMatrix(mat_A)
            print("input matrix : ")
            print(mat_A)
            eigenPair = performQRItertation(mat_A)
            print("Eigen Value and Eigen Vector Pair for the input matrix :")
            for eg_index in range(0,len(eigenPair[0])):
                print("eigen value : ",truncate(eigenPair[0][eg_index],3),",eigen vector : ",eigenPair[1][eg_index])
    if problemIndex == "problem2":
        gml_file = sys.argv[2]
        lclassico_data = nx.read_gml(gml_file)
        lclassico_mat = nx.adjacency_matrix(lclassico_data)
        lclassico_mat = lclassico_mat.todense().tolist()
        options = {
            'node_size': 30,
            'line_color': 'black',
            'linewidths': 0,
            'width': 0.1,
            'font_size':7 
            }
        indicatorVector = list()
        dim_mat = getDimensionOfMatrix(lclassico_mat)
        n_rows_ii = dim_mat[0]
        n_cols_ii = dim_mat[1]
        for cols_iter in range(0,n_cols_ii):
            indicatorVector.append(1)
        max_degree = n_cols_ii -1
        degree_vector = getMatrixVectorProduct(lclassico_mat, indicatorVector, n_rows_ii, n_cols_ii)
#         for degree_iter in range(0,len(degree_vector)):
#             degree_vector[degree_iter] /= float(max_degree)
        
        ''' reading player names'''
        playerNames = list(lclassico_data.nodes())
        degreeVectorDict = combineListToDict(playerNames, degree_vector)
        print("Degree Centrality of each node :")
        print(degreeVectorDict)
        
        plt.plot(degree_vector,marker='o')
        plt.title("Degree Centrality")
        plt.xlabel("Player Name")
        plt.ylabel("Degree Centrality")
        plt.grid()
#         x, y = zip(*degree_centrality(lclassico_data).items())
#         plt.plot(y,marker='o')
        plt.xticks(np.arange(0,len(playerNames)),playerNames,rotation='vertical')
        plt.tight_layout()
        plt.show()
        
        
        distances_nodes = all_pair_shortest_floyd(lclassico_mat)
        distance_vector = getMatrixVectorProduct(distances_nodes, indicatorVector, n_rows_ii, n_cols_ii)
        for distance_iter in range(0,len(distance_vector)):
            distance_vector[distance_iter] = float(max_degree)/(distance_vector[distance_iter])
        
        distance_vector_dict = combineListToDict(playerNames, distance_vector)
        print("Closeness Centrality of each node :")
        print(distance_vector_dict)
        
        plt.plot(distance_vector,marker='o')
        plt.title("Closeness Centrality")
        plt.xlabel("Player Name")
        plt.ylabel("Closeness Centrality")
        plt.grid()
#         x, y = zip(*closeness_centrality(lclassico_data).items())
#         plt.plot(y,marker='o')
        plt.xticks(np.arange(0,len(playerNames)), playerNames,rotation='vertical')
        plt.tight_layout()
        plt.show()
        
        
        betweeness_centrality = modBFS(lclassico_mat)
        betweenness_vector_dict = combineListToDict(playerNames,betweeness_centrality)
        print("Betweenness Centrality of each node :")
        print(betweenness_vector_dict)
        
        plt.plot(betweeness_centrality,marker='o')
        plt.title("Betweenness Centrality Un-normalized")
        plt.xlabel("Player Name")
        plt.ylabel("Betweenness Centrality")
        plt.grid()
#         x, y = zip(*betweenness_centrality(lclassico_data,normalized=False).items())
#         plt.plot(y,marker='o')
        plt.xticks(np.arange(0,len(playerNames)), playerNames,rotation='vertical')
        plt.tight_layout()
        plt.show()
       
       
       
        dominantEigenPair = getLargestEigenPair(lclassico_mat)
    
        dominantEigenVec = dominantEigenPair[1]
        dominantEigenValue =  dominantEigenPair[0]
        
        eigenVector_dict = combineListToDict(playerNames, dominantEigenVec)
        print("Eigenvector Centality of each node :")
        print(eigenVector_dict)
        
        plt.plot(dominantEigenVec,marker='o')
        plt.title("EigenVector Centrality")
        plt.xlabel("Player Name")
        plt.ylabel("EigenVector Centrality")
        plt.grid()
#         x, y = zip(*eigenvector_centrality_numpy(lclassico_data).items())
#         plt.plot(y,marker='o')
        plt.xticks(np.arange(0,len(playerNames)),playerNames,rotation='vertical')
        plt.tight_layout()
        plt.show()
        
        '''
        Sanity Code for Normalized Comparison of Centrality Measures   
        
        
        degree_vector = (degree_vector - np.min(degree_vector)) / float(np.max(degree_vector)-np.min(degree_vector))
        distance_vector = (distance_vector - np.min(distance_vector)) / float(np.max(distance_vector)-np.min(distance_vector))
        
        betweeness_centrality = (betweeness_centrality - np.min(betweeness_centrality)) / float(np.max(betweeness_centrality)-np.min(betweeness_centrality))
        
        dominantEigenVec = (dominantEigenVec - np.min(dominantEigenVec)) / float(np.max(dominantEigenVec)-np.min(dominantEigenVec))
        
        plt.plot(degree_vector,marker='o',label='Degree Centrality')
        plt.plot(distance_vector,marker='o',label='Closeness Centrality')
        plt.plot(betweeness_centrality,marker='o',label='Betweenness Centrality')
        plt.plot(dominantEigenVec,marker='o',label='EigenVector Centrality')
        plt.title("Normalized Centrality Measures For Comparison")
        plt.xlabel("Player Name")
        plt.ylabel("centrality measures")
        plt.grid()
        plt.xticks(np.arange(0,len(playerNames)), labels=playerNames,rotation='vertical')
        plt.legend()
        plt.tight_layout()
        plt.show()
        
    '''
        
    if problemIndex == "problem3":
        gml_file = sys.argv[2]
        lclassico_data = nx.read_gml(gml_file)#"./../L-Classico Network/lclassico.gml")

        lclassico_mat = nx.adjacency_matrix(lclassico_data)
        lclassico_mat = lclassico_mat.todense().tolist()
        
        indicatorVector = []
        dim_mat = getDimensionOfMatrix(lclassico_mat)
        n_rows_ii = dim_mat[0]
        n_cols_ii = dim_mat[1]
        for cols_iter in range(0,n_cols_ii):
            indicatorVector.append(1.0)
        degree_vector = getMatrixVectorProduct(lclassico_mat, indicatorVector, n_rows_ii, n_cols_ii)
        
        degree_mat = list()
        for r_i in range(0,n_rows_ii):
            r_list = []
            for c_i in range(0,n_cols_ii):
                if r_i == c_i:
                    r_list.append(float(degree_vector[r_i]))
                else:
                    r_list.append(0.0)
            degree_mat.append(r_list)
        
        print("spectral clustering initiated")
        lclassico_mat_norm = getLaplacianMatrix(degree_mat, lclassico_mat)
        
        eigenPairs = performQRItertation(lclassico_mat_norm)
      
        minEigenValue = 999999.9
        eigenValue_list = eigenPairs[0]
        eigenVector_list = eigenPairs[1]
        minIndex = 0
        second_smallest_eigIndex = 0
        for eig_index in range(0,len(eigenValue_list)):
            if eigenValue_list[eig_index] < minEigenValue:
                minEigenValue = eigenValue_list[eig_index] 
                second_smallest_eigIndex = minIndex
                minIndex = eig_index
                 
        
        print("spectral clustering finished")     
        colorMap = []
        fiedlerColorMap = []
        for i in eigenVector_list[second_smallest_eigIndex]:
            if(i >= 0):
                colorMap.append('cyan')
                fiedlerColorMap.append('red')
            else:
                colorMap.append('gold')
                fiedlerColorMap.append('blue')
        
        fiedlerVector = eigenVector_list[second_smallest_eigIndex]
         
        plt.figure(figsize=(7,5))
        plt.scatter(np.arange(0,len(eigenVector_list[second_smallest_eigIndex])),fiedlerVector,color=fiedlerColorMap)
        plt.title("Partitions of nodes from sign of fiedler vector components")
        red_patch = mpatches.Patch(color='red', label='Cluster-1')
        blue_patch = mpatches.Patch(color='blue', label='Cluster-2')
        plt.legend(handles=[red_patch,blue_patch])
        plt.grid()
        plt.tight_layout()
        plt.show()
        
        options = {
        'node_size': 30,
        'line_color': 'black',
        'linewidths': 0,
        'width': 0.1,
        'font_size':7 
        }
         
        plt.figure(figsize=(7,5))
        plt.title("Community Separartion From Fiedler Vector Component Clustering")
        red_patch = mpatches.Patch(color='cyan', label='Cluster-1')
        blue_patch = mpatches.Patch(color='gold', label='Cluster-2')
        plt.legend(handles=[red_patch,blue_patch])
        nx.draw(lclassico_data,node_color=colorMap,with_labels=True,hold=True,**options)
        plt.tight_layout()
        plt.show()
    if problemIndex == "problem4a":
        mnist_file = sys.argv[2]
        
#         mnist_file = "../MNIST-Dataset/mnist_train.csv"
        with open(mnist_file) as f:
            lines = f.readlines()
#         lines = [x.strip() for x in lines]
        
        lines = [line.split(',') for line in lines]
        lines = [ [float(x) for x in line] for line in lines ]
        
        mnist_mat = getTranspose(lines)
        
        labelColumn = mnist_mat.pop(0)
        
        mnist_mat = getTranspose(mnist_mat)
       
        
        print(" initiate svd")
        
        n_components_ = [2,5,10,20,50,100,200,500]
        
        rmse_list = []
        for d in n_components_:
            mat_u,mat_sigma,mat_V_T = getSVDFactors(mnist_mat,d)
            mnist_mat_d = np.matmul(mat_sigma, mat_V_T)
            mnist_mat_d = np.matmul(mat_u, mnist_mat_d)
            rmse_d = (np.array(mnist_mat).astype(float) - np.array(mnist_mat_d).astype(float))**2
            rmse_d = np.sum(rmse_d)
            rmse_d = np.sqrt(rmse_d)
            rmse_list.append(rmse_d)
        
        
        print("truncated svd d value : ",n_components_)
        print("corressponding rmse values : ")
        print(rmse_list)
        plt.scatter(n_components_,rmse_list)
        plt.plot(n_components_,rmse_list)
        plt.xlabel("No. of Components d")
        plt.ylabel("RMSE")
        plt.title("Reconstruction error for different d values in Truncated SVD")
        plt.grid()
        plt.show()
        
           
    if problemIndex == "problem4b":
        
        mnist_file = sys.argv[2]
#         mnist_file = "../MNIST-Dataset/mnist_train.csv"
        with open(mnist_file) as f:
            lines = f.readlines()
            
        lines = [line.split(',') for line in lines]
        lines = [ [float(x) for x in line] for line in lines ]
     
        mnist_mat = np.array(lines)
        
        labelColumn = np.copy(mnist_mat[:,0])
        
        mnist_mat = np.delete(mnist_mat,0,1)
        print("data loaded")
        print("dataset shape to work upon : ")
        print(mnist_mat.shape)
        
        print("t-SNE function invoked")
        projected_vec_2d = TSNE(n_components=2).fit_transform(mnist_mat)
   
        print(projected_vec_2d.shape)
         
        plt.scatter(projected_vec_2d[:,0],projected_vec_2d[:,1])
        plt.title("t-SNE projection of digit vectors on 2-D space")
        plt.show()
        
        
        
        
        
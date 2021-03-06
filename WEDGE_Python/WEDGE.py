import libNNLS as NNLS
import scipy.sparse as ss
from scipy.io import mmread
import  numpy as np
import pandas as pd
import os
import multiprocessing
from scipy.sparse.linalg import svds
import time

class WEDGE:
    def __init__(self,
                 input_path=None,
                 output_path='./WEDGE_res/',
                 output_format='',
                 normalization=True,
                 n_CPU=1,
                 options={},
                 recovery_fileName='WEDGE_recovery.csv'):
        """this function used to recovery scRNA-seq data    

        Args:
            input_path (str, optional): the dir of gene expression matrix. n_gene x n_cell. Defaults to None.
            output_path (str, optional): [description]. Defaults to './WEDGE_res/'.
            output_format (str, optional): [description]. 
            you can set output_format='' or output_format='W_H', Defaults output_format='' .
            normalization (bool, optional): we normalize the total expression of each cell to 10,000,
            and perform log-transform after adding a pseudocount of 1. . Defaults to True.
            n_CPU (int, optional): number of cpu. Defaults to 1.
            recovery_fileName (str, optional): [description]. Defaults to 'WEDGE_recovery.csv'.
            options (dictionary, optional): [description]. Defaults to {}.
            % options['n_rank'] : the rank of A_norm;
            % options['lambda']: the weights of  Zero elements;
            % options['max_iters']: the maximum number of iterations;
            % options['tolx'] the change value of object function when iter is convergence;
            % options['bound'] : the bound used to estimate rank, rank is smaller when bound is bigger. default  0.085;
            %there are two select for  lambda, default lambda = 'auto' 
            %1. setting postive number to lambda
            %2. lambda = 'auto';  %when r <= 0.25, lambda = 0.15;when r > 0.25,
            %lambda = r (r is the rate of no zero elements).
        """
        self.n_CPU = n_CPU
        self.input_path=input_path
        self.output_path=output_path
        self.output_format=output_format
        self.normalization=normalization
        self.options = options
        self.recovery_fileName=recovery_fileName
        self.n_rank = -1
        self.Lambda= 0.15
        self.tolX = 1e-5
        self.max_iters = 100
        self.bound = 0.085

        if 'n_rank' in options and type(options['n_rank']) == type(1):
            self.n_rank = options['n_rank']

        if 'lambda' in options:
            self.Lambda= options['lambda']

        if 'tolx' in options:
            self.tolX = options['tolx']

        if 'max_iters' in options:
            self.max_iters = options['max_iters']

        if 'bound' in options:
            self.bound = options['bound']
        
        if self.n_rank>0:
            self.n_pca =max(50, self.n_rank)
        else:
            self.n_pca =50 

    def WEDGE_recovery(self):
        print('It is reading')
        self.read_data()

        pass
        if self.normalization:
            print('It is normalizing')
            if type(self.Data)is np.ndarray:
                self.Data = ss.coo_matrix(self.Data)
            
            Data_ing = ss.csc_matrix(self.Data)
            self.Data = ss.coo_matrix(np.dot(Data_ing, ss.dia_matrix((10000/(Data_ing.sum(axis=0) + 1e-9), 0),
                                                       shape=(Data_ing.shape[1], Data_ing.shape[1]))))
            index_i = self.Data.row
            index_j = self.Data.col
            Data_v  = self.Data.data
            self.Data = ss.coo_matrix((np.log(Data_v +1), (index_i, index_j)))
    
            del index_i, index_j, Data_v, Data_ing

            if self.Lambda =='auto':
                nozeros_rate = np.sum(self.Data>0)/(np.shape(self.Data)[0]*np.shape(self.Data)[1])
                if nozeros_rate<= 0.25:
                    self.Lambda = 0.15
                else:
                    self.Lambda = nozeros_rate

            self.Data = self.Data.T
            n_rank_initial = self.n_rank
            self.WEDGE_initial(n_pca=self.n_pca);   
            if n_rank_initial < 0:
                self.n_rank =5* self.n_rank
                self.n_rank = min(self.n_rank,min(50,int(np.shape(self.Data)[1]*0.05)))    

        print('It is recovering')
        start = time.time()
        W, H = self.recoverWNMF()
        end = time.time()
        self.W = H.T
        self.H = W.T
        print('It has finished WEDGE Imputation;'+ " Total spend %f s"%(end-start))

        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path)

        if self.output_format=='W_H':
            np.savetxt(self.output_path+'W.csv', W, delimiter = ',') 
            np.savetxt(self.output_path+'H.csv', H, delimiter = ',')
            pd.DataFrame(data=self.geneName.tolist()).to_csv(self.output_path+'geneName.csv',
                                                             encoding='utf-8', index=False, header=False)
            pd.DataFrame(data=self.cellName.tolist()).to_csv(self.output_path+'cellName.csv',
                                                             encoding='utf-8', index=False, header=False)

        else:
            if len(self.geneName.shape)>1:
                self.geneName = self.geneName[:,0]
            
            A_recovery = np.dot(self.W, self.H)

            A_recovery = np.dot(self.H,
                                ss.dia_matrix((10000./(np.sum(A_recovery, axis=0) + 1e-9), 0),
                                              shape=(np.shape(self.H)[1],np.shape(self.H)[1])))
            A_recovery = np.dot(self.W, self.H)

            
            if self.input_path[-1]=='/':
                T = pd.DataFrame(A_recovery)
                T.columns = self.cellName
                T.index =self.geneName
            else:
                T = pd.DataFrame(A_recovery,index=self.geneName,columns = pd.Index.unique(self.cellName))
            T.to_csv(self.output_path+self.recovery_fileName, index=True, header=True)
    
        print('WEDGE finished')

    def read_data(self):
        """
        This function is used to read the gene expression matrix
        Note: this Function only used to read .csv, .tsv and .mtx
        """
        # read csv file
        if self.input_path[-3:]=='csv':
            sep_str = ','
        if self.input_path[-3:]=='tsv':
            sep_str = '\t'
        if self.input_path[-3:]=='csv' or self.input_path[-3:]=='tsv':
            self.Initial_Data = pd.read_csv(filepath_or_buffer=self.input_path,
                                        sep=sep_str,
                                        header=0,
                                        index_col=0)
            #Summing the expression of genes of the same name
            self.geneName = self.Initial_Data.index
            if np.size(self.geneName) > np.size(np.unique(self.geneName)):
                self.geneName = np.unique(self.geneName)
                self.Initial_Data = self.Initial_Data.groupby(level=0).sum()

            self.Data = ss.coo_matrix(self.Initial_Data.values)
            self.cellName = self.Initial_Data.columns
            return
        # read 10x data
        if self.input_path[-1]=='/':
            if os.path.exists(self.input_path+'/genes.tsv'):
                geneName = pd.read_csv(
                    filepath_or_buffer=self.input_path+'genes.tsv',
                    sep='\t',
                    header=None,
                    index_col=None)
                self.geneName = geneName.values
            else:
                print('Do not find the genes.tsv in: '+ self.input_path+'genes.tsv')
                return 

            if os.path.exists(self.input_path+'barcodes.tsv'):
                cellName = pd.read_csv(
                    filepath_or_buffer=self.input_path+'barcodes.tsv',
                    sep='\t',
                    header=None,
                    index_col=None)
                self.cellName = cellName.values[:,0]
            else:
                print('Do not find the barcodes.tsv in: '+ self.input_path+'barcodes.tsv')
                return

            if os.path.exists(self.input_path+'matrix.mtx'):
                self.Initial_Data = mmread(self.input_path+'matrix.mtx')

            else:
                print('Do not find the matrix.mtx in: '+ self.input_path +'matrix.mtx')
                return 
            
            _, c, b = np.unique(self.geneName[:,0],return_index= True, return_inverse=True)
            if len(c) < len(b):
                index_i = self.Initial_Data.row
                index_j = self.Initial_Data.col
                Data_v  = self.Initial_Data.data
                index_i_new = index_i
                for i in range(len(b)):
                    index_i_new[index_i==i] = b[i] 

                self.Data = ss.coo_matrix((Data_v, (index_i_new,index_j)))
                self.geneName = self.geneName[c,:]
            else:
                self.Data = self.Initial_Data
        else:
            print('the format of input file is .csv, .tsv and 10x path')

    def WEDGE_initial(self, n_pca=100, bound=0.085):
        """this function is used to calculate the initial value of SNMF  and estimate the rank of A

        Args:
            A (np.array()): n_gene * n_cell. Defaults to None.
            n_pca (int): the number of  singleValue will be calculated. Defaults to 100.
            bound (float): estimate the rank , rank is bigger when bound is smaller. Defaults to 0.085.

        Return:
            initial_W : the initial value of W in SNMF.
            n_rank : the rank of self.Data.

        """
        n_pca = min(n_pca, np.shape(self.Data)[1] - 1)
        W0, single_value, _ = svds(self.Data,
                                   k=n_pca,
                                   which='LM')

        single_value = single_value[range(np.shape(W0)[1] - 1, -1, -1)]
        W0 = W0[:, range(np.shape(W0)[1] - 1, -1, -1)]

        W0_P = np.maximum(W0, 0)
        W0_N = np.maximum(-W0, 0)
        n_P_big = np.sum(W0_P * W0_P, axis=0)
        n_N_big = np.sum(W0_N * W0_N, axis=0)
        if np.sum(n_P_big >= n_N_big) > 0:
            W0[:, n_P_big >= n_N_big] = W0_P[:, n_P_big >= n_N_big]

        if np.sum(n_P_big < n_N_big) > 0:
            W0[:, n_P_big < n_N_big] = W0_N[:, n_P_big < n_N_big]

        single_value = single_value[single_value > 0]
        n_svd = min(n_pca, len(single_value) - 1)
        single_value = single_value / single_value[0]
        latent_new_diff = single_value[0:n_svd] / single_value[1:n_svd + 1] - 1
        n_rank = 2
        for i in range(len(latent_new_diff) - 10):
            n_rank = i + 1
            if (latent_new_diff[i] >= bound) and all(
                    latent_new_diff[i + 1:11] < bound):
                break

        n_rank = max(n_rank, 3)
        self.n_rank = n_rank
        self.initial_W = W0

    def recoverWNMF(self):
        n = np.shape(self.Data)[1]
        W_in = self.initial_W[:,:self.n_rank].copy()
        H_in = np.zeros((self.n_rank ,n))
        A_F = np.sum(self.Data.data*self.Data.data)
        error_old = 10000000*A_F
        error_min = error_old      
        iter = 1
        count = 0
        W_old = W_in.copy()
        H_old = H_in.copy()
        count_rise = 0 # calculate the number of raising steps;
        # opt = optimset('TolX',1e-6);
        while count < 3 and iter<= self.max_iters:
            W_in, H_in, err = self.WNMF_inner(W_in, H_in, 1e-6)
            error_new = err/A_F
            
            if error_new > (error_old + 1e-5):
                count_rise = count_rise + 1
            else:
                count_rise = 0

            if error_new < error_min:
                error_min = error_new
                W_old = W_in
                H_old = H_in

            if count_rise>3:  
                W = W_old
                H = H_old
                return W, H

            if (abs(error_new - error_old)<= 1e-5):
                count = count +1                       
            else:
                count =0

            error_old = error_new
            # print("Iter:%d, error: %f"%(iter, error_old))
            iter = iter+1

            if iter == self.max_iters:
                    W = W_old
                    H = H_old
                    return W, H

        W = W_old
        H = H_old
        return W, H

    def WNMF_inner(self, W_in=None, H_in=None, TolX=1e-6):
        m, n = np.shape(self.Data)
        re_H = np.zeros(n) 
        # cores = multiprocessing.cpu_count()
        # pool = multiprocessing.Pool(processes=cores)
        error_ = 0.0
        Data_ing = ss.csc_matrix(self.Data)
        # results = []
        for i in range(np.shape(self.Data)[1]):                           
            C_cell = W_in.copy()
            Data_row_i = Data_ing.getcol(i).toarray()
            C_cell[Data_row_i[:,0]==0,] = self.Lambda*C_cell[Data_row_i[:,0]==0,]
            nnls_res = NNLS.lsqnonneg(C_cell, Data_row_i, 1e-6)    
            # results.append(pool.apply_async(NNLS.lsqnonneg, (C_cell, Data_row_i, 1e-6, )))         
            H_in[:,i] = nnls_res[1:]
            re_H[i] = nnls_res[0]
            error_ = error_ + nnls_res[0]
        '''
        pool.close()
        pool.join()
        print(np.array([r.get() for r in results]))
        nnls_res = [r.get() for r in results]
        H_in = nnls_res[1:,:]
        error_ = error_ + sum(nnls_res[0,:])
        # re_W = np.zeros(m) 
        '''
        Data_ing = ss.csr_matrix(self.Data)
        for i in range(np.shape(self.Data)[0]):
            C_gene = H_in.copy()
            Data_col_i = Data_ing.getrow(i).toarray()
            C_gene[:,Data_col_i[0,:]==0] = self.Lambda *H_in[:,Data_col_i[0,:]==0]
            nnls_res = NNLS.lsqnonneg(C_gene.T, Data_col_i.T, 1e-6)
            W_in[i,:] = nnls_res[1:]
            # re_W[i] = nnls_res[0]
            error_ = error_ + nnls_res[0]
        # return W_in, H_in, np.sum(re_H)+ np.sum(re_W)
        return W_in, H_in, error_

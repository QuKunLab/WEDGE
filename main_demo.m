clc;
clear;
close all;
input_path='./Data/Baron/Observed.csv'; %the dir of input file;
%'.csv', '.tsv' and '.mtx' format is ok. Genes in rows and cells in columns.
%If the file format is '.mtx',the input_path must end with '/', for example input_path= './Data/amount/'
%if the filr format is '.csv' or '.tsv',  Cell and gene names are mandatory.

output_format =[] ;%the saver format of recovery data. 
%if you set output_format = 'W_H', the output is W and H, the recovery data is equal to W*H;
%the default is output_format =[]; save recovery data as WEDGE_recovery.csv, rows is genes and columns is cells.

output_path = '';%the dir of output file, for example './res/'. 
%default output_path =[]; In this time if output_format = 'W_H', then output_path ='./res/'; otherwise, output_path ='./';

normalization = 1; 
%if normalization = 0, WEDGE do not normalize data.
%default, normalization = 1, WEDGE normalizes the total expression of each cell to 10,000 and performs log-transform after adding a pseudocount of 1. 

n_CPU = []; 
%if n_CPU = [], WEDGE will use all the CPU in this computer, otherwise,
%WEDGE use min(n_CPU, the number of CPU in this computer) CPU to impute

options.n_rank = []; 
%A positive integer, the rank of gene expression matrix; if options.n_rank = [], WEDGE will automatically estimate the rank of matrix
%otherwise, WEDGE will use the value set by the user.

options.lambda='';
%the weights of  Zero elements;
%setting a postive number to 0<=lambda<=1, default lambda = 0.15 

options.max_iters=[];% the maximum number of iterations for solving weight low-rank factorization;

options.tol=[]; %the error of object function when iter is convergence;

options.bound =[];%the minimum change value of single values when estimating the rank of gene expression matrix;

options.figure_state=0; %1 plot the change value of single values ; 0 don't plot the change value of single values . default 0;

[W,H]= WEDGE_recovery(input_path,output_path,output_format,normalization, n_CPU, options);
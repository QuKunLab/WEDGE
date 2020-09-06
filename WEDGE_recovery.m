function [W,H]= WEDGE_recovery(input_path,output_path,output_format,normalization, n_CPU, options,recovery_fileName)
    % this function used to recovery scRNA-seq data

    %input
    %path: the dir of gene expression matrix. n_gene x n_cell
    %normalization: we normalize the total expression of each cell to 10,000,
    %and perform log-transform after adding a pseudocount of 1. 
    %default, normalization = 1;
    %n_CPU,number of cpu, default use all the cpu on the computer

    % options.n_rank : the rank of A_norm;
    % options.lambda: the weights of  Zero elements;
    % options.max_iters: the maximum number of iterations;
    % options.tol£ºthe change value of object function when iter is convergence;
    % options.bound : the bound used to estimate rank, rank is smaller when bound is bigger. default  0.085;
    % options.figure_state: 1 plot single value ; 0 don't plot single value . default 0;
   
    %there are two select for  lambda, default lambda = 'auto' 
    %1. setting postive number to lambda
    %2. lambda = 'auto';  %when r <= 0.25, lambda = 0.15;when r > 0.25,
    %lambda = r (r is the rate of no zero elements).

if nargin <7
    recovery_fileName = 'WEDGE_recovery.csv';
end
%%
max_ncores = feature('numCores');
pool_state = gcp('nocreate');           

if ~exist('n_CPU', 'var') || isempty(n_CPU) || n_CPU <=0
    n_CPU = max_ncores;
end  
if max_ncores < n_CPU
  fprintf(['You requested %d CPUs,' ...
       'but the cluster "local" has the NumWorkers property set to allow a maximum of %d CPUs.\n' ...
       'WEDGE automatically sets Number of CPUs to %d .'], n_CPU,max_ncores,max_ncores);
   n_CPU = max_ncores;
end

if isempty(pool_state)  
       parpool(n_CPU);
else
   if n_CPU == pool_state.NumWorkers
       fprintf('Parallel pool on local has been run and number of CPUs is %d',n_CPU);
   else
       delete(gcp('nocreate')) 
       parpool(n_CPU);
   end               
end
 %%   
    if ~exist('normalization', 'var') || isempty(normalization) 
        normalization =1;
    end  
    
    if ~exist('output_path', 'var') || isempty(output_path) 
        if exist('output_format', 'var')&& strcmpi(output_format,'W_H')
            output_path ='./res/';
        else
            output_path ='./';
        end
    end   
    
    [n_rank, lambda, max_iters, tol, bound, figure_state] = set_options(options);
    %%
        sprintf('It is reading')
        [A , geneName, cellName] = read_data(input_path);
    %%
    if normalization ==1
        sprintf('It is normalizing')
        A_norm = A *sparse(1:size(A,2),1:size(A,2),10000./(sum(A) + 0.000000001));
        [index_i , index_j , A_value] = find(A_norm);
        A_norm = sparse(index_i ,index_j , log(A_value +1));
    else
        A_norm = A;
    end

    %%
        sprintf('It is recovering')
        begintime = tic;
        [W, H] =recovery(A_norm,n_rank, lambda, max_iters, tol, bound, figure_state); 
        sprintf('the time of recovery is %f s',toc(begintime))

    if exist(output_path,'dir')==0
       mkdir(output_path);
    end

    if strcmpi(output_format,'W_H')
        csvwrite([output_path,'W.csv'],W);
        csvwrite([output_path,'H.csv'],H);
        T = table(geneName);
        writetable(T, [output_path,'geneName.csv'],'WriteRowNames',0,'WriteVariableNames',0);
        T = table(cellName');
        writetable(T, [output_path,'cellName.csv'],'WriteRowNames',0,'WriteVariableNames',0);        
    else
        U = matlab.lang.makeUniqueStrings(cellName);
        A_recovery = W*H;
        n_col = size(A_recovery,2);
        A_recovery = A_recovery *sparse(1:n_col,1:n_col,10000./sum(A_recovery) + 0.000000001);
        try 
            T = array2table(A_recovery,'RowNames',geneName,'VariableNames',U);
        catch ME
            sprintf('Warning the cell name is change')
            strcat('WEDGE_',U)
            T = array2table(A_recovery,'RowNames',geneName,'VariableNames',U);
        end
        writetable(T, [output_path, recovery_fileName],'WriteRowNames',1,'WriteVariableNames',1);
    end
    
            sprintf('WEDGE finished')
end
function [W, H] =recovery(A_norm , n_rank , lambda ,step, error , bound,figure_state) 
% this function used to recovery single-cell data
%input

%A_norm: n_gene x n_cell
%n_rank: the rank of A_norm
%lambda: the weights of  Zero elements
%step: the maximum number of iterations
%error£ºConvergence error

%parameter of calculate initial value
%n_pca:  number of pca in SNMF_initial.  default 50 
%bound: the bound used to estimate rank£¬
%rank is smaller when bound is bigger. default  0.085
%figure_state: 1 plot single value ; 0 don't plot single value . default 1


%
%there are three selections for  lambda, default lambda = 0.15
%1. setting postive number to lambda
%2. lambda = 'Default' or 'default' , in this time lambda = 0.15
%3. lambda = 'auto'    ;we define r indicates the rate of no zero elements;
%when r <= 0.25, lambda = 0.15;when 0.25< r,  lambda = r;

if exist('n_rank', 'var') && ~isempty(n_rank)
    n_pca =max(50, n_rank);
else
    n_pca =50;
end 
if ~exist('bound', 'var') || isempty(bound) 
    bound = 0.085;
end

if ~exist('figure_state', 'var') || isempty(figure_state) 
    figure_state = 0;
end
    
A_norm = A_norm';
[initial_W,ReducedDimension] = SNMF_initial(A_norm, n_pca, bound, figure_state);   
if ~exist('n_rank', 'var') || isempty(n_rank)
     n_rank =5* ReducedDimension;
     n_rank = min(n_rank,min(50,floor(size(A_norm,2)*0.05)));     
end
    
if ~exist('lambda', 'var') || isempty(lambda) ||strcmpi(lambda,'default')
    lambda = 0.15;
end
  
if strcmpi(lambda,'auto')
        nozeros_rate = sum(A_norm(:)>0)/length(A_norm(:));
        if nozeros_rate<= 0.25
            lambda = 0.15;
        else
            lambda = nozeros_rate;
        end
end    

if ~exist('error', 'var') || isempty(error) 
    error = 1e-5;
end
    
if ~exist('step', 'var') || isempty(step) 
    step =100;
end
     
warning off; 

[W,H] =recoverWNMF(A_norm,initial_W ,n_rank,step, lambda,error);    

W_ing = W;
W = H';
H  = W_ing'; 

end

function [A , geneName,cellName] = read_data(path)
%This function is used to read the gene expression matrix
%Note: this Function only used to read .csv, .tsv and .mtx

%% read csv file
if  strcmpi(path(end-2:end),'csv')
    Data_no_change = importdata(path,',',1);
    A_source_head = Data_no_change.textdata;    
    % Summing the expression of genes of the same name
    cellName = A_source_head(1,2:end);
    [geneName,~,b] = unique(A_source_head(2:end,1),'stable');
    clear A_source_head
    if length(unique(b)) < length(b)
        [index_i, index_j, value_x] = find(Data_no_change.data);
        clear Data_no_change

        index_i_new=index_i;
        for i = 1: length(b)
            index_i_new(index_i==i) = b(i); 
        end
        A = sparse(index_i_new,index_j,value_x);

        clear index_i index_i_new index_j value_x
    else
        A = Data_no_change.data;
        A_row_sum = sum(A,2);
        A = A(A_row_sum>0,:);
        geneName = geneName(A_row_sum>0);
    end
    return;
end

%% read tsv file
if  strcmpi(path(end-2:end),'tsv')
    Data_no_change = importdata(path,'\t',1);
    A_source_head = Data_no_change.textdata;
    
    cellName = A_source_head(1,2:end);
    [geneName,~,b] = unique(A_source_head(2:end,1),'stable');
    clear A_source_head
    if length(unique(b)) < length(b)
        [index_i, index_j, value_x] = find(Data_no_change.data);
        clear Data_no_change

        index_i_new=index_i;
        for i = 1: length(b)
            index_i_new(index_i==i) = b(i); 
        end
        A = sparse(index_i_new,index_j,value_x);

        clear index_i index_i_new index_j value_x
    else
        A = Data_no_change.data;
        A_row_sum = sum(A,2);
        A = A(A_row_sum>0,:);
        geneName = geneName(A_row_sum>0);
    end
    return
end

%% read 10x data
try
    fid = fopen([path,'/','barcodes.tsv']);
    cellName = textscan(fid, '%s','delimiter', '\t');
    fclose(fid);
    cellName = cellName{1};

    fid = fopen([path,'/','genes.tsv']);
    geneName = textscan(fid, '%s %s','delimiter', '\t');
    geneName = geneName{1,1};
    fclose(fid);

    fileID = fopen([path,'/','matrix.mtx']);
    C = textscan(fileID,'%n %n %n','CommentStyle','%');
    fclose(fileID);
    A = sparse(C{1}(2:end),C{2}(2:end),C{3}(2:end),C{1}(1),C{2}(1));
    clear C;
catch
    sprintf('the format of input file is .csv, .tsv or 10x path')
    return;
end

 
[geneName, ~, b] = unique(geneName(1:end,1),'stable');

if length(unique(b)) < length(b)
    [index_i, index_j, value_x] = find(A);
    clear A

    index_i_new=index_i;
    for i = 1: length(b)
        index_i_new(index_i==i) = b(i); 
    end
    A = sparse(index_i_new,index_j,value_x);
    A_row_sum = sum(A,2);
    A = A(A_row_sum>0,:);
    geneName = geneName(A_row_sum>0);
end


end

function [initial_W,ReducedDimension] = SNMF_initial(A,n_pca,bound ,figure_state)
    %this function is used to calculate the initial value of SNMF  and
    %estimate the rank of A
    
    %input
    %A         : n_gene * n_cell 
    %n_pca  : the number of  singleValue will be calculated 
    %bound :  estimate the rank , rank is bigger when bound is smaller
    
    %output
    %initial_W : the initial value of W in SNMF  
    %ReducedDimension : the rank of A


    if ~exist('n_pca', 'var') || isempty(n_pca) 
        n_pca =100;
    end
    
    if ~exist('bound', 'var') || isempty(bound) 
        bound = 0.085;
    end
    
    if ~exist('figure_state', 'var') || isempty(figure_state) 
        figure_state = 1;
    end 
    
    [V,D,~] = svds(A ,n_pca);

    initial_W =V;
    latent_new = diag(D);
    latent_new = latent_new(latent_new > 0);

    initial_W_P = max(initial_W, 0);
    initial_W_N = max(-initial_W, 0);
    n_P_big = sum(initial_W_P.* initial_W_P);
    n_N_big = sum(initial_W_N.* initial_W_N);
    if sum(n_P_big >=n_N_big)>0
        initial_W(:,n_P_big >=n_N_big) = initial_W_P(:,n_P_big >=n_N_big);
    end
    if sum(n_P_big <n_N_big)>0
        initial_W(:,n_P_big <n_N_big) = initial_W_N(:,n_P_big <n_N_big);
    end

    n_svd = min(50,length(latent_new)-1);
    t = 1:n_svd;
    latent_new = latent_new/latent_new(1);
    %             rate = cumsum(latent_new)/sum(latent_new);
    latent_new_diff = latent_new(1:end-1)./ latent_new(2:end)-1;
    for i =1:length(latent_new_diff) -10
        ReducedDimension = i;
        if (latent_new_diff(i) >= bound) &&all(latent_new_diff(i+1 : i+10)<bound )
            break;
        end
    end
    ReducedDimension = max(ReducedDimension , 3);
    if figure_state
        figure
        hold on
        plot(t ,  latent_new_diff(1:n_svd),'xr');
        plot(t , repmat(bound , 1 ,n_svd),'k');
        hold off
    end
end

function [W,H] =recoverWNMF(A,initial_W , k,step,lambda,bound_error)
    [m , n] = size(A);           
    W_in = initial_W(:,1:k);%+ 0.000001 * rand(m ,k);
    H_in = zeros(k ,n);
    A_F = full(A(:)'*A(:));
    error_old = 10000000*A_F;
    error_min = error_old;            
    iter = 1;
    count = 0;
    W_old = W_in;
    H_old = H_in;
    count_rise = 0;%calculate the number of raising steps;
    opt = optimset('TolX',1e-6);             
    while count < 3  && iter<=step
            re_H = zeros(n,1); 
            parfor i = 1:n                                   
                C_cell = W_in;
                C_cell(full(A(:,i)==0),:) = lambda*C_cell(full(A(:,i)==0),:);
                b = C_cell' * A(:,i);
                C_cell = C_cell'*C_cell;
                [H_in(:,i),re_H(i,1)] = lsqnonneg(full(C_cell),full(b),opt);
            end
            re_W = zeros(m,1);
            parfor i = 1:m
                C_gene = H_in;
                C_gene(:,full(A(i,:)==0)) = lambda*C_gene(:,full(A(i,:)==0));
                b = C_gene*(A(i,:)');
                C_gene = C_gene*(C_gene');
                [W_in(i,:),re_W(i,1)] = lsqnonneg(full(C_gene),full(b),opt);
            end
            error_new = (sum(re_H) + sum(re_W))/A_F;
    %                     error_new = norm((A - W_in*H_in).*(A>0),'fro')^2/A_F;
            if error_new > (error_old + bound_error)
                count_rise = count_rise + 1;
            else
                count_rise = 0;
            end

            if error_new < error_min
                 error_min = error_new;
                 W_old = W_in;
                 H_old = H_in;
            end
            if count_rise>3    
                W = W_old;
                H = H_old;
                return;
            end

            if (abs(error_new - error_old)<= bound_error)
                count = count +1;                       
            else
                count =0;
            end
            error_old = error_new;

            iter = iter+1;
             pause(0.00001);
             if iter == step
                    W = W_old;
                    H = H_old;   
                    return;
             end
    end
    W = W_old;
    H = H_old; 
end

function [ n_rank, lambda, max_iters, tol, bound, figure_state] = set_options(options)
    % function [m,tol,max_iters] = set_options(options)
    % Set optionally defined user input parameters.
    % INPUTS:
    %  options: a matlab struct ('m','tol','max_iters').
    % OUTPUTS:
    %  m: the limited memory storage size.
    %  tol: the convergence tolerance criteria.
    %  max_iters: the maximum number of quasi-Newton iterations.
    %  display: true/false should iteration information be displayed?
    %  xhistory: true/false should the entire search history be stored?
    
    n_rank=[];
    lambda=[];
    max_iters=[];
    tol=[];
    bound =[];
    figure_state=[];
    if ( isfield(options, 'n_rank') )
      n_rank = options.n_rank;
    end
    if ( isfield(options, 'lambda') )
      lambda = options.lambda;
    end
    if ( isfield(options, 'tol') )
      tol = options.tol;
    end
    if ( isfield(options, 'max_iters') )
      max_iters = options.max_iters;
    end
    if ( isfield(options, 'bound') )
      bound = options.bound;
    end
    if ( isfield(options, 'figure_state') )
      figure_state = options.figure_state;
    end

end
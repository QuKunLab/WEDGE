function calculate_Pearson_Correlation
%%
clc
clear
Data_set = {'Baron' , 'Zeisel'}
index_data =2;
% Mark_gene = xlsread(".\Baron_old\Mark_gene.csv");
Mark_gene = xlsread(['.\',Data_set{index_data},'\Mark_gene.csv']);
Mark_gene(:,2) = Mark_gene(:,2);

    ID = importdata(['./',Data_set{index_data},'/Reference.csv'],',',1);
    ref = ID.data;

    ref = ref .*  repmat(10000./(sum(ref) + 0.000000001),size(ref,1) ,1); 
    ref = log(ref + 1);
    ref = ref(Mark_gene(:,2)>0,:);
    
    ID = importdata(['./',Data_set{index_data},'/','Observed.csv'],',',1);
    init = ID.data;
    init = init .*  repmat(10000./(sum(init) + 0.000000001),size(init,1) ,1); 
    init = log(init + 1);
 
    [cell(:,1) ,gene(:,1)]= coeff(ref,init(Mark_gene(:,2)>0,:));    
    
    ID = importdata(['./',Data_set{index_data},'/WEDGE_recovery.csv'],',',1);
    ID.data=  ID.data .*  repmat(10000./(sum( ID.data) + 0.000000001),size( ID.data,1) ,1);  
    [cell(:,2) ,gene(:,2)]= coeff(ref,ID.data(Mark_gene(:,2)>0,:));  

    
cell_corr = median(cell);
gene_corr =  median(gene);
    
sprintf('cell: meadian correlation of Obverved is %f, meadian correlation of WEDGE_recovery is %f',cell_corr(1),cell_corr(2))
sprintf('gene: meadian correlation of Obverved is %f, meadian correlation of WEDGE_recovery is %f',gene_corr(1),gene_corr(2))
    
    

end
function  [cell_wise ,gene_wise]= coeff(A,B)
    A_c = A - repmat(mean(A),size(A,1),1);
    B_c = B - repmat(mean(B),size(B,1),1);
    cell_wise = sum(A_c.*B_c)./(std(A).*std(B));
    cell_wise = cell_wise'/(size(A,1) - 1);

    A_g = A - repmat(mean(A,2),1,size(A,2));
    B_g = B - repmat(mean(B,2),1,size(B,2));
    gene_wise = sum(A_g.*B_g,2)./(std(A').*std(B'))';
    gene_wise = gene_wise/(size(A,2) - 1);
end
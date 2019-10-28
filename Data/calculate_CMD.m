function calculate_CMD
    %%
    clc
    clear
    Data_set = {'Baron' , 'Zeisel'}
    index_data =  1;
    Mark_gene = xlsread(['./',Data_set{index_data},'/Mark_gene.csv']);
    Mark_gene =  Mark_gene(:,2);
    
    ID = importdata(['./',Data_set{index_data},'/Reference.csv'],',',1);
    ID.data=  ID.data .*  repmat(10000./(sum( ID.data) + 0.000000001),size( ID.data,1) ,1);  
    ID.data = log( ID.data + 1);

    [ref_cell_wise ,ref_gene_wise]= coeff(ID.data(Mark_gene>0,:));  

    
    ID = importdata(['./',Data_set{index_data},'/','Observed.csv'],',',1);
    init = ID.data;
    init = init .*  repmat(10000./(sum(init) + 0.000000001),size(init,1) ,1); 
    init = log(init + 1);
    
    [cell_wise ,gene_wise]= coeff(init(Mark_gene>0,:));  
    cmd_cell(1) = calcuCMD(ref_cell_wise,cell_wise);
    cmd_gene(1) = calcuCMD(ref_gene_wise,gene_wise);    
    
    ID = importdata(['./',Data_set{index_data},'/WEDGE_recovery.csv'],',',1);
    ID.data=  ID.data .*  repmat(10000./(sum( ID.data) + 0.000000001),size( ID.data,1) ,1);  
    [cell_wise ,gene_wise]= coeff(ID.data(Mark_gene>0,:));  
    cmd_cell(2) = calcuCMD(ref_cell_wise,cell_wise);
    cmd_gene(2) = calcuCMD(ref_gene_wise,gene_wise);    

sprintf('cell: CMD of Obverved is %f, CMD of WEDGE_recovery is %f',cmd_cell(1),cmd_cell(2))

sprintf('gene: CMD of Obverved is %f, CMD of WEDGE_recovery is %f',cmd_gene(1),cmd_gene(2))

end


function  [cell_wise ,gene_wise]= coeff(A)
    cell_wise = corrcoef(A);
    gene_wise = corrcoef(A');
end

function cmd = calcuCMD(A,B)
    cmd = 1 - trace(A*B)/(norm(A,'fro')*norm(B,'fro'));
end
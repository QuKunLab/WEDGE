%%
%Jaccard Index
path = "./WEDGE/Data/figure4/";
T = readtable(path + "/Ref_Label.csv",'Delimiter',',','ReadVariableNames',true);
Ref=T.x;
T = readtable(path + "/WEDGE_Label.csv",'Delimiter',',','ReadVariableNames',true);
WEDGE=T.x;

WEDGE_m=[];
 Ref_m = [];
for i=1:length(Ref)
    WEDGE_m(i,1)=str2num(WEDGE{i})+1;
    Ref_m(i,1) = str2num(Ref{i})+1;
end

A = [];
for i =1:length(unique(Ref_m))
    for j = 1: length(unique(WEDGE_m))
        A(i,j) = sum(WEDGE_m==j&Ref_m==i)/sum(WEDGE_m==j|Ref_m==i);
    end
end

id_v =[];
for i = 1:size(A,1)
    id = find(A(i,:)==max(A(i,:)));
    id_v = [id_v,id];
end

id = [unique(id_v,"stable"),setdiff(1:size(A,2),unique(id_v,"stable"))];
cmap = colormap("parula");
h = heatmap(A(:,id),'Colormap',cmap);
C = {};
for i = 1:length(id)
    C{i} = num2str(id(i));
end
Cy = {};
for i = 1:size(A,1)
    Cy{i} = num2str(i);
end

h.XDisplayLabels = C;
h.YDisplayLabels = Cy;
grid off

%%
%plot ARI
figure
ARI = [0.70, 0.59, 0.58,0.56, 0.52, 0.43];
X = categorical({'WEDGE', 'ALRA', 'SAVERX','DCA', 'ENHANCE', 'MAGIC'});
X = reordercats(X,{'WEDGE', 'ALRA', 'SAVERX','DCA', 'ENHANCE', 'MAGIC'});   
b = bar(X,ARI',0.4,'FaceColor','flat');
ylim([0 0.75])

%%
%plot overlap
figure
X = categorical({'WEDGE', 'ALRA', 'SAVERX', 'ENAHNCE', 'DCA', 'MAGIC'});
X = reordercats(X,{'WEDGE', 'ALRA', 'SAVERX', 'ENAHNCE', 'DCA', 'MAGIC'});  
overlap = [0.993808  0.9597523 0.9287926 0.9845201 0.8482972 0.7182663];
b = bar(X,overlap,0.4,'FaceColor','flat');
ylim([0 1.1])



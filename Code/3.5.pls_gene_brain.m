%% PLS analysis for gene_brain
data_dir = '/Users/buxuan/BNU_Helab/data/nm_subtyping_adhd/analysis_v2/CT_subtypes_v2/pls/gene';
cortical_t = csvread([data_dir,'/deviation_stats/tval_subtype2_left.csv'],1,0);
gene = csvread([data_dir,'/left_DS01_expression.csv'],1,1);
out_dir = fullfile(data_dir,'/deviation_stats/output_gene/');

X = zscore(gene);
Y = zscore(cortical_t);
dim = 10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim,'CV',5);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);

%align PLS components with desired direction%
R1 = corr([XS(:,1),XS(:,2),XS(:,3)],Y(:,1));
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end
%% Spin test to assess the significance of PLS component variance explained ratios
surrogate_tmaps = csvread([data_dir,'/deviation_stats/nullspin_tmap5k_subtype2.csv'],1,0);

PCTVARrand = zeros(5000,dim);
Rsq = zeros(5000,1);
for j=1:5000
    disp(j);
    Yp=zscore(surrogate_tmaps(:,j));
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(zscore(gene),Yp,dim);
    PCTVARrand(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim); 
end 
p_single = zeros(1,dim);
for k=1:dim
    p_single(k)=length(find(PCTVARrand(:,k)>=PCTVAR(2,k)))/5000;
end
p_cum = length(find(Rsq>=Rsquared))/5000;
clear j Yp XLr YLr XSr YSr BETAr PCTVARr MSEr statsr temp
myStats=[PCTVAR; p_single];
csvwrite([out_dir,'PLS_stats2_spin.csv'],myStats);

%% save score,weights, and calculate correlation
csvwrite([out_dir,'PLS_scores2.csv'],XS(:,1));
csvwrite([out_dir,'PLS_gene_weights2.csv'],stats.W(:,1));

[r_corr,p_corr] = corr(XS(:,1),Y);
corr_real = r_corr;
% permutation 5000times for correlation
corr_surr = zeros(1,5000);
for j = 1:5000
    disp(j);
    corr_surr(1,j) = corr(XS(:,1),surrogate_tmaps(:,j));
end
p = zeros(1,1);
p = length(find(corr_surr(1,:)>r_corr))/5000;
corStats=[r_corr; p_corr];
csvwrite([out_dir,'PLS_cor2.csv'],corStats);

%% calculate corrected weight
gene_name = importdata([data_dir,'/gene_namesAHBA.txt']);
geneindex=1:size(gene,2);
bootnum=5000;
X = zscore(gene);
Y = zscore(cortical_t);
dim=1;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

[R1,p1]=corr(XS(:,1),cortical_t);
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=gene_name(x1);
geneindex1=geneindex(x1);

PLS1weights = zeros(length(gene_name),bootnum);

parfor i=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights(:,i) = newW;%store (ordered) weights from this bootstrap run
    
    %temp=stats.W(:,2);%extract PLS2 weights
    %newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    %if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
    %    newW=-1*newW;
    %end
    %PLS2weights(:,i) = newW; %store (ordered) weights from this bootstrap run    
end

PLS1sw = std(PLS1weights');
temp1=PLS1w./PLS1sw';
[Z1,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
fid1 = fopen([out_dir,'PLS_geneWeights_corrected_subtype2.csv'],'w');
for i=1:length(gene_name)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i},geneindex1(i), Z1(i));
end
fclose(fid1);
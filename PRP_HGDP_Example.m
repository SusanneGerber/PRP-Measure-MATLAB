%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matlab-Script from the paper
%%% S. Gerber and I. Horenko 
%%% "On computation of relations between categorical data 
%%% sets and application to genomics", 2016
%%%
%%% It demonstrates an application of the Posterior Relation Probability
%%% measure (PRP-measure) to an analysis of a part of the SNP data from  
%%%  
%%% L. Huang et.al. Genetic Epideiology, 2011, 35(8): 766-780
%%%
%%%
%%%     Copyright (C) 2016  (Illia Horenko and Susanne Gerber, 
%%%     http://icsweb.inf.unisi.ch/cms/index.php/people/20-illia-horenko.html) 
%%%     The full license and most recent version of the code can be found on GitHub
%%%
%%%     This program is free software: you can redistribute it and/or modify
%%%     it under the terms of the GNU General Public License as published by
%%%     the Free Software Foundation, either version 3 of the License, or
%%%     (at your option) any later version.
%%% 
%%%     This program is distributed in the hope that it will be useful,
%%%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%%%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%%     GNU General Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear all
close all

%%% Loads a subset of 300 SNPs for 1107 individuals from 
%%% L. Huang et.al. Genetic Epideiology, 2011, 35(8): 766-780
%%% M_{i,j}='+1' denotes homozigoticity wrt.the first allele of SNP j  
%%% for an individuum i;M_{i,j}='-1' denotes homozigoticity 
%%% wrt.the second allele of SNP j of individuum i; M_{i,j}='0' denotes heterozigoticity 
%%% of SNP j for individuum i  

load HGDP_Data_300SNPs_phased.mat

%%% Setting the significance level for the PRP- and for the
%%% cross-correlation-measures 

SignificanceLevel=0.05;

%%% Main program

[alphaX,N_alphaX] = GetCategoriesNames(reshape(M',1,size(M,1)*size(M,2)),'');
names=alphaX;
[T,N]=size(M);

clear PRP P EV C
C=zeros(1,N);PVAL=C;PRP=C;EV=C;

%%% Computation of the cross-correlation and the PRP-measure and
%%% identification of the significant relations for the both measures 

for i=1:N
    i
    for j=1:N
        [C(i,j),PVAL(i,j)]=corr(M(:,i),M(:,j));
        PRP(i,j) = PRP_measure_new(M(:,i)',M(:,j)',alphaX,N_alphaX,alphaX,N_alphaX,'BIC');
    end
end




%%% Setting font size and line width for the output, visualizing the
%%% results

FS=16;LW=3;

figure;imagesc(PRP);%imagesc((PRP>1-SignificanceLevel).*PRP);
caxis([0 1]);set(gca,'LineWidth',1);cmap = flipud(colormap('hot'));colormap(cmap);
cmap2=[cmap(size(cmap,1):-2:1,3) cmap(size(cmap,1):-2:1,2) cmap(size(cmap,1):-2:1,1);cmap(1:2:size(cmap,1),:)];
colorbar;xlabel('SNP j','FontSize',FS);ylabel('SNP i','FontSize',FS);
title('PRP-measure','FontSize',FS)

figure;imagesc((PVAL<SignificanceLevel).*C);
caxis([-1 1]);colormap(cmap2)
xlabel('SNP j','FontSize',FS);ylabel('SNP i','FontSize',FS);
title('(Cross-)Correlation measure','FontSize',FS)
   
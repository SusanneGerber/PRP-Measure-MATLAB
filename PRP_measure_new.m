function [PRP,P,c,WW,EV_measure,mIC] = PRP_measure_new(X,Y,alphaX,N_alphaX,alphaY,N_alphaY,IC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%     This code realizes a computation of the PRP relation measure between
%%%     a pair of categorical variables X and Y.
%%%
%%%     Copyright (C) 2016  (Illia Horenko and Susanne Gerber, 
%%%     http://icsweb.inf.unisi.ch/cms/index.php/people/20-illia-horenko.html) 
%%%     The full license and most recent version of the code can be found
%%%     on GitHub.
%%%
%%%     Citation: S. Gerber and I. Horenko 
%%%     "On computation of relations between categorical data 
%%%      sets and application to genomics", 2016
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
T=length(X);



N_cd=zeros(N_alphaY,N_alphaX);
LogL_cd=zeros(N_alphaY,N_alphaX);
N_ci=zeros(N_alphaY,1);
for state_i=1:N_alphaY
    for state_j=1:N_alphaX
        N_cd(state_i,state_j)=sum((Y==alphaY{state_i}).*(X==alphaX{state_j}));
    end
    N_ci(state_i)=sum((Y==alphaY{state_i}));
end


TTT=sum(sum(N_cd));
NN_ci=sum(N_ci);

Combinations=dec2bin([0:2^N_alphaY-1]');N_c=size(Combinations,1);
IC_cd=zeros(N_c,1);EV_measure=zeros(N_c,1);aa=zeros(N_c,N_alphaY);
for ind_m=1:N_c
    a=Combinations(ind_m,:);
    for i_a=1:N_alphaY
        aa(ind_m,i_a)=str2double(a(i_a));
    end
end

EV_measure=zeros(N_c,1);
LLL=EV_measure;
options=optimset('GradObj','on','Algorithm','sqp','Display','off');
for ind_m=1:N_c
    %ind_m
    [~,ii_ci]=find(aa(ind_m,:)==0);[~,ii_cd]=find(aa(ind_m,:)==1);
    N_par=0;
    A=[-eye(N_alphaY*N_alphaX)];b=zeros(N_alphaY*N_alphaX,1);
    for state_j=1:N_alphaX
        for state_i=1:N_alphaY
            Aeq1(state_j,state_j+N_alphaX*(state_i-1))=1;
        end
    end
    beq1=ones(N_alphaX,1);
    kkk=1;Aeq2=zeros(length(ii_ci)*(N_alphaX-1),N_alphaY*N_alphaX);
    beq2=zeros(length(ii_ci)*(N_alphaX-1),1);
    for state_i=ii_ci
        for state_j=1:N_alphaX-1
            Aeq2(kkk,state_j+N_alphaX*(state_i-1))=1;
            Aeq2(kkk,state_j+1+N_alphaX*(state_i-1))=-1;
            kkk=kkk+1;
        end
    end
    Aeq=[Aeq1;Aeq2];beq=[beq1;beq2];
    xxx0=ones(1,N_alphaY*N_alphaX)./(N_alphaY*N_alphaX);
    [xxx,fff,flag,output] =  fmincon(@(x)LogLikPairwiseCausality...
        (x,reshape(N_cd',N_alphaY*N_alphaX,1)')...
        ,xxx0,(A),(b),Aeq,beq,[],[],[],options);
    P(:,:,ind_m)=reshape(xxx,N_alphaX,N_alphaY)';
    EV_measure(ind_m)=0;
    for state_i=1:N_alphaY
        for state_j1=1:N_alphaX
            for state_j2=1:N_alphaX
                EV_measure(ind_m)=EV_measure(ind_m)+(P(state_i,state_j1,ind_m)...
                    -P(state_i,state_j2,ind_m))^2;
            end
        end
    end
    LLL(ind_m)=-fff;
    N_par(ind_m)=length(ii_ci)+N_alphaX*length(ii_cd);
    if strcmp(IC,'BIC'),
        IC_cd(ind_m)=-2*LLL(ind_m)+(N_par(ind_m))*(log(TTT)-log(2*pi));%+2*N_par_cd;%+2*N_par_cd*(N_par_cd+1)/(T-N_par_cd-1);%-log(2*pi));
    else
        IC_cd(ind_m)=-2*LLL(ind_m)+2*(N_par(ind_m))+2*(N_par(ind_m))*(N_par(ind_m)+1)/(TTT-N_par(ind_m)-1);%-log(2*pi));
    end
end
[mIC,ii]=min(IC_cd);
EV_measure=EV_measure(ii);
P=P(:,:,ii);
c=aa(ii,:);
WW(ii,1)=1/(1+sum(exp(-0.5*(IC_cd([1:ii-1 ii+1:N_c])-mIC))));WW([1:ii-1 ii+1:N_c],1)=WW(ii,1)*exp(-0.5*(IC_cd([1:ii-1 ii+1:N_c])-mIC));
PRP=1-WW(1,1);


end


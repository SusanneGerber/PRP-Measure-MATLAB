function [XX,label_new] = LogicalMerge(X,labelX,Y,labelY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%     This code realizes a computation of all non-redundant combinations
%%%     of the data X and Y (rows of X and Y contain 
%%%     single Boolean variable dimensions).
%%%
%%%     Copyright (C) 2016  (Illia Horenko and Susanne Gerber, 
%%%     http://icsweb.inf.unisi.ch/cms/index.php/people/20-illia-horenko.html) 
%%%     The full license and most recent version of the code can be found
%%%     on GitHub.
%%%
%%%     Citation: S. Gerber and I. Horenko 
%%%     "On computation of relations between categorical Omics data sets", 2016
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

NX=size(X,1);NY=size(Y,1);
kkk=1;XX=[];
for i=1:NX
    for j=i:NY
        if i>1
            if j>i
                sss=(X(i,:)).*(Y(j,:));
                dist=zeros(1,kkk-1);
                for iii=1:kkk-1
                    dist(iii)=norm(sss-XX(iii,:));
                end
                if and(sum(abs(diff(sss)))>0,min(dist)>0)
                    XX(kkk,:)=sss;
                    label_new{kkk}=strcat(labelX{i},'; ', labelY{j});
                    kkk=kkk+1;
                end
            end
        else
            XX(kkk,:)=Y(j,:);
            label_new{kkk}=strcat(labelY{j});
            kkk=kkk+1;
        end
    end
    disp([num2str(i/NX*100) ' %'])
end

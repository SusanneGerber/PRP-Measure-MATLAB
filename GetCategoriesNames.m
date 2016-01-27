function [alpha,N_alpha] = GetCategoriesNames_v2(X,empty)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%     This code obtaines all of the categories (alpha)
%%%     involved in the categorical variable X. N_alpha is the number of
%%%     involved categories
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

[T,M]=size(X);

k=1;alpha=[];
for t=1:T
    for n=1:M
        if and(~inset(X(t,n),alpha),~inset(X(t,n),empty))
            alpha{k}=X(t,n);
            k=k+1;
        end
    end
end
N_alpha=numel(alpha);

aa=sort(cell2mat(alpha));
for i=1:N_alpha
    alpha{i}=aa(i);
end


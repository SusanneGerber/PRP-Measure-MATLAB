function [flag,ind] = inset(a,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%     Copyright (C) 2016  (Illia Horenko and Susanne Gerber, 
%%%     http://icsweb.inf.unisi.ch/cms/index.php/people/20-illia-horenko.html) 
%%%     The full license and most recent version of the code can be found
%%%     on GitHub.
%%%
%%%     Citation: S. Gerber and I. Horenko 
%%%               "On computation of relations between categorical data 
%%%                sets and application to genomics", 2016
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
flag=0;
for i=1:numel(X)
    if a==X{i}
        flag=1;
        ind=i;
        break
    end
end
end


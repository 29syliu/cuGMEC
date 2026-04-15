function out = mbind_cell_col(f,C)
% f:: x -> [a]
% mbind_col:: [varargin] >>= f
out = f(C{1});
for i=2:length(C)
    out = [out;f(C{i})];
end
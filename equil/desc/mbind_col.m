function out = mbind_col(f,varargin)
% f:: x -> [a]
% mbind_col:: [varargin] >>= f
out = f(varargin{1});
for i=2:length(varargin)
    out = [out;f(varargin{i})];
end
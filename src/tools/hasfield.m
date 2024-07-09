function flag = hasfield(mfile, fieldname)
%HASFIELD Check if matfile has the result indexed by fieldname.
%
%   mfile: MAT-file object
%   fieldname: string
%
%   flag: boolean

% check if mfile has the entry indexed by fieldname
flag = ~isempty(who(mfile, fieldname));

function structure = rmfields(structure, fields)
%RMFIELDS Remove struct fields.
%   Do not throw error when some fields do not exist in struct.
%   
%   structure: struct
%   fields: cell
%
%   structure: struct


rm_inds = isfield(structure, fields);
structure = rmfield(structure, fields(rm_inds));

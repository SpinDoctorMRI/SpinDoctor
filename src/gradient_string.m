function [fullstring, ampstring, seqstring, dirstring] = gradient_string(amp, amptype, seq, g)
%GRADIENT_STRING Get string representation of gradient.
%
%   amp: [1 x 1]    - q-value or b-value.
%   amptype: string - "q" or "b"
%   seq: Sequence
%   g: [3 x 1]
%
%   fullstring: string - full gradient string
%   ampstring:  string - amplitude part
%   seqstring:  string - sequence part
%   dirstring:  string - direction part


ampstring = sprintf("%s%g", amptype, amp);
seqstring = sprintf("%s_d%g_D%g", class(seq), seq.delta, seq.Delta);
% seqstring = seq.string;
dirstring = sprintf("g[%.2f_%.2f_%.2f]", g);
fullstring = sprintf("%s_%s_%s", seqstring, ampstring, dirstring);

end


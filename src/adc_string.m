function [fullstring, seqstring, dirstring, cmptstring] = adc_string(icmpt, seq, g)
%ADC_STRING Get string representation of ADC iteration.
%
%   icmpt: int
%   seq: Sequence
%   g: [3 x 1]
%
%   fullstring: string - full ADC string
%   seqstring:  string - sequence part
%   dirstring:  string - direction part
%   cmptstring:  string - compartment part


cmptstring = sprintf("cmpt%d", icmpt);
seqstring = sprintf("%s_d%g_D%g", class(seq), seq.delta, seq.Delta);
% seqstring = seq.string;
dirstring = sprintf("g[%.2f_%.2f_%.2f]", g);
fullstring = sprintf("%s_%s_%s", seqstring, dirstring, cmptstring);

end


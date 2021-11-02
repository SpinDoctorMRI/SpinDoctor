% Author: Jan Simon
% Source: https://www.mathworks.com/matlabcentral/answers/3314-hash-function-for-matlab-struct#answer_4930

function H = DataHash(Data, n)
whichMEX  = which(['GetMD5.', mexext]);
if isempty(whichMEX)
   InstallMex('GetMD5.c');
end
whichMEX  = which(['GetMD5.', mexext]);
if isempty(whichMEX)
   error('GetMD5 installation failed. Please check it or replace DataHash with your favorite hash function.')
end

H = CoreHash(Data);
H = sprintf('%.2x', H);   % To hex string
if nargin == 2 && n>=1 && n<=32 
   H = H(1: round(n)); % Take only the first n characters
end

function H = CoreHash(Data)
if isstring(Data); Data = char(Data); end
% Consider the type of empty arrays:
S = [class(Data), sprintf('%d ', size(Data))];
H = GetMD5(typecast(uint16(S(:)), 'uint8'), 'bin', 'uint8');
if isa(Data, 'struct')
   n = numel(Data);
   if n == 1  % Scalar struct:
      F = sort(fieldnames(Data));  % ignore order of fields
      for iField = 1:length(F)
         H = bitxor(H, CoreHash(Data.(F{iField})));
      end
   else  % Struct array:
      for iS = 1:n
         H = bitxor(H, CoreHash(Data(iS)));
      end
   end
elseif isempty(Data)
   % No further actions needed
elseif isnumeric(Data)
   if isreal(Data)
      data = typecast(Data(:), 'uint8');
      H = bitxor(H, GetMD5(data, 'bin', 'uint8'));
   else
      data = typecast(real(Data(:)), 'uint8');
      H = bitxor(H, GetMD5(data, 'bin', 'uint8'));
      data = typecast(imag(Data(:)), 'uint8');
      H = bitxor(H, GetMD5(data, 'bin', 'uint8'));
   end
elseif ischar(Data)  % Silly TYPECAST cannot handle CHAR
   data = typecast(uint16(Data(:)), 'uint8');
   H = bitxor(H, GetMD5(data, 'bin', 'uint8'));
elseif iscell(Data)
   for iS = 1:numel(Data)
      H = bitxor(H, CoreHash(Data{iS}));
   end
elseif islogical(Data)
   data = typecast(double(Data(:)), 'uint8');
   H = bitxor(H, GetMD5(data, 'bin', 'uint8'));
elseif isa(Data, 'function_handle')
   H = bitxor(H, CoreHash(functions(Data)));
else
   warning(['Type of variable not considered: ', class(Data)]);
end

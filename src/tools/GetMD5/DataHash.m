% Author: Jan Simon
% Source: https://www.mathworks.com/matlabcentral/answers/3314-hash-function-for-matlab-struct#answer_4930

function H = DataHash(Data, n)
   % check if GetMD5 is compiled
   whichMEX  = which(['GetMD5.', mexext]);
   if isempty(whichMEX)
      InstallMex('GetMD5.c');
   end
   whichMEX  = which(['GetMD5.', mexext]);
   if isempty(whichMEX)
      error('GetMD5 installation failed. Please check it or replace DataHash with your favorite hash function.')
   end

   H = CoreHash(Data);
   if nargin == 2 && n>=1 && n<=32 
      H = H(1: round(n)); % Take only the first n characters
   end
end

function H = CoreHash(Data)
   H = '';
   if isstring(Data)
      Data = char(Data); 
   end
   % Consider the type of empty arrays:
   S = [class(Data), sprintf('%d ', size(Data))];
   H = append(H, GetMD5(S));
   if isa(Data, 'struct')
      n = numel(Data);
      if n == 1  % Scalar struct:
         F = sort(fieldnames(Data));  % ignore order of fields
         for iField = 1:length(F)
            H = append(H, CoreHash(Data.(F{iField})));
         end
      else  % Struct array:
         for iS = 1:n
            H = append(H, CoreHash(Data(iS)));
         end
      end
   elseif isempty(Data)
      % No further actions needed
   elseif isnumeric(Data)
      if isreal(Data)
         H = append(H, GetMD5(Data));
      else
         H = append(H, GetMD5(real(Data)), GetMD5(imag(Data)));
      end
   elseif ischar(Data)
      H = append(H, GetMD5(Data));
   elseif iscell(Data)
      for iS = 1:numel(Data)
         H = append(H, CoreHash(Data{iS}));
      end
   elseif islogical(Data)
      H = append(H, GetMD5(Data));
   elseif isa(Data, 'function_handle')
      H = append(H, GetMD5(func2str(Data)));
   elseif isa(Data, "Sequence")
      H = append(H, GetMD5(char(Data.string())));
   else
      warning(['Type of variable not considered: ', class(Data)]);
   end
   H = GetMD5(H);
end

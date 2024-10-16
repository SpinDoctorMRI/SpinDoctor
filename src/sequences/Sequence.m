classdef Sequence < AbsSequence
    %SEQUENCE Gradient sequence (time profile)
    %   This represents an abstract fixed direction gradient sequence.
    
    properties
        delta
        Delta
        TE
        t1
    end
    
    methods
        function obj = Sequence(delta, Delta,varargin)
            %SEQUENCE Construct an instance of this class.
            %   Here it is assumed that the sequence is parametrized by the two
            %   parameters `delta` and `Delta` only. Subclasses may have more
            %   parameters.
            assert(delta > 0 && Delta > 0);
            if Delta<delta
                error('Sequence: delta should be less than or equal to Delta.') 
            else
                obj.delta = delta;
                obj.Delta = Delta;
                if isempty(varargin)
                    obj.TE = obj.Delta + obj.delta;
                    obj.t1 = 0;
                elseif nargin == 3
                    obj.TE = varargin{1};
                    obj.t1 = 0;
                elseif nargin == 4
                    if isnumeric(varargin{2})
                        obj.t1 = varargin{1};
                        obj.TE = obj.t1 + delta +Delta + varargin{2};
                    elseif varargin{2} == "Symmetric" || varargin{2} == "symmetric" 
                        obj.TE = varargin{1};
                        obj.t1 = (obj.TE - Delta - delta)/2;
                    else
                        error("Error: invalid input into Sequence")
                    end
                else
                    error("Error: invalid input into Sequence")
                end
            end
        end
        
        function TE = echotime(obj)
            %ECHOTIME Get the echo time of the sequence.
            %   By default, this method returns `Delta + delta`.
            TE = obj.Delta + obj.delta;
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the time profile from `0` to `t`.
            %   Unless overwritten, it computes a numerical approximation.
            F = arrayfun(@(s) integral(@obj.call, 0, s, "AbsTol", 1e-6, ...
                "RelTol", 1e-3, "Waypoints", [obj.delta+obj.t1, obj.Delta+obj.t1,obj.Delta+obj.t1+obj.delta]), t);
        end

        
        function s = string(obj, simplified)
            %STRING Convert sequence to string.
            %   If there are other parameters than `delta` and `Delta`, this
            %   method should be overwritten.
            if nargin == 2 && simplified
                s = sprintf("%s_d%g_D%g", class(obj), obj.delta, obj.Delta);
            else
                s = sprintf("%s(delta=%g, Delta=%g)", class(obj), ...
                    obj.delta, obj.Delta);
            end
        end
        
    end
end

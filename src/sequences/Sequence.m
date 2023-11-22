classdef (Abstract) Sequence
    %SEQUENCE Gradient sequence (time profile)
    %   This represents an abstract gradient sequence.
    
    properties
        delta
        Delta
        TE
        t1
    end
    
    methods (Abstract)
        f = call(obj, t)
        %CALL Call the time profile at time `t`.
        
        t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of sequence.
        
        [timelist, interval_str, timeprofile_str] = intervals(obj)
        %INTERVALS Get intervals of the sequence.
        %   This function returns a list of important time steps (including
        %   start and stop), a list of strings representing the intervals
        %   between these time steps and a list of strings representing the
        %   behavior of the sequence on each of these intervals.
    end
    
    methods
        function obj = Sequence(delta, Delta,varargin)
            %SEQUENCE Construct an instance of this class.
            %   Here it is assumed that the sequence is parametrized by two
            %   required parameters `delta` and `Delta`. Optional parameters include:
            % (TE), (TE,"symmetric"), (t1,t2) where TE is the echotime, t1
            % a pause before the first pulse, t2 a pause before the second
            % pulse. Subclasses may have more parameters.
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
        
        function int = integral_F2(obj)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(obj, t)).
            %   Unless overwritten, it computes a numerical approximation.
            int = integral(@(s) obj.integral(s).^2, 0, obj.echotime, ...
            "AbsTol", 1e-6, "RelTol", 1e-3);
        end
        
        function b = bvalue_no_q(obj)
            %BVALUE_NO_Q Compute the time profile contribution to the b-value.
            b = obj.integral_F2;
        end
        
        function jn = J(obj, lambda, ninterval)
            %J Compute the quantity J(lambda) for the sequence.
            %   Unless overwritten, it computes a numerical approximation.
            if nargin == 2
                ninterval = 500;
            end
            ntime = ninterval + 1;
            time = linspace(0, obj.TE, ntime);
            dtime = obj.TE / ntime;
            
            fprintf("Use numerical integration to compute the the quantity J(lambda).\n");
            jn = @(t) integral(@(s) exp(lambda * (s - t)) .* obj.call(s), 0, t, ...
                    "Waypoints", [obj.delta+obj.t1, obj.Delta+obj.t1,obj.Delta+obj.t1+obj.delta]);
            jn = lambda * obj.integral(time) .* arrayfun(jn, time);
            jn = dtime * trapz(jn) / obj.integral_F2;
        end
        
        function s = string(obj, simplified)
            %STRING Convert sequence to string.
            %   If there are other parameters than `delta` and `Delta`, this
            %   method should be overwritten.
            if nargin == 2 && simplified
                s = sprintf("%s_d%g_D%g", class(obj), obj.delta, obj.Delta);
            else
                s = sprintf("%s(delta=%g, Delta=%g, TE=%g, t1=%g)", class(obj), ...
                    obj.delta, obj.Delta,obj.TE,obj.t1);
            end
        end
        
        function s = char(obj)
            %CHAR Convert sequence to character array.
            %   This is a wrapper for the STRING function.
            s = char(obj.string);
        end
    end
end

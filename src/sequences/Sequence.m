classdef (Abstract) Sequence
    %SEQUENCE Gradient sequence (time profile)
    %   This represents an abstract gradient sequence.
    
    properties
        delta
        Delta
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
        function obj = Sequence(delta, Delta)
            %SEQUENCE Construct an instance of this class.
            %   Here it is assumed that the sequence is parametrized by the two
            %   parameters `delta` and `Delta` only. Subclasses may have more
            %   parameters.
            if Delta<delta
                error('Sequence: delta should be less than or equal to Delta.') 
            else
                obj.delta = delta;
                obj.Delta = Delta;
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
                "RelTol", 1e-3, "Waypoints", [obj.delta, obj.Delta]), t);
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
            time = linspace(0, seq.echotime, ntime);
            dtime = seq.echotime / ntime;
            
            fprintf("Use numerical integration to compute the the quantity J(lambda).\n");
            jn = @(t) integral(@(s) exp(lambda * (s - t)) .* obj.call(s), 0, t, ...
                    "Waypoints", [obj.delta, obj.Delta]);
            jn = lambda * obj.integral(time) .* arrayfun(jn, time);
            jn = dtime * trapz(jn) / obj.integral_F2;
        end
        
        function s = string(obj)
            %STRING Convert sequence to string.
            %   If there are other parameters than `delta` and `Delta`, this
            %   method should be overwritten.
            % s = sprintf("%s(delta = %g, Delta = %g)", class(obj), obj.delta, obj.Delta);
            s = sprintf("%s(%g, %g)", class(obj), obj.delta, obj.Delta);
        end
        
        function s = char(obj)
            %CHAR Convert sequence to character array.
            %   This is a wrapper for the STRING function.
            s = char(obj.string);
        end
    end
end

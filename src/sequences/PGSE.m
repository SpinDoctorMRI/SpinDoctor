classdef PGSE < Sequence
    %PGSE Pulsed Gradient Spin Echo sequence.
    %   This sequence consists of two opposite pulses of duration `delta`
    %   separated by a pause of duration `Delta - delta`.
    
    methods
        function f = call(obj, t)
            %CALL Call the PGSE time profile at time t.
            %   The function is vectorized.
            f = (t < obj.delta) - (obj.Delta <= t);
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the PGSE sequence from 0 to t.
            %   An analytical expression is available.
            F = (t < obj.delta) .* t ...
                + (obj.delta <= t) .* obj.delta ...
                - (obj.Delta <= t) .* (t - obj.Delta);
        end
        
        function b = bvalue_no_q(obj)
            %BVALUE_NO_Q Compute the time profile contribution to the b-value.
            %   An analytical expression is available for the PGSE sequence.
            b = obj.delta^2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the PGSE sequence.
            t = obj.Delta - obj.delta / 3;
        end
    
        function [timelist, interval_str, timeprofile_str] = intervals(obj)
            %INTERVALS Get intervals of the sequence.
            %   This function returns a list of important time steps (including
            %   start and stop), a list of strings representing the intervals
            %   between these time steps and a list of strings representing the
            %   behavior of the sequence on each of these intervals.
            timelist = [0, obj.delta, obj.Delta, obj.Delta+obj.delta];
            interval_str = ["[0, delta]", "[delta, Delta]", "[Delta, Delta+delta]"];
            timeprofile_str = ["f(t) = 1 (constant)", "f(t) = 0 (constant)", "f(t) = -1 (constant)"];
            if obj.delta == obj.Delta
                % Remove unused interval
                timelist(3) = [];
                interval_str(2) = [];
                timeprofile_str(2) = [];
            end
        end
    end
end


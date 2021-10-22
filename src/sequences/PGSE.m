classdef PGSE < Sequence
    %PGSE Pulsed Gradient Spin Echo sequence.
    %   This sequence consists of two opposite pulses of duration `delta`
    %   separated by a pause of duration `Delta - delta`.
    
    methods
        function f = call(obj, t)
            %CALL Call the PGSE time profile at time t.
            %   The function is vectorized.
            f = (0 <= t & t < obj.delta) - (obj.Delta <= t & t <= obj.echotime);
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the PGSE sequence from 0 to t.
            %   An analytical expression is available.
            F = (0 <= t & t < obj.delta) .* t ...
                + (obj.delta <= t & t <= obj.echotime) .* obj.delta ...
                - (obj.Delta <= t & t <= obj.echotime) .* (t - obj.Delta);
        end
        
        function int = integral_F2(obj)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(obj, t)).
            %   An analytical expression is available for the PGSE sequence.
            int = obj.delta^2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time(obj)
            %DIFFUSION_TIME Get diffusion time of the PGSE sequence.
            t = obj.Delta - obj.delta / 3;
        end
        
        function t = diffusion_time_sta(obj)
            %DIFFUSION_TIME_STA Get STA diffusion time of the PGSE sequence.
            d = obj.delta;
            D = obj.Delta;
            out = (4 / 35) * ( ...
                    + (D + d)^(7 / 2) ...
                    + (D - d)^(7 / 2) ...
                    - 2 * D^(7 / 2) ...
                    - 2 * d^(7 / 2) ...
                    ) / (d^2 * (D - d/3));
            t = out^2;
        end

        function jn = J(obj, lambda)
            %J Compute the quantity J(lambda) for the sequence
            %   An analytical expression is available for the PGSE sequence
            d = obj.delta;
            D = obj.Delta;
            
            if lambda < 1e-7
                % Use Taylor expansion when lambda is close to 0 
                % to improve numerical stability
                jn = lambda ...
                    - lambda^2 * D^2 / (2 * (D - d/3)) ...
                    + lambda^3 * (10*D^3 + 5*D*d^2 - d^3) / (20 * (3*D - d)) ...
                    - lambda^4 * (D^4 + D^2*d^2) / (8 * (3*D - d)) ...
                    + lambda^5 * (21*D^5 + 35*D^3*d^2 + 7*D*d^4 - d^5) / (840 * (3*D - d));
            else
                jn = - 1 * ( ...
                    + exp(-lambda * (D + d)) ...
                    + exp(-lambda * (D - d)) ...
                    - 2 * exp(-lambda * d) ...
                    - 2 * exp(-lambda * D) ...
                    + 2 * (1 - lambda * d)) / ...
                    (lambda^2 * d^2 * (D - d/3));
            end
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

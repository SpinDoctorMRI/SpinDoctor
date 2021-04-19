classdef DoublePGSE < Sequence
    %DOUBLEPGSE Double Pulsed Gradient Spin Echo sequence.
    %   This sequence consists of two PGSE sequences without a pause between
    %   them. Each of the two PGSE sequences consists of two opposite pulses of
    %   duration delta, separated by a pause of duration Delta - delta.
    
    methods
        function TE = echotime(obj)
            %ECHOTIME Get the echo time of the DoublePGSE sequence.
            %   The echo time is twice the normal echo time.
            TE = 2 * obj.Delta + 2 * obj.delta;
        end
        
        function f = call(obj, t)
            %CALL Call the DoublePGSE time profile at time t.
            %   The function is vectorized.
            d = obj.delta;
            D = obj.Delta;
            f = (t < d) ...
                - (D <= t & t < D + d) ...
                + (D + d <= t & t < D + 2 * d) ...
                - (2 * D + d <= t);
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the time profile from 0 to t.
            %   For the DoublePGSE sequence, an analytical expression is
            %   available.
            d = obj.delta;
            D = obj.Delta;
            F = (t < d) .* t ...
                + (d <= t & t < D + d) .* d ...
                - (D <= t & t < D + d) .* (t - D) ...
                + (D + d <= t & t < D + 2 * d) .* (t - (D + d)) ...
                + (D + 2 * d <= t) .* d ...
                - (2 * D + d <= t) .* (t - (2 * D + d));
        end
        
        function b = bvalue_no_q(obj)
            %BVALUE_NO_Q Compute the time profile contribution to the b-value.
            %   An analytical expression is available for the DoublePGSE
            %   sequence.
            b = 2 * obj.delta^2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the DoublePGSE sequence.
            t = 2 * (obj.Delta - obj.delta / 3);
        end
        
        function [timelist, interval_str, timeprofile_str] = intervals(obj)
            %INTERVALS Get intervals of the sequence.
            %   This function returns a list of important time steps (including
            %   start and stop), a list of strings representing the intervals
            %   between these time steps and a list of strings representing the
            %   behavior of the sequence on each of these intervals.
            d = obj.delta;
            D = obj.Delta;
            timelist = [0 d D D+d D+2*d 2*D+d 2*(D+d)];
            interval_str = ["[0, delta]", "[delta, Delta]", "[Delta, Delta+delta]", ...
                "[Delta+delta, Delta+2*delta]", "[Delta+2*delta, 2*Delta+delta]", ...
                "[2*Delta+delta, 2*(Delta+delta)]"];
            timeprofile_str = repmat(["f(t) = 1 (constant)", "f(t) = 0 (constant)", "f(t) = -1 (constant)"], 1, 2);
            if d == D
                % Remove unused intervals
                timelist([3 6]) = [];
                interval_str([2 5]) = [];
                timeprofile_str([2 5]) = [];
            end
        end
    end
end


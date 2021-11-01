classdef DoublePGSE < Sequence
    %DOUBLEPGSE Double Pulsed Gradient Spin Echo sequence.
    %   This sequence consists of two PGSE sequences without a pause between
    %   them. Each of the two PGSE sequences consists of two opposite pulses of
    %   duration delta, separated by a pause of duration Delta - delta.
    
    properties
        tpause
    end
    
    methods
        function obj = DoublePGSE(delta, Delta, tpause)
            %DoublePGSE Construct an instance of the DoublePGSE class.
            %   `tpause` is the pause between the two PGSE sequences.
            obj@Sequence(delta, Delta);
            if nargin < nargin(@DoublePGSE)
                tpause = 0;
            end
            assert(tpause >= 0);
            obj.tpause = tpause;
        end
        
        function TE = echotime(obj)
            %ECHOTIME Get the echo time of the DoublePGSE sequence.
            %   The echo time is twice the normal echo time.
            TE = 2 * obj.Delta + 2 * obj.delta + obj.tpause;
        end
        
        function f = call(obj, t)
            %CALL Call the DoublePGSE time profile at time t.
            %   The function is vectorized.
            d = obj.delta;
            D = obj.Delta;
            p = obj.tpause;
            f = (0 <= t & t < d) ...
                - (D <= t & t < D + d) ...
                + (p + D + d <= t & t < p + D + 2 * d) ...
                - (p + 2 * D + d <= t & t <= obj.echotime);
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the time profile from 0 to t.
            %   For the DoublePGSE sequence, an analytical expression is
            %   available.
            d = obj.delta;
            D = obj.Delta;
            p = obj.tpause;
            F = (0 <= t & t < d) .* t ...
                + (d <= t & t < D + d) .* d ...
                - (D <= t & t < D + d) .* (t - D) ...
                + (p + D + d <= t & t < p + D + 2 * d) .* (t - (p + D + d)) ...
                + (p + D + 2 * d <= t & t <= obj.echotime) .* d ...
                - (p + 2 * D + d <= t & t <= obj.echotime) .* (t - (p + 2 * D + d));
        end
        
        function int = integral_F2(obj)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(obj, t)).
            %   An analytical expression is available for the DoublePGSE sequence.
            int = 2 * obj.delta^2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time(obj)
            %DIFFUSION_TIME Get diffusion time of the DoublePGSE sequence.
            t = 2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time_sta(obj)
            %DIFFUSION_TIME_STA Get STA diffusion time of the DoublePGSE sequence.
            d = obj.delta;
            D = obj.Delta;
            tm = obj.delta + obj.tpause;
            out = (2 / 35) * ( ...
                        + (2*D + tm + d)^(7 / 2) ...
                        + (2*D + tm - d)^(7 / 2) ...
                        + (tm + d)^(7 / 2) ...
                        + (tm - d)^(7 / 2) ...
                        - 2 * (2*D + tm)^(7 / 2) ...
                        - 2 * (D + tm + d)^(7 / 2) ...
                        - 2 * (D + tm - d)^(7 / 2) ...
                        + 2 * (D + d)^(7 / 2) ...
                        + 2 * (D - d)^(7 / 2) ...
                        - 2 * tm^(7 / 2) ...
                        + 4 * (D + tm)^(7 / 2) ...
                        - 4 * D^(7 / 2) ...
                        - 4 * d^(7 / 2) ...
                    ) / (d^2 * (D - d/3));
            t = out^2;
        end
        
        function jn = J(obj, lambda)
            %J Compute the quantity J(lambda) for the sequence
            %   An analytical expression is available for the DoublePGSE sequence
            d = obj.delta;
            D = obj.Delta;
            tm = obj.tpause + obj.delta;
            
            if lambda < 1e-7
                % Use Taylor expansion when lambda is close to 0 
                % to improve numerical stability
                jn = lambda ...
                    - lambda^2 * 3 * D^2 / (3*D - d) ...
                    + lambda^3 * (40*D^3 + 5*D*d^2 - d^3 + 30*D^2*tm) / (20 * (3*D - d)) ...
                    - lambda^4 * (4*D^4 + D^2*d^2 + 6*D^3*tm + 3*D^2*tm^2) / (4 * (3*D - d)) ...
                    + lambda^5 * (336*D^5 + 140*D^3*d^2 + 7*D*d^4 - d^5 + 735*D^4*tm ...
                        + 105*D^2*d^2*tm + 630*D^3*tm^2 + 210*D^2*tm^3) / (840 * (3*D - d));
            else
                jn = - 1 * ( ...
                    + exp(-lambda * (2*D + tm - d)) ...
                    + exp(-lambda * (2*D + tm + d)) ...
                    + exp(-lambda * (tm - d)) ...
                    + exp(-lambda * (tm + d)) ...
                    - 2 * exp(-lambda * (D + tm - d)) ...
                    - 2 * exp(-lambda * (D + tm + d)) ...
                    - 2 * exp(-lambda * (2*D + tm)) ...
                    + 2 * exp(-lambda * (D - d)) ...
                    + 2 * exp(-lambda * (D + d)) ...
                    - 2 * exp(-lambda * tm) ...
                    + 4 ...
                    + 4 * exp(-lambda * (D+tm)) ...
                    - 4 * exp(-lambda * D) ...
                    - 4 * exp(-lambda * d) ...
                    - 4 * lambda * d ) / ...
                    (2 * lambda^2 * d^2 * (D - d/3));
            end
        end

        function [timelist, interval_str, timeprofile_str] = intervals(obj)
            %INTERVALS Get intervals of the sequence.
            %   This function returns a list of important time steps (including
            %   start and stop), a list of strings representing the intervals
            %   between these time steps and a list of strings representing the
            %   behavior of the sequence on each of these intervals.
            d = obj.delta;
            D = obj.Delta;
            p = obj.tpause;
            timelist = [0 d D D+d p+D+d p+D+2*d p+2*D+d p+2*(D+d)];
            interval_str = ["[0, delta]", ...
                "[delta, Delta]", ...
                "[Delta, Delta+delta]", ...
                "[Delta+delta, tpause+Delta+delta]", ...
                "[tpause+Delta+delta, tpause+Delta+2*delta]", ...
                "[tpause+Delta+2*delta, tpause+2*Delta+delta]", ...
                "[tpause+2*Delta+delta, tpause+2*(Delta+delta)]"];
            
            timeprofile_str = repmat("f(t) = %d (constant)", 1, 7);
            timeprofile_str = arrayfun(@(s, i) sprintf(s, i), timeprofile_str, ...
                [1 0 -1 0 1 0 -1]);
            % Remove unused intervals
            if d == D
                timelist(7) = [];
                interval_str(6) = [];
                timeprofile_str(6) = [];
            end
            if p == 0
                timelist(5) = [];
                interval_str(4) = [];
                timeprofile_str(4) = [];
            end
            if d == D
                timelist(3) = [];
                interval_str(2) = [];
                timeprofile_str(2) = [];
            end
        end
        
        function s = string(obj)
            %STRING Convert sequence to string.
            if nargin == 2 && simplified
                s = sprintf("%s_d%g_D%g_tm%g", class(obj), ...
                    obj.delta, obj.Delta, obj.delta+obj.tpause);
            else
                s = sprintf("%s(delta=%g, Delta=%g, tm=%g)", class(obj), ...
                    obj.delta, obj.Delta, obj.delta+obj.tpause);
            end
        end
    end
end

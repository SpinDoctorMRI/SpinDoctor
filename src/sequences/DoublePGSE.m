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
            f = (t < d) ...
                - (D <= t & t < D + d) ...
                + (p + D + d <= t & t < p + D + 2 * d) ...
                - (p + 2 * D + d <= t);
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the time profile from 0 to t.
            %   For the DoublePGSE sequence, an analytical expression is
            %   available.
            d = obj.delta;
            D = obj.Delta;
            p = obj.tpause;
            F = (t < d) .* t ...
                + (d <= t & t < D + d) .* d ...
                - (D <= t & t < D + d) .* (t - D) ...
                + (p + D + d <= t & t < p + D + 2 * d) .* (t - (p + D + d)) ...
                + (p + D + 2 * d <= t) .* d ...
                - (p + 2 * D + d <= t) .* (t - (p + 2 * D + d));
        end
        
        function b = bvalue_no_q(obj)
            %BVALUE_NO_Q Compute the time profile contribution to the b-value.
            %   An analytical expression is available for the DoublePGSE
            %   sequence.
            b = 2 * obj.delta^2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the DoublePGSE sequence.
            t = obj.tpause + 2 * (obj.Delta - obj.delta / 3);
        end
        
        function t = diffusion_time_sta(obj)
            %DIFFUSION_TIME_STA Get STA diffusion time of the DoublePGSE sequence.
            tm = obj.delta + obj.tpause;
            out = (2/35).*((2.*obj.Delta+tm-obj.delta).^(7/2)-2.*(obj.Delta+tm-obj.delta).^(7/2)+...
                         (obj.delta+2.*obj.Delta+tm).^(7/2)+2.*(obj.Delta-obj.delta).^(7/2)+...
                         (tm-obj.delta).^(7/2)-2.*(2.*obj.Delta+tm).^(7/2)-2*...
                         (obj.delta+obj.Delta+tm).^(7/2)+2*(obj.delta+obj.Delta).^(7/2)+4*...
                         (obj.Delta+tm).^(7/2)+(obj.delta+tm).^(7/2)-2*tm.^(7/2)-4*...
                         obj.delta.^(7/2)-4*obj.Delta.^(7/2))./(obj.delta.^2.*(obj.Delta-obj.delta/3));
            t = out.^2;
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
    end
end


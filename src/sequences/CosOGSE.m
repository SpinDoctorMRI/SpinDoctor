classdef CosOGSE < Sequence
    %COSOGSE Oscillating Gradient Spin Echo.
    %   This sequence consists of two cos-pulses of duration delta, separated by
    %   a pause of duration Delta - delta, with nperiod periods per pulse.

    properties
        nperiod
    end

    methods
        function obj = CosOGSE(delta, Delta, nperiod)
            %CosOGSE Construct an instance of this class.
            %   The constructor stores the parameters.
            obj@Sequence(delta, Delta)
            obj.nperiod = nperiod;
        end

        function f = call(obj, t)
            %CALL Call the CosOGSE time profile at time t.
            %   The function is vectorized.
            f = (t < obj.delta) .* cos(2 * pi * obj.nperiod / obj.delta * t) ...
                - (obj.Delta <= t) .* cos(2 * pi * obj.nperiod / obj.delta * (t - obj.Delta));
        end

        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the CosOGSE sequence from 0 to t.
            %   An analytical expression is available.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            F = ((t < d) .* sin(2 * pi * n / d * t) ...
                + (d <= t) .* sin(2 * pi * n) ...
                - (D <= t) .* sin(2 * pi * n / d * (t - D))) ...
                * d / (2 * pi * n);
        end

        function b = bvalue_no_q(obj)
            %BVALUE_NO_Q Compute the time profile contribution to the b-value.
            %   An analytical expression is available for the CosOGSE sequence.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            omega = 2 * pi * n / d;
            b = -(cos(omega * d)^2 * omega * D ...
                - omega * (D + d) ...
                - cos(omega * d) * sin(omega * d) ...
                + 2 * sin(omega * d)) / omega^3;
        end

        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the PGSE sequence.
            t = 3 / 8 * obj.nperiod / obj.delta;
        end

        function [timelist, interval_str, timeprofile_str] = intervals(obj)
            %INTERVALS Get intervals of the sequence.
            %   This function returns a list of important time steps (including
            %   start and stop), a list of strings representing the intervals
            %   between these time steps and a list of strings representing the
            %   behavior of the sequence on each of these intervals.
            timelist = [0, obj.delta, obj.Delta, obj.Delta+obj.delta];
            interval_str = ["[0, delta]", "[delta, Delta]", "[Delta, Delta+delta]"];
            funcname = sprintf("cos(2*pi*%d/delta*t)", obj.nperiod);
            timeprofile_str = ["f(t) = " + funcname, "f(t) = 0 (constant)", ...
                "f(t) = -" + funcname];
            if obj.delta == obj.Delta
                % Remove unused interval
                timelist(3) = [];
                interval_str(2) = [];
                timeprofile_str(2) = [];
            end
        end

        function s = seq2str(obj)
            %SEQ2STR Convert sequence to string.
            s = sprintf("%s(delta=%g, Delta=%g, nperiod=%g)", class(obj), ...
                obj.delta, obj.Delta, obj.nperiod);
        end
    end
end

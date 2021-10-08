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
            if (nperiod == round(nperiod)) && nperiod>0
                obj.nperiod = nperiod;
            else
                error('CosOGSE (Sequence): nperiod must be a positive integer.')
            end
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
        %DIFFUSION_TIME Get diffusion time of the CosOGSE sequence.
            t = 1 / 8 * obj.delta / obj.nperiod;
        end
        
        function t = diffusion_time_sta(obj)
        %DIFFUSION_TIME_STA Get STA diffusion time of the CosOGSE sequence.
        %Matlab's fresnel functions requires Symbolic Math Toolbox. We use
        %a free version of fresnel functions instead.
            n = obj.nperiod;
            S = sin(2*pi*n*obj.Delta./obj.delta);
            C = cos(2*pi*n*obj.Delta./obj.delta);
            FS1 = fresnelS(2*sqrt(n));
            FS2 = fresnelS(2*sqrt(obj.Delta*n./obj.delta));
            FS3 = fresnelS(2*sqrt((obj.Delta+obj.delta)*n./obj.delta));
            FS4 = fresnelS(2*sqrt((obj.Delta-obj.delta)*n./obj.delta));
            FC1 = fresnelC(2*sqrt(n));
            FC2 = fresnelC(2*sqrt(obj.Delta*n./obj.delta));
            FC3 = fresnelC(2*sqrt((obj.Delta+obj.delta)*n./obj.delta));
            FC4 = fresnelC(2*sqrt((obj.Delta-obj.delta)*n./obj.delta));

            out = (3/4)*( (C.*((8.*FC2-4.*FC3-4.*FC4).*obj.Delta./(8*sqrt(obj.delta)) + (-4.*FC3+4.*FC4).*sqrt(obj.delta)/8)...
                + S.*((8.*FS2-4.*FS3-4.*FS4).*obj.Delta./(8*sqrt(obj.delta)) + (-4.*FS3+4.*FS4).*sqrt(obj.delta)/8) + FC1.*sqrt(obj.delta))/sqrt(n)...
                + (C.*sqrt(obj.delta).*(6.*FS2-3.*FS3-3.*FS4)/8 + S.*sqrt(obj.delta).*(-6.*FC2+3.*FC3+3.*FC4)/8 + 3*FS1.*sqrt(obj.delta)/4)./(n^(3/2)*pi) );
            t = out.^2;
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

        function s = string(obj)
            %STRING Convert sequence to string.
            s = sprintf("%s(delta=%g, Delta=%g, nperiod=%g)", class(obj), ...
                obj.delta, obj.Delta, obj.nperiod);
        end
    end
end

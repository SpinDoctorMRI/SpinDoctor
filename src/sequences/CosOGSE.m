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
            f = (0 <= t & t < obj.delta) .* cos(2 * pi * obj.nperiod / obj.delta * t) ...
                - (obj.Delta <= t & t <= obj.echotime) .* cos(2 * pi * obj.nperiod / obj.delta * (t - obj.Delta));
        end

        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the CosOGSE sequence from 0 to t.
            %   An analytical expression is available.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            F = ((0 <= t & t < d) .* sin(2 * pi * n / d * t) ...
                + (d <= t & t <= obj.echotime) .* sin(2 * pi * n) ...
                - (D <= t & t <= obj.echotime) .* sin(2 * pi * n / d * (t - D))) ...
                * d / (2 * pi * n);
        end
        
        function int = integral_F2(obj)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(obj, t)).
            %   An analytical expression is available for the CosOGSE sequence.
            int = obj.delta^3 / (4 * pi^2 * obj.nperiod^2);
        end

        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the CosOGSE sequence.
            t = 1 / 8 * obj.delta / obj.nperiod;
        end
        
        function t = diffusion_time_sta(obj)
            %DIFFUSION_TIME_STA Get STA diffusion time of the CosOGSE sequence.
            %   Matlab's fresnel functions requires Symbolic Math Toolbox. We use
            %   a free version of fresnel functions instead.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            S = sin(2*pi * n * D / d);
            C = cos(2*pi * n * D / d);
            FS1 = fresnelS(2 * sqrt(n));
            FS2 = fresnelS(2 * sqrt(D * n / d));
            FS3 = fresnelS(2 * sqrt((D + d) * n / d));
            FS4 = fresnelS(2 * sqrt((D - d) * n / d));
            FC1 = fresnelC(2 * sqrt(n));
            FC2 = fresnelC(2 * sqrt(D * n / d));
            FC3 = fresnelC(2 * sqrt((D + d) * n / d));
            FC4 = fresnelC(2 * sqrt((D - d) * n / d));
            
            out = 3 / 8 / sqrt(n * d) * ( ...
                + 2 * D * FC2 * C ...
                + d * FC4 * C ...
                - d * FC3 * C ...
                - D * FC3 * C ...
                - D * FC4 * C ...
                + 2 * D * FS2 * S ...
                + d * FS4 * S ...
                - d * FS3 * S ...
                - D * FS3 * S ...
                - D * FS4 * S ...
                + 2 * d * FC1 ...
                ) + ...
                9 * sqrt(d) / (32*pi) / n^(3/2) * ( ...
                + 2 * FS2 * C ...
                - FS3 * C ...
                - FS4 * C ...
                - 2 * FC2 * S ...
                + FC3 * S ...
                + FC4 * S ...
                + 2 * FS1 ...
                );
            t = out^2;
        end
        
        function jn = J(obj, lambda)
            %J Compute the quantity J(lambda) for the sequence
            %   An analytical expression is available for the CosOGSE sequence
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            
            if lambda < 1e-7
                % Use Taylor expansion when lambda is close to 0 
                % to improve numerical stability
                jn = lambda ...
                    - lambda^3 * d^2 * 3 / (4 * (n*pi)^2) ...
                    + lambda^5 * (12*D*n^2*pi^2*d^3 + 15*d^4 - 4*n^2*pi^2*d^4) / (48 * n^4 * pi^4) ...
                    - lambda^6 * D^2 * d^3 / (8 * (n*pi)^2);
            else
                jn = 4 * n^2 * pi^2 * lambda * ( ...
                    - exp(-lambda * (D+d)) * lambda * d ...
                    - exp(-lambda * (D-d)) * lambda * d ...
                    + 2 * exp(-lambda * D) * lambda * d ...
                    + 2 * exp(-lambda * d) * lambda * d ...
                    + 4 * n^2 * pi^2 ...
                    + (lambda * d - 2) * lambda * d ) ...
                    / (4*n^2*pi^2 + lambda^2*d^2)^2;
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

        function s = string(obj, simplified)
            %STRING Convert sequence to string.
            if nargin == 2 && simplified
                s = sprintf("%s_d%g_D%g_n%g", class(obj), ...
                    obj.delta, obj.Delta, obj.nperiod);
            else
                s = sprintf("%s(delta=%g, Delta=%g, nperiod=%g)", class(obj), ...
                    obj.delta, obj.Delta, obj.nperiod);
            end
        end
    end
end

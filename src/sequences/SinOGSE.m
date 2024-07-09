classdef SinOGSE < Sequence
    %SINOGSE Oscillating Gradient Spin Echo.
    %   This sequence consists of two sin-pulses of duration delta, separated by
    %   a pause of duration Delta - delta, with nperiod periods per pulse. 
    %   An initial pause and a fixed echo time can be set through optional
    %   Sequence parameters.
    
    properties
        nperiod
    end
    
    methods
        function obj = SinOGSE(delta, Delta, nperiod,varargin)
            %SINOGSE Construct an instance of this class.
            %   The constructor stores the parameters.
            obj@Sequence(delta, Delta,varargin{:})
            if (nperiod == round(nperiod)) && nperiod>0
                obj.nperiod = nperiod;
            else
                error('SinOGSE (Sequence): nperiod must be a positive integer.')
            end
        end
        
        function f = call(obj, t)
            %CALL Call the SinOGSE time profile at time t.
            %   The function is vectorized.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            t = t-obj.t1;
            f = (0 <= t & t < d) .* sin(2 * pi * n / d * t) ...
                - (D <= t & t <= obj.echotime) .* sin(2 * pi * n / d * (t - D));
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the SinOGSE sequence from 0 to t.
            %   An analytical expression is available.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            t = t-obj.t1;
            F = ((0 <= t & t < d) .* (1 - cos(2 * pi * n / d * t)) ...
                + (d <= t & t <= obj.echotime) .* (1 - cos(2 * pi * n / d * d)) ...
                - (D <= t & t <= obj.echotime) .* (1 - cos(2 * pi * n / d * (t - D)))) ...
                * d /(2 * pi * n);
        end
        
        function int = integral_F2(obj)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(obj, t)).
            %   An analytical expression is available for the SinOGSE sequence.
            int = 3 * obj.delta^3 / (4 * pi^2 * obj.nperiod^2);
        end
        
        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the SinOGSE sequence.
            t = 3 / 8 * obj.delta / obj.nperiod;
        end
        
        function t = diffusion_time_sta(obj)
            %DIFFUSION_TIME_STA Get STA diffusion time of the SinOGSE sequence.
            %Matlab's fresnel functions requires Symbolic Math Toolbox. We use
            %a free version of fresnel functions instead.
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
            out = ( ...
                    + 64*pi*n^(3/2) * D^(3/2) ...
                    + 64*pi*n^(3/2) * d^(3/2) ...
                    + 42*d^(3/2) * FS2 * C ...
                    + 42*d^(3/2) * FS1 ...
                    - 42*d^(3/2) * FC2 * S ...
                    + 32*pi*n^(3/2) * sqrt(D - d) * d ...
                    - 32*pi*n^(3/2) * sqrt(D + d) * D...
                    - 32*pi*n^(3/2) * sqrt(D + d) * d ...
                    - 32*pi*n^(3/2) * sqrt(D - d) * D ...
                    + 24*pi*n*D*sqrt(d) * FC2 * C ...
                    + 24*pi*n*D*sqrt(d) * FS2 * S ...
                    + 24*pi*n*d^(3/2) * FC1...
                    + 21*d^(3/2) * FC3 * S ...
                    + 21*d^(3/2) * FC4 * S ...
                    - 21*d^(3/2) * FS3 * C ...
                    - 21*d^(3/2) * FS4 * C ...
                    + 12*pi*n*d^(3/2) * FS4 * S ...
                    + 12*pi*n*d^(3/2) * FC4 * C...
                    - 12*pi*n*d^(3/2) * FC3 * C ...
                    - 12*pi*n*d^(3/2) * FS3 * S ...
                    - 12*pi*n*D*sqrt(d) * FC4 * C ...
                    - 12*pi*n*D*sqrt(d) * FC3 * C ...
                    - 12*pi*n*D*sqrt(d) * FS3 * S ...
                    - 12*pi*n*D*sqrt(d) * FS4 * S ...
                    ) / (96*pi*n^(3/2) * d);
            t = out^2;
        end

        function jn = J(obj, lambda)
            %J Compute the quantity J(lambda) for the sequence
            %   An analytical expression is available for the SinOGSE sequence
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            if lambda < 1e-7
                % Use Taylor expansion when lambda is close to 0 
                % to improve numerical stability
                jn = lambda ...
                    + lambda^3 * d * (4*n^2*pi^2*d ...
                        - 12*D*n^2*pi^2 - 15*d) ...
                        / (36 * n^2 * pi^2) ...
                    + lambda^4 * d * D^2 / 6 ...
                    + lambda^5 * d * (120*D*n^2*pi^2*d^2 + 4*n^4*pi^4*d^3 ...
                        - 20*D*n^4*pi^4*d^2 + 105*d^3 ...
                        - 40*n^2*pi^2*d^3 - 40*D^3*n^4*pi^4) ...
                        / (720 * n^4 * pi^4) ...
                    + lambda^6 * ((D^4*d^2 + D^2*d^4)/24 ...
                        - (D^2*d^4) / (4*n^2*pi^2)) ...
                        / (3*d);
            else
                jn = 4 * pi^2 * n^2 * ( ...
                    (lambda*d)^3 + 4*n^2*pi^2*( ...
                        + exp(-lambda * (D + d)) ...
                        + exp(-lambda * (D - d)) ...
                        - 2 * exp(-lambda * d) ...
                        - 2 * exp(-lambda * D) ...
                        + lambda * d ...
                        + 2) ) ...
                    / (3 * d * (4*n^2*pi^2 + lambda^2*d^2)^2);
            end
        
        end

        function [timelist, interval_str, timeprofile_str] = intervals(obj)
            %INTERVALS Get intervals of the sequence.
            %   This function returns a list of important time steps (including
            %   start and stop), a list of strings representing the intervals
            %   between these time steps and a list of strings representing the
            %   behavior of the sequence on each of these intervals.
            timelist = [obj.t1 ,obj.t1 + obj.delta, obj.t1 + obj.Delta, obj.t1 + obj.Delta+obj.delta];
            interval_str = ["[t1,t1 + delta]", "[t1 + delta,t1 + Delta]", "[t1+Delta, t1+Delta+delta]"];
            funcname = sprintf("sin(2*pi*%d/delta*t)", obj.nperiod);
            timeprofile_str = ["f(t) = " + funcname, "f(t) = 0 (constant)", ...
                "f(t) = -" + funcname];
            if obj.delta == obj.Delta
                % Remove unused interval
                timelist(3) = [];
                interval_str(2) = [];
                timeprofile_str(2) = [];
            end
            if obj.TE > obj.t1 + obj.Delta+obj.delta
                timelist = [timelist, obj.TE];
                interval_str = [interval_str,"[t1+Delta+delta,TE]"];
                timeprofile_str = [timeprofile_str,"f(t) = 0 (constant)"];
            end
            if obj.t1 > 0
                timelist = [0,timelist];
                interval_str = ["[0,t1]",interval_str];
                timeprofile_str = ["f(t) = 0 (constant)",timeprofile_str];
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
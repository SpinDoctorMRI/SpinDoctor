classdef SinOGSE < Sequence
    %SINOGSE Oscillating Gradient Spin Echo.
    %   This sequence consists of two sin-pulses of duration delta, separated by
    %   a pause of duration Delta - delta, with nperiod periods per pulse.
    
    properties
        nperiod
    end
    
    methods
        function obj = SinOGSE(delta, Delta, nperiod)
            %SINOGSE Construct an instance of this class.
            %   The constructor stores the parameters.
            obj@Sequence(delta, Delta)
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
            f = (t < d) .* sin(2 * pi * n / d * t) ...
                - (D <= t) .* sin(2 * pi * n / d * (t - D));
        end
        
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the SinOGSE sequence from 0 to t.
            %   An analytical expression is available.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            F = ((t < d) .* (1 - cos(2 * pi * n / d * t)) ...
                + (d <= t) .* (1 - cos(2 * pi * n / d * d)) ...
                - (D <= t) .* (1 - cos(2 * pi * n / d * (t - D)))) ...
                * d /(2 * pi * n);
        end
        
        function b = bvalue_no_q(obj)
            %BVALUE_NO_Q Compute the time profile contribution to the b-value.
            %   An analytical expression is available for the SinOGSE sequence.
            d = obj.delta;
            D = obj.Delta;
            n = obj.nperiod;
            omega = 2 * pi * n / d;
            b = (cos(omega*d)^2 * omega * D ...
                - 2 * cos(omega * d) * omega * (D - d) ...
                + omega * (D + d) ...
                - cos(omega * d) * sin(omega * d) ...
                - 2 * sin(omega * d)) / omega^3;
        end
        
        function t = diffusion_time(obj)
        %DIFFUSION_TIME Get diffusion time of the SinOGSE sequence.
            t = 3 / 8 * obj.delta / obj.nperiod;
        end
        
        function t = diffusion_time_sta(obj)
            %DIFFUSION_TIME_STA Get STA diffusion time of the SinOGSE sequence.
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
            out = (-32*pi*n^(3/2).*sqrt(obj.Delta+obj.delta).*obj.delta ...
                + 32*pi*n^(3/2).*sqrt(obj.Delta-obj.delta).*obj.delta + 24*pi*n*obj.delta.^(3/2).*FC1...
                - 32*pi*n^(3/2).*sqrt(obj.Delta-obj.delta).*obj.Delta - 32*pi*n^(3/2).*sqrt(obj.Delta+obj.delta).*obj.Delta...
                - 12.*obj.delta.^(3/2).*FC3.*n.*pi.*C - 12.*obj.delta.^(3/2).*FS3.*n.*pi.*S...
                + 12.*obj.delta.^(3/2).*FS4.*n.*pi.*S + 12.*obj.delta.^(3/2).*FC4.*n.*pi.*C...
                + 42.*obj.delta.^(3/2).*FS2.*C - 21*obj.delta.^(3/2).*FS3.*C...
                - 21.*obj.delta.^(3/2).*FS4.*C + 21.*obj.delta.^(3/2).*FC4.*S...
                - 42.*obj.delta.^(3/2).*FC2.*S + 21.*obj.delta.^(3/2).*FC3.*S...
                + 64.*obj.delta.^(3/2).*pi.*n.^(3/2) + 64.*obj.Delta.^(3/2).*pi.*n^(3/2)...
                + 42.*obj.delta.^(3/2).*FS1 - 12.*FC4.*n.*obj.Delta.*pi.*sqrt(obj.delta).*C...
                + 24.*FC2.*n.*obj.Delta.*pi.*sqrt(obj.delta).*C - 12.*FC3.*n.*obj.Delta.*pi.*sqrt(obj.delta).*C...
                - 12.*FS3.*n.*obj.Delta.*pi.*sqrt(obj.delta).*S + 24.*FS2.*n.*obj.Delta.*pi.*sqrt(obj.delta).*S...
                - 12.*FS4.*n.*obj.Delta.*pi.*sqrt(obj.delta).*S)./(96*pi*n^(3/2).*obj.delta);
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
            funcname = sprintf("sin(2*pi*%d/delta*t)", obj.nperiod);
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


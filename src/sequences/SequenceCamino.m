classdef SequenceCamino  < AbsSequence
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        dt
        K
        g
        echotime
        b_tensor
        name
    end

    methods
        function obj = SequenceCamino(K,dt,g,name)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.K = K;
            obj.dt =dt;  % Conversion from s to micro s
            obj.g = g;
            obj.echotime = K*dt;
            obj.name = name;
        end

        function f = vec_call(obj, t)
        %CALL Call the time profile at time `t`. Vectorised routine
            find_imin = @(s) min([max([0,floor(s./obj.dt)]), obj.K - 1 ]);
            imin = 1+arrayfun(find_imin,t);
            find_imax = @(s) max([min([1+floor(s./obj.dt),obj.K - 1]), 0]);
            imax = 1+arrayfun(find_imax,t);
            f = (obj.g(:,imin)+obj.g(:,imax))/2;
        end
        
        function b = b_value(seq)
            b = seq.b_tensor(1) + 2*seq.b_tensor(2);
        end

        function f = call(obj, t)
        %CALL Call the time profile at time `t`. Vectorised routine to
        %compute magnitude of sequence.
            f = vecnorm(obj.vec_call(t),2,1);
        end

        function t = diffusion_time(obj)
            % Should be specified by the user.
            t = obj.K*obj.dt/2;
        end
        function [timelist, interval_str, timeprofile_str] = intervals(obj)
            timelist = obj.dt*(0:(obj.K-1));
            t_1 = string(timelist(1:end-1)); t_2 = string(timelist(2:end));
            interval_str = sprintf("[%s,%s]!",[t_1;t_2]);
            interval_str = split(interval_str, "!");
            interval_str = interval_str(1:end-1);
            timeprofile_str = repmat("",obj.K - 1,1);
        end
        function F = integral(obj, t)
            %INTEGRAL Compute the integral of the time profile from `0` to
            %`t`. Uses step function approximation
            find_imin = @(s) min([max([0,floor(s./obj.dt)]), obj.K - 1 ]);
            imin= arrayfun(find_imin,t);
            fun = @(i) obj.dt/2 + obj.dt*(0:i-1);
            delta = (t - obj.dt.*(imin)).*(imin>=0 & imin < obj.K);
            s = arrayfun(fun,imin,'UniformOutput',false);
            F = zeros(3,length(t));
            for i = 1: length(t)
                F(:,i) = sum(obj.call(s{i}),2)*obj.dt + delta(i)*(obj.call(t(i)));
            end
        end

        function jn = J(obj, lambda, ninterval)
            %J Compute the quantity J(lambda) for the sequence.
            %   Unless overwritten, it computes a numerical approximation.
            if nargin == 2
                ninterval = 500;
            end
            ntime = ninterval + 1;
            time = linspace(0, obj.echotime, ntime);
            dtime = obj.echotime / ntime;
            fprintf("Use numerical integration to compute the the quantity J(lambda).\n");
            jn = @(t) integral(@(s) exp(lambda * (s - t)) .* obj.call(s), 0, t);
            jn = lambda * obj.integral(time) .* arrayfun(jn, time);
            jn = dtime * trapz(jn) / obj.integral_F2;
        end
        function int = integral_F2(obj)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(obj, t)).
            %   Unless overwritten, it computes a numerical approximation.
            dt = obj.dt;
            t = dt*(0:obj.K);
            F = obj.integral(t);
            int = 0;
            for i = 1:obj.K
                g = obj.vec_call(mean(t(i) + t(i+1)));
                int = int + norm(F(:,i))^2*dt + dot(F(:,i),g)*dt^2 + (norm(g)^2)*(dt^3)/3;
            end
        end

        function plot(obj)
            % figure;
            hold on;
            t = obj.dt*(0:obj.K -1);
            g_x = obj.g(1,:); g_y = obj.g(2,:); g_z = obj.g(3,:);
            plot(t,g_x,'DisplayName','$g_x$');
            plot(t,g_y,'DisplayName','$g_y$');
            plot(t,g_z,'DisplayName','$g_z$');
            xlabel("t / \mu" + "s")
            ylabel("g / mT/m");
            legend(Interpreter='latex');
            hold off
        end
        function s = string(obj,vargin)
            s = obj.name;
        end

    end
end
classdef SequenceCamino  < AbsSequence
    %SEQUENCECAMINO A direction varying sequence.
    %   This is is obtained from a camino .scheme file.

    properties
        dt
        K
        g
        echotime
        b_tensor
        name
    end

    methods
        function seq = SequenceCamino(K,dt,g,name)
            %SequenceCamino Construct an instance of this class
            seq.K = K;
            seq.dt =dt;  
            seq.g = g;
            seq.echotime = K*dt;
            seq.name = name;
        end

        function f = vec_call(seq, t)
        %VEC_CALL Call the time profile at time `t`. Vectorised routine
            find_imin = @(s) min([max([0,floor(s./seq.dt)]), seq.K - 1 ]);
            imin = 1+arrayfun(find_imin,t);
            find_imax = @(s) max([min([1+floor(s./seq.dt),seq.K - 1]), 0]);
            imax = 1+arrayfun(find_imax,t);
            f = (seq.g(:,imin)+seq.g(:,imax))/2;
        end
        
        function b = b_value(seq)
            b = seq.b_tensor(1) + 2*seq.b_tensor(2);
        end

        function f = call(seq, t)
        %CALL compute the magnitude of the sequence at time `t`. Vectorised routine to
            f = vecnorm(seq.vec_call(t),2,1);
        end

        function t = diffusion_time(seq)
            %DIFFUSION_TIME note not necessarily well-defined for a general
            %sequence
            t = seq.K*seq.dt/2;
        end
        function [timelist, interval_str, timeprofile_str] = intervals(seq)
            %INTERVALS Extract time intervals for sequence
            timelist = seq.dt*(0:(seq.K-1));
            t_1 = string(timelist(1:end-1)); t_2 = string(timelist(2:end));
            interval_str = sprintf("[%s,%s]!",[t_1;t_2]);
            interval_str = split(interval_str, "!");
            interval_str = interval_str(1:end-1);
            timeprofile_str = repmat("",seq.K - 1,1);
        end
        function F = integral(seq, t)
            %INTEGRAL Compute the integral of the time profile from `0` to
            %`t`. Uses step function approximation
            find_imin = @(s) min([max([0,floor(s./seq.dt)]), seq.K - 1 ]);
            imin= arrayfun(find_imin,t);
            fun = @(i) seq.dt/2 + seq.dt*(0:i-1);
            delta = (t - seq.dt.*(imin)).*(imin>=0 & imin < seq.K);
            s = arrayfun(fun,imin,'UniformOutput',false);
            F = zeros(3,length(t));
            for i = 1: length(t)
                F(:,i) = sum(seq.call(s{i}),2)*seq.dt + delta(i)*(seq.call(t(i)));
            end
        end

        function jn = J(seq, lambda, ninterval)
            %J Compute the quantity J(lambda) for the sequence.
            %   Unless overwritten, it computes a numerical approximation.
            if nargin == 2
                ninterval = 500;
            end
            ntime = ninterval + 1;
            time = linspace(0, seq.echotime, ntime);
            dtime = seq.echotime / ntime;
            fprintf("Use numerical integration to compute the the quantity J(lambda).\n");
            jn = @(t) integral(@(s) exp(lambda * (s - t)) .* seq.call(s), 0, t);
            jn = lambda * seq.integral(time) .* arrayfun(jn, time);
            jn = dtime * trapz(jn) / seq.integral_F2;
        end
        function int = integral_F2(seq)
            %INTEGRAL_F2 Compute the temporal integration of F^2 (F = integral(seq, t)).
            %   Unless overwritten, it computes a numerical approximation.
            dt = seq.dt;
            t = dt*(0:seq.K);
            F = seq.integral(t);
            int = 0;
            for i = 1:seq.K
                g = seq.vec_call(mean(t(i) + t(i+1)));
                int = int + norm(F(:,i))^2*dt + dot(F(:,i),g)*dt^2 + (norm(g)^2)*(dt^3)/3;
            end
        end

        function B = get_B_tensor(seq,gamma)
            G = 0.5*(seq.g(:,1:end-1) +seq.g(:,2:end));
            delta_t = seq.dt;
            G = G';
            q = gamma*delta_t*cumsum(G); 
            q  = [zeros(1,3); q];
            B =zeros(3,3);
            for i = 1:seq.K-1
                for k = 1:3
                    for l = 1:k
                        Q = gamma^2*delta_t^3*G(i,k)*G(i,l)/3;
                        B(k,l) = B(k,l) + Q;
                        C = delta_t*q(i,k)*q(i,l);
                        L = gamma*delta_t^2*(q(i,k)*G(i,l)+q(i,l)*G(i,k))/2;
                        B(k,l) = B(k,l) + C + L;
                    end
                end
            end
            B = B + B'.*(1-eye(3));
            
        end


        function plot(seq)
            %PLOT plot the sequence.
            hold on;
            t = seq.dt*(0:seq.K -1);
            g_x = seq.g(1,:); g_y = seq.g(2,:); g_z = seq.g(3,:);
            plot(t,g_x,'DisplayName','$g_x$');
            plot(t,g_y,'DisplayName','$g_y$');
            plot(t,g_z,'DisplayName','$g_z$');
            xlabel("t / \mu" + "s")
            ylabel("g / mT/m");
            legend(Interpreter='latex');
            hold off
        end
        function s = string(seq,vargin)
            s = seq.name;
        end

    end
end
function smpl = ignpoiScalar(mu)
% IGNPOI generates a Poisson random deviate.
%
%  Discussion:
%
%    This procedure generates a single random deviate from a Poisson
%    distribution with given mean.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    01 April 2013
%
%  Author:
%
%    Original FORTRAN77 version by Barry Brown, James Lovato.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Joachim Ahrens, Ulrich Dieter,
%    Computer Generation of Poisson Deviates
%    From Modified Normal Distributions,
%    ACM Transactions on Mathematical Software,
%    Volume 8, Number 2, June 1982, pages 163-179.
%
%  Parameters:
%
%    Input, real MU, the mean of the Poisson distribution
%    from which a random deviate is to be generated.
%
%    Output, integer smpl, a random deviate from
%    the distribution.

%-------------------------------------------------------------------------%
% Modified: 08/06/14
%   -the point of this version is that it operates on scalars only and has
%   no subfunctions and hence can be used by arrayfun and, therefore, a GPU
%   by JGM
%-------------------------------------------------------------------------%


smpl = NaN;
a0 = -0.5;
a1 =  0.3333333;
a2 = -0.2500068;
a3 =  0.2000118;
a4 = -0.1661269;
a5 =  0.1421878;
a6 = -0.1384794;
a7 =  0.1250060;


% Start new table and calculate P0.
if (mu < 10.0)

    
    % getSamplesForSmallMu
    p = exp(-mu);
    q = p;
    p0 = p;
    
    % getPoisCumProb ("Uniform sample for inversion method.")
    while true
        
        u = rand;
        smpl = 0;
        
        if u <= p0
            return;
        end
        
        %  ("Creation of new Poisson probabilities.")
        for k = 1:35
            p = p*mu/k;
            q = q + p;
            if u <= q
                smpl = k;
                return
            end
        end
        
    end
    
else
    
    % getSamplesForBigMu
    s = sqrt(mu);
    g = mu + s*randn;
    
    if 0.0 <= g
        
        % getSamplesForPosNormRnds
        l = floor(mu - 1.1484);
        smpl = floor(g);
        
        %  ("Immediate acceptance if large enough.")
        if l <= smpl, return; end
        
        
        % getSamplesForMediateAcceptance ("Squeeze acceptance.")
        difmuk = mu - smpl;
        u = rand;
        d = 6.0*mu*mu;
        
        if difmuk*difmuk*difmuk <= d*u, return; end
        
        
        % getSamplesViaPxPyAcceptance ("Preparation for steps P and Q.")
        omega = 0.3989423/s;
        b1 = 0.04166667/mu;
        b2 = 0.3*b1*b1;
        c3 = 0.1428571*b1*b2;
        c2 = b2 - 15.0*c3;
        c1 = b1 - 6.0*b2 + 45.0*c3;
        c0 = 1.0 - b1 + 3.0*b2 - 15.0*c3;
        c = 0.1069/mu;
        
        if smpl < 10
            
            % getPxPyForSmallSamples
            px = -mu;
            if smpl==0, fct = 1.0;
            elseif smpl==1, fct = 1.0;
            elseif smpl==2, fct = 2.0;
            elseif smpl==3, fct = 6.0;
            elseif smpl==4, fct = 24.0;
            elseif smpl==5, fct = 120.0;
            elseif smpl==6, fct = 720.0;
            elseif smpl==7, fct = 5040.0;
            elseif smpl==8, fct = 40320.0;
            elseif smpl==9, fct = 362880.0;
            else fct=NaN; %%% to satisfy compiler
            end
            py = mu^smpl/fct;   %%% fact(smpl+1);
            
        else
            
            % getPxPyForLargeSamples
            del = 0.8333333e-01/smpl;
            del = del - 4.8*del*del*del;
            v = difmuk/smpl;
            
            if 0.25 < abs(v)
                % getPxForBigAbsv
                px = smpl*log(1.0 + v) - difmuk - del;
            else
                % getPxForSmallAbsv
                px = smpl*v*v*(((((((a7*v + a6 )*v + a5)*v + a4)*v +...
                    a3)*v + a2)*v + a1)*v + a0) - del;
            end
            
            py = 0.3989423/sqrt(smpl);
            
        end
        
        x = (0.5 - difmuk)/s;
        xx = x*x;
        fx = -0.5*xx;
        fy = omega*(((c3*xx + c2)*xx + c1)*xx + c0);
        
        if (fy - u*fy <= py*exp(px - fx)), return; end
        
        
    end
    
    
    % getSamplesViaExpoSamples ("Exponential sample.")
    while true
        
        e = -log(rand);
        u = 2.0*rand - 1.0;
        t = 1.8 + e*sign(u);
        
        if t <= -0.6744, continue; end
        
        % getNewSamplesViaEtc
        smpl = floor(mu + s*t);
        difmuk = mu - smpl;
        
        % getSamplesViaPxPyAcceptance ("Calculation of PX, PY, FX, FY.")
        omega = 0.3989423/s;
        b1 = 0.04166667/mu;
        b2 = 0.3*b1*b1;
        c3 = 0.1428571*b1*b2;
        c2 = b2 - 15.0*c3;
        c1 = b1 - 6.0*b2 + 45.0*c3;
        c0 = 1.0 - b1 + 3.0*b2 - 15.0*c3;
        c = 0.1069/mu;
        
        if smpl < 10
            
            % getPxPyForSmallSamples
            px = -mu;
            if smpl==0, fct = 1.0;
            elseif smpl==1, fct = 1.0;
            elseif smpl==2, fct = 2.0;
            elseif smpl==3, fct = 6.0;
            elseif smpl==4, fct = 24.0;
            elseif smpl==5, fct = 120.0;
            elseif smpl==6, fct = 720.0;
            elseif smpl==7, fct = 5040.0;
            elseif smpl==8, fct = 40320.0;
            elseif smpl==9, fct = 362880.0;
            else fct = NaN; %%% to satisfy compiler
            end
            py = mu ^ smpl / fct;  % fact(smpl+1);
            
        else
            
            % getPxPyForLargeSamples
            del = 0.8333333e-01/smpl;
            del = del - 4.8*del*del*del;
            v = difmuk/smpl;
            
            if 0.25 < abs(v)
                % getPxForBigAbsv
                px = smpl*log(1.0 + v) - difmuk - del;
            else
                % getPxForSmallAbsv
                px = smpl*v*v*(((((((a7*v + a6 )*v + a5)*v + a4)*v +...
                    a3)*v + a2)*v + a1)*v + a0) - del;
            end
            
            py = 0.3989423/sqrt(smpl);
            
        end
        
        
        x = (0.5 - difmuk)/s;
        xx = x*x;
        fx = -0.5*xx;
        fy = omega*(((c3*xx + c2)*xx + c1)*xx + c0);
        
        if (c*abs(u) <= py*exp(px + e) - fy*exp(fx + e)), return; end
        
    end
    
end

end

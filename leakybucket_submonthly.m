%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% LEAKY BUCKET WITH SUBSTEPPING %%%%%%%%%%%%%%%%%%%%%%%%
function [M,potEv,ndl,cdays] = leakybucket_submonthly(syear,eyear,phi,T,P,E,...
    Mmax,Mmin,alph,m_th,mu_th,rootd,M0,pet_model)
% leackybucket_submonthly.m - Simulate soil moisture; substeps within monthly timesteps
% to better capture nonlinearities and improve moisture estimates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [M,potEv,ndl,cdays] = leakybucket_submonthly(syear,eyear,phi,T,P,...
%                   Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
%    outputs simulated soil moisture and potential evapotranspiration.
%
% Inputs:
%   syear = start year of simulation.
%   eyear = end year of simulation.
%   phi = latitude of site (in degrees N)
%   T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
%   P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
%   E = (12 x Nyrs) matrix of ordered accumulated monthly potential evapotranspiration (in mm) if  using Penman-Monteith, Priestley-Taylor, or Hargreaves (overwritten with Thornthwaite otherwise)
%   Mmax = scalar maximum soil moisture held by the soil (in v/v)
%   Mmin = scalar minimum soil moisture (for error-catching) (in v/v)
%   alph = scalar runoff parameter 1 (in inverse months)
%   m_th = scalar runoff parameter 3 (unitless)
%   mu_th = scalar runoff parameter 2 (unitless)
%   rootd = scalar root/"bucket" depth (in mm)
%   M0 = initial value for previous month's soil moisture at t = 1 (in v/v)
%
% Outputs:
%   M = soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
%   potEv = potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
%
% SETW+ N. Graham and K. Georgakakos 2011

% modified by Nick G. and K. Georgakakos - to sub-step the monthly steps. Also this version has added
% soil moisture initial conditions for restarts, or spin-up.  Hands back monthly soil moisture
% and summer soil moisture as well - see varargout.  Nick G. 2011/06
% 
% Modified by MPD, 2019-2021: allow input of PET from other methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iyear = syear:eyear;
nyrs = length(iyear);
% Storage for growth response output variables (size [12 x Nyears]):
M   = NaN(12,nyrs);
potEv = NaN(12,nyrs);

% ADDED BY NICK
if(M0 < 0.)
    M0=200/rootd;
end

% Compute normalized daylength (neglecting small difference in calculation for leap-years)
latr = phi*pi/180;  % change to radians
ndays = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdays = cumsum(ndays);
sd = asin(sin(pi*23.5/180) * sin(pi * (((1:365) - 80)/180)))';   % solar declination
y = -tan(ones(365,1).* latr) .* tan(sd);
if ~isempty(find(y>=1,1))
    y(y>=1) = 1;
end
if ~isempty(find(y<=-1,1))
    y(y<=-1) = -1;
end
hdl = acos(y);
dtsi = (hdl.* sin(ones(365,1).*latr).*sin(sd))+(cos(ones(365,1).*latr).*cos(sd).*sin(hdl));
ndl=dtsi./max(dtsi); % normalized day length

% calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
jday = cdays(1:12) +.5*ndays(2:13);
m_star = 1-tand(phi)*tand(23.439*cos(jday*pi/182.625));
mmm = NaN*ones(1,12);
for mo = 1:12
    if m_star(mo) < 0
        mmm(mo) = 0;
    elseif m_star(mo) >0 && m_star(mo) < 2
        mmm(mo) = m_star(mo);
    elseif m_star(mo) > 2
        mmm(mo) = 2;
    end
end
nhrs = 24*acosd(1-mmm)/180; % the number of hours in the day in the middle of the month
L = (ndays(2:13)/30).*(nhrs/12); % mean daylength in each month.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -- year cycle -- %%%%
% syear = start (first) year of simulation
% eyear = end (last) year of simulation
% cyear = year the model is currently working on
% iyear = index of simulation year

for cyear=1:length(iyear)      % begin cycling over years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 1:12  % begin cycling over months in a year
        %%%%% Compute potential evapotranspiration for current month after Thornthwaite:
        if strcmp(pet_model,'Th')
            if T(t,cyear) < 0
                Ep = 0;
            elseif T(t,cyear)>=0 && T(t,cyear) < 26.5;
                istar = (T(:,cyear)/5); istar(istar<0) = 0;
                I = sum(istar.^1.514);
                a = (6.75e-7)*I^3 - (7.71e-5)*I^2 + (1.79e-2)*I + .49;
                Ep = 16*L(t)*(10*T(t,cyear)/I)^a;
            elseif T(t,cyear) >= 26.5;
                Ep = -415.85 + 32.25*T(t,cyear) - .43* T(t,cyear)^2;
            end
            potEv(t,cyear) = Ep;
        else
            Ep = E(t,cyear); % if not using Thornthwaite, use the PET provided in the function input
            potEv(t,cyear) = Ep;
        end
        %%%%% Now calculate soil moisture according to the CPC Leaky Bucket model
        %%%%% (see J. Huang et al, 1996). Set n-steps according to 2 mm increments
        %%%%% have to update alpha and Ep as well - 2 mm increments came from
        %%%%% testing by K. Georgakakos, but one could use 5 or more, with less "accurate" results.
        %%%%% Stepping is necessary because the parametization is linearized around init condition.
        %%%%%%%%%%%%%%%%%
        dp = 2.0; % mm of precip per increment
        nstep = floor(P(t,cyear)/dp)+1; % number of sub-monthly substeps
        Pinc = P(t,cyear)/nstep; % precip per substep
        alphinc = alph/nstep; % runoff rate per substep time interval
        Epinc = Ep/nstep; % potential evapotrans per substep.
        %%%%%%%%%%%%%%%%%
        % handling for sm_init
        if (t > 1)
            M0=M(t-1,cyear);
        elseif (t == 1 && cyear > 1)
            M0=M(12,cyear-1);
        end
        sm0=M0;
        
        for istep=1:nstep
            % evapotranspiration:
            Etrans = Epinc*sm0*rootd/(Mmax*rootd);
            % groundwater loss via percolation:
            G = mu_th*alphinc/(1+mu_th)*sm0*rootd;
            % runoff; contributions from surface flow (1st term) and subsurface (2nd term)
            R = Pinc*(sm0*rootd/(Mmax*rootd))^m_th + (alphinc/(1+mu_th))*sm0*rootd;
            dWdt = Pinc - Etrans - R - G;
            sm1 = sm0 + dWdt/rootd;
            %
            sm0=max(sm1,Mmin);
            sm0=min(sm0,Mmax);
        end
        M(t,cyear) = sm0;
        % error-catching:
        if M(t,cyear) <= Mmin; M(t,cyear) = Mmin; end;
        if M(t,cyear) >= Mmax; M(t,cyear) = Mmax; end;
        if isnan(M(t,cyear))==1; M(t,cyear) = Mmin; end;
    end % end month (t) cycle
end % end year cycle

end

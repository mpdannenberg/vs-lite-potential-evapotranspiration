function [trw,varargout] = VSLite_v3(syear,eyear,phi,T1,T2,M1,M2,K1,K2,Tmin,Tmax,Tdmean,P,coast,elev,varargin)
% VSLite_v2_3.m - Simulate tree ring width index given monthly climate inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Usage:
%    trw = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P)
%    gives just simulated tree ring as ouput.
%
%   [trw,gT,gM,gD,gE,M] = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P))
%    also includes growth response to temperature, growth response to soil
%    moisture, growth response to VPD, scaled insolation index, and soil moisture estimate in outputs.
%
% Basic Inputs:
%   syear = start year of simulation.
%   eyear = end year of simulation.
%   phi = latitude of site (in degrees N)
%   T1 = scalar temperature threshold below which temp. growth response is zero (in deg. C)
%   T2 = scalar temperature threshold above which temp. growth response is one (in deg. C)
%   M1 = scalar soil moisture threshold below which moist. growth response is zero (in v/v)
%   M2 = scalar soil moisture threshold above which moist. growth response is one (in v/v)
%   K1 = scalar VPD threshold below which VPD growth response is one (in kPa)
%   K2 = scalar VPD threshold above which VPD growth response is zero (in kPa)
%     (Note that optimal growth response parameters T1, T2, M1, M2 may be estimated
%      using code estimate_vslite_params_v2_3.m also freely available at
%      the NOAA NCDC Paleoclimatology software library.)
%   Tmin = (12 x Nyrs) matrix of ordered mean minimum monthly temperatures (in degees C)
%   Tmax = (12 x Nyrs) matrix of ordered mean maximum monthly temperatures (in degees C)
%   Tdmean = (12 x Nyrs) matrix of ordered mean monthly dewpoint temperature (in degEes C)
%   P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
%   coast = coastal indicator (1=coastline, 0=interior) for P-M estimation
%   elev = elevation of site (in m)
%
% Advanced Inputs (must be specified as property/value pairs):
%   'lbparams':  Parameters of the Leaky Bucket model of soil moisture.
%                These may be specified in an 8 x 1 vector in the following
%                order (otherwise the default values are read in):
%                   Mmax: scalar maximum soil moisture content (in v/v),
%                     default value is 0.76
%                   Mmin: scalar minimum soil moisture (in v/v), default
%                     value is 0.01
%                   alph: scalar runoff parameter 1 (in inverse months),
%                     default value is 0.093
%                   m_th: scalar runoff parameter 3 (unitless), default
%                     value is 4.886
%                   mu_th: scalar runoff parameter 2 (unitless), default
%                     value is 5.80
%                   rootd: scalar root/"bucket" depth (in mm), default
%                     value is 1000
%                   M0: initial value for previous month's soil moisture at
%                     t = 1 (in v/v), default value is 0.2
%                   substep: logical 1 or 0; perform monthly substepping in
%                     leaky bucket (1) or not (0)? Default value is 0.
%   'intwindow': Integration window. Which months' growth responses should
%                be intregrated to compute the annual ring-width index?
%                Specified as a 2 x 1 vector of integer values. Both
%                elements are given in integer number of months since January
%                (July) 1st of the current year in the Northern (Southern)
%                hemisphere, and specify the beginning and end of the integration
%                window, respectively. Defaults is [1 ; 12] (eg. integrate
%                response to climate over the corresponding calendar year,
%                assuming location is in the northern hemisphere).
%   'hydroclim': Value is a single character either taking value ['P'] or ['M'].
%                If ['M'], then 9th input is interpreted as an estimate of
%                soil moisture content (in v/v) rather than as precipitation.
%                Model default is to read in precipitation and use the CPC's
%                Leaky Bucket model of hydrology to estimate soil moisture,
%                however if soil moisture observations or a more sophisticated
%                estimate of moisture accounting for snow-related processes
%                is available, then using these data directly are recommended
%                (and will also speed up code).
%   'pet_model': Potential evapotranspiration model to use for leaky bucket
%                model. Options are: 'Th'=Thornthwaite (DEFAULT),
%                'PM'=Penman-Monteith, 'PT'=Priestley-Taylor, and
%                'Hg'=Hargreaves.
%   'vpd_scalar':logical 1 or 0; include VPD scalar in addition to soil
%                moisture scalar? no = 0 (default), yes = 1.
%   'gT_0':User-input temperature scalar (12 x Nyrs) that overrides the gT calculation.
%   'gM_0':User-input soil moisture scalar (12 x Nyrs) that overrides the gM calculation.
%   'gE_0':User-input energy scalar (12 x 1) that overrides the gE calculation.
%   'gD_0':User-input VPD scalar (12 x Nyrs) that overrides the gD calculation.
%
% For more detailed documentation, see:
% 1) Tolwinski-Ward et al., An efficient forward model of the climate
% controls on interannual variation in tree-ring width, Climate Dynamics (2011)
% DOI: 10.1007/s00382-010-0945-5
%
% 2) Tolwinski-Ward et al., Erratum to: An efficient forward model of the climate
% controls on interannual variation in tree-ring width, Climate Dynamics (2011)
% DOI: 10.1007/s00382-011-1062-9
%
% 3) Tolwinski-Ward et al., Bayesian parameter estimation and
% interpretation for an intermediate model of tree-ring width, Clim. Past
% (2013), DOI: 10.5194/cp-9-1-2013
%
% 4) Documentation available with the model at http://www.ncdc.noaa.gov/paleo/softlib/softlib.html
% 
% 5) Dannenberg, M. P. (2021), Modeling tree radial growth in a warming 
% climate: Where, when, and how much do potential evapotranspiration models
% matter?, Environmental Research Letters.
%
% Revision History
% v0.1 - Original coding at monthly timestep from full daily timestep model (SETW, 4/09)
% v1.0 - Changed soil moisture module to the CPC Leaky Bucket model (SETW, 5/09)
% v1.1 - No upper parametric bounds for gT, gW as in full model; no density module (SETW, 9/09)
% v1.2 - Added adjustable integration window parameters (SETW, 1/10)
% v2.0 - Minor debugging for Octave compatibility, final version for publication (SETW, 10/10)
% v2.1 - Error in evapotranspiration calculation corrected (SETW, 7/11)
% v2.2 - Add switch to allow for monthly sub-stepping in soil moisture computation (SETW, N.Graham, K.Georgakaos, 9/11)
% v2.3 - Add switch to allow moisture M to be given as input rather than estimated
%        from T and P; add variable input options and improve commenting (SETW, 7/13)
% v3.0 - Add additional inputs and switch for Penman-Monteith and VPD
%        response function (MPD, 3/19)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iyear = syear:eyear;
nyrs = length(syear:eyear);
T = (Tmin + Tmax) / 2; % estimate monthly mean temperature from minimum and maximum temperatures
[E, ~, D] = fao_pm(Tmax, Tmin, Tdmean, phi, elev, coast, iyear); % Estimate potential evapotranspiration (E) and VPD (D) from climate inputs using FAO Penman-Monteith method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in advanced inputs if user-specified; else read in parameter defaults:
if nargin > 15 
    %%%% First fill parameter values in with defaults: %%%%%
    % Parameters of the Leaky Bucket model:
    Mmax = 0.76;
    Mmin = 0.01;
    alph = 0.093;
    m_th = 4.886;
    mu_th = 5.80;
    rootd = 1000;
    M0 = 0.2;
    substep = 0;
    % Integration window parameters:
    I_0 = 1;
    I_f = 12;
    % Hydroclimate variable:
    hydroclim = 'P';
    % PET model
    pet_model = 'Th';
    % VPD scalar
    vpd_scalar = 0;
    for i = 1:length(varargin)/2
        namein = varargin{2*(i-1)+1};
        valin = varargin{2*i};
        switch namein
            case 'lbparams'
                Mmax = valin(1);
                Mmin = valin(2);
                alph = valin(3);
                m_th = valin(4);
                mu_th = valin(5);
                rootd = valin(6);
                M0 = valin(7);
                substep = valin(8);
            case 'intwindow'
                I_0 = valin(1);
                I_f = valin(2);
            case 'hydroclim'
                hydroclim = valin;
            case 'pet_model'
                pet_model = valin;
            case 'vpd_scalar'
                vpd_scalar = valin;
            case 'gT_0'
                gT_0 = valin;
            case 'gM_0'
                gM_0 = valin;
            case 'gE_0'
                gE_0 = valin;
            case 'gD_0'
                gD_0 = valin;
        end
    end
else % otherwise, read in defaults:
    % Parameters of the Leaky Bucket model:
    Mmax = 0.76;
    Mmin = 0.01;
    alph = 0.093;
    m_th = 4.886;
    mu_th = 5.80;
    rootd = 1000;
    M0 = 0.2;
    substep = 0;
    % Integration window parameters:
    I_0 = 1;
    I_f = 12;
    % Hydroclimate variable:
    hydroclim = 'P';
    % PET model
    pet_model = 'Th';
    % VPD scalar
    vpd_scalar = 0;
end
%%% Pre-allocate storage for outputs: %%%%
Gr = NaN(12,nyrs); % Growth index
gT = NaN(12,nyrs); % Temperature response
gM = NaN(12,nyrs); % Soil moisture response
gD = NaN(12,nyrs); % VPD response
%%% Optional switch for PET model (if pet_model = 'Th', PET is overwritten
%%% with Thornthwaite method in the leakybucket functions)
if strcmp(pet_model,'Hg')
    E = hargreaves(Tmax, Tmin, phi, iyear); % Overwrite P-M with Hargreaves if desired
elseif strcmp(pet_model,'PT')
    E = priestley_taylor(Tmax, Tmin, Tdmean, phi, elev, coast, iyear); % Overwrite P-M with Priestley-Taylor if desired
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in or estimate soil moisture:
if strcmp(hydroclim,'M')
    % Read in soil moisture:
    M = P;
else
    % Compute soil moisture:
    if substep == 1
        [M, potEv] = leakybucket_submonthly(syear,eyear,phi,T,P,E,Mmax,Mmin,alph,m_th,mu_th,rootd,M0,pet_model);
    elseif substep == 0
        [M, potEv] = leakybucket_monthly(syear,eyear,phi,T,P,E,Mmax,Mmin,alph,m_th,mu_th,rootd,M0,pet_model);
    elseif substep ~=1 && substep ~= 0
        disp('''substep'' must either be set to 1 or 0.');
        return
    end
end
% Compute gE, the scaled monthly proxy for insolation:
gE = Compute_gE(phi);
if exist('gE_0','var')
    gE = gE_0;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now compute growth responses to climate, and simulate proxy:
%%%%%%%%%%%%%%%%
% syear = start (first) year of simulation
% eyear = end (last) year of simulation
% cyear = year the model is currently working on
% iyear = index of simulation year
% Compute monthly growth response to T & M, and overall growth response G:
for cyear=1:length(iyear)      % begin cycling over years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 1:12  % begin cycling over months in a year
        %%% Calculate Growth Response functions gT(t) and gM(t)
        % First, temperature growth response:
        x = T(t,cyear);
        if (x < T1)
            gT(t,cyear) = 0;
        elseif (x >= T1) & (x <= T2)
            gT(t,cyear) = (x - T1)/(T2 - T1);
        elseif (x >= T2)
            gT(t,cyear) = 1;
        end
        if exist('gT_0','var')
            gT(t,cyear) = gT_0(t,cyear);
        end
        % Next, Soil moisture growth response:
        x = M(t,cyear);
        if (x < M1)
            gM(t,cyear) = 0;
        elseif (x >= M1) & (x <= M2)
            gM(t,cyear) = (x - M1)/(M2 - M1);
        elseif (x >= M2)
            gM(t,cyear) = 1;
        end
        if exist('gM_0','var')
            gM(t,cyear) = gM_0(t,cyear);
        end
        % Next, VPD growth response:
        if vpd_scalar==1
            x = D(t,cyear);
            if (x < K1)
                gD(t,cyear) = 1;
            elseif (x >= K1) & (x <= K2)
                gD(t,cyear) = (K2 - x)/(K2 - K1);
            elseif (x >= K2)
                gD(t,cyear) = 0;
            end
        end
        if exist('gD_0','var')
            gD(t,cyear) = gD_0(t,cyear);
        end
    end % end month (t) cycle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute overall growth rate:
    if vpd_scalar==0
        Gr(:,cyear) = gE.*min(gT(:,cyear),gM(:,cyear));
    elseif vpd_scalar==1
        Gr(:,cyear) = gE.*min(gT(:,cyear),gM(:,cyear).*gD(:,cyear)); % Multiply VPD and M scalars (under assumption that they reinforce each other)
    end
end % end year cycle
%%%%%%%%%%%%%% Compute proxy quantity from growth responses %%%%%%%%%%%%%%%
width = NaN*ones(length(syear:eyear),1);
if phi>0 % if site is in the Northern Hemisphere:
    if I_0<0 % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        % use average of growth data across modeled years to estimate first year's growth due
        % to previous year:
        width(1) = sum(Gr(1:endmo,1)) + sum(mean(Gr(startmo:12,:),2));
        for cyear = 2:length(syear:eyear)
            width(cyear) = sum(Gr(startmo:12,cyear-1)) + sum(Gr(1:endmo,cyear));
        end
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        for cyear = 1:length(syear:eyear)
            width(cyear) = sum(Gr(startmo:endmo,cyear));
        end
    end
elseif phi<0 % if site is in the Southern Hemisphere:
    % (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo = 7+I_0; % (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo = I_f-6; % (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
    for cyear = 1:length(syear:eyear)-1
        width(cyear) = sum(Gr(startmo:12,cyear)) + sum(Gr(1:endmo,cyear+1));
    end
    % use average of growth data across modeled years to estimate last year's growth due
    % to the next year:
    width(length(syear:eyear)) = sum(Gr(startmo:12,length(syear:eyear)))+...
        sum(mean(Gr(1:endmo,:),2));
end
%
trw = ((width-mean(width))/std(width))'; % proxy series is standardized width.
%
if nargout >=1
    varargout(1) = {gT};
    varargout(2) = {gM};
    varargout(3) = {gD};
    varargout(4) = {gE};
    varargout(5) = {M};
    varargout(6) = {potEv};
    varargout(7) = {width};
    varargout(8) = {mean(width)};
    varargout(9) = {std(width)};
end
%
end


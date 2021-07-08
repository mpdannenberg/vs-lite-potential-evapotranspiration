%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SCALED DAYLENGTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gE] = Compute_gE(phi)
% Just what it sounds like... computes just gE from latitude a la VS-Lite,
% but without all the other stuff.
%
% Usage: gE = Compute_gE(phi)
%
% SETW 3/8/13

%
gE = NaN(12,1);
%
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
%nhrs = 24*acosd(1-mmm)/180; % the number of hours in the day in the middle of the month
%L = (ndays(2:13)/30).*(nhrs/12);
%
for t = 1:12
    gE(t) = mean(ndl(cdays(t)+1:cdays(t+1),1));
end
%%%%%%%%%%%%%%%
end

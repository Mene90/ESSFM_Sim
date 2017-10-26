function Hf=myfilter(ftype,f,bw,ord)

% CONSTANTS

      r4p2r2=2.61312592975275;      % =sqrt(2*(2+sqrt(2))); % used in Butterworth 4th order
      b1=3.86370330515627315;
      b2=7.4641016151377546;
      b3=9.1416201726856413;
      b4=b2;
      b5=b1;
      Bb=0.3863;
      d0=945;
      d1=945;
      d2=420;
      d3=105;
      d4=15;

x = f(:)/bw;     % frequency normalized to the bandwidth
ftype = lower(ftype);
switch ftype
    
    case 'movavg'
% Short term integrator

      Hf = sinc(x);
      
    case 'gauss' 
% Gaussian

      Hf = exp(-0.5*log(2).*x.*x);
      
    case 'gauss_off'
% Gaussian with offset

      Hf = exp(-0.5*log(2).*(x-ord/bw).*(x-ord/bw));
      
    case 'butt2'
% Butterworth 2nd order

      Hf = 1./(1-x.*x + i*sqrt(2)*x);
      
    case 'butt4'
% Butterworth 4th order

      x2 = x.*x;
      umx2 = 1-x2;
      Hf = 1./(umx2.*umx2-sqrt(2)*x2 + i*r4p2r2.*x.*umx2);
      
    case 'butt6'
% Butterworth 6th order

      x2 = x.*x;
      x3 = x2.*x;
      x4 = x3.*x;
      x5 = x4.*x;
      x6 = x5.*x;
      Hf = 1./(1.-b2*x2+b4*x4-x6 + i*(b1*x-b3*x3+b5*x5));
      
    case 'ideal'
% Ideal filter

      Hf = (abs(x) <= 0.5);
    
    case 'bessel5'
% Bessel 5th order

      om  = 2*pi*x*Bb;
      om2 = om.*om;
      om3 = om2.*om;
      om4 = om3.*om;
      om5 = om4.*om;

      pre = d0-d2*om2+d4*om4;
      pim = d1*om-d3*om3+om5;
      Hf = d0./(pre + i*pim);
      
     case 'rc1'
% RC filter 1st order == Butterworth 1st order

      Hf = 1./(1 + i*x);
      
     case 'rc2'
% RC filter 2nd order

      Hf = 1./(1 + i*sqrt(sqrt(2)-1)*x).^2;     % |Hf|^2=0.5 for x=+-1
      
    case 'supergauss'
% Super-Gaussian of order ORD        
      if nargin ~= 4
          error('missing superGauss order');
      end
      
      Hf = exp(-0.5*log(2).*x.^(2*ord));      
    
  otherwise
      error('the filter ftype does not exist.');
end


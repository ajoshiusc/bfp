function [c, colors] = bipolarcmapB(m, cRange, method, colors, method_param)
% Returns an approx. bipolar colormap centered at zero (in cRange). Zero
% maps to black.
%
% method is type of interpolation to be used. It can be
%     ['sigmoid'], 'linear', 'quadratic', or 'cubic'
%
% method_param defines the parameter for method.
%
% colors is the color combination to be used for +ve and -ve part. It can be
%     ['rg'] - red for -ve, green for +ve
%  or  'rb' - red for -ve, blue for +ve
%  or  'bg',  'br',  'gr', or 'gb'
% 
% When a parmeter is not specified, it gets the parameter from the current
% figure. When there is no figure, matlab creates one. 
%
% Usage: 
%   c = bipolarcmapB(m, cRange, method, colors, method_param);
%   c = bipolarcmapB(m, cRange, method, colors)
%   c = bipolarcmapB(m, cRange, method)
%   c = bipolarcmapB(m, cRange)
%   c = bipolarcmapB(method, colors)
%   c = bipolarcmapB(method)
%   c = bipolarcmapB(colors)
%   c = bipolarcmapB
%
% NOTE: Any other usage may result into weird result or error. 
%

if nargin < 1
   m = size(get(gcf,'colormap'),1);
   cRange = get(gca, 'clim');
   method = 'sigmoid';
   colors = 'rg';
   
elseif nargin==1 && length(m)==2
   colors = m;
   m = size(get(gcf,'colormap'),1);
   cRange = get(gca, 'clim');
   method = 'sigmoid';
   
elseif nargin==1 && length(m)>2
   method = m;
   m = size(get(gcf,'colormap'),1);
   cRange = get(gca, 'clim');
   colors = 'rg';
   
elseif nargin==2 && isnumeric(m)
   method = 'sigmoid';
   colors = 'rg';
   
elseif nargin==2 && ischar(m)
   method = m;
   colors = cRange;
   m = size(get(gcf,'colormap'),1);
   cRange = get(gca, 'clim');
   
elseif nargin==3
   colors = 'rg';
end

if nargin<=4 
   method_param = [];
end

if cRange(2)<cRange(1)
   error('HIGH must be greater than LOW in display range');
elseif prod(cRange)>=0
   error('LOW should be a -ve number & HIGH should be +ve number')
end

% deal with small inputs in a consistent way
if m < 3
   if m == 0
      c = zeros(0,3);
   elseif m == 1
      c = zeros(1,3);
   else
      c = [0 1 0; 1 0 0];
   end
   return;
end


m_len1 = nextOdd(m*max(abs(cRange)) / (cRange(2)-cRange(1)) ) - 1; % length of longer side
m_len2 = nextOdd(m*min(abs(cRange))/(cRange(2)-cRange(1))) - 1;
c1 = rgCMap(m_len1, method, method_param); % length (2*m_len1 + 1)
Ix = m_len1 + 1; % location of zero row in c1

% get symmetric colormap of correct size
if abs(cRange(1))<cRange(2)
   cB = c1( (Ix-m_len2):end, :);
else
   cB = c1(1:(Ix+m_len2), :);  
end
clear c1

cB = sum(cB,2);
[~, Ix] = sort(cB, 'ascend');
Ix = Ix(1); % location of zero row in cB
c = zeros(length(cB),3);

switch colors 
   case 'rg'
      p1=1; n1=2;
   case 'rb'
      p1=1; n1=3; 
   case 'gr'
      p1=2; n1=1;
   case 'gb'
      p1=2; n1=3;
   case 'br'
      p1=3; n1=1;
   case 'bg'
      p1=3; n1=2;
   otherwise
      error('unknown color combination!')
end

c(1:Ix,p1) = cB(1:Ix);
c(Ix+1:end,n1) = cB(Ix+1:end);

end


function outint = nextOdd(n)

if mod(round(n),2)==0
   outint = round(n)+1;
else
   outint = round(n);
end
end


function p = rgCMap(m, method, mparam)
% Returns reg-green colormap of size (2*m + 1) using 'method'.
% Similar output as redgreencmap() from bioinformatics toolbox, but this
% one has more options.

interpCol = (1/(m):1/(m):1);

switch lower(method)
   case 'tanh'
      if isempty(mparam)
         mparam = pi;
      end
      interpCol = tanh(mparam*(interpCol));
      
   case 'erf'
      if isempty(mparam)
         mparam = 1.85;
      end
      interpCol = erf(mparam*(interpCol));
      
   case 'sigmoid' % x/sqrt(1+x^2)
      if isempty(mparam)
         mparam = 2;
      end
      interpCol = (mparam*interpCol);
      interpCol = interpCol./sqrt(1+interpCol.^2);
      
   case 'linear'
      % Don't need to do anything
      
   case 'quadratic'
      interpCol = (interpCol).^(1/2);
      
   case 'cubic'
      interpCol = (interpCol).^(1/3);
      
   otherwise
      error('Unknown method used.')
end

red = [zeros(size(interpCol)) 0 interpCol];
green = fliplr(red);
blue = zeros(size(red));

p = [red', green', blue'];

end





function c = bipolarcmapW(varargin)
% Returns an approx. bipolar colormap centered at zero (in cRange). Zero
% maps to white.
%
% method is type of interpolation to be used. It can be
%     ['sigmoid'], 'linear', 'quadratic', or 'cubic'
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
%   c = bipolarcmapW(m, cRange, method, colors) 
%   c = bipolarcmapW(m, cRange, method)
%   c = bipolarcmapW(m, cRange)
%   c = bipolarcmapW(method, colors)
%   c = bipolarcmapW(method)
%   c = bipolarcmapW(colors)
%   c = bipolarcmapW
%
% NOTE: Any other usage may result into weird result or error. 
%


if nargin<1
   [cB, colors] = bipolarcmapB;
elseif nargin<2
   [cB, colors] = bipolarcmapB(varargin{1});
elseif nargin<3
   [cB, colors] = bipolarcmapB(varargin{1}, varargin{2});
elseif nargin<4
   [cB, colors] = bipolarcmapB(varargin{1}, varargin{2}, varargin{3});
elseif nargin==4
   [cB, colors] = bipolarcmapB(varargin{1}, varargin{2}, varargin{3}, varargin{4});
elseif nargin==5
   [cB, colors] = bipolarcmapB(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
else
   error('Too many input parameters.')
end

cB = sum(cB,2);
[~, Ix] = sort(cB, 'ascend');
Ix = Ix(1);
c = zeros(length(cB),3);

switch colors 
   case 'rg'
      p1=1; p2=2; p3=3;
      n1=2; n2=1; n3=3;
   case 'rb'
      p1=1; p2=2; p3=3;
      n1=3; n2=2; n3=1;
   case 'gr'
      p1=2; p2=1; p3=3;
      n1=1; n2=2; n3=3;
   case 'gb'
      p1=2; p2=1; p3=3;
      n1=3; n2=2; n3=1;
   case 'br'
      p1=3; p2=1; p3=2;
      n1=1; n2=2; n3=3;
   case 'bg'
      p1=3; p2=1; p3=2;
      n1=2; n2=1; n3=3;
end

c(1:Ix,p1) = 1;
c(1:Ix,p2) = 1-cB(1:Ix);
c(1:Ix,p3) = 1-cB(1:Ix);

c(Ix+1:end,n1) = 1;
c(Ix+1:end,n2) = 1-cB(Ix+1:end);
c(Ix+1:end,n3) = 1-cB(Ix+1:end);

end

% 
% BDP BrainSuite Diffusion Pipeline
% 
% Copyright (C) 2019 The Regents of the University of California and
% the University of Southern California
% 
% Created by Chitresh Bhushan, Divya Varadarajan, Justin P. Haldar, Anand A. Joshi,
%            David W. Shattuck, and Richard M. Leahy
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; version 2.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
% USA.
% 


function [x_out] = gradient_descent_adjust_step(funfcn, x_init, options)
% Gradient descent optimization with heuristic line-search (by reducing step_size). All options must
% be explicitly defined. options.step_size_scale_iter controls how step size is decreased. The
% optimizer returns a solution when any of options.MaxFunEvals, options.TolX, or options.TolFun is
% satisfied. 
%
% This is a conservative approach to backtracking line search and is substantially slower to converge
% to solution, as it requires heavy computation of gradients. However, it may be more stable in
% converging to the local valley (i.e. wont get out of the current valley easily). Experiment! Also
% see gradient_descent_backtrack().
%
%   funfcn - function handle which returns [fval grad]
%   x_init - Initial estimate of unknown parameters
%
%   options.step_size : Amplitude scaling of returned gradient. A value of 1.0 will use same
%                       magnitude of gradient as returned by funfcn. This step_size is decreased by
%                       a factor of options.step_size_scale_iter when a smaller fval (returned by
%                       funfcn) is not found with current step_size. In this case the optimizer
%                       backtracks and tries to find more suitable solution with smaller steps. 
%
%   options.step_size_scale_iter : See above
%   options.Display : Shows details of optimization at each step when set to 'iter'. Otherwise no
%                     display.
%   options.GradObj : Must be set to 'on'. 
%   options.MaxFunEvals : Maximum number of function evaluations allowed, a positive integer.
%   options.TolX : A positive scalar. Minimum change allowed in solution for optimizer to continue.
%   options.TolFun : A positive scalar. Minimum change allowed in fval for optimizer to continue.
%

if options.step_size_scale_iter>=1
   options.step_size_scale_iter = 0.7;
end

step_size = options.step_size;
funcCount = 0;
TolX = Inf;
X = x_init;
x_out = x_init;

fval = Inf;
fval_old = Inf;
TolFun = Inf;
timeID_start = tic;
exitflag = 0;

if strcmp(options.Display,'iter')
    disp('   Func-Grad-count         f(x)            fold()-f()     Step-size        Time Taken');
end


timeID_iter = tic;
[fval, grad_X] = feval(funfcn, X); % compute fval for x_init
funcCount = funcCount+1;

X_old = X;
fval_old = fval;
grad_Xold = grad_X;

if strcmp(options.Display,'iter')
   s=sprintf('  %5.0f -      %20.12g   %13.12g   %13.6g            %5.2f', funcCount, fval_old, TolFun, NaN, toc(timeID_iter)); disp(s);
end

if strcmpi(options.GradObj,'on')
   while ~exitflag
      timeID_iter = tic;
      
      % take a step and compute fval and grad for new X
      X = X_old - (step_size .* grad_Xold);
      [fval, grad_X] = feval(funfcn, X);
      
      funcCount = funcCount+1;
      TolX = norm(X_old - X);
      TolFun = fval_old - fval;
      
      if strcmp(options.Display,'iter')
         if TolFun>0, accp_str='Y'; else accp_str='N'; end
         s=sprintf('  %5.0f %c      %20.12g   %13.12g    %13.6g            %5.2f', funcCount, accp_str, fval, TolFun, TolX, toc(timeID_iter));
         disp(s);
      end
      
      if TolFun>0 % current X has lower fval, accept current X
         x_out = X;
         X_old = X;
         fval_old = fval;
         grad_Xold = grad_X;
         
      else % adjust stepsize and/or request more accurate gradient computation
         step_size = step_size * options.step_size_scale_iter;
         % request for more accurate gradient computation - NOT implemented 
      end
      
      
      if(funcCount > options.MaxFunEvals)
         exitflag = 2;
      elseif (TolX < options.TolX)
         exitflag = 3;
      elseif TolFun>0 && (TolFun < options.TolFun)
         exitflag = 4;
      end     
   end
   
else
   error('Function funfcn must provide gradient and options.GradObj must be ON.')
end

if strcmp(options.Display,'iter')
   if exitflag == 1,  disp('Cannot find minimum along gradient direction'); end
   if exitflag == 2,  disp('Function evaluations exceeded MaxFunEvals'); end
   if exitflag == 3,  disp('TolX is smaller than specified TolX'); end
   if exitflag == 4,  disp('TolFun is smaller than specified TolFun'); end
   s = sprintf('Total time taken = %10.2f sec.', toc(timeID_start)); disp(s);
end


end


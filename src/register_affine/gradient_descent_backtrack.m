function [x_out] = gradient_descent_backtrack(funfcn, x_init, options)
% Gradient descent optimization with backtracking search. All options must be explicitly defined.
% Also see backtracking_linesearch().
%
%   funfcn - function handle which returns [fval grad]
%   x_init - Initial estimate of unknown parameters
%
%   options.step_size : Initial step_size for backtracking. Should be a large number. 
%   options.step_size_scale_iter : backtracking tau. step_size is scacled by this value on each
%                                  iteration.
%   options.backtrack_c: Minmum fraction required for Armijo criteria.
%   options.Display : Shows details of optimization at each step when set to 'iter'. Otherwise no
%                     display.
%   options.GradObj : Must be set to 'on'. 
%   options.MaxFunEvals : Maximum number of function evaluations allowed, a positive integer.
%   options.TolX : A positive scalar. Minimum change allowed in solution for optimizer to continue.
%

if options.step_size_scale_iter>=1
   options.step_size_scale_iter = 0.5;
end

BTs_params.Display = options.Display;
BTs_params.init_alpha = options.step_size;
BTs_params.c = options.backtrack_c;
BTs_params.maxiters = options.MaxFunEvals/2; 
BTs_params.tau = options.step_size_scale_iter;
BTs_params.tolX = options.TolX;

funcCount = 0;
X = x_init;
x_out = x_init;
TolFun = Inf;
timeID_start = tic;
exitflag = 0;

if strcmp(options.Display,'iter')
   disp(' '); % new line
   disp('   Func-Grad-count         f(x)            f_old()-f()     Step-size        Time Taken');
end


% compute fval for x_init
timeID_iter = tic;
[fval, grad_X] = feval(funfcn, X);
funcCount = funcCount+1;

X_old = X;
fval_old = fval;
grad_Xold = grad_X;

if strcmp(options.Display,'iter')
   s = sprintf('  %5.0f -      %20.12g   %13.12g   %13.6g            %5.2f', funcCount, fval_old, TolFun, NaN, toc(timeID_iter)); 
   disp(s); drawnow();
end

if strcmpi(options.GradObj,'on')
   while ~exitflag
      timeID_iter = tic;
      
      % Line search and compute fval and grad for new X
      X = backtracking_linesearch(funfcn, X_old, fval_old, grad_Xold, -1*grad_Xold, BTs_params);
      [fval, grad_X] = feval(funfcn, X);
      
      funcCount = funcCount+1;
      TolX = norm(X_old - X);
      TolFun = fval_old - fval;
      
      if strcmp(options.Display,'iter')
         if TolFun>0, accp_str='Y'; else accp_str='N'; end
         s = sprintf('  %5.0f %c      %20.12g   %13.12g    %13.6g            %5.2f', funcCount, accp_str, fval, TolFun, TolX, toc(timeID_iter));
         disp(s); drawnow();
      end
      
      if TolFun>0 % current X has lower fval, accept current X
         x_out = X;
         X_old = X;
         fval_old = fval;
         grad_Xold = grad_X;
      end
      
      if(funcCount > options.MaxFunEvals)
         exitflag = 2;
      elseif (TolX < options.TolX)
         exitflag = 3;
      end
   end
   
else
   error('Function funfcn must provide gradient and options.GradObj must be ON.')
end

if strcmp(options.Display,'iter')
   if exitflag == 2,  disp('Function evaluations exceeded MaxFunEvals'); end
   if exitflag == 3,  disp('TolX is smaller than specified TolX'); end
   s = sprintf('Total time taken = %10.2f sec.', toc(timeID_start)); 
   disp(s); drawnow();
end


end


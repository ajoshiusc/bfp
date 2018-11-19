function [X_out, alpha, cntr] = backtracking_linesearch(funfcn, X, fval, gradX, search_dir, params)
% Backtracking with Armijo-condition. No computation of gradients.
%
%  funfcn - function handle which returns function evaluation
%  X           - current estimate of solution
%  fval, gradX - value & derivative of function at current estimate, i.e. [fval, gradX] = feval(funfcn, X)
%  search_dir  - search direction (NOT a unit vector)
%  params      - matlab structure, all field are required
%      params.Display    - 'iter', 'off', 'final'
%      params.init_alpha - Initial step size estimate; Should be large (but absolute value depends
%                          on the specific cost function)
%      params.c          - should be in (0,1); typically 1e-4; see algorithm below
%      params.maxiters   - max number of backtracking steps. should be ~100 unless feval is expensive
%      params.tau        - should be in (0,1); typically 0.5;  see algorithm below
%      params.tolX       - Minimum tolerance on change in X to continue with iterations
%
% Objective
%   If P is search direction, then find \alpha (stepsize) such that
%       f(X + \alpha P) <= f(X) + (\alpha * c * m)
%   where, m = <P, f'(X)> and c \in (0,1). That is, we want that the achieved reduction in
%   f(X + \alpha P) be at least a fixed fraction 'c' of the reduction promised by the tangent at X
%   (or that by the first-oder Taylor approximation of f() at X).
%
% Alogrithm:
%    0. set c \in (0,1), tau \in (0,1)
%    1. set t = -c * m; m = P'*search_dir
%    2. while f(X + \alpha P) > f(X) - \alpha * t
%          \alpha = tau * \alpha
%       end_while
%
%    3. X_out = X + \alpha P
%


if strcmp(params.Display,'iter'), timeID_start = tic(); end
m = search_dir' * gradX; % must be -ve
s_mag = norm(search_dir);

if m>=0
   X_out = X;
   exitflag = 2;
   
else
   t = -1 * params.c * m;
   alpha = params.init_alpha;
   cntr = 1;
   
   if strcmp(params.Display,'iter')
      fprintf('\t\t\t Cntr  alpha      f_new()      f()-f_new()   Time-taken\n');
      timeID_iter = tic();
   end
   
   fval_new = feval(funfcn, X + alpha*search_dir);
   
   while (fval_new > (fval - alpha*t)) && (cntr<params.maxiters) && (s_mag*alpha > params.tolX)
      
      if strcmp(params.Display,'iter')
         fprintf('\t\t\t   %d  %1.5f  %13.12g  %13.12g  %5.2f sec\n', cntr, alpha, fval_new, (fval-fval_new), toc(timeID_iter));
         timeID_iter = tic();
      end
      
      alpha = params.tau * alpha;
      fval_new = feval(funfcn, X + alpha*search_dir);
      cntr = cntr+1;
   end
   if strcmp(params.Display,'iter')
      fprintf('\t\t\t   %d  %1.5f  %13.12g  %13.12g  %5.2f sec\n', cntr, alpha, fval_new, (fval-fval_new), toc(timeID_iter));
   end
   
   
   if fval_new>=fval || (fval-fval_new)<=eps()
      X_out = X; % could not find lower fval
      exitflag = 2;
   else
      X_out = X + alpha*search_dir;
      if cntr>=params.maxiters
         exitflag = 1;
      elseif (s_mag*alpha < params.tolX)
         exitflag = 4;
      else
         exitflag = 3;
      end
   end
end

if strcmp(params.Display,'iter') || strcmp(params.Display,'final')
   if exitflag == 1,  fprintf('\t\t\t Number of backtracking iterations exceeded maxiters.'); end
   if exitflag == 2,  fprintf('\t\t\t Could not find a lower fval along search direction!'); end   
   if exitflag == 3,  fprintf(['\t\t\t Armijo condition satistisfied in ' num2str(cntr) ' iterations.']); end
   if exitflag == 4,  fprintf('\t\t\t Minimum tolX reached.'); end
   fprintf(' Total time taken = %f sec\n\n', toc(timeID_start));
end

end


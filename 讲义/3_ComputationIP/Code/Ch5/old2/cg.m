  function [x,term_code] = cg(x0,b,params,a_mult,m_invert);

%  function [x,term_code] = cg(x0,b,params,a_mult,m_invert);
%  function [x,term_code] = cg(x0,b,params,a_mult);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Use (preconditioned) CG iteration to minimize the quadratic functional
%       J(x) = b'*x + x'*A*x/2, 
%  where A is symmetric positive definite (SPD). Equivalently, solve
%  the linear system A*x = -b.
%
%  Input Variables:
%    x0       --- vector containing the initial guess.
%    b        --- vector containing right hand side of system A*x = -b.
%    a_mult   --- text string containing name of the MATLAB function
%       which implements multiplication by A. To compute y=A*x, call 
%       y = feval(a_mult,x,params). A must be SPD.
%    m_invert --- text string containing the name of the MATLAB
%       function which implements the preconditioner M. To solve M*x=y,
%       call x = feval(m_invert,u,params). M must be SPD.
%    params   --- MATLAB structure array containing CG parameters
%       and information used by a_mult and m_invert.
%    params.max_cg_iter      Maximimum number of CG iterations.
%    params.cg_step_tol      Stop CG when ||x_k+1 - x_k|| < step_tol.
%    params.grad_tol_factor  Stop CG when ||g_k|| < grad_tol_factor*||g_0||.
%    params.cg_io_flag       Output CG info if ioflag = 1.
%    params.cg_figure_no     Figure number for CG output.
%
%  Output Variables:
%    x         --- Approximate solution obtained by CG.
%    term_code --- CG termination code.
%    term_code = 1    Maximum CG iterations reached.
%    term_code = 2    CG step norm stopping tolerance reached.
%    term_code = 3    CG gradient norm stopping tolerance reached.
%    term_code = 4    Negative curvature detected. This suggests that
%                     either A or M is not SPD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  max_cg_iter     = params.max_cg_iter;
  step_tol        = params.cg_step_tol;
  grad_tol_factor = params.grad_tol_factor;
  io_flag         = params.cg_io_flag;
  cg_fig          = params.cg_figure_no;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CG initialization.                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if nargin == 5
    precond_flag = 1;   %  Use preconditioning.
  elseif nargin == 4
    precond_flag = 0;   %  No preconditioning.
  else
    fprintf('***** Incorrect no. of input arguments to cg.m.\n');
    return
  end

  x = x0;
  Ax = feval(a_mult,x,params);
  g = Ax - b;                       %  Initial gradient g = A*x0 - b.
  if precond_flag
    z = feval(m_invert,g,params);   %  Solve M*z = g.
  else
    z = g;                          %  No preconditioning.
  end
  d = -z;                           %  Initial search direction.
  delta = g(:)'*z(:);               %  delta = g' * M^{-1} * g.
  grad_tol = grad_tol_factor * sqrt(delta);  %  gradient stopping tolerance.
  stepnormvec = [];
  gradnormvec = sqrt(delta);
  term_code = 0;
  cg_iter = 0;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CG iteration.                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  while term_code == 0
    cg_iter = cg_iter + 1;
    
    Ad = feval(a_mult,d,params);    %  Compute A*d.
    dAd = d(:)'*Ad(:);              %  Compute d'*A*d.
    
    if dAd <= 0
      term_code = 1;                %  Negative curvature detected.
      fprintf('***** Negative curvature detected in CG.\n');
      x = [];
      return
    end

    tau = delta / dAd;              %  Line search parameter.
    x = x + tau*d;                  %  Update approximate solution.
    g = g + tau*Ad;                 %  Update gradient g = b + A*x.
    if precond_flag
      z = feval(m_invert,g,params); %  Solve M*z = g.
    else
      z = g;                        %  No preconditioning.
    end
    delta_new = g(:)'*z(:);         %  delta = g' * M^{-1} * g.
    my_beta = delta_new / delta;    %  Note: beta is a MATLAB function.
    d = -z + my_beta*d;             %  Update CG search direction.
    delta = delta_new;

    %  Compute and store CG convergence information.
    
    snorm = abs(tau)*norm(d(:));
    gnorm = sqrt(delta);
    stepnormvec = [stepnormvec; snorm];
    gradnormvec = [gradnormvec; gnorm];
    
    %  Plot CG convergence information.
    
    if ~isempty(cg_fig)
      figure(cg_fig)
      
      indx = [1:max(size(gradnormvec))]';
      subplot(221)
      if grad_tol > 0
        gconst = grad_tol * ones(size(indx));
        semilogy(indx,gradnormvec,'o-', indx,gconst,'-')
      else
	semilogy(indx,gradnormvec,'o-')
      end
      title('CG Gradient Norm')
      xlabel('CG iterate')
      
      if ~isempty(stepnormvec)
	subplot(222)
	indx = [1:max(size(stepnormvec))]';
	if step_tol > 0
          sconst = step_tol * ones(size(indx));
          semilogy(indx,stepnormvec,'o-', indx,sconst,'-')
	else
	  semilogy(indx,stepnormvec,'o-')
	end
	title('Norm of CG Step')
	xlabel('CG iterate')
      end
      drawnow
    end

    %  Check stopping criteria.
    
    if cg_iter >= max_cg_iter
      term_code = 1;
      fprintf('***** Max CG iterations exceeded.\n');
      return
    elseif snorm <= step_tol
      term_code = 2;
      fprintf('***** Step norm stopping tolerance reached in CG.\n');
      return
    elseif gnorm <= grad_tol
      term_code = 3;
      fprintf('***** Gradient norm stopping tolerance reached in CG.\n');
      return
    end
    
  end %(while term_code == 0)

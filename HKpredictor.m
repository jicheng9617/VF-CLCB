function [pre_hf,mse_hf,pre_lf,mse_lf] = HKpredictor(x,kriging_hf,kriging_lf)
nx=size(x,1);
if nx>1
    [p_lf,s2_lf]=predictor(x,kriging_lf);
    [p_hf,s2_hf]=predictorh(x,kriging_hf,kriging_lf);
else
    [p_lf,~,s2_lf]=predictor(x,kriging_lf);
    [temp_p,temp_s2]=predictorh([x;x],kriging_hf,kriging_lf);
    p_hf=temp_p(1,:); s2_hf=temp_s2(1,:);
end
pre_lf=p_lf;
mse_lf=s2_lf;
pre_hf=p_hf;
mse_hf=s2_hf;
end


function  [y, or1, or2, dmse] = predictorh(x, dmodel,dmodell)
%PREDICTOR  Predictor for y(x) using the given DACE model.
%
% Call:   y = predictor(x, dmodel)
%         [y, or] = predictor(x, dmodel)
%         [y, dy, mse] = predictor(x, dmodel) 
%         [y, dy, mse, dmse] = predictor(x, dmodel) 
%
% Input
% x      : trial design sites with n dimensions.  
%          For mx trial sites x:
%          If mx = 1, then both a row and a column vector is accepted,
%          otherwise, x must be an mx*n matrix with the sites stored
%          rowwise.
% dmodel : Struct with DACE model; see DACEFIT
%
% Output
% y    : predicted response at x.
% or   : If mx = 1, then or = gradient vector/Jacobian matrix of predictor
%        otherwise, or is an vector with mx rows containing the estimated
%                   mean squared error of the predictor
% Three or four results are allowed only when mx = 1,
% dy   : Gradient of predictor; column vector with  n elements
% mse  : Estimated mean squared error of the predictor;
% dmse : Gradient vector/Jacobian matrix of mse

% hbn@imm.dtu.dk
% Last update August 26, 2002
 
  or1 = NaN;   or2 = NaN;  dmse = NaN;  % Default return values
  if  isnan(dmodel.beta)
    y = NaN;   
    error('DMODEL has not been found')
  end

  [m n] = size(dmodel.S);  % number of design sites and number of dimensions
  sx = size(x);            % number of trial sites and their dimension
  if  min(sx) == 1 & n > 1 % Single trial point 
    nx = max(sx);
    if  nx == n 
      mx = 1;  x = x(:).';
    end
  else
    mx = sx(1);  nx = sx(2);
  end
  if  nx ~= n
    error(sprintf('Dimension of trial sites should be %d',n))
  end
  
  % Normalize trial sites 
  x1=x;
  x = (x - repmat(dmodel.Ssc(1,:),mx,1)) ./ repmat(dmodel.Ssc(2,:),mx,1);
  q = size(dmodel.Ysc,2);  % number of response functions
  y = zeros(mx,q);         % initialize result
  
  if  mx == 1  % one site only
    dx = repmat(x,m,1) - dmodel.S;  % distances to design sites
    if  nargout > 1                 % gradient/Jacobian wanted
      [f df] = feval(dmodel.regr, x1,dmodell);
      f = (f - repmat(dmodel.Ysc(1,:),mx,1))./ repmat(dmodel.Ysc(2,:),mx,1);
      [r dr] = feval(dmodel.corr, dmodel.theta, dx);
      % Scaled Jacobian
      dy = (df * dmodel.beta).' + dmodel.gamma * dr;
      % Unscaled Jacobian
      or1 = dy .* repmat(dmodel.Ysc(2, :)', 1, nx) ./ repmat(dmodel.Ssc(2,:), q, 1);
      if q == 1
        % Gradient as a column vector
        or1 = or1';
      end
      if  nargout > 2  % MSE wanted
        
        rt = dmodel.C \ r;
        u = dmodel.Ft.' * rt - f.';
        v = dmodel.G \ u;
        or2 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q);
        
        if  nargout > 3  % gradient/Jacobian of MSE wanted
          % Scaled gradient as a row vector
          Gv = dmodel.G' \ v;
          g = (dmodel.Ft * Gv - rt)' * (dmodel.C \ dr) - (df * Gv)';
          % Unscaled Jacobian
          dmse = repmat(2 * dmodel.sigma2',1,nx) .* repmat(g ./ dmodel.Ssc(2,:),q,1);
          if q == 1
            % Gradient as a column vector
            dmse = dmse';
          end
        end
        
      end
      
    else  % predictor only
      f = feval(dmodel.regr, x1,dmodell);
      r = feval(dmodel.corr, dmodel.theta, dx);
    end
    
    % Scaled predictor
    sy = f * dmodel.beta + (dmodel.gamma*r).';
    % Predictor
    y = (dmodel.Ysc(1,:) + dmodel.Ysc(2,:) .* sy)';
    
  else  % several trial sites
    % Get distances to design sites  
    dx = zeros(mx*m,n);  kk = 1:m;
    for  k = 1 : mx
      dx(kk,:) = repmat(x(k,:),m,1) - dmodel.S;
      kk = kk + m;
    end
    % Get regression function and correlation
    f = feval(dmodel.regr, x1,dmodell);
    f = (f - repmat(dmodel.Ysc(1,:),mx,1))./ repmat(dmodel.Ysc(2,:),mx,1);
    r = feval(dmodel.corr, dmodel.theta, dx);
    r = reshape(r, m, mx);
    
    % Scaled predictor 
    sy = f * dmodel.beta + (dmodel.gamma * r).';
    % Predictor
    y = repmat(dmodel.Ysc(1,:),mx,1) + repmat(dmodel.Ysc(2,:),mx,1) .* sy;
    
    if  nargout > 1   % MSE wanted
      rt = dmodel.C \ r;
      u = dmodel.G \ (dmodel.Ft.' * rt - f.');
      or1 = repmat(dmodel.sigma2,mx,1) .* repmat((1 + colsum(u.^2) - colsum(rt.^2))',1,q);
      if  nargout > 2
        disp('WARNING from PREDICTOR.  Only  y  and  or1=mse  are computed')
      end
    end
    
  end % of several sites
end
  
% >>>>>>>>>>>>>>>>   Auxiliary function  ====================

function  s = colsum(x)
% Columnwise sum of elements in  x
if  size(x,1) == 1,  s = x; 
else,                s = sum(x);  end
end


function  f = regpoly_c(S1,dmodel)

f = predictor(S1,dmodel);
% if  nargout > 1
%   [df] = predictor(S1,dmodel);
end
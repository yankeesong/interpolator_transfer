rm(list = ls())

require(pacman)
p_load(glmnet,matrixStats,Matrix,MASS,purrr,expm,ggplot2,nleqslv)


limit_model_shift = function(n1, n2, p, lambda, sigma, sigma_shift, sigma_beta){
  # This function calculates the limiting risk of the joint interpolator under model shift.
  # We assume two tasks case, where 1 represents source task and n2 represents target task
  # Generalization to multiple tasks from this code is possible.
  # Author: Yanke Song (ysong@g.harvard.edu)

  # Inputs:
  # n1, n2, p: dimensionality parameters
  # sigma, sigma_shift, sigma_beta: variances for noise, signal shift and target signal
  # lambda: regularization parameter.
  #   if lambda == 0, computes the interpolator risk
  #   if lambda > 0, computes the ridge estimator risk
  n = n1 + n2
  alpha = n2 / n1
  s = lambda * (1 + alpha)
  gamma = p / n
  gamma1 = p / n1
  gamma2 = p / n2
  if (lambda > 0) # Ridge
  {
    # ########################### Ridge ######################################
    # Equation for 3
    equation = function(x)
    {
      f = s*gamma1*x^2 + (1+alpha-gamma1+s)*x - 1
      return (f)
    }
    result = uniroot(equation, c(0,1/s))
    
    f1 = result$root
    f2 = f1^2*(1+gamma1*f1)^2 / ((1+gamma1*f1)^2-gamma1*f1^2)
    f3 = alpha / (1+gamma1*f1) + s
    
    
    B2_inter_limit = (1-2*f1*f3+f2*f3^2) / (1-gamma1*f2*(f3-s)/(1+gamma1*f1))
    
    #### From HMRT
    result = limit_HMRT(gamma, c(1), isotropic=TRUE, lambda=lambda)
    B1_inter_limit = result[1]
    V_inter_limit = result[2]
    
    risk_inter_limit = sigma_beta^2 * B1_inter_limit + 
      sigma_shift^2 * B2_inter_limit + sigma^2 * V_inter_limit
    
  } else # Interpolator
  {
    if (n == p)
    {
      stop("n and p cannot be equal, otherwise the interpolator risk limit is infinity")
    } else if (n < p)
    {
      ##################### Interpolator ##################################
      result = limit_HMRT(gamma, c(1), isotropic=TRUE)
      B1_inter_limit = result[1]
      V_inter_limit = result[2]
      ######### New results
      
      B2_inter_limit = n1*(p-n1) / p / (p-n1-n2)
      
      risk_inter_limit = sigma_beta^2 * B1_inter_limit + 
        sigma_shift^2 * B2_inter_limit + sigma^2 * V_inter_limit
      
    } else if (n > p)
    {
      V_inter_limit = p / (n-p)
      B_inter_limit = (n1^2*(n-p)+p*n1*n2) / (n^2*(n-p))
      risk_inter_limit = sigma^2 * V_inter_limit + sigma_shift^2 * B_inter_limit
    }
  }
  
  return (risk_inter_limit)
}

limit_HMRT = function(gamma, vs, isotropic = FALSE, lambda=0){
  # This is a helper function for results from HMRT paper:
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9481183/
  if (isotropic)
  {
    if (lambda == 0)
    {
      return (c(1-1/gamma, 1/(gamma-1)))
    } else
    {
      sq_helper = sqrt((1-gamma+lambda)^2+4*gamma*lambda)
      m = (1-gamma+lambda-sq_helper) / (-2*gamma*lambda)
      mp = (gamma*(sq_helper+lambda-2) - sq_helper + gamma^2 + lambda + 1) / (2*gamma*lambda^2*sq_helper)
      return (c(lambda^2*mp, gamma*(m-lambda*mp)))
    }
  }
  
  
  nv = length(vs)
  equation = function(x)
  {
    f = 1 - 1/gamma
    for (iv in 1:nv)
    {
      f = f - 1/nv / (1+x*gamma*vs[iv])
    }
    return (f)
  }
  
  result = uniroot(equation, c(0,10))
  c0 = result$root
  
  V_inter_limit_nume = 0
  V_inter_limit_denom = 0
  for (iv in 1:nv)
  {
    V_inter_limit_nume = V_inter_limit_nume + 1/nv * vs[iv]^2 / (1+c0*gamma*vs[iv])^2
    V_inter_limit_denom = V_inter_limit_denom + 1/nv * vs[iv] / (1+c0*gamma*vs[iv])^2
  }
  V_inter_limit = gamma * c0 * V_inter_limit_nume / V_inter_limit_denom
  B1_inter_limit = (1 + V_inter_limit) * V_inter_limit_denom
  
  return (c(B1_inter_limit, V_inter_limit))
}



limit_design_shift = function(n1, n2, p, lambda, sigma, sigma_beta, v1s, v2s){
  # This function calculates the limiting risk of the joint interpolator under design shift.
  # We assume two tasks case, where 1 represents source task and n2 represents target task
  # Generalization to multiple tasks from this code is possible.
  # Author: Yanke Song (ysong@g.harvard.edu)
  
  # Inputs:
  # n1, n2, p: dimensionality parameters
  # sigma, sigma_beta: variances for noise and target signal
  # v1s, v2s: the list of eigenvalues for source and target covariance matrices
  # lambda: regularization parameter.
  #   if lambda == 0, computes the interpolator risk
  #   if lambda > 0, computes the ridge estimator risk
  n = n1 + n2
  nv = length(v1s)
  if (lambda > 0) # Ridge
  {
    var_ridge_list_limit = function(lambda)
    {
      eqs = function(a)
      {
        a1 = a[1]
        a2 = a[2]
        f1 = a1 + a2 - 1
        f2 = a1 - n1/n
        for (iv in 1:nv)
        {
          f1 = f1 + 1/n * p/nv * (v1s[iv]*a1+v2s[iv]*a2)/(v1s[iv]*a1+v2s[iv]*a2+lambda)
          f2 = f2 + 1/n * p/nv * (v1s[iv]*a1)/(v1s[iv]*a1+v2s[iv]*a2+lambda)
        }
        
        c(f1,f2)
      }
      start = c(0,0)
      solution = nleqslv(start, eqs)
      if (solution$termcd > 1) {print('Warning: solution has not converged')}
      alpha1 = solution$x[1]
      alpha2 = solution$x[2]
      
      eqsp = function(ap)
      {
        a1p = ap[1]
        a2p = ap[2]
        f1 = a1p + a2p
        f2 = a1p
        for (iv in 1:nv)
        {
          f1 = f1 + 1/n * p/nv * (v1s[iv]*(lambda*a1p-alpha1)+v2s[iv]*(lambda*a2p-alpha2))/(v1s[iv]*alpha1+v2s[iv]*alpha2+lambda)^2
          f2 = f2 + 1/n * p/nv * (v1s[iv]*(lambda*a1p-alpha1)+v1s[iv]*v2s[iv]*(a1p*alpha2-a2p*alpha1))/(v1s[iv]*alpha1+v2s[iv]*alpha2+lambda)^2
        }
        c(f1,f2)
      }
      start = c(0,0)
      solution = nleqslv(start, eqsp)
      if (solution$termcd > 1) {print('Warning: solution has not converged')}
      alpha1p = solution$x[1]
      print(alpha1)
      print(alpha1/lambda)
      print(alpha1p)
      print(alpha1p-alpha1/lambda)
      alpha2p = solution$x[2]
      
      var_inter_limit = 0
      for (iv in 1:nv)
      {
        var_inter_limit = var_inter_limit + 1/n * p/nv * (v1s[iv]*v2s[iv]*(alpha1-lambda*alpha1p) + v2s[iv]^2*(alpha2-lambda*alpha2p)) / (v1s[iv]*alpha1+v2s[iv]*alpha2+lambda)^2
      }
      # print(var_inter_limit)
      return (var_inter_limit)
    }
    
    bias_ridge_list_limit = function(eta)
    {
      ve2s = 1 + eta * v2s
      
      eqs = function(a)
      {
        a1 = a[1]
        a2 = a[2]
        f1 = a1 + a2 - 1
        f2 = a1 - n1/n
        for (iv in 1:nv)
        {
          f1 = f1 + 1/n * p/nv * (v1s[iv]*a1+v2s[iv]*a2)/(v1s[iv]*a1+v2s[iv]*a2+lambda)
          f2 = f2 + 1/n * p/nv * (v1s[iv]*a1)/(v1s[iv]*a1+v2s[iv]*a2+lambda)
        }
        c(f1,f2)
      }
      start = c(0,0)
      solution = nleqslv(start, eqs)
      if (solution$termcd > 1) {print('Warning: solution has not converged')}
      alpha1 = solution$x[1]
      alpha2 = solution$x[2]
      
      eqsp = function(ap)
      {
        a1p = ap[1]
        a2p = ap[2]
        f1 = a1p + a2p
        f2 = a1p
        for (iv in 1:nv)
        {
          f1 = f1 + lambda/n * p/nv * (v1s[iv]*(a1p-alpha1*v2s[iv])+v2s[iv]*(a2p-alpha2*v2s[iv]))/(v1s[iv]*alpha1+v2s[iv]*alpha2+lambda)^2
          f2 = f2 + 1/n * p/nv * (v1s[iv]*v2s[iv]*(a1p*alpha2-alpha1*a2p)+v1s[iv]*lambda*(a1p-alpha1*v2s[iv]))/(v1s[iv]*alpha1+v2s[iv]*alpha2+ve2s[iv]*lambda)^2
        }
        c(f1,f2)
      }
      start = c(0,0)
      solution = nleqslv(start, eqsp)
      if (solution$termcd > 1) {print('Warning: solution has not converged')}
      alpha1p = solution$x[1]
      alpha2p = solution$x[2]
      
      
      bias_inter_limit = 0
      for (iv in 1:nv)
      {
        bias_inter_limit = bias_inter_limit + lambda * 1/nv * (v1s[iv]*alpha1p+v2s[iv]*(alpha2p+lambda))/(v1s[iv]*alpha1+v2s[iv]*alpha2+lambda)^2
      }
      return (bias_inter_limit)
    }
    
    risk_inter_limit = sigma_beta^2 * bias_ridge_list_limit(0) + sigma^2 * var_ridge_list_limit(lambda)
  } else # Interpolator
  {
    if (n == p)
    {
      stop("n and p cannot be equal, otherwise the interpolator risk limit is infinity")
    } else if (n < p)
    {
      var_ridge_list_limit = function(lambda)
      {

        eqs_inter = function(a)
        {
          a1 = a[1]
          a2 = a[2]
          f1 = -1
          f2 = -n1/n
          for (iv in 1:nv)
          {
            f1 = f1 + 1/n * p/nv * (v1s[iv]*a1+v2s[iv]*a2)/(v1s[iv]*a1+v2s[iv]*a2+1)
            f2 = f2 + 1/n * p/nv * (v1s[iv]*a1)/(v1s[iv]*a1+v2s[iv]*a2+1)
          }
          c(f1,f2)
        }
        start = c(0,0)
        solution = nleqslv(start, eqs_inter)
        if (solution$termcd > 1) {print('Warning: solution has not converged')}
        alpha1_inter = solution$x[1]
        alpha2_inter = solution$x[2]
        
        eqsp_inter = function(ap)
        {
          a1p = ap[1]
          a2p = ap[2]
          f1 = alpha1_inter + alpha2_inter
          f2 = alpha1_inter
          for (iv in 1:nv)
          {
            f1 = f1 + 1/n * p/nv * (v1s[iv]*a1p+v2s[iv]*a2p)/(v1s[iv]*alpha1_inter+v2s[iv]*alpha2_inter+1)^2
            f2 = f2 + 1/n * p/nv * (v1s[iv]*v2s[iv]*(a1p*alpha2_inter-a2p*alpha1_inter)+v1s[iv]*a1p)/(v1s[iv]*alpha1_inter+v2s[iv]*alpha2_inter+1)^2
          }
          c(f1,f2)
        }
        start = c(0,0)
        solution = nleqslv(start, eqsp_inter)
        if (solution$termcd > 1) {print('Warning: solution has not converged')}
        alpha1p_inter = solution$x[1]
        alpha2p_inter = solution$x[2]
        
        var_inter_limit = 0
        for (iv in 1:nv)
        {
          var_inter_limit = var_inter_limit - 1/n * p/nv * v2s[iv]*(alpha1p_inter*v1s[iv]+alpha2p_inter*v2s[iv]) / (alpha1_inter*v1s[iv]+alpha2_inter*v2s[iv]+1)^2
        }
        return (var_inter_limit)
      }
      
      bias_ridge_list_limit = function(eta)
      {
        ve2s = 1 + eta * v2s
        
        eqs_inter = function(a)
        {
          a1 = a[1]
          a2 = a[2]
          f1 = -1
          f2 = -n1/n
          for (iv in 1:nv)
          {
            f1 = f1 + 1/n * p/nv * (v1s[iv]*a1+v2s[iv]*a2)/(v1s[iv]*a1+v2s[iv]*a2+ve2s[iv])
            f2 = f2 + 1/n * p/nv * (v1s[iv]*a1)/(v1s[iv]*a1+v2s[iv]*a2+ve2s[iv])
          }
          c(f1,f2)
        }
        start = c(0,0)
        solution = nleqslv(start, eqs_inter)
        if (solution$termcd > 1) {print('Warning: solution has not converged')}
        alpha1_inter = solution$x[1]
        alpha2_inter = solution$x[2]
        
        eqsp_inter = function(ap)
        {
          a1p = ap[1]
          a2p = ap[2]
          f1 = 0
          f2 = 0
          for (iv in 1:nv)
          {
            f1 = f1 + 1/n * p/nv * (v1s[iv]*(a1p*ve2s[iv]-alpha1_inter*v2s[iv])+v2s[iv]*(a2p*ve2s[iv]-alpha2_inter*v2s[iv]))/(v1s[iv]*alpha1_inter+v2s[iv]*alpha2_inter+ve2s[iv])^2
            f2 = f2 + 1/n * p/nv * (v1s[iv]*v2s[iv]*(a1p*alpha2_inter-alpha1_inter*a2p)+v1s[iv]*(a1p*ve2s[iv]-alpha1_inter*v2s[iv]))/(v1s[iv]*alpha1_inter+v2s[iv]*alpha2_inter+ve2s[iv])^2
          }
          c(f1,f2)
          c(f1,f2)
        }
        start = c(0,0)
        solution = nleqslv(start, eqsp_inter)
        if (solution$termcd > 1) {print('Warning: solution has not converged')}
        alpha1p_inter = solution$x[1]
        alpha2p_inter = solution$x[2]
        
        # Here both alpha_inter and alphap_inter are scaled by 1/lambda
        bias_inter_limit = 0
        for (iv in 1:nv)
        {
          bias_inter_limit = bias_inter_limit + 1/nv * (v1s[iv]*alpha1p_inter+v2s[iv]*(alpha2p_inter+1))/(v1s[iv]*alpha1_inter+v2s[iv]*alpha2_inter+ve2s[iv])^2
        }
        return (bias_inter_limit)
      }
      
      risk_inter_limit = sigma_beta^2 * bias_ridge_list_limit(0) + sigma^2 * var_ridge_list_limit(lambda)
    } else if (n > p)
    {
      eqs = function(a)
      {
        a1 = a[1]
        a2 = a[2]
        f1 = a1 + a2 - 1 + p/n
        f2 = a1 - n1/n
        for (iv in 1:nv)
        {
          f2 = f2 + 1/n * p/nv * (v1s[iv]*a1)/(v1s[iv]*a1+v2s[iv]*a2)
        }
        c(f1,f2)
      }
      start = c(0.1,0.1)
      solution = nleqslv(start, eqs)
      if (solution$termcd > 1) {print('Warning: solution has not converged')}
      alpha1 = solution$x[1]
      alpha2 = solution$x[2]
      
      var_over_limit = 0
      for (iv in 1:nv)
      {
        var_over_limit = var_over_limit + 1/n * p/nv * v2s[iv] / (v1s[iv]*alpha1+v2s[iv]*alpha2)
      }
      risk_inter_limit = sigma^2 * var_over_limit
    }
  }
}


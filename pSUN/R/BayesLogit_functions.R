TRUNC = 0.64
cutoff = 1 / TRUNC;

a.coef <- function(n,x)
{
  if ( x>TRUNC )
    pi * (n+0.5) * exp( -(n+0.5)^2*pi^2*x/2 )
  else
    (2/pi/x)^1.5 * pi * (n+0.5) * exp( -2*(n+0.5)^2/x )
}

mass.texpon <- function(Z)
{
  x = TRUNC;
  fz = pi^2 / 8 + Z^2 / 2;
  b = sqrt(1.0 / x) * (x * Z - 1);
  a = -1.0 * sqrt(1.0 / x) * (x * Z + 1);

  x0 = log(fz) + fz * TRUNC;
  xb = x0 - Z + pnorm(b, log.p=TRUE);
  xa = x0 + Z + pnorm(a, log.p=TRUE);

  qdivp = 4 / pi * ( exp(xb) + exp(xa) );

  1.0 / (1.0 + qdivp);
}


rpg.devroye.1 <- function(Z)
{
  Z = abs(Z) * 0.5;

  ## PG(1,z) = 1/4 J*(1,Z/2)
  fz = pi^2 / 8 + Z^2 / 2;
  ## p = (0.5 * pi) * exp( -1.0 * fz * TRUNC) / fz;
  ## q = 2 * exp(-1.0 * Z) * pigauss(TRUNC, 1.0/Z, 1.0);

  num.trials = 0;
  total.iter = 0;

  rtigauss <- function(Z, R=TRUNC)
{
  Z = abs(Z);
  mu = 1/Z;
  X = R + 1;
  if (mu > R) {
    alpha = 0.0;
    while (runif(1) > alpha) {
      ## X = R + 1
      ## while (X > R) {
      ##     X = 1.0 / rgamma(1, 0.5, rate=0.5);
      ## }
      E = rexp(2)
      while ( E[1]^2 > 2 * E[2] / R) {
        E = rexp(2)
      }
      X = R / (1 + R*E[1])^2
      alpha = exp(-0.5 * Z^2 * X);
    }
  }
  else {
    while (X > R) {
      lambda = 1.0;
      Y = rnorm(1)^2;
      X = mu + 0.5 * mu^2 / lambda * Y -
        0.5 * mu / lambda * sqrt(4 * mu * lambda * Y + (mu * Y)^2);
      if ( runif(1) > mu / (mu + X) ) {
        X = mu^2 / X;
      }
    }
  }
  X;
}

  while (TRUE)
    {
      num.trials = num.trials + 1;

      if ( runif(1) < mass.texpon(Z) ) {
        ## Truncated Exponential
        X = TRUNC + rexp(1) / fz
      }
      else {
        ## Truncated Inverse Normal
        X = rtigauss(Z)
      }

      ## C = cosh(Z) * exp( -0.5 * Z^2 * X )

      ## Don't need to multiply everything by C, since it cancels in inequality.
      S = a.coef(0,X)
      Y = runif(1)*S
      n = 0

      while (TRUE)
        {
          n = n + 1
          total.iter = total.iter + 1;
          if ( n %% 2 == 1 )
            {
              S = S - a.coef(n,X)
              if ( Y<=S ) break
            }
          else
            {
              S = S + a.coef(n,X)
              if ( Y>S ) break
            }
        }

      if ( Y<=S ) break
    }

  ## 0.25 * X
  list("x"=0.25 * X, "n"=num.trials, "total.iter"=total.iter)
}

rpg.devroye.R <- function(num=1, h=1, z=0.0)
{
  n = h
    
  z = array(z, num);
  n = array(n, num);

  total.trials = 0;

  x = rep(0, num);
  for (i in 1:num) {
    x[i] = 0;
    for (j in 1:n[i]) {
      ## x[i] = x[i] + rpg.devroye.1(z[i])
      temp = rpg.devroye.1(z[i]);
      x[i] = x[i] + temp$x;
      total.trials = total.trials + temp$n;
    }
  }
  ## list("x"=x, "rate"=sum(n)/total.trials)
  x
}
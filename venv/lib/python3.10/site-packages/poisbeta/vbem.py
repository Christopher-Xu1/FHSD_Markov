import numpy as np
import scipy.optimize as so
import scipy.special as sp

def poisson_beta_elbo(theta, x, mu, alpha, beta):
  """Return the evidence lower bound

  theta - [ln k_r, ln k_on, ln k_off]
  x - array-like [n,]
  mu - array-like [n,]
  alpha - array-like [n,]
  beta - array-like [n,]

  """
  kr, kon, koff = np.exp(theta)
  return ((x + kon - alpha) * (sp.digamma(alpha) - sp.digamma(alpha + beta))
          + (mu + koff - beta) * (sp.digamma(beta) - sp.digamma(alpha + beta))
          + (x + mu) * np.log(kr) - kr - mu * np.log(mu) + mu - sp.gammaln(x + 1)
          - sp.betaln(kon, koff) + sp.betaln(alpha, beta)).sum()

def poisson_beta_delbo_dkon(ln_kon, ln_koff, alpha, beta):
  """Return the partial derivative of ELBO wrt kon"""
  return (sp.digamma(alpha)
          - sp.digamma(alpha - beta)
          - sp.digamma(np.exp(ln_kon))
          + sp.digamma(np.exp(ln_kon) + np.exp(ln_koff))).sum()

def poisson_beta_delbo_dkoff(ln_koff, ln_kon, alpha, beta):
  """Return the partial derivative of ELBO wrt koff"""
  return (sp.digamma(beta)
          - sp.digamma(alpha - beta)
          - sp.digamma(np.exp(ln_koff))
          + sp.digamma(np.exp(ln_kon) + np.exp(ln_koff))).sum()

def fit_poisson_beta_vbem(x, init=None, atol=1e-8, max_iters=1000, verbose=False):
  """Return kr, kon, koff

  init - [ln k_r, ln k_on, ln k_off]

  """
  if init is None:
    init = np.log(fit_poisson_beta_moment(x))
  theta = init
  mu = np.zeros(x.shape)
  alpha = np.ones(x.shape)
  beta = np.ones(x.shape)

  obj = -np.inf
  for t in range(max_iters):
    alpha = x + np.exp(theta[1])
    beta = mu + np.exp(theta[2])
    mu = np.exp(sp.digamma(beta) - sp.digamma(alpha + beta) + theta[0])
    theta[0] = np.log((x + mu).mean())
    opt_kon = so.newton(poisson_beta_delbo_dkon, x0=theta[1], args=(theta[2], alpha, beta))
    if not opt_kon.success:
      raise RuntimeError(f'k_on update failed: {opt_kon.message}')
    theta[1] = opt_kon.x
    opt_koff = so.root_scalar(poisson_beta_delbo_dkoff, x0=theta[2], args=(theta[1], alpha, beta))
    if not opt_koff.success:
      raise RuntimeError(f'k_on update failed: {opt_koff.message}')
    theta[1] = opt_koff.x
    update = poisson_beta_elbo(theta, x, mu, alpha, beta)
    if verbose:
      print(f'Epoch {t}: {update}')
    if abs(obj - update) < atol:
      return theta
    else:
      obj = update
  raise RuntimeError('failed to converge')

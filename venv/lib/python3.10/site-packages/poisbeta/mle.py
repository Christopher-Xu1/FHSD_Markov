import poisbeta
import numpy as np
import scipy.optimize as so
import scipy.special as sp
import scipy.stats as st

def poisson_beta_logpmf(theta, x, s, n_points=50):
  kr, kon, koff = np.exp(theta) + 1e-8
  # Important: Gauss-Jacobi quadrature computes the integral over t ∈ [-1, 1],
  # but we want the integral over p ∈ [0, 1]
  t, w = sp.roots_jacobi(n=n_points, alpha=koff - 1, beta=kon - 1)
  # (n_points, 1)
  p = ((1 + t) / 2).reshape(-1, 1)
  # (1, n_points) @ (n_points, n)
  px = w.reshape(1, -1) @ st.poisson(mu=s * kr * p).pmf(x.reshape(1, -1))
  # Important: extra 1/2 comes from u-substitution
  return np.log(px) - sp.betaln(kon, koff) - (kon + koff - 1) * np.log(2)

def poisson_beta_pmf(theta, x, s, n_points=50):
  return np.exp(poisson_beta_logpmf(theta, x, s, n_points))

def poisson_beta_neg_llik(theta, x, s, n_points=50):
  """Return the negative log likelihood of the data

  theta - [ln k_r, ln k_on, ln k_off]
  n_points - number of points used in numerical integration

  """
  return -poisson_beta_logpmf(theta, x, s, n_points).sum()

def fit_poisson_beta_mle(x, s=None, init=None, max_iters=1000, n_points=50):
  """Return ln k_r, ln k_on, ln k_off

  x - array-like [n,]
  init - [ln k_r, ln k_on, ln k_off]
  n_points - number of points used in numerical integration

  """
  if init is None:
    try:
      init = poisbeta.fit_poisson_beta_moment(x)
    except RuntimeError:
      init = np.zeros(3)
  if s is None:
    s = 1
  opt = so.minimize(poisson_beta_neg_llik, x0=init, args=(x, s, n_points),
                    method='Nelder-Mead', options={'maxiter': max_iters})
  if not opt.success:
    raise RuntimeError(f'failed to converge: {opt.message}')
  return opt.x

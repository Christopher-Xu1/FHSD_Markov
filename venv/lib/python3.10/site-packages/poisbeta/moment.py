import numpy as np

def fit_poisson_beta_moment(x, **kwargs):
  """Return ln kr, ln kon, ln koff

  Estimate kr, kon, koff using the first three exponential moments (Peccoud &
  Ycart 1995).

  x - array-like [n,]

  """
  moments = np.array([1, x.mean(), (x * (x - 1)).mean(), (x * (x - 1) * (x - 2)).mean()])
  ratios = moments[1:] / moments[:-1]
  kr = (2 * ratios[0] * ratios[2] - ratios[0] * ratios[1] - ratios[1] * ratios[2]) / (ratios[0] - 2 * ratios[1] + ratios[2])
  kon = (2 * ratios[0] * (ratios[2] - ratios[1])) / (ratios[0] * ratios[1] - 2 * ratios[0] * ratios[2] + ratios[1] * ratios[2])
  koff = (2 * (ratios[2] - ratios[1]) * (ratios[0] - ratios[2]) * (ratios[1] - ratios[0])) / ((ratios[0] * ratios[1] - 2 * ratios[0] * ratios[2] + ratios[1] * ratios[2]) * (ratios[0] - 2 * ratios[1] + ratios[2]))
  result = np.array([kr, kon, koff])
  if not (np.isfinite(result).all() and (result > 0).all()):
    raise RuntimeError('moment estimation failed')
  return np.log(result)

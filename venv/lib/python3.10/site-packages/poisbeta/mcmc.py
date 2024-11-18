import numpy as np
import poisbeta
import scipy.special as sp
import scipy.stats as st

def slice_sample(f, init, width=0.1, max_steps=100, bounds=None, **kwargs):
  """Return samples from the density proportional to exp(f)

  f - log density (up to a constant)
  init - initial value
  width - typical width of f
  max_steps - maximum number of times to step out
  bounds - support of f
  kwargs - additional arguments to f

  """
  if bounds is None:
    bounds = (-np.inf, np.inf)
  # Auxiliary variable defines the slice
  z = f(init, **kwargs) - np.random.exponential()
  # Step out
  left = init - width * np.random.uniform()
  right = left + width
  left_steps = int(np.random.uniform() * max_steps)
  for _ in range(left_steps):
    if left < bounds[0] or z < f(left, **kwargs):
      break
    left -= width
  left = np.clip(left, *bounds)
  for _ in range(max_steps - left_steps):
    if right > bounds[1] or z < f(right, **kwargs):
      break
    right += width
  right = np.clip(right, *bounds)
  # Step in
  while right > left:
    proposal = left + np.random.uniform() * (right - left)
    if z < f(proposal, **kwargs):
      return proposal
    elif proposal < init:
      left = proposal
    else:
      right = proposal
  raise RuntimeError('failed to find an acceptable sample')

def _safe_log(x, eps=1e-8):
  return np.log(x + eps)

def _cond_logpdf_p(p, x, kr, kon, koff):
  return (x + kon - 1) * _safe_log(p) + (koff - 1) * _safe_log(1 - p) - kr * p

def _cond_logpdf_kr(kr, x, p, ar, br):
  return ((x + ar - 1) * _safe_log(kr) - (p + br) * kr).sum()

def _cond_logpdf_kon(kon, p, koff, aon, bon):
  return (kon * _safe_log(p) + (aon - 1) * _safe_log(kon)
          - bon * kon - sp.gammaln(kon) + sp.gammaln(kon + koff)).sum()

def _cond_logpdf_koff(koff, p, kon, aoff, boff):
  return (koff * _safe_log(1 - p) + (aoff - 1) * _safe_log(koff)
          - boff * koff - sp.gammaln(koff) + sp.gammaln(kon + koff)).sum()

def fit_poisson_beta_mcmc(x, n_samples, ar, br, aon, bon, aoff, boff, max_steps=100, verbose=False, trace=False):
  """Return samples from the posterior p(kon, koff, kr | x)

  x - counts (n,)
  n_samples - number of samples to draw
  ar, br - prior parameters kr ~ Gamma(ar, br)
  aon, bon - prior parameters kon ~ Gamma(aon, bon)
  aoff, boff - prior parameters koff ~ Gamma(aoff, boff)
  max_steps - maximum number of steps in slice sampling
  verbose - if True, print log joint after each sample

  """
  # Important: these are fixed
  Fr = st.gamma(a=ar, scale=1 / br)
  Fon = st.gamma(a=aon, scale=1 / bon)
  Foff = st.gamma(a=aoff, scale=1 / boff)
  # Initialize at the moment estimate
  try:
    kr, kon, koff = np.exp(poisbeta.fit_poisson_beta_moment(x))
  except:
    kr = Fr.rvs(size=1)
    kon = Fon.rvs(size=1)
    koff = Foff.rvs(size=1)
  p = np.full(x.shape, 0.5)
  samples = []
  log_joint_trace = []
  for t in range(n_samples):
    samples.append((kr, kon, koff))
    for i in range(x.shape[0]):
      p[i] = slice_sample(_cond_logpdf_p, init=p[i], x=x[i], kr=kr, kon=kon, koff=koff, bounds=[1e-8, 1 - 1e-8])
    kr = slice_sample(_cond_logpdf_kr, init=kr, x=x, p=p, ar=ar, br=br, bounds=[1e-8, np.inf])
    kon = slice_sample(_cond_logpdf_kon, init=kon, p=p, koff=koff, aon=aon, bon=bon, bounds=[1e-8, np.inf])
    koff = slice_sample(_cond_logpdf_koff, init=koff, p=p, kon=kon, aoff=aoff, boff=boff, bounds=[1e-8, np.inf])
    if trace or verbose:
      log_joint = (st.poisson(mu=kr * p).logpmf(x).sum()
                   + st.beta(a=kon, b=koff).logpdf(p).sum()
                   + Fr.logpdf(kr)
                   + Fon.logpdf(kon)
                   + Foff.logpdf(koff))
      if trace:
        log_joint_trace.append(log_joint)
      if verbose:
        print(f'sample {t}: {log_joint}')
  return np.array(samples).reshape(-1, 3), log_joint_trace

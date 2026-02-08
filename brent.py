import numpy as np

def zero(f, a, b, t, eps=None):
    """
    Brent's method for root-finding, a direct translation of the ALGOL 60
    routine from Brent (1973). The original description from Brent's book
    is reproduced below.

    Procedure `zero` returns a zero `x` of the function `f` in the interval
    [a, b], to within a tolerance `6*eps*abs(x) + 2*t`, where `eps` is the
    relative machine precision and `t` is a positive tolerance. The procedure
    assumes that `f(a)` and `f(b)` have different signs.

    Parameters
    ----------
    f: Callable
        Function for which to find a root. Should take a single float
        argument and return a float value.
    a, b: float
        Initial bracket: two values on either side of the root, for which
        f(a)*f(b) < 0.
    t: float
        Absolute tolerance for the result. Should be positive.
    eps: float
        Relative tolerance. If `None` (the default), machine epsilon is used.

    Returns
    -------
    zero: Estimated value of the root.
    """
    fa = f(a)
    fb = f(b)

    # label: int
    c, fc = a, fa
    d = b - a
    e = d

    # label: ext
    if np.abs(fc) < np.abs(fb):
        a, fa = b, fb
        b, fb = c, fc
        c, fc = a, fa

    tol = 2*eps*np.abs(b) + t
    m = 0.5*(c - b)

    if np.abs(m) > tol and fb != 0:
        # See if a bisection is forced
        if np.abs(e) < tol and np.abs(fa) <= np.abs(fb):
            e = m
            d = e
        else:
            s = fb/fa
            if a == c:
                # Linear interpolation
                p = 2*m*s
                q = 1 - s
            else:
                # Inverse quadratic interpolation
                q = fa/fc
                r = fb/fc
                p = s*(2*m*q*(q - r) - (b - a)*(r - 1))
                q = (q - 1)*(r - 1)*(s - 1)
            if p > 0:
                q = -q
            else:
                p = -p
            s = e
            e = d
            if 2*p < 3*m*q - np.abs(tol*q) and p < np.abs(0.5*s*q):
                d = p/q
            else:
                e = m
                d = e
        a, fa = b, fb
        b += d if np.abs(d) > tol else (tol if m > 0 else -tol)
        fb = f(b)
        if (fb > 0) == (fc > 0):
            # go to: int
            pass
        else:
            # go to: ext
            pass

    return zero
        
    

import math
import sys
from collections.abc import Callable

macheps = sys.float_info.epsilon

def zero(
    f: Callable[[float], float],
    a: float,
    b: float,
    t: float,
    eps: float | None = None
) -> float:
    """
    Brent's method for root-finding, a direct translation of the ALGOL 60
    procedure from Richard P. Brent's book "Algorithms for Minimization Without
    Derivatives" (Prentice-Hall, 1973). The original description from Brent's
    book is reproduced below.

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
    x: float
        Estimated value of the root.
    """
    if eps is None:
        eps = macheps

    fa = f(a)
    fb = f(b)

    c, fc = a, fa
    d = b - a
    e = d

    while True:
        if abs(fc) < abs(fb):
            a, fa = b, fb
            b, fb = c, fc
            c, fc = a, fa

        tol = 2*eps*abs(b) + t
        m = 0.5*(c - b)

        if abs(m) <= tol or fb == 0:
            break

        # See if a bisection is forced
        if abs(e) < tol and abs(fa) <= abs(fb):
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
            if 2*p < 3*m*q - abs(tol*q) and p < abs(0.5*s*q):
                d = p/q
            else:
                e = m
                d = e
        a, fa = b, fb
        b += d if abs(d) > tol else (tol if m > 0 else -tol)
        fb = f(b)

        if (fb > 0) == (fc > 0):
            c, fc = a, fa
            d = b - a
            e = d

    return b

def localmin(
    f: Callable[[float], float],
    a: float,
    b: float,
    x: float,
    t: float,
    eps: float | None = None,
) -> float:
    """
    Brent's method for minimization, a direct translation of the ALGOL 60
    procedure from Richard P. Brent's book "Algorithms for Minimization Without
    Derivatives" (Prentice-Hall, 1973). The original description from Brent's
    book is reproduced below.

    If the function `f` is defined on the interval (a, b), then `localmin`
    finds an approximation `x` to the point at which `f` attains its maximum
    (or the appropriate limit point) and returns the value of `x`.
    `t` and `eps` define a tolerance `tol = eps*abs(x) + t`, and `f` is never
    evaluated at two points closer together than `tol`. If `f` is δ-unimodal
    for some δ < `tol`, then `x` approximates the global minimum of `f` with
    an error less than `3*tol`. If `f` is not δ-unimodal on (a, b), then `x`
    may approximate a local, but non-global, minimum. `eps` should be no
    smaller than `2*macheps`, and preferably not much less than `sqrt(macheps)`,
    where `macheps` is the relative machine precision. `t` should be positive.

    The method used is a combination of golden section search and successive
    parabolic interpolation. Convergence is never much slower than for a
    Fibonacci search. If `f` has a continuous second derivative which is
    positive at the minimum (not at `a` or `b`) then, ignoring rounding errors,
    convergence is superlinear, and usually the order is at least 1.3247...

    Parameters
    ----------
    f: Callable
        Function for which to find a root. Should take a single float
        argument and return a float value.
    a, b: float
        Initial bracket: two values on either side of the minimum.
    x: float
        Initial estimate of the minimum. Should satisfy f(x) < f(a) and
        f(x) < f(b).
    t: float
        Absolute tolerance for the result. Should be positive.
    eps: float
        Relative tolerance. If `None`, a default value is calculated based as
        the square root of the machine precision.

    Returns
    -------
    x: float
        Estimated value of the argument at which `f` attains its minimum.
    """
    if eps is None:
        eps = math.sqrt(macheps)

    c = (3 - math.sqrt(5))/2
    x = a + c*(b - a)
    e = 0
    fx = f(x)
    w, fw = x, fx
    v, fv = w, fw

    while True:
        m = 0.5*(a + b)
        tol = eps*abs(x) + t
        t2 = 2*tol

        if abs(x - m) <= t2 - 0.5*(b - a):
            break

        p = 0
        q = 0
        r = 0
        if abs(e) > tol:
            # Fit parabola
            r = (x - w)*(fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v)*q - (x - w)*r
            q = 2*(q - r)
            if q > 0:
                p = -p
            else:
                q = -q
            r = e
            e = d
        if abs(p) < abs(0.5*q*r) and p < q*(a - x) and p < q*(b - x):
            # A "parabolic interpolation" step
            d = p/q
            u = x + d
            # f must not be evaluated too close to a or b
            if u - a < t2 and b - u < t2:
                d = tol if x < m else -tol
        else:
            # A "golden section" step
            e = (b if x < m else a) - x
            d = c*e

        # f must not be evaluated too close to x
        u = x + (d if abs(d) >= tol else (tol if d > 0 else -tol))
        fu = f(u)

        # Update a, b, v, w, and x
        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu

    return x

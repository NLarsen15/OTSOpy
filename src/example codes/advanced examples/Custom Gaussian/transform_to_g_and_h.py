import numpy as np
from pymagglobal.utils import lm2i, i2lm_l


def transform_to_g_and_h(coeffs, lmax_out=None):
    """ Separate a set of coefficients in the pymagglobal "standard" order
    into g and h coefficients. Note that for compatibility with OTSO, this
    will append a degree and order 0 to the list.

    Parameters
    ----------
    coeffs : array
        The coefficients to be transformed

    Returns
    -------
    ns : array
        The coefficient degrees
    ms : array
        The coefficient orders
    coeffs_gh : array
        The coefficients as g and h lists. coeffs_gh[0] are the g coefficients,
        coeffs_gh[1] the h coefficients
    """
    lmax = i2lm_l(len(coeffs) - 1)
    if lmax_out is None:
        lmax_out = lmax

    ns = [0]
    ms = [0]

    coeffs_gh = np.zeros((2, (lmax_out * (lmax_out + 3)) // 2 + 1))

    idx = 1
    for ell in range(1, lmax+1):
        for emm in range(ell+1):
            g_idx = lm2i(ell, emm)
            ns.append(ell)
            ms.append(emm)

            coeffs_gh[0, idx] = coeffs[g_idx]

            if emm != 0:
                h_idx = lm2i(ell, -emm)
                coeffs_gh[1, idx] = coeffs[h_idx]

            idx += 1

    return np.array(ns, dtype=int), np.array(ms, dtype=int), coeffs_gh


if __name__ == '__main__':
    from pymagglobal import Model

    ggf = Model('GGF100k')

    coeffs = ggf.coefficients(-38_000)
    print(coeffs)
    ns, ms, gh = get_g_and_h_list(coeffs)
    print(ns)
    print(ms)
    print(gh)

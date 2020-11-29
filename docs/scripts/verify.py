# Independent re-computation of the 2020 MATLAB results (numpy/scipy), to check the code's math.
import sys, os
import numpy as np
from scipy.optimize import least_squares

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data")
sw = np.loadtxt(os.path.join(D, "swaminathan.txt"))          # rows: T, a, c
ju = np.loadtxt(os.path.join(D, "julian.txt"))               # cols: T, a, c
an = np.loadtxt(os.path.join(D, "anthracene.txt"))           # cols: P, a, b, c, beta

def aicc(rss, n, k):
    return n * np.log(rss / n) + 2 * k * n / (n - k - 1)

def poly_aicc(x, y, deg):
    X = np.vander(x, deg + 1, increasing=True)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    rss = float(np.sum((y - X @ beta) ** 2))
    return aicc(rss, len(y), deg + 2), beta

def weights(a):
    a = np.asarray(a); d = a - a.min(); w = np.exp(-d / 2); return w / w.sum()

for name, x, y, P in [("swaminathan a", sw[0], sw[1], 3), ("swaminathan c", sw[0], sw[2], 3),
                      ("julian a", ju[:, 0], ju[:, 1], 3), ("julian c", ju[:, 0], ju[:, 2], 3),
                      ("anthracene a", an[:, 0], an[:, 1], 10)]:
    vals = [poly_aicc(x, y, d)[0] for d in range(1, P + 1)]
    print(name, "AICc:", " ".join("%.4f" % v for v in vals))

# Einstein / Debye on Swaminathan a
T, a = sw[0], sw[1]
ein = lambda p, T: p[0] + p[1] / (np.exp(p[2] / T) - 1)
p2020 = np.array([5.56491, 0.0373932, 133.132])
print("Einstein AICc from 2020 params:", "%.4f" % aicc(np.sum((ein(p2020, T) - a) ** 2), 7, 4))
fit = least_squares(lambda p: ein(p, T) - a, [5.640, 9.4 * 280e-5, 280], method="lm")
print("Einstein refit:", fit.x, "AICc %.4f" % aicc(2 * fit.cost, 7, 4))

c = [3953.632, -800.6087, 85.07724, -4.43282, 0.0946173, 682.0012, 143.155337, 15.121491]
def deb(p, T):
    z = p[2] / T
    num = c[0] + c[1] * z + c[2] * z**2 + c[3] * z**3 + c[4] * z**4
    den = c[0] + c[5] * z + c[6] * z**2 + c[7] * z**3 + z**4
    return p[0] + p[1] * T * num / den
d2020 = np.array([5.56596, 0.000357386, 280])
print("Debye AICc from 2020 params:", "%.4f" % aicc(np.sum((deb(d2020, T) - a) ** 2), 7, 4))
fit = least_squares(lambda p: deb(p, T) - a, [5.640, 9.4 * 280e-5, 280], method="lm")
print("Debye refit:", fit.x, "AICc %.4f" % aicc(2 * fit.cost, 7, 4), fit.message)

print("weights (p1,p2,p3,E,D) 2020 AICc:", weights([-73.4674, -74.8968, -34.2499, -76.65, -68.3644]).round(4))

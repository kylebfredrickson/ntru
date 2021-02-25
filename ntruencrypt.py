from math import gcd
import numpy as np
from random import SystemRandom
from sympy import isprime

from polynomial import Polynomial

def uniform_sample_T(N, d_1, d_2):
    assert(N >= d_1 + d_2)

    random_gen = SystemRandom()

    coeffs = np.zeros(N + 1)
    indices = [i for i in range(N + 1)]

    for i in range(d_1):
        idx = random_gen.randrange(0, len(indices))
        coeffs[indices[idx]] = 1
        indices.remove(indices[idx])

    for i in range(d_2):
        idx = random_gen.randrange(0, len(indices))
        coeffs[indices[idx]] = -1
        indices.remove(indices[idx])

    return Polynomial(coeffs)

class NTRUEncrypt:
    def __init__(self, N, p, q, d = None):
        assert(isprime(N))
        assert(gcd(N, q) == 1)
        assert(gcd(p, q) == 1)
        assert(p < q)

        self.N = N
        poly_mod_coeffs = np.zeros(N + 1)
        poly_mod_coeffs[N] = 1
        poly_mod_coeffs[0] = -1
        self.poly_mod = Polynomial(poly_mod_coeffs)

        self.p = p
        self.q = q
        if d != None:
            self.d = d
        else:
            self.d = round(N / 3)

    def key_gen(self, f = None, g = None):
        if f is None:
            f = uniform_sample_T(self.N, self.d + 1, self.d)
        if g is None:
            g = uniform_sample_T(self.N, self.d, self.d)

        while True:
            try:
                f_inv_p = f.inverse(self.poly_mod, self.p)
                f_inv_q = f.inverse(self.poly_mod, self.q)

                self.f = f
                self.g = g
                self.f_inv_p = f_inv_p
                self.f_inv_q = f_inv_q
                self.h = (f_inv_q.__mul__(g, self.q)).__mod__(self.poly_mod, self.q)
                break
            except:
                print("Inverses of f do not exist. Choosing a new key.")
                f = uniform_sample_T(self.N, self.d + 1, self.d)

    def set_pub_key(self, h):
        self.h = h

    def get_pub_key(self):
        if self.h is None:
            raise ValueError("Key not set.")

        return self.h

    def encrypt(self, m, r = None):
        if self.h is None:
            raise ValueError("Public key not set.")

        if r is None:
            r = uniform_sample_T(self.N, self.d, self.d)

        blind = Polynomial(self.p).__mul__(self.h.__mul__(r, self.q), self.q)
        return (blind.__add__(m, self.q)).__mod__(self.poly_mod, self.q)

    def decrypt(self, e):
        if self.f is None or self.f_inv_p is None:
            raise ValueError("Private key not set.")

        a = self.f.__mul__(e, self.q).__mod__(self.poly_mod, self.q).center_lift(self.q)
        a.reduce(self.p)
        m = self.f_inv_p.__mul__(a, self.p).__mod__(self.poly_mod, self.p).center_lift(self.p)
        return m


def main():
    dec = NTRUEncrypt(251, 3, 257, 5)
    enc = NTRUEncrypt(251, 3, 257, 5)

    f_coeffs = np.zeros(252)
    f_coeffs[251] = 1
    f_coeffs[249] = 1
    f_coeffs[211] = 1
    f_coeffs[188] = 1
    f_coeffs[134] = -1
    f_coeffs[107] = 1
    f_coeffs[57] = 1
    f_coeffs[46] = -1
    f_coeffs[44] = -1
    f_coeffs[8] = -1
    f_coeffs[2] = -1
    f = Polynomial(f_coeffs)

    g_coeffs = np.zeros(252)
    g_coeffs[229] = 1
    g_coeffs[192] = -1
    g_coeffs[181] = 1
    g_coeffs[128] = -1
    g_coeffs[103] = 1
    g_coeffs[92] = -1
    g_coeffs[88] = 1
    g_coeffs[74] = -1
    g_coeffs[33] = 1
    g_coeffs[29] = -1
    g = Polynomial(g_coeffs)

    dec.key_gen(f, g)

    enc.set_pub_key(dec.get_pub_key())
    # m = Polynomial([1, -1, 1, 1, 0, -1])
    m = Polynomial([1, -1, 1, 1, 0, -1])
    # r = Polynomial([-1, 1, 0, 0, 0, -1, 1])
    r = Polynomial([1,1,1,1])
    e = enc.encrypt(m, r)

    m1 = dec.decrypt(m)

    print("Keys: ")
    print("f:", dec.f)
    print("g:", dec.g)
    print("p:", dec.f_inv_p)
    print("q:", dec.f_inv_q)
    print("h:", dec.h)

    print()

    print("Encryption")
    print("m:", m)
    print("r:", r)
    print("e:", e)
    print("dec(m):", m1)

    print()

    print("m = dec(m):", m == m1)

if __name__ == "__main__":
    main()

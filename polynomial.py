import numpy as np

def inv(a, mod):
    if mod == None:
        return 1 / a
    elif a == int(a):
        try:
            return pow(int(a), -1, mod)
        except:
            raise ValueError("Possible division by zero")
    else:
        raise TypeError("Tried to take find inverse of a float mod " + str(mod))

class Polynomial:
    def __init__(self, coefficients = None):
        if coefficients is None:
            self.coefficients = np.zeros(1)
        elif type(coefficients) is int or type(coefficients) is float:
            self.coefficients = np.array([coefficients])
        else:
            self.coefficients = np.array(coefficients)

    ### Boolean ###

    def __eq__(self, other):
        if type(other) is Polynomial:
            self._trim()
            other._trim()
            return np.array_equal(self.coefficients, other.coefficients)
        else:
            raise TypeError("Bad type in __eq__")

    ### Addition ###

    def __add__(self, other, mod = None):
        if type(other) is Polynomial:
            return self._polynomial_add(other, mod)
        else:
            raise TypeError("Bad type in __add__")

    def __sub__(self, other, mod = None):
        if type(other) is Polynomial:
            return self.__add__(Polynomial(-1 * other.coefficients), mod)
        else:
            raise TypeError("Bad type in __sub__")

    def _polynomial_add(self, other, mod = None):
        if len(self.coefficients) < len(other.coefficients):
            padded_coeffs = np.append(self.coefficients,
                                      np.zeros(abs(len(self.coefficients) - len(other.coefficients))))
            coeffs = padded_coeffs + other.coefficients
        else:
            padded_coeffs = np.append(other.coefficients,
                                      np.zeros(abs(len(self.coefficients) - len(other.coefficients))))
            coeffs = padded_coeffs + self.coefficients

        if mod != None:
            coeffs %= mod
            coeffs = np.asarray(coeffs, dtype=int)

        p = Polynomial(coeffs)
        p._trim()
        return p

    ### Multiplication ###

    def __mul__(self, other, mod = None):
        if type(other) is Polynomial:
            return self._polynomial_mul(other, mod)
        else:
            raise TypeError("Bad type in __mul__")

    def _polynomial_mul(self, other, mod):
        acc = Polynomial([0])
        for i in range(len(self.coefficients)):
            if (self.coefficients[i] == 0):
                continue
            p = Polynomial(np.append(np.zeros(i), self.coefficients[i] * other.coefficients))
            acc = acc.__add__(p, mod)

        acc._trim()
        return acc

    ### Division ###

    def __floordiv__(self, other, mod = None):
        if not type(other) is Polynomial:
            raise TypeError("Bad type in __truediv__")

        return self._divide_with_remainder(other, mod)[0]

    def __mod__(self, other, mod = None):
        if type(other) is Polynomial:
            return self._divide_with_remainder(other, mod)[1]
        else:
            raise TypeError("Bad type in __mod__")

    def _divide_with_remainder(self,other, mod = None):
        quotient = Polynomial(0)
        remainder = Polynomial(self.coefficients)
        while remainder.degree() >= other.degree() and not remainder == Polynomial(0):
            deg = remainder.degree() - other.degree()
            p_coeffs = np.zeros(deg + 1)
            p_coeffs[deg] = remainder.coefficients[remainder.degree()] * inv(other.coefficients[other.degree()], mod)
            p = Polynomial(p_coeffs)
            quotient = quotient.__add__(p, mod)
            remainder = remainder.__sub__(p.__mul__(other, mod), mod)
        return (quotient, remainder)

    ### Inverse ###

    # See https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
    def inverse(self, other, mod = None):
        if not type(other) is Polynomial:
            raise TypeError("Bad type in inverse")

        t, new_t = Polynomial(0), Polynomial(1)
        r, new_r = Polynomial(other.coefficients), Polynomial(self.coefficients)

        while new_r != Polynomial(0):
            quotient = r.__floordiv__(new_r, mod)
            (t, new_t) = new_t, t.__sub__(quotient.__mul__(new_t, mod), mod)
            (r, new_r) = new_r, r.__sub__(quotient.__mul__(new_r, mod), mod)

        if r.degree() > 0:
            raise ValueError("Inverse does not exist")

        return t.__mul__(Polynomial(1).__floordiv__(r, mod), mod)

    ### Misc. ###

    def center_lift(self, q):
        center_lift_coeffs = self.coefficients % q
        center_lift_coeffs = [c if c <= q / 2 else c - q for c in center_lift_coeffs]
        return Polynomial(center_lift_coeffs)

    def reduce(self, q):
        self.coefficients = self.coefficients % q

    def degree(self):
        nonzero = np.flatnonzero(self.coefficients)
        if (len(nonzero) != 0):
            return int(nonzero[len(nonzero) - 1])
        else:
            return 0

    def _trim(self):
        nonzero = np.flatnonzero(self.coefficients)
        if (len(nonzero) != 0):
            last = int(nonzero[len(nonzero) - 1])
            self.coefficients = self.coefficients[0: last + 1]
        else:
            self.coefficients = np.zeros(1)

    def __str__(self):
        coeff_vars = [str(self.coefficients[i]) + "x^" + str(i) for i in range(len(self.coefficients)) if self.coefficients[i] != 0]
        coeff_vars.reverse()
        return " + ".join(coeff_vars)

def main():
    a = Polynomial([1, 1])
    print("Basic Polynomial Ops")
    print(str(Polynomial(0)) + " + " + str(Polynomial(1)) + " = " + str(Polynomial(0) + Polynomial(1)))
    print(str(a) + " + " + str(a) + " = " + str(a + a))
    print(str(a) + " - " + str(a) + " = " + str(a - a))
    print(str(a) + " * " + str(a) + " = " + str(a * a))
    print(str(a) + " // " + str(a) + " = " + str(a // a))
    print(str(a) + " % " + str(a) + " = " + str(a % a))

    print()

    print("Mod Polynomial Ops (same signature for all ops)")
    b = Polynomial([3, 5])
    print(str(b) + " + " + str(b) + " mod 3 = " + str(b.__add__(b, 3)))

    print()

    a = Polynomial([-1, 0, 1, 1, -1, 0, 1])
    b = Polynomial([-1, 0, 0, 0, 0, 0, 0, 1])
    print("Inverses w.r.t Polynomials")
    print("Inv: " + str(a) + " w.r.t. " + str(b) + " mod 3 = " + str(a.inverse(b, 41)))

if (__name__ == "__main__"):
    main()


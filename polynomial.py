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

    def degree(self):
        nonzero = np.flatnonzero(self.coefficients)
        if (len(nonzero) != 0):
            return int(nonzero[len(nonzero) - 1])
        else:
            return 0

    def __eq__(self, other):
        self._trim()
        other._trim()
        return np.array_equal(self.coefficients, other.coefficients)

    def __add__(self, other, mod = None):
        if type(other) is Polynomial:
            return self._polynomial_add(other, mod)
        elif type(other) is int or type(other) is float:
            return self._polynomial_add(Polynomial([other]), mod)
        else:
            raise TypeError("Bad type in __add__")

    def __sub__(self, other, mod = None):
        if type(other) is Polynomial:
            return self.__add__(Polynomial(-1 * other.coefficients), mod)
        elif type(other) is int or type(other) is float:
            return self.__add__(-1 * other, mod)
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
        p = Polynomial(coeffs)
        p._trim()
        return p

    def __mul__(self, other, mod = None):
        if type(other) is Polynomial:
            return self._polynomial_mul(other, mod)
        elif type(other) is int or type(other) is float:
            coeffs = other * self.coefficients
            if mod != None:
                coeffs %= mod
            p = Polynomial(coeffs)
            p._trim()
            return p
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

    def __truediv__(self, other, mod = None):
        if type(other) is Polynomial:
            return self._divide_with_remainder(other, mod)[0]
        elif type(other) is int or type(other) is float:
            return self._divide_with_remainder(Polynomial(other), mod)[0]
        else:
            print(type(other))
            raise TypeError("Bad type in __truediv__")

    def __mod__(self, other, mod = None):
        if other is Polynomial:
            return self._divide_with_remainder(other, mod)[1]
        elif type(other) is int or type(other) is float:
            return self._divide_with_remainder(Polynomial(other), mod)[1]
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

    def _trim(self):
        nonzero = np.flatnonzero(self.coefficients)
        if (len(nonzero) != 0):
            last = int(nonzero[len(nonzero) - 1])
            self.coefficients = self.coefficients[0: last + 1]
        else:
            self.coefficients = np.zeros(1)

    def __str__(self):
        coeff_vars = [str(self.coefficients[i]) + "x^" + str(i) for i in range(len(self.coefficients))]
        coeff_vars.reverse()
        return " + ".join(coeff_vars)

def main():
    a = Polynomial([1,0, 2])
    b = Polynomial([1,2, 1])

if (__name__ == "__main__"):
    main()


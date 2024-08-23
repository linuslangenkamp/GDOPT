import math

class Dual:
    def __init__(self, real, dual=0.0):
        self.real = real  # The real part of the dual number
        self.dual = dual  # The dual part of the dual number

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.real + other.real, self.dual + other.dual)
        else:
            return Dual(self.real + other, self.dual)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Dual):
            return Dual(self.real - other.real, self.dual - other.dual)
        else:
            return Dual(self.real - other, self.dual)

    def __rsub__(self, other):
        return Dual(other - self.real, -self.dual)

    def __mul__(self, other):
        if isinstance(other, Dual):
            return Dual(self.real * other.real, self.real * other.dual + self.dual * other.real)
        else:
            return Dual(self.real * other, self.dual * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Dual):
            real_part = self.real / other.real
            dual_part = (self.dual * other.real - self.real * other.dual) / (other.real ** 2)
            return Dual(real_part, dual_part)
        else:
            return Dual(self.real / other, self.dual / other)

    def __rtruediv__(self, other):
        real_part = other / self.real
        dual_part = -other * self.dual / (self.real ** 2)
        return Dual(real_part, dual_part)

    def __pow__(self, power):
        if isinstance(power, Dual):
            real_part = self.real ** power.real
            dual_part = real_part * (power.real * self.dual / self.real + power.dual * math.log(self.real))
            return Dual(real_part, dual_part)
        else:
            return Dual(self.real ** power, power * self.real ** (power - 1) * self.dual)

    def __repr__(self):
        return f"Dual(real={self.real}, dual={self.dual})"


# Overriding math functions to handle Dual numbers

def sin(x):
    if isinstance(x, Dual):
        return Dual(math.sin(x.real), x.dual * math.cos(x.real))
    else:
        return math.sin(x)

def cos(x):
    if isinstance(x, Dual):
        return Dual(math.cos(x.real), -x.dual * math.sin(x.real))
    else:
        return math.cos(x)

def exp(x):
    if isinstance(x, Dual):
        return Dual(math.exp(x.real), x.dual * math.exp(x.real))
    else:
        return math.exp(x)

def log(x):
    if isinstance(x, Dual):
        return Dual(math.log(x.real), x.dual / x.real)
    else:
        return math.log(x)

# Example usage:

def f(x):
    return exp(x)

x = Dual(2.0, 1.0)  # (2.0 + ε) where ε^2 = 0
result = f(x)

print(f"f(x) at x = 2 is: {result.real}")
print(f"f'(x) at x = 2 is: {result.dual}")

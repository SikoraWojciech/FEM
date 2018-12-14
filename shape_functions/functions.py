

# shape functions
def N1(ksi, eta):
    return (1/4)*(1 - ksi)*(1 - eta)


def N2(ksi, eta):
    return (1/4)*(1 + ksi)*(1 - eta)


def N3(ksi, eta):
    return (1/4)*(1 + ksi)*(1 + eta)


def N4(ksi, eta):
    return (1/4)*(1 - ksi)*(1 + eta)


# ksi derivatives
def dN1_dksi(eta):
    return (-1/4)*(1 - eta)


def dN2_dksi(eta):
    return (1/4)*(1 - eta)


def dN3_dksi(eta):
    return (1/4)*(1 + eta)


def dN4_dksi(eta):
    return (-1/4)*(1 + eta)


# eta derrivatives
def dN1_deta(ksi):
    return (-1/4)*(1 - ksi)


def dN2_deta(ksi):
    return (-1/4)*(1 + ksi)


def dN3_deta(ksi):
    return (1/4)*(1 + ksi)


def dN4_deta(ksi):
    return (1/4)*(1 - ksi)


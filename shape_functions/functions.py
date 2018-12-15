

# shape functions
def N(i, ksi, eta):
    if i == 0:
        return (1 / 4) * (1 - ksi) * (1 - eta)
    elif i == 1:
        return (1 / 4) * (1 + ksi) * (1 - eta)
    elif i == 2:
        return (1 / 4) * (1 + ksi) * (1 + eta)
    elif i == 3:
        return (1 / 4) * (1 - ksi) * (1 + eta)
    else:
        return 0


# ksi derivatives
def dN_dksi(i, eta):
    if i == 0:
        return (-1 / 4) * (1 - eta)
    elif i == 1:
        return (1 / 4) * (1 - eta)
    elif i == 2:
        return (1 / 4) * (1 + eta)
    elif i == 3:
        return (-1 / 4) * (1 + eta)
    else:
        return 0


# eta derrivatives
def dN_deta(i, ksi):
    if i == 0:
        return (-1 / 4) * (1 - ksi)
    elif i == 1:
        return (-1 / 4) * (1 + ksi)
    elif i == 2:
        return (1 / 4) * (1 + ksi)
    elif i == 3:
        return (1 / 4) * (1 - ksi)
    else:
        return 0

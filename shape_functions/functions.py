

# shape functions
def N(index, ksi, eta):
    if index == 0:
        return (1 / 4) * (1 - ksi) * (1 - eta)
    elif index == 1:
        return (1 / 4) * (1 + ksi) * (1 - eta)
    elif index == 2:
        return (1 / 4) * (1 + ksi) * (1 + eta)
    elif index == 3:
        return (1 / 4) * (1 - ksi) * (1 + eta)
    else:
        return 0


# ksi derivatives
def dN_dksi(index, eta):
    if index == 0:
        return (-1 / 4) * (1 - eta)
    elif index == 1:
        return (1 / 4) * (1 - eta)
    elif index == 2:
        return (1 / 4) * (1 + eta)
    elif index == 3:
        return (-1 / 4) * (1 + eta)
    else:
        return 0


# eta derrivatives
def dN_deta(index, ksi):
    if index == 0:
        return (-1 / 4) * (1 - ksi)
    elif index == 1:
        return (-1 / 4) * (1 + ksi)
    elif index == 2:
        return (1 / 4) * (1 + ksi)
    elif index == 3:
        return (1 / 4) * (1 - ksi)
    else:
        return 0

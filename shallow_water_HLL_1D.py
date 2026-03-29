import numpy as np
import matplotlib.pyplot as plt

g = 9.81

def primitive(U):
    h = U[0].copy()
    hu = U[1].copy()
    u = np.zeros_like(h)
    mask = h > 1e-12
    u[mask] = hu[mask] / h[mask]
    return h, u

def hll_flux(UL, UR):
    hL, huL = UL
    hR, huR = UR

    uL = huL / hL if hL > 1e-12 else 0.0
    uR = huR / hR if hR > 1e-12 else 0.0
    cL = np.sqrt(g * max(hL, 0.0))
    cR = np.sqrt(g * max(hR, 0.0))

    sL = min(uL - cL, uR - cR)
    sR = max(uL + cL, uR + cR)

    FL = np.array([huL, huL * uL + 0.5 * g * hL**2])
    FR = np.array([huR, huR * uR + 0.5 * g * hR**2])

    if sL >= 0:
        return FL
    if sR <= 0:
        return FR
    return (sR * FL - sL * FR + sL * sR * (UR - UL)) / (sR - sL)

def solve_shallow_water_1d(L=1.0, nx=300, T=0.08, cfl=0.45, save_prefix="sw_hll_1d"):
    x = np.linspace(0.0, L, nx)
    dx = x[1] - x[0]

    h0 = 1.0 + 0.2 * np.exp(-400.0 * (x - 0.5) ** 2)
    u0 = np.zeros_like(x)

    U = np.zeros((2, nx))
    U[0] = h0
    U[1] = h0 * u0

    t = 0.0
    while t < T:
        h, u = primitive(U)
        amax = np.max(np.abs(u) + np.sqrt(g * h))
        dt = cfl * dx / max(amax, 1e-8)
        if t + dt > T:
            dt = T - t

        fluxes = np.zeros((2, nx - 1))
        for i in range(nx - 1):
            fluxes[:, i] = hll_flux(U[:, i], U[:, i + 1])

        U_new = U.copy()
        U_new[:, 1:-1] = U[:, 1:-1] - (dt / dx) * (fluxes[:, 1:] - fluxes[:, :-1])

        # transmissive boundaries
        U_new[:, 0] = U_new[:, 1]
        U_new[:, -1] = U_new[:, -2]

        # positivity fix
        U_new[0] = np.maximum(U_new[0], 1e-8)

        U = U_new
        t += dt

    h, u = primitive(U)

    plt.figure(figsize=(8, 4))
    plt.plot(x, h)
    plt.xlabel("x")
    plt.ylabel("h")
    plt.title("1D Shallow Water: Water Height")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{save_prefix}_height.png", dpi=150)
    plt.close()

    plt.figure(figsize=(8, 4))
    plt.plot(x, u)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.title("1D Shallow Water: Velocity")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{save_prefix}_velocity.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    solve_shallow_water_1d()

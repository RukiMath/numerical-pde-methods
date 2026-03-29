import numpy as np
import matplotlib.pyplot as plt

g = 9.81

def primitive(U):
    h = U[0]
    hu = U[1]
    hv = U[2]
    u = np.zeros_like(h)
    v = np.zeros_like(h)
    mask = h > 1e-12
    u[mask] = hu[mask] / h[mask]
    v[mask] = hv[mask] / h[mask]
    return h, u, v

def solve_shallow_water_2d(Lx=1.0, Ly=1.0, nx=80, ny=80, T=0.02, cfl=0.2, save_path="sw_2d.png"):
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, Ly, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    X, Y = np.meshgrid(x, y, indexing="ij")

    h0 = 1.0 + 0.1 * np.exp(-80 * ((X - 0.5) ** 2 + (Y - 0.5) ** 2))
    u0 = np.zeros_like(h0)
    v0 = np.zeros_like(h0)

    U = np.zeros((3, nx, ny))
    U[0] = h0
    U[1] = h0 * u0
    U[2] = h0 * v0

    t = 0.0
    while t < T:
        h, u, v = primitive(U)
        wave = np.max(np.abs(u) + np.sqrt(g * h) + np.abs(v) + np.sqrt(g * h))
        dt = cfl * min(dx, dy) / max(wave, 1e-8)
        if t + dt > T:
            dt = T - t

        h, u, v = primitive(U)

        F1 = U[1]
        F2 = U[1] * u + 0.5 * g * h**2
        F3 = U[1] * v

        G1 = U[2]
        G2 = U[2] * u
        G3 = U[2] * v + 0.5 * g * h**2

        U_new = U.copy()
        U_new[:, 1:-1, 1:-1] = (
            U[:, 1:-1, 1:-1]
            - (dt / dx) * (np.stack([F1, F2, F3], axis=0)[:, 1:-1, 1:-1] - np.stack([F1, F2, F3], axis=0)[:, :-2, 1:-1])
            - (dt / dy) * (np.stack([G1, G2, G3], axis=0)[:, 1:-1, 1:-1] - np.stack([G1, G2, G3], axis=0)[:, 1:-1, :-2])
        )

        # copy boundaries
        U_new[:, 0, :] = U_new[:, 1, :]
        U_new[:, -1, :] = U_new[:, -2, :]
        U_new[:, :, 0] = U_new[:, :, 1]
        U_new[:, :, -1] = U_new[:, :, -2]
        U_new[0] = np.maximum(U_new[0], 1e-8)

        U = U_new
        t += dt

    plt.figure(figsize=(7, 6))
    cs = plt.contourf(X, Y, U[0], levels=30)
    plt.colorbar(cs, label="Water height")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("2D Shallow Water")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()

if __name__ == "__main__":
    solve_shallow_water_2d()

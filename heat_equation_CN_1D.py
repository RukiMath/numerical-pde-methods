import numpy as np
import matplotlib.pyplot as plt

def solve_heat_cn_1d(L=1.0, T=0.2, alpha=1.0, nx=50, nt=200, save_path="heat_cn_1d.png"):
    dx = L / (nx - 1)
    dt = T / nt
    r = alpha * dt / dx**2

    x = np.linspace(0.0, L, nx)
    u = np.sin(np.pi * x)
    u[0] = 0.0
    u[-1] = 0.0

    n = nx - 2  # interior points

    A = np.diag((1 + r) * np.ones(n))
    A += np.diag((-r / 2) * np.ones(n - 1), 1)
    A += np.diag((-r / 2) * np.ones(n - 1), -1)

    B = np.diag((1 - r) * np.ones(n))
    B += np.diag((r / 2) * np.ones(n - 1), 1)
    B += np.diag((r / 2) * np.ones(n - 1), -1)

    for _ in range(nt):
        rhs = B @ u[1:-1]
        u[1:-1] = np.linalg.solve(A, rhs)

    u_exact = np.exp(-np.pi**2 * alpha * T) * np.sin(np.pi * x)

    plt.figure(figsize=(8, 5))
    plt.plot(x, u, label="Crank–Nicolson")
    plt.plot(x, u_exact, "--", label="Exact")
    plt.xlabel("x")
    plt.ylabel("u(x,T)")
    plt.title("1D Heat Equation")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()

if __name__ == "__main__":
    solve_heat_cn_1d()

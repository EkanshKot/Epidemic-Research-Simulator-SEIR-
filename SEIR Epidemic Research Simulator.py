# -*- coding: utf-8 -*-
"""
Created on Tue Feb  11 18:47:56 2024

@author: ekansh
"""

import numpy as np
import tkinter as tk
from tkinter import ttk
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

N = 1_000_000
y0 = [N - 30, 20, 10, 0]
days = 160
t = np.linspace(0, days, 400)

lockdown_start, lockdown_end = 30, 90
vaccination_start = 60
healthcare_capacity = 0.02

def seir(t, y, beta, sigma, gamma, lockdown, v_rate):
    S, E, I, R = y
    b = beta * (1 - lockdown if lockdown_start <= t <= lockdown_end else 1)
    v = v_rate if t >= vaccination_start else 0
    return [
        -b * S * I / N - v * S,
        b * S * I / N - sigma * E,
        sigma * E - gamma * I,
        gamma * I + v * S
    ]

def solve(beta, sigma, gamma, lockdown, v_rate):
    sol = solve_ivp(
        seir,
        (0, days),
        y0,
        t_eval=t,
        args=(beta, sigma, gamma, lockdown, v_rate)
    )
    return sol.y

root = tk.Tk()
root.title("SEIR Epidemic Research Simulator")

root.geometry("1400x800")

left = ttk.Frame(root)
left.pack(side=tk.LEFT, fill=tk.Y, padx=10)

right = ttk.Frame(root)
right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

def slider(text, frm, to, init):
    v = tk.DoubleVar(value=init)
    ttk.Label(left, text=text).pack(anchor="w")
    ttk.Scale(left, from_=frm, to=to, variable=v, orient=tk.HORIZONTAL).pack(fill="x", pady=2)
    return v

beta_v = slider("β – Transmission Rate (contact intensity)", 0.1, 0.6, 0.35)
sigma_v = slider("σ – Incubation Rate (1 / latent period)", 0.05, 0.5, 1/5.2)
gamma_v = slider("γ – Recovery Rate (1 / infectious period)", 0.05, 0.5, 1/10)
lock_v = slider("Lockdown Strength (reduction in β)", 0.0, 0.9, 0.6)
vax_v = slider("Vaccination Rate (fraction per day)", 0.0, 0.01, 0.002)

stochastic_v = tk.BooleanVar()
ttk.Checkbutton(
    left,
    text="Show stochastic uncertainty",
    variable=stochastic_v
).pack(pady=6)

fig, axs = plt.subplots(2, 2, figsize=(10, 8))
canvas = FigureCanvasTkAgg(fig, master=right)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

info = tk.Text(left, height=18, wrap=tk.WORD)
info.pack(fill=tk.X, pady=10)

def update():
    for ax in axs.flat:
        ax.clear()

    beta, sigma, gamma = beta_v.get(), sigma_v.get(), gamma_v.get()
    lockdown, v_rate = lock_v.get(), vax_v.get()

    S, E, I, R = solve(beta, sigma, gamma, lockdown, v_rate)

    axs[0,0].plot(t, S/N, label="Susceptible (S)")
    axs[0,0].plot(t, E/N, label="Exposed (E)")
    axs[0,0].plot(t, I/N, label="Infectious (I)")
    axs[0,0].plot(t, R/N, label="Recovered (R)")
    axs[0,0].set_title("SEIR Compartment Dynamics")
    axs[0,0].set_xlabel("Time (days)")
    axs[0,0].set_ylabel("Population Fraction")
    axs[0,0].legend()
    axs[0,0].grid(True)

    axs[0,1].plot(t, I/N, color="red")
    axs[0,1].axhline(healthcare_capacity, linestyle="--")
    axs[0,1].set_title("Infected Population vs Healthcare Capacity")
    axs[0,1].set_xlabel("Time (days)")
    axs[0,1].set_ylabel("Infected Fraction")
    axs[0,1].grid(True)

    axs[1,0].plot(S/N, I/N)
    axs[1,0].set_title("Phase Space: Susceptible vs Infected")
    axs[1,0].set_xlabel("Susceptible Fraction (S/N)")
    axs[1,0].set_ylabel("Infected Fraction (I/N)")
    axs[1,0].grid(True)

    if stochastic_v.get():
        dt = t[1] - t[0]
        for _ in range(15):
            y = np.zeros((4, len(t)))
            y[:,0] = y0
            for i in range(1, len(t)):
                S0,E0,I0,R0 = y[:,i-1]
                b = beta * np.random.normal(1,0.05)
                dS = -b*S0*I0/N
                dE = b*S0*I0/N - sigma*E0
                dI = sigma*E0 - gamma*I0
                dR = gamma*I0
                y[:,i] = np.maximum(y[:,i-1] + np.array([dS,dE,dI,dR])*dt,0)
            axs[1,1].plot(t, y[2]/N, alpha=0.3)
    else:
        axs[1,1].plot(t, I/N)

    axs[1,1].set_title("Stochastic Epidemic Trajectories")
    axs[1,1].set_xlabel("Time (days)")
    axs[1,1].set_ylabel("Infected Fraction")
    axs[1,1].grid(True)

    peak = np.max(I)/N
    t_peak = t[np.argmax(I)]
    total = R[-1]/N
    R0 = beta/gamma

    info.delete("1.0", tk.END)
    info.insert(tk.END,
        f"INTERPRETATION PANEL\n\n"
        f"β (Transmission Rate): {beta:.3f}\n"
        f"σ (Incubation Rate): {sigma:.3f}\n"
        f"γ (Recovery Rate): {gamma:.3f}\n"
        f"Basic Reproduction Number R₀ = β/γ = {R0:.2f}\n\n"
        f"Peak infected fraction: {peak:.2%}\n"
        f"Day of peak infection: {t_peak:.1f}\n"
        f"Total infected fraction: {total:.2%}\n\n"
        f"Graph meanings:\n"
        f"• Top-left: Flow of population through SEIR compartments\n"
        f"• Top-right: Healthcare stress relative to capacity\n"
        f"• Bottom-left: Stability and epidemic trajectory (phase space)\n"
        f"• Bottom-right: Uncertainty due to stochastic effects\n"
    )

    fig.tight_layout()
    canvas.draw_idle()

for v in [beta_v, sigma_v, gamma_v, lock_v, vax_v]:
    v.trace_add("write", lambda *args: update())

update()
root.mainloop()

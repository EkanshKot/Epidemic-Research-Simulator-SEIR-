#  Epidemic Research Simulator&#x20;

## Before you read:

I built this at the  beginning of 2025,My main goal was to explore how disease dynamics evolve under different interventions.

Some of the key features i added: 

1\)To representing uncertainty inherent in real epidemics, stochastic simulations are included. Random perturbations are introduced to the transmission rate, producing multiple plausible outbreak trajectories. This illustrates variability beyond deterministic predictions.

2\)The system of nonlinear differential equations is solved numerically using the Runge–Kutta method via SciPy’s solve\_ivp function. This allows accurate approximation of solutions over time when closed-form solutions are not available. 



Also While building this i learnt How reducing transmission rates or increasing recovery rates can flatten infection curves and reduce healthcare burden.



Also I believe that this project can serve as a learning tool for understanding conditional models and also a research oriented platform for experimenting with intervention strategies.



Now here is a formal read me:





## Overview

The SEIR Epidemic Research Simulator is an interactive epidemiological modeling tool developed in Python. It simulates the spread of infectious diseases such as COVID-19 using a compartmental SEIR model, enhanced with:

- Lockdown interventions
- Vaccination dynamics
- Healthcare capacity constraints
- Stochastic uncertainty modeling
- Real-time parameter tuning through a graphical interface

The simulator is designed for demonstration purposes, bridging mathematical epidemiology with computational modeling.

## Features

- Deterministic SEIR compartmental modeling
- Interactive Tkinter GUI with real-time updates
- Adjustable epidemiological parameters
- Lockdown and vaccination policy simulation
- Healthcare capacity stress visualization
- Stochastic epidemic trajectory generation
- Phase-space analysis for system stability
- Quantitative interpretation panel

## Mathematical Model

The population is divided into four compartments:

- **S (Susceptible)** – individuals who can become infected
- **E (Exposed)** – infected but not yet infectious
- **I (Infectious)** – actively transmitting the disease
- **R (Recovered)** – immune or removed individuals

Total population is conserved: S + E + I + R = N

### Governing Equations

```
dS/dt = −βSI/N − vS
dE/dt = βSI/N − σE
dI/dt = σE − γI
dR/dt = γI + vS
```

Where:

- β = transmission rate
- σ = incubation rate
- γ = recovery rate
- v = vaccination rate

## Epidemiological Interpretation

- **Basic Reproduction Number (R₀):** R₀ = β / γ; determines whether an outbreak grows (R₀ > 1) or dies out (R₀ < 1)
- **Peak Infection Load:** Maximum fraction of infected individuals over time
- **Healthcare Stress:** Infection curve compared against system capacity threshold
- **Phase-Space Plot (S vs I):** Shows nonlinear dynamics and epidemic stability
- **Stochastic Trajectories:** Capture uncertainty and variability seen in real epidemics

## Graph Descriptions

- **SEIR Compartment Dynamics:** Tracks population movement across S, E, I, R compartments
- **Infected vs Healthcare Capacity:** Indicates risk of healthcare system overload
- **Phase Space (S vs I):** Reveals epidemic trajectory and equilibrium behavior
- **Stochastic Epidemic Paths:** Demonstrates randomness in disease spread

##

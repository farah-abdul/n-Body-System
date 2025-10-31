# n-Body-System
This project models three-body and N-body orbital systems using numerical integration techniques to simulate gravitational interactions and orbital trajectories.   Developed as a group project as part of Introduction to Computational Fluid Dynamics at The Ohio State University

---

## Purpose
To numerically solve Newton’s equations of motion for multiple gravitational bodies and visualize resulting trajectories in 2D and 3D.  
The project demonstrates how integration accuracy and step size influence orbital stability and energy conservation in complex dynamical systems.

---

##  Features
- N-body gravitational simulation for arbitrary mass systems  
- Implemented **Runge–Kutta 4** solver  
- Energy and momentum conservation checks for validation  
- Configurable initial conditions (masses, positions, velocities)  
- 2D and 3D trajectory visualization using MATLAB plotting functions  

---

> Example output below shows a simulated three-body system illustrating gravitational perturbations and orbital resonance patterns.

![Three Body Trajectory](plots/three_body_trajectory.png)

---

## 🧮 Code Overview
| File | Description |
|------|--------------|
| `n_body_sim.m` | Main simulation function for N-body propagation |
| `N_body_script.m` | Main script for n-body trajectory |
| `earths_orbit.m` | mpdel's earth's 2D trajectory as it orbits the sun |
| `AE5615_P1_C.m` | script that determines the initial velocity that minimizes the distance between Earth's initial and final position after one simulated year |
| `Figure8Orbit_Partf.m` |Script that creates a figure-8 orbit between 3 celestial trajectories |


---

### Team Project Disclosure
This repository includes original MATLAB code developed by me and fellow teammates for AE5615 coursework.  

---

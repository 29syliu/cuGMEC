# cuGMEC

<p align="center">
  <b>cuGMEC: A High-Performance Code for Gyrokinetic-MHD Hybrid Simulation on GPUs with CUDA C++</b>
</p>

<p align="center">
  <a href="https://doi.org/10.1016/j.cpc.2026.110249">
    <img src="https://img.shields.io/badge/CPC-10.1016%2Fj.cpc.2026.110249-blue" alt="CPC DOI">
  </a>
  <img src="https://img.shields.io/badge/CUDA-76B900?style=flat&logo=nvidia&logoColor=white" alt="CUDA">
  <img src="https://img.shields.io/badge/C++-00599C?style=flat&logo=cplusplus&logoColor=white" alt="C++">
  <img src="https://img.shields.io/badge/License-GPLv3-blue" alt="GPLv3">
</p>

---

## 📌 Overview

**cuGMEC** is the fully reconstructed GPU version of the **GMEC** (*Gyrokinetic-MHD Energetic-particle Code*) family:

* **P. Y. Jiang, Z. Y. Liu, S. Y. Liu, J. Bao, and G. Y. Fu**<br>
  *Development of a gyrokinetic-MHD energetic particle simulation code. I. MHD version*<br>
  **Physics of Plasmas 31, 073904 (2024)**<br>
  DOI: [10.1063/5.0203252](https://doi.org/10.1063/5.0203252)

* **Z. Y. Liu, P. Y. Jiang, S. Y. Liu, L. L. Zhang, and G. Y. Fu**<br>
  *Development of a gyrokinetic-MHD energetic particle simulation code. II. Linear simulations of Alfvén eigenmodes driven by energetic particles*<br>
  **Physics of Plasmas 31, 073905 (2024)**<br>
  DOI: [10.1063/5.0206762](https://doi.org/10.1063/5.0206762)

* **S. Y. Liu, P. Y. Jiang, and G. Y. Fu**<br>
  *cuGMEC: A High-Performance Code for Gyrokinetic-MHD Hybrid Simulation on GPUs with CUDA C++*<br>
  **Computer Physics Communications (2026)**<br>
  DOI: [10.1016/j.cpc.2026.110249](https://doi.org/10.1016/j.cpc.2026.110249)

---

## 🧲 Physical Model

cuGMEC solves a nonlinear gyrokinetic-MHD hybrid model with fluid electrons and gyrokinetic ions, including both thermal ions and energetic particles.

The code can be used to study energetic-particle-driven Alfvén eigenmodes, such as TAE, RSAE, and BAE. It can also be applied to drift-wave and electromagnetic microinstabilities, such as ITG and KBM.

For details, please refer to our CPC paper.

### MHD Component

The MHD component includes:

* Gyrokinetic vorticity equation
* Gyrokinetic Poisson equation
* Parallel Ampère’s law
* Generalized Ohm’s law
* Electron continuity equation
* Electron isothermal condition

### PIC Component

Thermal ions and energetic particles are modeled gyrokinetically and advanced using the delta-f particle-in-cell method.

The MHD component and the PIC component are coupled through the perturbed pressure terms in the gyrokinetic vorticity equation.

---

## 🧮 Numerical Methods

cuGMEC uses:

* Five-point central finite-difference scheme for spatial discretization
* Fourth-order Runge-Kutta method for time integration
* Shifted metric coordinates
* Delta-f method for particle simulation

For details, please refer to our CPC paper.

---

## ⚡ Performance

The performance reported in our CPC paper submitted on November 6, 2025 is no longer the current best performance of cuGMEC. Further optimizations of both the MHD and PIC components have been carried out.

A comprehensive performance analysis, with direct comparisons to the published version, is under construction.

---

## 💻 How to Use

Under construction.

---

## 📄 License

cuGMEC is released under the **GNU General Public License v3.0**.

---

## ✉️ Contact

Using cuGMEC generally requires nontrivial preprocessing and postprocessing, including equilibrium preparation, metric conversion, input-file generation, and output analysis. If you are interested in using cuGMEC, it is strongly recommended to contact the developer first.

For questions, suggestions, or collaboration, please contact:

* [s.y.liu@zju.edu.cn](mailto:s.y.liu@zju.edu.cn)
* [29.sy.liu@gmail.com](mailto:29.sy.liu@gmail.com)

---

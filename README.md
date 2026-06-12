# cuGMEC

<p align="center">
  <b>cuGMEC: A High-Performance Code for Gyrokinetic-MHD Hybrid Simulation on GPUs with CUDA C++</b>
</p>

<p align="center">
  <a href="https://doi.org/10.1016/j.cpc.2026.110249"><img src="https://img.shields.io/badge/CPC-10.1016%2Fj.cpc.2026.110249-blue" alt="CPC DOI"></a>
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
  **Computer Physics Communications 327, 110249 (2026)**<br>
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

Thermal ions and energetic particles are modeled gyrokinetically and advanced using the delta-f particle-in-cell method. PIC-computed pressures enter the gyrokinetic vorticity equation through the pressure-curvature term.

---

## 🧮 Numerical Methods

cuGMEC uses:

* Shifted metric coordinates
* Delta-f method for particle simulation
* Five-point central finite-difference scheme for spatial discretization
* Fourth-order Runge-Kutta method for time integration

For details, please refer to our CPC paper.

---

## ⚡ Performance

The performance reported in our CPC paper no longer represents cuGMEC's best performance, as both the MHD and PIC components have since been optimized. The updated single-GPU speed tests below are based on the ITPA TAE benchmark case, using a 256x64x16 grid and dt = 0.025 Alfvén time (major radius divided by Alfvén velocity), and we believe the results are quite fast among codes of the same type. The MHD component uses double precision; the PIC component uses 400 keV fast ions with a maximum velocity of about 1.2 times the Alfvén velocity and is tested in both double and float precision. Note that particle deposition is performed in double precision in both double- and float-precision PIC runs. GYRO denotes the number of gyro-average points, and P/G denotes the average number of particles per grid. For convenience, all GPUs are tested using the same gridDim and blockDim, so the timings may differ very slightly from the absolute optimum.

Click the triangle below for results.

<details>
  <summary>NVIDIA GeForce RTX 4090 D</summary>

<p style="margin-top: 16px; margin-bottom: 2px;">The MHD per-step time is 14.6ms. The PIC per-step times are shown in the table.</p>

<table border="0" cellspacing="0" cellpadding="0" bgcolor="#ffffff" style="background-color: #ffffff; border: none; border-collapse: collapse;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff; border: none;">
    <td valign="top" bgcolor="#ffffff" style="background-color: #ffffff; border: none; padding-right: 36px;">
<table bgcolor="#ffffff" style="background-color: #ffffff;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td colspan="5" align="center" bgcolor="#ffffff" style="background-color: #ffffff;">double</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" width="82" height="42" bgcolor="#ffffff" style="background-color: #ffffff; width: 82px; min-width: 82px; max-width: 82px; height: 42px; padding: 0; line-height: 0;">
      <svg width="82" height="42" viewBox="0 0 82 42" xmlns="http://www.w3.org/2000/svg" style="display:block;">
        <line x1="0" y1="0" x2="74" y2="42" stroke="#d0d7de" stroke-width="1"/>
        <text x="74" y="17" text-anchor="end" font-size="13">GYRO</text>
        <text x="16" y="37" text-anchor="start" font-size="13">P/G</text>
      </svg>
    </td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">0</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">16</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">32</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">30.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">102ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">165ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">296ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">64</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">60.1ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">203ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">329ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">584ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">128</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">119ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">403ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">656ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">1.16s</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">256</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">237ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">806ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">1.33s</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">2.31s</td>
  </tr>
</table>
    </td>
    <td valign="top" bgcolor="#ffffff" style="background-color: #ffffff; border: none; padding-left: 0px;">
<table bgcolor="#ffffff" style="background-color: #ffffff;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td colspan="5" align="center" bgcolor="#ffffff" style="background-color: #ffffff;">float</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" width="82" height="42" bgcolor="#ffffff" style="background-color: #ffffff; width: 82px; min-width: 82px; max-width: 82px; height: 42px; padding: 0; line-height: 0;">
      <svg width="82" height="42" viewBox="0 0 82 42" xmlns="http://www.w3.org/2000/svg" style="display:block;">
        <line x1="0" y1="0" x2="74" y2="42" stroke="#d0d7de" stroke-width="1"/>
        <text x="74" y="17" text-anchor="end" font-size="13">GYRO</text>
        <text x="16" y="37" text-anchor="start" font-size="13">P/G</text>
      </svg>
    </td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">0</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">16</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">32</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">2.82ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8.96ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">15.0ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">28.5ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">64</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4.27ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">15.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">26.1ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">49.7ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">128</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">7.50ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">28.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">48.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">93.1ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">256</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">14.0ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">54.9ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">94.0ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">178ms</td>
  </tr>
</table>
    </td>
  </tr>
</table>

</details>

<details>
  <summary>NVIDIA A800-SXM4-80GB</summary>

<p style="margin-top: 16px; margin-bottom: 2px;">The MHD per-step time is 13.8ms. The PIC per-step times are shown in the table.</p>

<table border="0" cellspacing="0" cellpadding="0" bgcolor="#ffffff" style="background-color: #ffffff; border: none; border-collapse: collapse;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff; border: none;">
    <td valign="top" bgcolor="#ffffff" style="background-color: #ffffff; border: none; padding-right: 36px;">
<table bgcolor="#ffffff" style="background-color: #ffffff;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td colspan="5" align="center" bgcolor="#ffffff" style="background-color: #ffffff;">double</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" width="82" height="42" bgcolor="#ffffff" style="background-color: #ffffff; width: 82px; min-width: 82px; max-width: 82px; height: 42px; padding: 0; line-height: 0;">
      <svg width="82" height="42" viewBox="0 0 82 42" xmlns="http://www.w3.org/2000/svg" style="display:block;">
        <line x1="0" y1="0" x2="74" y2="42" stroke="#d0d7de" stroke-width="1"/>
        <text x="74" y="17" text-anchor="end" font-size="13">GYRO</text>
        <text x="16" y="37" text-anchor="start" font-size="13">P/G</text>
      </svg>
    </td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">0</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">16</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">32</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">7.26ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">26.3ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">44.4ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">78.4ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">64</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">13.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">48.9ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">83.8ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">148ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">128</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">26.2ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">94.2ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">162ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">288ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">256</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">51.2ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">184ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">319ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">566ms</td>
  </tr>
</table>
    </td>
    <td valign="top" bgcolor="#ffffff" style="background-color: #ffffff; border: none; padding-left: 0px;">
<table bgcolor="#ffffff" style="background-color: #ffffff;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td colspan="5" align="center" bgcolor="#ffffff" style="background-color: #ffffff;">float</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" width="82" height="42" bgcolor="#ffffff" style="background-color: #ffffff; width: 82px; min-width: 82px; max-width: 82px; height: 42px; padding: 0; line-height: 0;">
      <svg width="82" height="42" viewBox="0 0 82 42" xmlns="http://www.w3.org/2000/svg" style="display:block;">
        <line x1="0" y1="0" x2="74" y2="42" stroke="#d0d7de" stroke-width="1"/>
        <text x="74" y="17" text-anchor="end" font-size="13">GYRO</text>
        <text x="16" y="37" text-anchor="start" font-size="13">P/G</text>
      </svg>
    </td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">0</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">16</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">32</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">3.72ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">13.9ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">23.1ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">42.0ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">64</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">6.65ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">25.4ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">41.7ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">75.8ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">128</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">12.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">48.6ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">79.8ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">143ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">256</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">24.4ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">95.7ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">161ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">281ms</td>
  </tr>
</table>
    </td>
  </tr>
</table>

</details>

<details>
  <summary>NVIDIA RTX PRO 6000 Blackwell Server Edition</summary>

<p style="margin-top: 16px; margin-bottom: 2px;">The MHD per-step time is 9.76ms. The PIC per-step times are shown in the table.</p>

<table border="0" cellspacing="0" cellpadding="0" bgcolor="#ffffff" style="background-color: #ffffff; border: none; border-collapse: collapse;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff; border: none;">
    <td valign="top" bgcolor="#ffffff" style="background-color: #ffffff; border: none; padding-right: 36px;">
<table bgcolor="#ffffff" style="background-color: #ffffff;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td colspan="5" align="center" bgcolor="#ffffff" style="background-color: #ffffff;">double</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" width="82" height="42" bgcolor="#ffffff" style="background-color: #ffffff; width: 82px; min-width: 82px; max-width: 82px; height: 42px; padding: 0; line-height: 0;">
      <svg width="82" height="42" viewBox="0 0 82 42" xmlns="http://www.w3.org/2000/svg" style="display:block;">
        <line x1="0" y1="0" x2="74" y2="42" stroke="#d0d7de" stroke-width="1"/>
        <text x="74" y="17" text-anchor="end" font-size="13">GYRO</text>
        <text x="16" y="37" text-anchor="start" font-size="13">P/G</text>
      </svg>
    </td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">0</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">16</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">32</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">22.4ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">76.5ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">125ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">222ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">64</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">41.0ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">140ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">229ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">407ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">128</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">81.6ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">280ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">457ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">809ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">256</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">162ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">556ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">906ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">1.61s</td>
  </tr>
</table>
    </td>
    <td valign="top" bgcolor="#ffffff" style="background-color: #ffffff; border: none; padding-left: 0px;">
<table bgcolor="#ffffff" style="background-color: #ffffff;">
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td colspan="5" align="center" bgcolor="#ffffff" style="background-color: #ffffff;">float</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" width="82" height="42" bgcolor="#ffffff" style="background-color: #ffffff; width: 82px; min-width: 82px; max-width: 82px; height: 42px; padding: 0; line-height: 0;">
      <svg width="82" height="42" viewBox="0 0 82 42" xmlns="http://www.w3.org/2000/svg" style="display:block;">
        <line x1="0" y1="0" x2="74" y2="42" stroke="#d0d7de" stroke-width="1"/>
        <text x="74" y="17" text-anchor="end" font-size="13">GYRO</text>
        <text x="16" y="37" text-anchor="start" font-size="13">P/G</text>
      </svg>
    </td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">0</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">8</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">16</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">32</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">1.68ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">6.24ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">9.71ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">25.4ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">64</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">2.65ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">10.4ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">17.0ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">43.6ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">128</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">4.09ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">18.7ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">31.0ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">66.3ms</td>
  </tr>
  <tr bgcolor="#ffffff" style="background-color: #ffffff;">
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">256</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">7.90ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">36.6ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">60.3ms</td>
    <td align="center" bgcolor="#ffffff" style="background-color: #ffffff;">127ms</td>
  </tr>
</table>
    </td>
  </tr>
</table>

</details>

The total runtime of a simulation can be estimated as (MHD per-step time + PIC per-step time) × the total number of time steps. For other numerical parameters or multi-GPU runs, the runtime can be estimated by simple scaling. For details, please refer to our CPC paper.

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

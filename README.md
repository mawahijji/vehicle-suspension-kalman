# vehicle-suspension-kalman
MATLAB implementation of a Discrete Kalman Filter from scratch
# Vehicle Suspension State Estimator (Kalman Filter)

## 📌 Project Overview
This project implements a **Discrete Kalman Filter from scratch** (without using pre-built sensor fusion libraries) to estimate hidden vehicle states—specifically suspension deflection and vertical wheel velocity—using noisy accelerometer data.

The algorithm is based on a **Quarter-Car Model** and was developed to demonstrate the mathematical implementation of state-space estimation ($P = APA^T + Q$) in MATLAB.

## 🚀 Key Features
* **From-Scratch Implementation:** Manually programmed the Prediction and Correction steps to demonstrate understanding of linear algebra and control theory.
* **Physics Modeling:** Derived State-Space matrices (A, B, C, D) for a 4-DOF Quarter-Car system.
* **Sensitivity Analysis:** Tuned Process Noise ($Q$) and Measurement Noise ($R$) covariance matrices to balance responsiveness vs. noise rejection.

## 📊 Results
* **Accuracy:** Achieved a **6.57x reduction in RMSE** compared to raw sensor measurements.
* **Performance:** Validated zero-phase lag tracking of rapid transient road disturbances (bumps).
* **Verification:** Cross-validated results against a Simulink Model-Based Design.

## 🛠️ Tech Stack
* **Language:** MATLAB (R2023b+)
* **Concepts:** State Estimation, Linear Algebra, Control Systems

## 📂 File Structure
* `main.m`: The core algorithm (Physics setup, Filter loop, Plotting).

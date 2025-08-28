# Time-series-analysis  

Project1 – Detection of periodic signals in AGN red noise light curves: empirical tests on the Auto-Correlation Function and Phase Dispersion Minimization. Fortran codes are used to search for periodicities in the time series from large data set. Python codes are used for Data visualization.  

---

## Summary  

1. **Data Simulation**  
   - AGN-like light curves are simulated using the **Timmer & Koenig (1995)** algorithm.  
   - This method generates stochastic time series by randomizing Fourier phases and amplitudes for a chosen PSD shape.  
   - The resulting flux distributions follow Gaussian statistics.  
   - Multiple realizations are produced via **Monte Carlo simulations**, with different sampling strategies:  
     - **samp1** – 250 days, avg Δt ≈ 2 days  
     - **samp2** – 2 years, avg Δt ≈ 7 days  
     - **samp3** – 10 years, LSST uniform cadence (~3-day sampling)  
     - **samp4** – 10 years, yearly Sun gaps (40% missing coverage)  

   ![Data Simulation Figure](time_series_sims.png)   

2. **Validation**  
   - Testing **ACF (Auto-Correlation Function)** and **PDM (Phase Dispersion Minimization)** techniques for period searching.  
   - Null hypothesis testing is performed to distinguish genuine periodic signals from stochastic red noise variability.  

3. **Performance Testing**  
   - Empirical tests on the sensitivity and reliability of ACF and PDM in detecting periodic signals.  

---

## Repository Structure  

- **Fortran_codes/** – Periodicity detection codes (ACF & PDM implementations).  
- **Python_codes/** – Data visualization and plotting scripts.  
- **publications/** – Related publication(s).  
- **time_series_sims.png** – Example figure of simulated light curves.  
- **README.md** – Project overview (this file).  

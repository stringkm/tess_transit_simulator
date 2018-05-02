# tess_transit_simulator
Determine the physical parameter space over which TESS can detect an exoplanet

## Prototype Overall Plan for Repo

### Wrapper Function
* Use argparse to read in the parameters of the star and planet that you want
* Calls associated functions to simulate planet transit

### Simulate base transit curve
* Uses planet and star parameters to return base flux curve

### Add noise
* Applies distance modulus
* Adds statistical scatter based on flux

### Determine mag
* Convolves flux curve with TESS filter over course of transit
* multiplies flux density units over area of star
* Applies dimming due to extinction

### Determine S/N
* Estimates SNR from scatter
* Saves to file

### Plotter
* Plots SNR curves in parameter phase space

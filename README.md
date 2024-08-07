# Wind field correlation function inversion 

![Screenshot from 2024-06-26 09-08-30](https://github.com/jvierine/cfi/assets/16102494/c6e2d40f-88c5-4bec-bd5d-961cafd9df25)

A set of Python scripts to allow estimating second order correlation functions of wind from one dimensional projections of wind measured at random points in space and time. The main application is studies of turbulence in the mesosphere-lower thermosphere when using multi-static meteor radar measurements of specular meteor trail echo Doppler shifts. 

# Code citation

Please use as you wish. Reuse, copy, modify etc. I would appreciate if you cite my paper if you find the code or the method useful:

Vierinen, Juha, et al. "Observing Mesospheric Turbulence With Specular Meteor Radars: A Novel Method for Estimating Second‐Order Statistics of Wind Velocity." Earth and Space Science 6.7 (2019): 1171-1195.

https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2019EA000570

Don't hesitate to contact me if you are interested in collaborating.

# Other papers on the topic

Poblet FL, Vierinen J, Avsarkisov V, Conte JF, Charuvil Asokan H, Jacobi C, Chau JL. Horizontal Correlation Functions of Wind Fluctuations in the Mesosphere and Lower Thermosphere. Journal of Geophysical Research: Atmospheres. 2023 Mar 27;128(6):e2022JD038092.

Poblet FL, Chau JL, Conte JF, Avsarkisov V, Vierinen J, Charuvil Asokan H. Horizontal wavenumber spectra of vertical vorticity and horizontal divergence of mesoscale dynamics in the mesosphere and lower thermosphere using multistatic specular meteor radar observations. Earth and Space Science. 2022 Sep;9(9):e2021EA002201.

Charuvil Asokan H, Chau JL, Marino R, Vierinen J, Vargas F, Urco JM, Clahsen M, Jacobi C. Frequency spectra of horizontal winds in the mesosphere and lower thermosphere region from multistatic specular meteor radar observations during the SIMONe 2018 campaign. Earth, Planets and Space. 2022 May 11;74(1):69.

Vierinen, J., Chau, J. L., Charuvil, H., Urco, J. M., Clahsen, M., Avsarkisov, V., ... & Volz, R. (2019). Observing mesospheric turbulence with specular meteor radars: A novel method for estimating second‐order statistics of wind velocity. Earth and Space Science, 6(7), 1171-1195.

# Usage 

calculate correlation functions:
> mpirun -np 8 python3 hor_acf_daily.py -c norway_2022.ini
plot:
> python3 plot_hor_acf_hall.py -c norway_2022.ini



# Synthetic test for Radon based receiver function preprocessing 
![logo](https://user-images.githubusercontent.com/97296586/148636639-2e19b033-0dc5-41ca-bd5d-b8a727240a3f.png)

Quan Zhang Dec/30/2021

This code is a reproducible code for Fig.2 and Fig.4 in Zhang, 2022 (in prep.).

This code relys on the open-source package of jsonlab (https://github.com/fangq/jsonlab.git), SeismicLab (http://seismic-lab.physics.ualberta.ca/), and processRFmatlab (https://github.com/iwbailey/processRFmatlab.git). Before you start running this code, you need to install all of them and modify the path to where you have just installed.

Simulation is conducted with the 2D finite-difference code of SOFI2D (https://git.scc.kit.edu/GPIAG-Software/SOFI2D.git). High-resolution Radon transform is based on this paper «Chen, 2018, GEO, Automatic velocity analysis using high-resolution hyperbolic Radon transform».

![compare](https://user-images.githubusercontent.com/97296586/148637570-5abce4bc-36dd-48dc-a8cf-dae22e2cfded.png)
![compare_rf](https://user-images.githubusercontent.com/97296586/148637573-2944c9c1-90d8-42f0-8300-c4a206b3115d.png)

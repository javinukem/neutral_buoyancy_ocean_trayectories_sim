# neutral_buoyancy_ocean_trayectories_sim
Python code that computes neutral buoyancy trajectories for imaginary closed horizontal trajectories through the ocean. Bachelor Thesis done in 2024 at Imperial College with Juan Alfárez under Czaja Arnaud's supervision.

The WOCE dataset required to use this repository's code is available at:
[Download dataset](https://icdc.cen.uni-hamburg.de/thredds/catalog/ftpthredds/woce/catalog.html?dataset=ftpthreddsscan%2Fwoce%2Fwghc_params.nc)

## Abstract
The helix phenomena are one of several advection processes which determine how properties, like heat, salt or nutrients, are vertically transported in the ocean. This phenomenon is a result of studying trajectories that suffer a 0 net force. This happens when the buoyancy force, which results from being in a fluid, and the gravitational force, are equal, and therefore no energy needs to be spent in the movement. In doing so for a closed trajectory for the ocean (in the horizontal plane) it is found that the start and endpoint aren’t the same, there’s a height difference, a pitch. As a result, one can’t define a surface for the whole ocean in which every trajectory suffers a 0 net force, but it can be approximated, a technique that’s used by oceanographers nowadays. 

In this project, we developed a code in Python that computes these helices and their pitch for imaginary closed horizontal trajectories. To do so we make use of a data set of real ocean data points (WOCE) and interpolate between them to get the desired data, as well as the theory of neutral trajectories in the ocean studied by McDougall and Jackett in 1988. Using our code we’ve been able of validating the neutral surface approximation. We also checked that our code reproduces real physical phenomena like the Antarctic Polar Front and the Gulf Stream. A lot of work needs to be done to improve this code, but we are happy with this first attempt on it, as it seems to be already a useful tool to determine the difference in height of neutral trajectories and its relationship with thermodynamical properties
## Author

Javier Palau Alegria
- Contact email: javier.palau.alegria@gmail.com
- [LinkedIn](https://www.linkedin.com/in/javier-palau-alegria-215444237/)

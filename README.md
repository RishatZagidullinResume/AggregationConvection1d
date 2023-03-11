# Code for modelling spatially-heterogeneous coagulation in 1D with monomer source.

## Code description

Monomer source lets us find stationary solutions

* Parameters:
`h` - particle size step, 
`dt` - time step, 
`J` - source intencity, 
`TIME` - total number of iterations, 
`MOD` - periodicity of dumping data,
`size` - total number of particle sizes, 
`max_x` - total number of space steps, 
`dx` - space step, 
`vel_coefs` - advection coefficients array for different particle sizes, 
`dif_coefs` - diffusion coefficients array for different particle sizes, 
`initial_layer` - array with the initial condition (zero by default).

* Monomer source is incorporated into the boundary condition (see functions **AdvDifCoag1d::iteration** and **AdvDifCoag1d::fillMatrix** in __solver.cpp__)

* Data is dumped every `MOD` time iterations. `size` values are written after a newline. As a result we get a `size` x `TIME/MOD * max_x` matrix. Reading the data can be done in python using `.reshape` (see file __graph_output.ipynb__).

* There are 5 launch examples:
1. res_diff - No advection, constant diffusion, constant coagulation.
2. res_adv - Constant advection, constant diffusion, constant coagulation.
3. res_full - Variable advection, variable diffusion, coagulation kernel with $\frac{4}{3}$ homogeneity.
4. res_ballistic - Variable advection, variable diffusion, ballistic coagulation.
5. res_baikal - based on res_ballistic for comparison with experimental data from Baikal lake.

To launch the needed configuration comment and uncomment specified parts of the program.

* Command line arguments. For example, if you launch with `solver.exe output.txt 0 I.txt`, __output.txt__ file will appear where data will be dumped, __output.txt__ will appear where data on ozone concentrations are dumped for comparision. If you launch with `solver.exe output.txt 5` then it's assumed that __output.txt__ already exists and it has 6 iterations of data written. Then 5 iterations are skipped and 6th iteration is loaded as an initial condition.

## Demo

Based on [this publication](https://iopscience.iop.org/article/10.1088/1751-8121/ac711a/meta).

* Thanks to the stationary solution it is possible to derive analytical formulas. The presented numerical results show good agreement with analytical predictions:

![Fig1.jpg](/Fig1.jpg)

![Fig4.jpg](/Fig4.jpg)

* There is an attempt to compare the results with experimental data. To do that the model is coupled with a set of ODEs:

$$\frac{dC_{oz}}{dt} = I - \alpha C_{oz} -\beta \sum_{k=1}^\infty C_k k^{2/3} C_{oz}$$

$$\frac{dq}{dt} = \gamma C_{oz} (q_{max}-q) \Theta((q_{max}-q)) - \xi q$$

$$K_{ij} \rightarrow K_{ij}(1+q)$$

where $C_{oz}$ - ozone concentration, $I$ - source, $\alpha$, $\beta$, $\gamma$,  $\xi$, $q_{max}$ - model parameters, $q$ - reaction parameter, $\Theta$ - Heaviside function. The ODEs have a multiplication effect on the coagulation kernel. Here is the current result (work in progress):

![Fig5.png](/Fig5.png)



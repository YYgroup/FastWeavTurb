# Description

`FastWeavTurb` is used for the fast construction of multi-scale coherent vortices to synthesize turbulence.
The method generates multi-scale Burgers vortex tubes based on stochastic centerlines constructed by the fractional Brownian bridge.
The largest and smallest vortices govern the integral and dissipation scales, respectively, enabling precise control over the energy spectrum across the energy-containing, inertial, and dissipation ranges.
The adjustable vortex density allows tailored intermittency in velocity statistics. 
Remarkably, the woven turbulence achieves an extremely low computational cost, proportional to the total number of grid points.
This approach can simultaneously reproduce Gaussian velocity distributions, Kolmogorov's five-thirds scaling for the energy spectrum, as well as intertwined coherent vortices and intermittency akin to real turbulence.
This approach not only bridges the gap between structural and statistical turbulence modelling, but also offers an efficient tool for applications requiring realistic turbulent fields.

The code is written in FORTRAN and supports MPI parallel computing. 
If you are interested in using the code for your own research, please contact zishuo_han@stu.pku.edu.cn and yyg@pku.edu.cn for more details.

# Environment

The fortran code can be run on the `Tianhe-xy` supercomputer. The environment is set by the following commands: 

```
module load intel/oneapi2024.2_impi
export I_MPI_PMI_LIBRARY='/usr/lib/libpmi.so'
```

# Setup

Run `FBB_with_parameters.m` to generate vortex centerlines (`centerline_input.dat`) and parameters (`parameters_input.dat`) as inputs for the construction program of multi-scale vortex tubes. The default settings correspond to the case  "WT2".

# Run 

```
make
sh job.sh
```

# Postprocess

- to show the isosurfaces of vorticity magnitude, use the `./output/field_wh_norm.py`
- to show the statistics, see `./output/stat` ...

## References

Z. Han, W. Shen, and Y. Yang, Fast synthesis of turbulence with multi-scale coherent vortices, arXiv preprint, arXiv:2506.01022 (2025). [arXiv](https://arxiv.org/abs/2506.01022)

W. Shen, J. Yao, and Y. Yang, "Designing turbulence with entangled vortices", Proc. Natl. Acad. Sci. U.S.A., 121, e2405351121 (2024). [PNAS](https://www.pnas.org/doi/abs/10.1073/pnas.2405351121)

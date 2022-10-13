# README #

This is a modified (discrete dislocation dynamics) DDD framework based on ParaDiS-v2.5.1.
Please refer to the [manual](./docs/ParaDiSUsersGuide-V2.5.1-Pub.pdf) for the basic usage of ParaDiS.
The cited literature for the basic theoretical framework of DDD is
```latex
@Article{Arsenlis-Enablingstrainhardening-MaSiMSaE-2007,
  author    = {Arsenlis, A. and Cai, W. and Tang, M. and Rhee, M. and Oppelstrup, T. and Hommes, G. and Pierce, T. G. and Bulatov, V. V.},
  journal   = {Modelling and Simulation in Materials Science and Engineering},
  title     = {Enabling strain hardening simulations with dislocation dynamics},
  year      = {2007},
  issn      = {0965-0393},
  month     = {jul},
  number    = {6},
  pages     = {553--595},
  volume    = {15},
  abstract  = {Numerical algorithms for discrete dislocation dynamics simulations are investigated for the 
  purpose of enabling strain hardening simulations of single crystals on massively parallel computers. The 
  algorithms investigated include the ##IMG## [http://ej.iop.org/images/0965-0393/15/6/001/msmse240243in001.gif]
  {{\mathcal O}(N)} calculation of forces, the equations of motion, time integration, adaptive mesh 
  refinement, the treatment of dislocation core reactions and the dynamic distribution of data and work on 
  parallel computers. A simulation integrating all these algorithmic elements using the Parallel Dislocation 
  Simulator (ParaDiS) code is performed to understand their behaviour in concert and to evaluate the overall 
  numerical performance of dislocation dynamics simulations and their ability to accumulate percent of 
  plastic strain.},
  doi       = {10.1088/0965-0393/15/6/001},
  publisher = {{IOP} Publishing},
  type      = {Journal Article},
}

```

If you want to use the interface polycrystalline framework of this paper, please cite the following literature
```latex
@Article{Wei-Effectstwinboundary-AM-2019,
  author    = {Wei, DeAn and Zaiser, Michael and Feng, Zhiqiang and Kang, Guozheng and Fan, Haidong and Zhang, Xu},
  journal   = {Acta Materialia},
  title     = {Effects of twin boundary orientation on plasticity of bicrystalline copper micropillars:
  A discrete dislocation dynamics simulation study},
  year      = {2019},
  issn      = {1359-6454},
  month     = sep,
  pages     = {289--296},
  volume    = {176},
  abstract  = {Twin boundaries (TBs) constitute a special type of symmetric grain boundary (GB). TBs 
  influence plastic deformation in a complex manner. They not only act as dislocation obstacles but can also 
  accommodate twinning dislocations (TDs) whose motion enables twin boundary migration as an alternative 
  deformation mechanism. Exploiting this dual effect offers interesting perspectives in view of designing 
  materials that combine high strength and ductility. In the present work, we propose a model for 
  dislocation-TB interactions in the framework of discrete dislocation dynamics (DDD) simulations, which we 
  use to investigate the mechanical properties of bicrystalline copper micropillars containing a TB. We 
  systematically investigate how the compressive response depends on the orientation of the TB with respect 
  to the micropillar cross-section. The simulations show significant strengthening effects for TB 
  orientation angles less than 45°, where the interaction between mixed dislocations and TBs plays an 
  important role. For angles larger than 45°, the interaction between screw dislocations and TB dominates, 
  leading to weak strengthening. At 45°, dislocations on the dominant slip systems glide parallel to the TB 
  and no strengthening is observed.},
  doi       = {10.1016/j.actamat.2019.07.007},
  keywords  = {Discrete dislocation dynamics simulation Bicrystals Twin boundary Micropillar},
  publisher = {Elsevier {BV}},
  type      = {Journal Article},
}

```
### What is this repository for? ###

1. Dislocation - twin boundary model
2. Dislocation velocity kernel
3. Random internal stress field for HEAs research

### Dependency Libraries ###
1. MPI (necessary)
3. [PETSCs](https://petsc.org/release/) for velocity kernel
4. [voro++](https://math.lbl.gov/voro++/) for polycrystal
5. [GSL](https://www.gnu.org/software/gsl/) for data processing

### Install ###
* ` mkdir build && cd build && make -j4 && make install `

### Todo list ###
1. Merge the branche dev
2. Couple finite element method
3. Fast Fourier transform
4. thermal activation mechanism of cross-slip shift
5. thermal activation mechanism of dislocation-shaped kernels
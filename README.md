# SuperBand
The SuperBand database is an open access energy structure database for superconductors that have been synthesized experimentally before 2024-04.  It provides tools for data-mining and machine learning techniques. The universal features provided on our web https://www.superband.work/ interface facilitate the design of novel superconductors with a wide-range of applications. We describe the database and the algorithm to generate it in our paper [1] https://arxiv.org/abs/2409.09419. 

# Using the SuperBand
The SuperBand website http://www.superband.work/ offers a visual interface for exploring lattice structures, electronic band structures, and Fermi surface plots, providing an intuitive means of understanding the electronic properties of superconductors.

Our goal is to provide energy bands informations of all superconductors. Including
-1 Basic material informations, including Chemical Formula, Tc, reference and so on. .
-2 Material structure information in CIF format.
-3 Band Structure and DOS diagram.
-4 Fermi surface diagrams of (001) and (100) planes.
-5 3D Fermi surface map.
-6 Normalized band data for Machine Learning.

Note that in the github version of this dataset, in directory `superband/`, we only give four examples in SuperBand. 

# Using the  training set in SuperBand.

SuperBand training set has band structure data for 2474 superconductors.[2] Each material has 18 electronic bands around to the Fermi surface. Each band is mapped onto a 32 × 32 × 32 grid, yielding band data with dimensions of 18 × 32 × 32 × 32. 

SuperBand training set can be download from https://www.scidb.cn/en/file?id=c2487ff51a6230ba2ff55f0b5a8bfe8f

# Prerequisites

This code was developed and tested with and for Linux. Most likely it will throw errors for other OS. To install the Python packages we used miniconda 3.9 and jupyter-notebbok. 

Some package should be installed to run the choosen program, including ifermi, pymatgen, atomate, h5py, torch, vit_pytorch, skimage, BoltzTraP2, fireworks. 

# Installation

1. Download the superband repository into the current directory
   ```sh
   git clone https://github.com/ljcj007/superband.git
   ```
2. Use jupyter-notebbok to open the superband fold.

3. Test each `.ipynb`

# License
The SuperBand database is subject to the Creative Commons Attribution 4.0 License, implying that the content may be copied, distributed, transmitted, and adapted, without obtaining specific permission from the repository owner, provided proper attribution is given to the repository owner. All software in this repository is subject to the MIT license. See `LICENSE.md` for more information.

# Origin of data

We are grateful to the provider of different databases which have made the SuperBand possible:

- 1. The superconductor data is mainly sourced by SuperCon database, which is provided at [3] under a CC BY 4.0 license. 

- 2. The crystal structures are freely accessible and provided by the Materials Project database[4] and Open Quantum Materials Database (OQMD)[5] under a CC BY 4.0 license.

- 3. We applied the 3DSC methodology [6] to handle chemical formulas, including the definitions for exact matching, similarities, doping, and unmatched cases.

-4. High-throughput DFT calculations are facilitated by the Atomate open-source package [7], with parameter settings derived from the MIT High-Throughput Project [8]. For workflow automation, we employ the FireWorks package [9], which efficiently manages the task flow for structure optimization, static calculations, self-consistent field (SCF) calculations, and band structure determinations at high-symmetry points in reciprocal space.

-5. To extract relevant data for analysis, including band structure and DOS, we utilize the Pymatgen package [10]. For Fermi surface generation, analysis, and visualization, the Ifermi package [11] is used, enabling detailed examination of electronic properties crucial for understanding superconducting mechanisms.

## References
1. T. Zhang, C. Suo, Y. Wu, X. Xu, Y. Liu, D.-X. Yao, and J. Li, Superband: an Electronic-band and Fermi surface structure database of superconductor (2024), arXiv:2409.09419 https://arxiv.org/abs/2409.09419

2. J. Li, W. Fang, S. Jin, T. Zhang, Y. Wu, X. Xu, Y. Liu, and D.-X. Yao, A deep learning approach to search for superconductors from electronic bands (2024), arXiv:2409.07721 https://arxiv.org/abs/2409.07721

3. C. for Basic Research on Materials, Mdr supercon datasheet ver.240322 (2024). https://mdr.nims.go.jp/concern/datasets/m900p111d?locale=en

4. A. Jain, S. P. Ong, G. Hautier, W. Chen, W. D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, and K. A. Persson, Commentary: The Materials Project: A materials genome approach to accelerating materials innovation, APL Materials 1, 011002 (2013). https://materialsproject.org

5. S. Kirklin, J. E. Saal, B. Meredig, A. Thompson, J. W. Doak, M. Aykol, S. Ruhl, and C. Wolverton, The Open Quantum Materials Database (OQMD): assessing the accuracy of DFT formation energies, npj Computational Materials 1, 1 (2015). https://www.oqmd.org/

6. T. Sommer, R. Willa, J. Schmalian, and P. Friederich, 3DSC - a dataset of superconductors including crystal structures, Scientific Data 10, 1 (2023). https://github.com/aimat-lab/3DSC

7. K. Mathew, J. H. Montoya, A. Faghaninia, S. Dwarakanath, M. Aykol, H. Tang, I.-h. Chu, T. Smidt, B. Bocklund, M. Horton, J. Dagdelen, B. Wood, Z.-K. Liu, J. Neaton, S. P. Ong, K. Persson, and A. Jain, Atomate: A high-level interface to generate, execute, and analyze computational materials science
workflows, Computational Materials Science 139, 140 (2017). https://atomate.org/

8. A. Jain, G. Hautier, C. J. Moore, S. Ping Ong, C. C. Fischer, T. Mueller, K. A. Persson, and G. Ceder, A highthroughput infrastructure for density functional theory calculations, Computational Materials Science 50, 2295 (2011). https://doi.org/10.1016/j.commatsci.2011.02.023

9. A. Jain, S. P. Ong, W. Chen, B. Medasani, X. Qu, M. Kocher, M. Brafman, G. Petretto, G.-M. Rignanese,
G. Hautier, D. Gunter, and K. A. Persson, Fireworks: a dynamic workflow system designed for high-throughput applications, Concurrency and Computation: Practice and Experience 27, 5037 (2015). https://materialsproject.github.io/fireworks/

10. S. P. Ong, S. Cholia, A. Jain, M. Brafman, D. Gunter, G. Ceder, and K. A. Persson, The Materials Application Programming Interface (API): A simple, flexible and efficient API for materials data based on REpresentational State Transfer (REST) principles, Computational Materials Science 97, 209 (2015).
https://doi.org/10.1016/j.commatsci.2014.10.037

11. A. M. Ganose, A. Searle, A. Jain, and S. M. Griffin, IFermi: A python library for Fermi surface generation
and analysis, Journal of Open Source Software 6, 3089 (2021). https://fermisurfaces.github.io/IFermi/index.html




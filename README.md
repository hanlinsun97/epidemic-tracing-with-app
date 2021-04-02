# Epidemic tracing with app
Code for message-passing algorithm and Monte Carlo simulation on random networks and real-world networks.

This repository contains:

- 'GC_from_generator.jl': Implementation of the (averaged) message passing algorithm and Monte Carlo simulation on random poisson and scale-free networks.
- 'GC_from_dataset.jl': Implementation of the (averaged) message passing algorithm on real datasets (Friendship networks from the music streaming site Deezer in the countries of Romania, Hungary and Croatia). The original data can be found at https://github.com/benedekrozemberczki/GEMSEC.
- 'GC_MP_from_generator.c': Message passing algorithm on random Poisson network. Implemented in C. Output is a txt file containing three columns `$p$` (Probability of retaining a link), `$\rho$` (probability of adopting the app) and `$R$` (fraction of node in the giant component).
- 'GC_MC_from_generator.c': Monte Carlo simulation on random Poisson network, implemented in C. 
- 'GC_MP_from_dataset.c': Message Passing algorithm on real datasets, implemented in C. 
- 'GC_MC_from_dataset.c': Monte Carlo simulation on real datasets, implemented in C. 
- 'REF9_HR_edges.txt': Data of Croatia
- 'REF9_HU_edges.txt': Data of Hungary
- 'REF9_RO_edges.txt': Data of Romania

# How to use
The Julia code is implemented in Julia 1.6. Packages used are `DelimitedFiles` and `Plots`.

# Citing
If you find the code useful in your research, please cite the following paper:

```latex
@article{PhysRevResearch.3.L012014,
  title = {Message-passing approach to epidemic tracing and mitigation with apps},
  author = {Bianconi, Ginestra and Sun, Hanlin and Rapisardi, Giacomo and Arenas, Alex},
  journal = {Phys. Rev. Research},
  volume = {3},
  issue = {1},
  pages = {L012014},
  numpages = {6},
  year = {2021},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevResearch.3.L012014},
  url = {https://link.aps.org/doi/10.1103/PhysRevResearch.3.L012014}
}
```
# License
This code can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  
This program is distributed ny the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

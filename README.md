# Epidemic tracing with app
Code for message-passing algorithm and Monte Carlo simulation on random networks and real-world networks.

This repository contains:

- 'GC_from_generator.jl': Implementation of the (averaged) message passing algorithm and Monte Carlo simulation on random poisson and scale-free networks.
- 'GC_from_dataset.jl': Implementation of the (averaged) message passing algorithm on real dataset (Friendship networks from the music streaming site Deezer in the countries of Romania, Hungary and Croatia). The original data can be found at https://github.com/benedekrozemberczki/GEMSEC.

# How to use
The code is implemented in Julia 1.6. Packages used are DelimitedFiles and Plots.

# Citing
If you find the code useful in your research, please cite the following paper:

'''latex
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
'''

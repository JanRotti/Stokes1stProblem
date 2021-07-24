# Lattice Boltzmann Algorithm 
Implementation of Stokes first problem utilizing LBA with Bhatnagar-Gross-Krook model.

Rectangular fluid domain defined by equidistant grid with fixed length and height.

Boundary conditions:
  - Top: specular reflection
  - Bottom: bounce back with fixed wall velocity
  - Left+Right: periodict boundary with no overlap

Mesh:
![Mesh](/img/Mesh.png)
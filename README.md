# Coral Cover Dynamic Toy Model

This is a toy model to represent the evolution of the coral cover of a reef over the years. Here the corals are divided into three size classes $s = 1,2,3$ (representing small, medium and large), according to their diameter. Corals with diameter $d$ between $d_0 = 0$ and $d_1$ are small; corals with diameter between $d_1$ and $d_2$ are medium; and corals with diameter between $d_2$ and $d_3$ are large. Within each size class the corals are supposed to have an uniform distribution of diameters.

At each time step all corals grow in size proportionally to the available space. As a consequence, at each time step part of the corals change from a size class to another. At each time step all small corals either die or go to the medium size class, and we get a bundle of new small corals (called settlers). The number of settlers and how much the corals grow every time step are proportional to the k-area (which is the total area available for corals to live). For each size class there is a background mortality. The coral cover correspondent to each size class at each time step is given by:

- Small size class:
    - Settlers cover (currently this is a constant number weighted by the k-area, but in the future it should also be proportional to the number of medium and large coral covers)
- Medium size class:
    - Fraction of small cover that migrated to the medium size class after growing
    - Fraction of medium cover not dead that didn't migrate to the large size class, after growing
- Large size class
    - Fraction of medium cover that migrated to the large size class after growing
    - Fraction of large cover not dead, after growing

In the next section I present how to calculate the cover for the corals within each size class.

## Coral cover  calculation

We treat each coral as a perfect circumference, so the area $A(d)$ of a coral with diameter $d$ is given by:

$A(d) = \frac{\pi d^2}{4}$

If there are $N_s > 0$ corals within a size class $s$, then the density of corals within that size class is given by:

$\lambda_s = \frac{N_s}{\Delta d_s}$

Where $\Delta d_s = d_s - d_{s-1}$.

The distribution of areas per diameter for the size class $s$ is given, then, by:

$\Gamma_s(d) = A(d) \lambda_s = \frac{\pi d^2}{4} \frac{N_s}{\Delta d_s}$

<!-- ![image](https://latex.codecogs.com/gif.image?\int^{\infty}_{0}) -->

<!-- ![image](https://github.com/Zapiano/Coral-Cover-Dynamic-Toy-Model/assets/8040719/6e2601ea-e574-443a-95e4-5d79efc72b36) -->

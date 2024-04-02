# Coral Cover Dynamic Toy Model

This is a toy model to represent a Reef's coral cover change in time over the years. Here the corals are divided into three size classes $s = 1,2,3$ (representing small, medium and large), according to their diameter. Small corals have a diameter $d$ between $d_0 = 0$ and $d_1$; medium ones have a diameter between $d_1$ and $d_2$; and large ones have a diameter between $d_2$ and $d_3$. Within each size class, the corals are supposed to have a uniform diameter distribution. We treat each coral as a perfect circumference, so that the area $a(d)$ of a coral with diameter $d$ is given by:

$a(d) = \frac{\pi d^2}{4} \qquad\qquad\qquad\qquad (1)$

The diameter increase at each time step is represented by a parameter $\Delta d = l \cdot k$, where $l$ is called *linear extension* (the base growth of that coral' size class) and $k$ is the *k-area*, which is the total area available for corals to live. The summed change in diameter of all corals within all size classes causes an increase in the total coral cover and also makes a fraction of the corals within small and medium size classes to move to the next size classes.

## The Dynamics

In this model what changes over time is an area (coral cover cover), and not the number of individuals. The maximum available area is given by $k$ and we start with a certain amount of area within each size class. There are three phases within a time step iteration: a mortality phase; a growth phase; and a settlement phase. In the mortality phase, a fraction of each size class area is removed. In the growth phase the corals grow in diameter and, as a consequence, part of that area changes to another size class (what can cause the area of each size class to decrease in the end), resulting in a total area change. All small corals either change to the medium size class or die in that phase. In the settlement phase we get a bundle of new small corals (called settlers) that is proportional to $k$.

The coral cover correspondent to each size class at the end of each time step is given by:

- Small size class:
    - Settlers cover (currently this is a constant weighted by $k$, but in the future it should also be proportional to the number of medium and large coral covers)
- Medium size class:
    - Fraction of the survival small cover that grew and migrated to the medium size
    - Fraction of the survival medium cover that grew and didn't migrate to the large size class
- Large size class
    - Fraction of the survival medium cover that grew and migrated to the large size class
    - Fraction of the survival large cover that grew

In the next section I present how to calculate the cover for the corals within each size class.

## Coral cover  calculation

Suppose there are $N_s(t) > 0$ corals within a size class $s$, at time step $t$, and that they are uniformly distributed among the diameter space, then the density of corals per diameter within that size class is given by:

$\lambda_s(t) = \frac{N_s(t)}{\Delta d_s} \qquad\qquad\qquad\qquad (2)$

Where $\Delta d_s = d_s - d_{s-1}$ is the bin size of that size class. Note that, although the number of corals change in time, the bin sizes don't.

The distribution of area per diameter for the size class $s$, at the time step $t$, is given, then, by:

$\Gamma_s(t; d) = a(d) \lambda_s(t) \\
= \frac{\pi d^2}{4} \frac{N_s(t)}{\Delta d_s} \qquad\qquad\qquad\qquad\quad (3)$

One thing to note here is that this distribution has dimension of `Area/diameter`. The cover $C_s^{d_i, d_f}(t)$ correspondent to the corals within size class $s$ whose diameter is between $d_i$ and $d_f$, where $d_{s-1} \leq d_i \leq d_f \leq d_s$ is given by:

$C_s^{d_i, d_f}(t) = \int_{d_i}^{d_f} \Gamma_s(t;x) dx\\
= \int_{d_i}^{d_f} \frac{\pi x^2}{4} \frac{N_s(t)}{\Delta d_s} dx\\
= \frac{\pi N_s(t)}{4 \Delta d_s} \int_{d_i}^{d_f} x^2dx \Rightarrow \\$

$C_s^{d_i, d_f}(t) = \frac{N_s(t)}{\Delta d_s} \cdot \frac{\pi}{12} (d_f^3 - d_i^3) \qquad\qquad\qquad\qquad (4)$

> Quick note: I am using $x$ as integration variable instead of $d$ to avoid writing $dd$ ("infinitesimal delta d").

Since the cover $C_s(t) \equiv C_s^{d_{s-1}, d_s}$ of size class $s$ at time $t$ is known, we can use this equation with integration limits going from $d_{s-1}$ to $d_{s}$ to find an expression for $N_s(t)$:

$C_s(t) = \frac{N_s(t)}{\Delta d_s} \cdot \frac{\pi}{12} (d_{s-1}^3 - d_{s}^3) \Rightarrow$

$\boxed{\frac{N_s(t)}{\Delta d_s} = C_s(t) \cdot \frac{12}{\pi} \cdot \frac{1}{(d_{s-1}^3 - d_{s}^3)}} \qquad\qquad\qquad\qquad (5)$

And using eq. (5) in (4) we have:

$C_s^{d_i, d_f}(t) = C_s(t) \cdot \frac{12}{\pi} \cdot \frac{1}{(d_{s-1}^3 - d_{s}^3)} \cdot \frac{\pi }{12} (d_f^3 - d_i^3) \Rightarrow$

$\boxed{C_s^{d_i, d_f}(t) = C_s(t)  \cdot \frac{d_f^3 - d_i^3}{
d_{s-1}^3 - d_{s}^3}}  \qquad\qquad\qquad\qquad  (6)$

Now we can use eq. (6) to find the cover correspondent to the corals with diameter constrained between $d_i$ and $d_f$. Next we are going to use this equation to find how much the cover changes in a time step and which fraction of that cover goes to the next size class.

<!-- ![image](https://latex.codecogs.com/gif.image?\int^{\infty}_{0}) -->

<!-- ![image](https://github.com/Zapiano/Coral-Cover-Dynamic-Toy-Model/assets/8040719/6e2601ea-e574-443a-95e4-5d79efc72b36) -->

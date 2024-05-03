# Quick Start

The following is for running the model for 75 timesteps.

```
import DynamicCoralCoverModel: blocks_model

cover = blocks_model.run_model(75)
```

# Coral Cover Dynamic Toy Model

This model to represents a Reef's coral cover change in time over the years. It allows for the representation of distinct functional groups, indexed $f=1,2,...$. A single coral's area is approximated by the area of a circumference with diameter $x$:

$$
a(x) = \frac{\pi x^2}{4} \tag{1}
$$

Hence, when we talk about the size of a coral it means its diameter. The sum of the area of all individual corals within a particular functional group corresponds to that functional group's coral cover.

At each time step, a fraction of the corals die and the survivals grow, and we have a batch of new corals (settlers) due to reproduction. The background mortality and growth rates depend on the size of the corals and functional group. To represent this difference, we split the space of possible diameters for each functional group $f$ into $s$ subspaces $\sigma_{f,1},...,\sigma_{f,s}$ so all corals with diameter within the same subspace have the same ecological behavior. A subspace $\sigma_{f,s}$ is an interval (in the diameter space) constrained by a lower bound $\sigma_{f,s}^-$ and a upper bound $\sigma_{f,s}^+$ or, in other words, $\sigma_{f,s} = [\sigma_{f,s}^-, \sigma_{f,s}^+]$. As a consequence, $\sigma_{f,s}^+ = \sigma_{f,s+1}^-$.

For each functional group and diameter subspace we have a base value for the increment in the diameter size, called *linear extension* and represented by $l_{f,s}$. The increase in diameter for a single coral, at each time step, is given by

$$
D_{f,s}(t) = k(t) \cdot l_{f,s} \tag{2}
$$

where

$$
k(t) = 1 - C(t)/k_{max} \tag{3}
$$

is the available space at time $t$. $k_{max}$ is the carrying capacity and $C(t)$ is the total coral cover at time $t$ considering all functional groups and diameters subspaces.

In the next section we will learn how to calculate the cover for the corals within each diameter subspace so we can build the dynamic equations.

## Coral cover  calculation

Suppose there are $N_{f,s}(t) > 0$ corals of functional group $f$ within some diameter subspace $s$, at time step $t$, and that they are uniformly distributed along that subspace. The coral density per diameter is given by:

$$
\lambda_{f,s}(t) = \frac{N_{f,s}(t)}{\Delta \sigma_{f,s}} \tag{4}
$$

Where $\Delta \sigma_{f,s} = \Delta \sigma_{f,s}^+ - \Delta \sigma_{f,s}^-$ is the bin size of the diameter subspace $\sigma_{f,s}$.

The distribution of area per diameter for the functional group $f$ and subspace $s$, at the time step $t$, is given, by:

$$
\begin{align*}
\Gamma_{f,s}(t; x) =& a(x) \lambda_{f,s}(t) \\
=& \frac{\pi x^2}{4} \frac{N_{f,s}(t)}{\Delta \sigma_{f,s}} \tag{5}
\end{align*}
$$

One thing to note here is that this distribution has dimension of `Area/diameter`[^1].

[^1]: Since the area and diameter can be measured in `m^2`/`m`, `cm^2`/`cm`, etc, the actual dimension here should be based on whatever units you use to measure those. But I decided to leave it this way to emphasize that `Area` is measured in the "physical" space, whereas `diameter` is measured in the diameter space.

The cover $C_{f,s}(t;a,b)$ correspondent to the corals from the functional group $f$ within the interval $[a, b]$, where $[a,b] \subset \sigma_{f,s}$ is given by:

$$
\begin{aligned}
C_{f,s}(t; a,b) =& \int_{a}^{b} \Gamma_{f,s}(t; x) dx \\
=& \int_{a}^{b} \frac{\pi x^2}{4} \lambda_{f,s}(t) dx \\
=& \lambda_{f,s}(t)\frac{\pi }{4} \int_{a}^{b} x^2dx \Rightarrow \\
\end{aligned}
$$

$$
\boxed{C_{f,s}(t; a,b) = \lambda_{f,s}(t) \frac{\pi}{12} (b^3 - a^3)} \tag{6}
$$

Since the cover $C_{f,s}(0) \equiv C_{f,s}(0; \sigma_{f,s}^-, \sigma_{f,s}^+)$ of all corals from functional group $f$ with diameter within the diameter subspace $s$ at time $0$ is known, we can use eq. (6) to find an expression for the initial densities $\lambda_{f,s}(0)$:

$$
C_{f,s}(0) = \lambda_{f,s}(0) \cdot \frac{\pi}{12} \left[(\sigma_{f,s}^+)^3 - (\sigma_{f,s}^-)^3\right] \Rightarrow
$$

$$
\boxed{\lambda_{f,s}(0) = C_{f,s}(0) \cdot \frac{12}{\pi} \cdot \frac{1}{(\sigma_{f,s}^+)^3 - (\sigma_{f,s}^-)^3}} \tag{7}
$$

<!--
And using eq. (5) in (4) we have:

$$
C_s(t; d_i, d_f) = C_s(t) \cdot \frac{12}{\pi} \cdot \frac{1}{(d_{s+1}^3 - d_{s}^3)} \cdot \frac{\pi }{12} (d_f^3 - d_i^3) \Rightarrow
$$

$$
\boxed{C_s(t; d_i, d_f) = C_s(t)  \cdot \frac{d_f^3 - d_i^3}{d_{s+1}^3 - d_{s}^3}}  \tag{8}
$$

Eq. (8) gives cover correspondent to the corals with diameter between $d_i$ and $d_f$. Note that, when $d_i = d_{s}$ and $d_f = d_{s+1}$, we have $C_s(t; d_{s}, d_{s+1}) = C_s(t)$ as expected; and when $d_i = d_f$ the cover goes to zero. Next we are going to use this equation to find the change in cover for each size class at each time step due to growth and size class changing.
-->

> **_NOTE:_** This documentation from this point up to the end still need to be updated. Don't trust anything you read from now on.

## The Dynamics

In this model what changes over time is an area (coral cover cover), and not the number of individuals. There are three phases within a time step iteration: a *mortality phase*; a *growth phase*; and a *settlement phase*.

In the *mortality phase*, a fraction $m_s(t)$ of each size class $s$ is removed. The mortality affects all corals within a size class equally, so it **doesn't change the size class bin sizes/limits**, it just changes the coral density. Therefore, if we want to find the cover between $d_i$ and $d_f$ **after** the mortality phase, we would just have to multiply eq. (8) by $(1 - m_s)$.

In the *growth phase* the corals grow in diameter and, as a consequence, part of the cover changes to another size class, what can cause a size class cover to decrease while the total cover increases. All small corals either change to the medium size class or die in that phase.

In the *settlement phase* we get a bundle of new small corals (called settlers):

$$
\zeta(t) = r \cdot k(t) \tag{9}
$$

 where $r$ is a constant - although in the future this should be proportional to the number of medium and large covers.

Before we build the cover dynamic equations, it's important to note that the small and large size classes have a special behavior. If we add more size classes, all of them, except the smallest and largest ones, would behave exactly like the medium size we present here. Next we build each size class equation separately before putting it all together, highlighting inside the box the meaning of each term in that equation.

### Small size class
All small size class corals either die or migrate to the next step, so the only coral remaining in that size class at the end of the iteration are the settlers. The dynamics is given by:

> - $\zeta(t)$ - Settlers cover (currently this is a constant weighted by $k$, but )

$$
C_1(t+1) = \zeta(t) \tag{10.1}
$$

### Medium size class
There are two factors that affect this size class: internal grow (growth of the corals within size class 2 that don't migrate to the large size class) and migration of the small size class after growing. The cover after the internal grow is given by:

$$
\lambda_2 \int_{d_2}^{d_3 - D_2(t)}  \frac{\pi \Big(x + D_2(t)\Big)^2}{4} dx
$$

Which is the base of eq. (6) but with diameter given by $x + D_2(t)$ instead of $x$, since their diameter has increased. We can do a variable change $x + D_2(t) = y$ to rewrite this equation as:

$$
\lambda_2 \int_{d_2 + D_2(t)}^{d_3}  \frac{\pi y^2}{4} dy
$$

But now is easy to see that this is the same equation that lead us to eq. (8), only with intervals $[d_2 + D_2(t), d_3]$, or, in other words, $C_2(t; d_2 + D_2(t), d_3)$.

By the same logic, we can write the fraction of small cover that grew and migrated to the medium size class as:

$$
\lambda_1 \int_{d_2 - D_1(t)}^{d_2}  \frac{\pi \Big(x + D_1(t)\Big)^2}{4} dx
$$

Which, with the variable change $x + D_1(t) = z$ becomes:

$$
\lambda_1 \int_{d_2}^{d_2 + D_1(t)}  \frac{\pi z^2}{4} dz
$$

Which, again, can be written using eq. (8) as $C_1(t; d_2, d_2 + D_1(t))$. So, to summarize, we have:

>- $C_2(t; d_2 + D_2(t), d_3) \cdot (1 - m_2)$ - Fraction of the medium cover survivals that grew and didn't migrate to the large size class
>- $C_1(t; d_2, d_2 + D_1(t)) \cdot (1 - m_1)$ - Fraction of the small cover survivals that grew and migrated to the medium size

$$
\begin{align*}
C_2(t+1) =& C_2(t; d_2 + D_2(t), d_3) \cdot (1 - m_2) + C_1(t; d_2, d_2 + D_1(t)) \cdot (1 - m_1) \\
=& C_2(t) \cdot (1 - m_2) \cdot \frac{d_3^3 - (d_2 + D_2(t))^3}{d_{3}^3 - d_{2}^3} + C_1(t) \cdot (1 - m_1) \cdot \frac{(d_2 + D_1(t))^3 - d_2^3}{d_{2}^3 - d_{1}^3} \tag{10.2}
\end{align*}
$$

### Large size class
The corals within the large size class have all the same area and are not allowed to grow (otherwise we would have infinite growth). To get the cover correspondent to the corals that go from the medium size class to the large, multiply the number of survival corals within size class 2 that are inside the migration range by the area of a coral with maximum diameter $d_{max}$:

$$
\begin{align*}
&N_2(t) \cdot (1- m_2) \cdot \frac{d_2 - (d_2 - D_2(t))}{\Delta d_2} \cdot \frac{\pi}{4}d_{max}^2 \\
&= C_2(t) \cdot (1- m_2) \cdot \frac{12}{\pi} \frac{1}{d_{2}^3 - d_{1}^3} \Delta d_2 \cdot \frac{D_2(t)}{\Delta d_2} \cdot \frac{\pi}{4} d_{max}^2 \\
&= C_2(t) \cdot (1- m_2) \cdot \frac{3 D_2(t) \cdot d^2_{max}}{d_{2}^3 - d_{1}^3}
\end{align*}
$$

The dynamics is given by:

> - $C_3(t) \cdot (1 - m_3)$ - Fraction of large cover that survived
> - $C_2(t) \cdot (1 - m_2) \cdot\frac{3 D_2(t) \cdot d^2_{max}}{d_{3}^3 - d_{2}^3}$ - Fraction of medium cover survivals that grew and migrated to the large size class


$$
C_3(t+1) = C_3(t) \cdot (1 - m_3) + C_2(t) \cdot (1 - m_2) \cdot \frac{3 D_2(t) \cdot d^2_{max}}{(d_{3}^3 - d_{2}^3)} \tag{10.3}
$$

### Putting it all together

The final cover dynamics equations are, then:

$$
\begin{cases}
C_1(t+1) = \zeta(t)\\
C_2(t+1) = C_2(t) \cdot (1 - m_2) \cdot \frac{d_3^3 - (d_2 + D_2(t))^3}{d_{3}^3 - d_{2}^3} + C_1(t) \cdot (1 - m_1) \cdot \frac{(d_2 + D_1(t))^3 - d_2^3}{d_{2}^3 - d_{1}^3} \\
C_3(t+1) = C_3(t) \cdot (1 - m_3) + C_2(t) \cdot (1 - m_2) \cdot \frac{3 D_2(t) \cdot d^2_{max}}{d_{3}^3 - d_{2}^3} \tag{11}
\end{cases}
$$

A one important things to note is the, when $k(t) = 0$ (meaning that there is no available space), we have $Z(t) = 0$ and $D_s(t) = 0$ for all $s$. In other words, no settlers and no diameter growth. As a consequence, the equations above become:

$$
\begin{cases}
C_1(t+1) = 0\\
C_2(t+1) = C_2(t) \cdot (1 - m_2)\\
C_3(t+1) = C_3(t) \cdot (1 - m_3)
\end{cases}
$$

That means that we won't have any new coral entering the system and also no cover increase, only background mortality.

## Results

The results below are from the most recent version of this model with two species and 4 size classes. All the logic explained in this document remains the same.

![Total coral cover](/figures/cover_species_sizeclasses.png)

![Coral cover by size class](/figures/cover_species_tot.png)

# Coral Cover Dynamic Toy Model

This is a toy model to represent a Reef's coral cover change in time over the years. Here the corals are divided into three size classes $s = 1,2,3$ (representing small, medium and large), according to their diameter. Small corals have a diameter $d$ between $d_0 = 0$ and $d_1$; medium ones have a diameter between $d_1$ and $d_2$; and large ones have a diameter between $d_2$ and $d_3$. Within each size class, the corals are supposed to have a uniform diameter distribution. We treat each coral as a perfect circumference, so that the area $a(d)$ of a coral with diameter $d$ is given by:

$a(d) = \frac{\pi d^2}{4} \qquad\qquad\qquad\qquad (1)$

The diameter increase at each time step is represented by a parameter $D_s(t) = l_s \cdot k(t)$, where $l_s$ is called *linear extension* of the size class $s$ (the base growth in diameter for corals within size class $s$) and $k(t) = k_{max} - C(t)$ is the *k-area*, which is the available area for corals to live at time $t$. Here, $k_max$ is the maximum possible cover and $C(t)$ is the total coral cover at time $t$.

The change in diameter at each time step causes an increase in the total coral cover and also makes a fraction of the corals within small and medium size classes to move to medium and large size classes. In the next section we will learn how to calculate the cover for the corals within each size class so we can build the cover dynamic equations.

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

$\boxed{\frac{N_s(t)}{\Delta d_s} = C_s(t) \cdot \frac{12}{\pi} \cdot \frac{1}{(d_{s}^3 - d_{s-1}^3)}} \qquad\qquad\qquad\qquad (5)$

And using eq. (5) in (4) we have:

$C_s^{d_i, d_f}(t) = C_s(t) \cdot \frac{12}{\pi} \cdot \frac{1}{(d_{s}^3 - d_{s-1}^3)} \cdot \frac{\pi }{12} (d_f^3 - d_i^3) \Rightarrow$

$\boxed{C_s^{d_i, d_f}(t) = C_s(t)  \cdot \frac{d_f^3 - d_i^3}{
d_{s}^3 - d_{s-1}^3}}  \qquad\qquad\qquad\qquad  (6)$

Now we can use eq. (6) to find the cover correspondent to the corals with diameter between $d_i$ and $d_f$. Note that, when $d_i = d_{s-1}$ and $d_f = d_s$, we have $C_s^{d_{s-1}, d_s}(t) = C_s(t)$ as expected; and when $d_i = d_f$ the cover goes to zero.

Next we are going to use equation (6) to find the change in cover for each size class at each time step due to growth and size class changing.

## The Dynamics

In this model what changes over time is an area (coral cover cover), and not the number of individuals. There are three phases within a time step iteration: a *mortality phase*; a *growth phase*; and a *settlement phase*.

In the *mortality phase*, a fraction $m_s(t)$ of each size class $s$ is removed. The mortality affects all corals within a size class equally, so it **doesn't change the size class bin sizes/limits**, it just changes the coral density. Therefore, if we want to find the cover between $d_i$ and $d_f$ **after** the mortality phase, we would just have to multiply eq. (6) by $(1 - m_s)$.

In the *growth phase* the corals grow in diameter and, as a consequence, part of the cover changes to another size class, what can cause a size class cover to decrease while the total cover increases. All small corals either change to the medium size class or die in that phase.

In the *settlement phase* we get a bundle $\Zeta(t) = \zeta \cdot k(t)$ of new small corals (called settlers) where $\zeta$ is a constant - although in the future this should be proportional to the number of medium and large covers.

Before we build the cover dynamic equations, it's important to note that the small and large size classes have a special behavior. If we add more size classes, all of them, except the smallest and largest ones, would behave exactly like the medium size we present here. Next we build each size class equation separately before putting it all together, highlighting inside the box the meaning of each term in that equation.

### Small size class
All small size class corals either die or migrate to the next step, so the only coral remaining in that size class at the end of the iteration are the settlers. The dynamics is given by:

> - $\Zeta(t)$ - Settlers cover (currently this is a constant weighted by $k$, but )

\
$C_1(t+1) = \Zeta(t)$

### Medium size class
The dynamics is given by:

> - $C_1^{d_1-D_1(t), d_1}(t) \cdot (1 - m_1)$ - Fraction of the small cover survivals that grew and migrated to the medium size
> - $C_2^{d_1, d_2 - D_2(t)}(t) \cdot (1 - m_2)$ - Fraction of the medium cover survivals that grew and didn't migrate to the large size class

\
$C_2(t+1) = C_1^{d_1-D_1(t), d_1}(t) \cdot (1 - m_1) + C_2^{d_1, d_2 - D_2(t)}(t) \cdot (1 - m_2)\\
= C_1(t) \cdot \frac{d_{d_1}^3 - d_{d_1-D_1(t)}^3}{d_{1}^3 - d_{0}^3} \cdot (1 - m_1) + C_2(t) \cdot \frac{d_{d_2 - D_2(t)}^3 - d_{d_1}^3}{d_{2}^3 - d_{1}^3} \cdot (1 - m_2)\\$

### Large size class
The corals within the large size class have all the same area and are not allowed to grow (otherwise we would have infinite growth). To get the cover correspondent to the corals that go from the medium size class to the large, multiply the number of survival corals within size class 2 that are inside the migration range by the area of a coral with maximum diameter $d_{max}$:

$N_2(t) \cdot (1- m_2) \cdot \frac{d_2 - (d_2 - D_2(t))}{\Delta d_2} \cdot \frac{\pi}{4}d_{max}^2 \\
= C_2(t) \cdot (1- m_2) \cdot \frac{12}{\pi} \frac{1}{(d_{2}^3 - d_{1}^3)} \Delta d_2 \cdot \frac{D_2(t)}{\Delta d_2} \cdot \frac{\pi}{4} d_{max}^2\\
= C_2(t) \cdot (1- m_2) \cdot \frac{3 D_2(t) \cdot d^2_{max}}{(d_{2}^3 - d_{1}^3)} \qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad (7)$

The dynamics is given by:

> - $C_2(t) \cdot (1 - m_2) \cdot\frac{3 D_2(t) \cdot d^2_{max}}{(d_{2}^3 - d_{1}^3)}$ - Fraction of medium cover survivals that grew and migrated to the large size class
> - $C_3(t) \cdot (1 - m_3)$ - Fraction of large cover survivals that grew

\
$C_3(t+1) = C_2(t) \cdot (1 - m_2) \cdot \frac{3 D_2(t) \cdot d^2_{max}}{(d_{2}^3 - d_{1}^3)} + C_3(t)(1 - m_3)\\$

### Putting it all together

The final cover dynamics equations are, then:

$
    C_1(t+1) = \Zeta(t)\\
    C_2(t+1) = C_1(t) \cdot \frac{d_{d_1}^3 - d_{d_1-D_1(t)}^3}{d_{1}^3 - d_{0}^3} \cdot (1 - m_1) + C_2(t) \cdot \frac{d_{d_2 - D_2(t)}^3 - d_{d_1}^3}{d_{2}^3 - d_{1}^3} \cdot (1 - m_2) \\
    C_3(t+1) = C_2(t) \cdot (1 - m_2) \cdot \frac{3 D_2(t) \cdot d^2_{max}}{(d_{2}^3 - d_{1}^3)} + C_3(t)(1 - m_3)
$

A one important things to note is the, when $k(t) = 0$ (meaning that there is no available space), we have $Z(t) = 0$ and $D_s(t) = 0$ for all $s$. In other words, no settlers and no diameter growth. As a consequence, the equations above become:

$
    C_1(t+1) = 0\\
    C_2(t+1) = C_2(t) \cdot (1 - m_2)\\
    C_3(t+1) = C_3(t) \cdot (1 - m_3)
$

That means that we won't have any new coral entering the system and also no cover increase, only background mortality.

## Results

[To Do]

<!-- ![image](https://latex.codecogs.com/gif.image?\int^{\infty}_{0}) -->

<!-- ![image](https://github.com/Zapiano/Coral-Cover-Dynamic-Toy-Model/assets/8040719/6e2601ea-e574-443a-95e4-5d79efc72b36) -->

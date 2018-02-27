---
  indent: true
---
# Derivation and estimation of a fluorescence calibration factor
Griffin Chure - Feb. 20, 2018

## Derivation
### Relating cellular intensity to fluorophore copy number
  Imagine that we have a cell that has a fixed number of fluorescent proteins.
As this cell divides, these proteins will be Binomially partitioned into the
two daughter cells with some probability $p$ dictating the proteins' likelihood
to end up in one cell over the other. By examining the difference in
fluorescence between the two daughters, we can determine just how bright a
single molecule is on average.

We begin by positing that the fluorescence is conserved from mother to
daughters, where we assume that production and degradation rates of the
fluorophore are negligible. Mathematically, we can state this as
$$
I_\text{tot} = I_1 + I_2
$$ {#eq:fluo_cons}
where $I_\text{tot}$ is the total fluorescence of the mother cell and $I_1$ and
$I_2$ are the intensities of the two daughter cells. If we assume that the
measured fluorescence arises only from the fluorophores, we say that the
total intensity of the mother cell should be proportional to the number of
fluorescent proteins $N_\text{tot}$,
$$
I_\text{tot} = \alpha N_\text{tot}.
$${#eq:antot}
This can be similarly stated for each individual daughter cell,
$$
I_1 = \alpha N_1\,\, ; \,\, I_2 = \alpha N_2
$${#eq:in1}
which follows the assumption that the protein copy number is conserved as well,

$$
N_\text{tot} = N_1 + N_2
$${#eq:ntotn1n2}
Upon a division event, the proteins tumbling about the cytoplasm  are partitioned
between the two daughter cells. We can write the probability distribution $g(n)$
of finding $n$ proteins in one cell or another as

$$
g(n | N_\text{tot}, p) = {N_\text{tot}! \over n! (N_\text{tot}-n)!}p^n(1-p)^{N_\text{tot} - n},
$${#eq:binomial}
where $p$ is the probability of a protein being partitioned into one cell versus the other.
Using our assumption in [@eq:ntotn1n2] that the total protein copy number is conserved,
we should be able to predict the average square fluctuation in intensity between two
daughter cells,

$$
\langle (I_1 - I_2)^2 = \langle (2I_1 - I_\text{tot})^2 \rangle.
$${#eq:simple_fluct}
Using [@eq:in1] and expanding the right-hand side of [@eq:simple_fluct], we arrive at

$$
\langle (I_1 - I_2)^2 \rangle = 4\alpha^2\langle n^2 \rangle - 4\alpha^2\langle n \rangle N_\text{tot} + \alpha^2N_\text{tot}^2.
$${#eq:expanded_fluct}
To beat this into a more manageable form, we must know the first two moments of
the Binomial distribution. While we likely have these committed to memory, it's
a useful exercise to derive them explicitly.

### The moments of the Binomial distribution
Odds are that when you have a probability distribution in your hands, you will
want to know its moments, at least up to the second order. One of the most
useful tricks to have up your sleeve for this type of calculation is the
moment generating function. For a generic probability distribution $g(x)$,
the moment generating function is defined as
$$
M_x(t) = \mathbb{E}e^{xt}
$${#eq:def_mgf}
Where $\mathbb{E}$ is the expectation value and $t \in \mathbb{R}$. For the
Binomial distribution ([@eq:binomial]), we can compute the moment generating
function as
$$
\mathbb{E}e^{nt} = \sum\limits_{i=0}^{N_\text{tot}}e^{nt}{N_\text{tot}! \over n_i!(N_\text{tot} - n_i)!}p^{n}(1-p)^{N_\text{tot} - n}.
$${#eq:binom_mgf_exp}
Applying the Binomial theorem to [@eq:binom_mgf_exp] allows for a beautiful
simplification to
$$
\mathbb{E}e^{nt} = (1 - p + pe^t)^{N_\text{tot}} = M_x(t),
$${#eq:binom_mgf}
which is the moment generating function. In general, you can find whichever
moment you please by differentiating this function and evaluating at $t=0$,
$$
m_n = {d^nM_x \over dt^n} \biggr\rvert_{t = 0}.
$${#eq:moment_n}
With this handy trick, we can compute the first moment of the Binomial distribution
by computing the derivative,
$$ \langle n \rangle = {d\over dt}(1 - p + pe^t)^{N_\text{tot}}\biggr\rvert_{t=0} =N_\text{tot}pe^t( 1 - p + pe^t)^{N_\text{tot} - 1},
$${#eq:first_moment_deriv}
and evaluating at $t=0$,

$$
\langle n \rangle = N_\text{tot}p.
$$
This shouldn't be too surprising as this is what we would have said as a gut
reaction.
Similarly, we can compute the second moment $\langle n^2 \rangle$ by differentiating,
$$
\langle n^2 \rangle = {d^2 M_n\over dt^2} = N_\text{tot}pe^t(1 - p + pe^t)^{N_\text{tot} - 2}(1 - p + N_\text{tot}pe^t),
$${#eq:second_moment_deriv}
and evaluating at $t=0$ to arrive at
$$
\langle n^2 \rangle = N_\text{tot}p(1 - p).
$${#eq:binom_second_moment}

### Tying it all together
Now that we have the moments of the Binomial distribution in hand, we can return
to our primary goal. Thus far, we have always kept $p$ around in our calculations.
In the situation where we are completely ignorant regarding the localization
of our protein of interest or have reasons to believe that the partitioning is
actively driven by some process, it is fair to keep it around as a parameter
in our model. However, we can simplify our calculations here by assuming that
we know the partitioning of our protein of interest is fair (or at least very
close to fair), meaning that $p= 1/2$. We can use this assumption to condense
our moments to
$$
\langle n \rangle = {N_\text{tot} \over 2}
$${#eq:simp_mom_one}
and
$$
\langle n^2 \rangle = {2N_\text{tot} - N_\text{tot} \over 4}.
$${#eq:simp_mom_two}

We can now plug these moments into [@eq:expanded_fluct] generating
$$
\langle (I_1 - I_2)^2 \rangle = 4\alpha^2{2N_\text{tot} - N_\text{tot} \over 4} - 4\alpha^2 {N_\text{tot}^2 \over 2} + \alpha^2N_\text{tot}^2.
$${#eq:plugged_fluct}
Some simplification delivers us to a beautiful result,
$$
\langle(I_1 - I_2)^2 \rangle = \alpha^2N_\text{tot} = \alpha I_\text{tot}.
$${#eq:golden_rule}

This incredibly simple prediction allows us to turn fluorescence counts typically
reported as arbitrary units (a.u.) to copy number, which has an obvious connection
to the physics of the system. To show the power of this system, this process of
partitioning proteins through a series of cell divisions is shown in [@fig:dilution_sim]. Here, we have examined the partitioning of proteins
from ten copies per cell to one thousand. Each simulation of this parititioning
was repeated one hundred times, allowing us to compute the average fluctuation,
which is what [@eq:golden_rule] predicts. These averages are shown in [@fig:dilution_sim]
as red points. The blue line represents the prediction of [@eq:golden_rule].

![**Simulated protein partitioning with perfect measurement**.
 Intensities were calculated by multiplying the number of partitioned proteins by a set
calibration factor $\alpha = 150$ a.u. per molecule. Black points correspond to
individual simulations, red points are the averages over all simulations, and
the blue line is the prediction given in
[@eq:golden_rule].)](../figs/simulated_dilution_simple.pdf){#fig:dilution_sim}


## Practical estimation of $\alpha$
With [@eq:golden_rule] at our disposal, we must develop some scheme for estimating the value of $\alpha$ from a given data set. Unlike in our simulations, we cannot break down each cell into bins of a single copy number. One approach is to break the data set up into discrete bins of a certain number of events, compute the necessary statistics, and then perform a linear regression for $\alpha$ on the binned data. While this is a completely valid approach, it requires an arbitrary decision of how many events to choose per bin. To ensure the binning scheme is valid, you must perform the parameter estimation over a range of bin widths and choose one in which the estimation appears to converge on one value.

Rather than binning, we have taken a Bayesian approach in which each individual division is treated as an individual experiment. This approach completely removes the requirement of binning in exchange for more complexity.

### A Bayesian approach
We would like to estimate $\alpha$, the calibration factor, given a set of sister cell intensities $[I_1, I_2]$. For the following derivation, we will neglect the fact that we have a set of $I_1$ and $I_2$ values and will consider only a single division event. We can begin our derivation by stating Bayes' theorem,
$$
g(\alpha\vert I_1, I_2) = {f(I_1, I_2, \vert \alpha)g(\alpha) \over f(I_1, I_2)},
$${#eq:bayes_rule_complete}
where $f$ represents probability distributions over observed variables and $g$
represents those over unobserved parameters. As we are not particularly
interested in the proper normalization of the posterior distribution, we can
simply neglect the denominator of [@eq:bayes_rule_complete] and say
$$
g(\alpha\vert I_1, I_2) \propto f(I_1, I_2\vert\alpha)g(\alpha) = g(\alpha\vert I_1, I_2).
$${#eq:bayes_rule}

As the value of $I_1$ is dependent on $I_2$ through [@eq:fluo_cons], the likelihood can be broken down to
$$
f(I_1, I_2\vert \alpha) = f(I_1\vert I_2, \alpha)g(I_2),
$${#eq:broken_like}
where we have introduced $g(I_2)$ as another prior. To be completely fair and uninformative, we can assume that these parameters can take on any value, so long as it is positive. Mathematically, this means that they are constant and can be dropped from our formulation of the posterior. This makes our posterior distribution pretty convenient to handle,
$$
g(\alpha \vert I_1, I_2) = f(I_1\vert \alpha, I_2).
$${#eq:simplified_posterior}
Since we are assuming that fluorescence is conserved, we must say that the total protein copy number is conserved, as is stated in [@eq:ntotn1n2]. We can use this assumption to our advantage and express our likelihood  in terms of protein copy number $N$. This can be achieved through change of variables,
$$
f(I_1\vert\alpha, I_2) = f(N_1\vert \alpha, I_2)\biggl\lvert {dN_1 \over dI_1}\biggr\rvert = {1 \over \alpha}f(N_1\vert \alpha, I_2).
$$
Piecing the rest of the likelihood together is now a slice of cake. As we stated in [@eq:binomial], we expect the proteins to be partitioned between daughter cells in a Binomial manner. If we assume that the partitioning is fair (i.e. $p=1/2$), we can write our likelihood as
$$
f(I_1\vert \alpha, I_2) = {1 \over \alpha} {N_\text{tot}! \over N_1!(N_\text{tot} - N_1)}2^{-N_\text{tot}}.
$${#eq:binomial_revisited}
Since we've included the prefactor necessary to change variables, we are free to use [@eq:antot] to rewrite $N_1$ and $N_\text{tot}$ as
$$
f(I_1 \vert \alpha I_2) = {1 \over \alpha}{{I_1 + I_2 \over \alpha}! \over {I_1 \over \alpha}!{I_2 \over \alpha}!}2^{-(I_1 + I_2) / \alpha}.
$$
However, it doesn't make much sense to compute factorials of the intensity since it is a continuous quantity. We can make the approximation
$$
n! \approx n\Gamma(n) = \Gamma(n+1),
$${#eq:gamma_approx}
to write our likelihood as
$$
f(I_1\vert \alpha, I_2) = {1 \over \alpha}{\Gamma\left({I_1 + I_2 \over \alpha} + 1\right) \over \Gamma\left({I_1 \over \alpha} + 1\right)\Gamma\left({I_2 \over \alpha} + 1\right)}2^{-(I_1 + I_2) / \alpha}.
$${#eq:single_div_post}
As we have neglected the prior distributions for $I_2$ and $\alpha$, [@eq:single_div_post] happens to be our posterior distribution for a single division event. To generalize this to a set of $k$ divisions, we
can say
$$
g(\alpha\vert [I_1, I_2]) = {1 \over \alpha^k}\prod\limits_{i=0}^k{\Gamma\left({I_{1,i} + I_{2,i} \over \alpha} + 1\right) \over \Gamma\left({I_{1,i} \over \alpha} + 1\right)\Gamma\left({I_{2,i} \over \alpha} + 1 \right)}{2^{-(I_{1,i} + I_{2, i} / \alpha)}}.
$${#eq:posterior}

This result allows us to estimate the best-fit value for $\alpha$ without relying on any binning procedure. As this is a distribution containing only one parameter, it is trivial to apply to even large data sets through your favorite
optimization procedure. [@fig:param_estimation] shows the posterior distribution evaluated over a range of $\alpha$ values for the data shown in [@fig:dilution_sim]. The red dashed line is an approximation of the posterior as a Gaussian, allowing us to use the standard deviation as a measurement of the statistical error.

![**Posterior probability distribution for $\alpha$**. (A) The
  normalized posterior probability distribution $g(\alpha\, \vert\, [I_1, I_2])$
  from Eq. [@eq:posterior] for data shown in Fig. \ref{fig:dilution_sim}.
  The dark line and shaded region represent the actual value of the posterior
  distribution. The dashed red line represents a Gaussian approximation of the
  posterior. The true value of $\alpha$ is shown as a purple vertical line.
  The best-fit value for $\alpha$ in this data set is $149 \pm 1$
  a.u. per molecule.](../figs/alpha_simple_minimization.pdf){#fig:param_estimation}

---
  title: Sensitivity of fluorescent fluctuations to systematic noise
  author: Griffin Chure
  date: \today
  indent: true
  bibliography: mwc_growth.bib
---


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


## Practical estimation of $\alpha${#sec:bayesian_estimation}
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
optimization procedure.

![**Posterior probability distribution for $\alpha$**. (A) The
  normalized posterior probability distribution $g(\alpha\, \vert\, [I_1, I_2])$
  from Eq. [@eq:posterior] for data shown in Fig. \ref{fig:dilution_sim}.
  The dark line and shaded region represent the actual value of the posterior
  distribution. The dashed red line represents a Gaussian approximation of the
  posterior. The true value of $\alpha$ is shown as a purple vertical line.
  The best-fit value for $\alpha$ in this data set is $149 \pm 1$
  a.u. per molecule.](../figs/alpha_simple_minimization.pdf){#fig:param_estimation}

### Confirming the method
While this all makes sense on paper, it's a wise idea to confirm that it works when applied to
several datasets and perturbed in as many ways you can think of. In the sections that follow, we
will turn to the simulation shown in [@fig:dilution_sim] and test various aspects of this approach.

To test the validity of [@eq:posterior], we took the simulated data from [@fig:dilution_sim] and
found the most-likely value of the calibration factor through minimization. [@Fig:param_estimation]
shows the full posterior distribution evaluated over a broad range of calibration factor values,
revealing a unimodal distribution peaked around the seeded value of $\alpha$, 150 a.u. per molecule.
To get a measure of the statistical error in the determination of alpha, we can approximate
the posterior as a Gaussian distribution by computing the Hessian of the posterior. The Gaussian
approximation of this posterior is shown in [@fig:param_estimation] as a dashed red line. The final
best estimate for the calibration factor from the data shown in [@fig:dilution_sim] is 149
$\pm$  1 a.u., which agrees nicely with the true seeded value of 150 a.u..


To confirm that this method works for experimental data (and all of the small errors that have thus
far been neglected), we applied this method of estimation to data found in the literature where the
calibration factor had been estimated through linear regression. [@Fig:brewster_method_agreement] (A)
shows data taken from Brewster et al. [-@Brewster2014a] where they used this method to count
the number of fluorescently labeled transcription factors.  They determined the best-fit calibration
factor by binning the data arbitrarily and performing linear regression on the log of the
fluctuations, yielding a calibration factor of 156 $\pm$ 8 a.u. per molecule. Applying [@eq:posterior]
to this raw data results in a best-fit measurement of 169 $\pm$ 10 a.u. per molecule. The posterior
distribution for this data set is shown in [@fig:brewster_method_agreement] (B) where the vertical
line indicates the reported value. While the calibration factor determined via [@eq:posterior] is
different than that reported by Brewster et al. [-@Brewster2014a], it is rather close. The author's
estimation of this parameter showed to be dependent on the binning size used, varying by around 10 a.u. per molecule.


![**Comparison of estimation methods on a standard data set.** (A) Data from Brewster et al. [-@Brewster2014a] .
Black points are individual cell division events. Red points correspond to the mean value of 54
separate events. These binned points were used for linear regression to estimate the calibration
factor. [@Eq:golden_rule] evaluated using their calibration factor of 156 a.u. per molecule
is shown as a red line. [@Eq:golden_rule] evaluated using the most likely value of $\alpha$ defined
in [@eq:posterior] is shown as a dashed blue line. (B) The full posterior distribution for [@eq:posterior]
is shown in gray. The black vertical line represents the best-fit value for the calibratoin factor
presented in Brewster et al. [-@Brewster2014a].](../figs/brewster_method_agreement.pdf){#fig:brewster_method_agreement}




![**Calibration factor estimation under three different illumination sources.**
Representative examples of experimental measurements using the three different
sources of excitation illumination. Top-left and top-right correspond to experiments performed
with an Hg lamp source and a LED, respectivel. The bottom left corner was collected using a 589 nm laser source. The bottom right panel shows the fold-change in gene expression measurements from all three experiments shown in the other panels.](../figs/light_source_experiments.pdf){#fig:example_experiments}

## Including various flavors of error

Up to this point, we've imagined a scenario in which the protein partitioning
is completely random, proteins are infinitely stable, and there is no error
in the measurement of the cellular intensity. Reality, how, is often more lemon
than lemonade. Our actual measurements of these experiments will be rife with
error ranging from the stability of the light source to the biological
peculiarities. In this section, we will thoroughly dissect several models of
systematic measurement error. We will rely on numerical simulations of all
of the coming models to get a handle on the type of error. The procedure for
this approach is diagrammed in [@fig:simulation_schematic]. We will quantify
the error in the measurement of $\alpha$ when [@eq:posterior] is applied to
noisy data.

![**Pipeline for numerical calculation of estimation error.** The following
steps are implemented in the simulations for the forthcoming noise models.
First, an appropriate range of noise is defined. For each value in this
range, the complete simulation of the experiment is performed $N_{sim}$ times
as is described in [@fig:dilution_sim], this time including the error. [@Eq:posterior]
is then minimized given this simulated data and a best-fit estimate of $\alpha$
is obtained. This value is compared with the true (seeded) value of $\alpha$
and is stored in a vector for later analysis. This cycle is repeated for
every value defined in the original noise range.
 ](../figs/simulation_scheme.pdf){#fig:simulation_schematic}

### Simple measurement error

Perhaps the easiest type of noise to incorporate is the error in measurement
due to random shot noise of the camera. This type of error (at least for the
purposes of this work) can be considered to be independent of space, time,
and all features of the cell. As this is independent, we can say that measured
cellular intensity is
$$
I = \alpha N + \epsilon,
$${#eq:model1_ian}
where $\epsilon$ is normally distributed with zero mean and a variance $\sigma_1$.
This type of measurement error has been thoroughly dissected using Bayesian
methods similar to those discussed previously [@Rosenfeld2006]. Implementation
of this inference method is computationally costly and chock-full of caveats
that are beyond the task of this writing.

There can be multiple sources for $\epsilon$ entering this model. For example,
this could be the photon counting error in the camera, an autofluorescence
distribution which is not a delta function, or even blinking fluorophores. The
measurement error through the camera is the easiest to access experimentally.

[@Fig:shot_noise_measurement] shows a measurement of random measurement error
for three different sources of illumination. To make this measurement, a homogeneously
fluorescent slide was imaged for a single exposure. A 100-by-100 pixel region
directly in the center of the imaged was then selected and a distribution of
the pixel intensities was generated. The center of the image was chosen as
this should be the brightest region and the intensity distribution can be
approximated to be uniform. As there is vignetting with all three sources
of illumination, it would not be fair to use the total fluorescence of the
entire image. Representative images are shown in [@fig:shot_noise_measurement] (A).
The distribution of pixel intensities (normalized and centered at zero) are shown in
panel (B). All three sources of illumination appear to be relatively tightly distributed
about the mean with variation ranging from one to two percent variation, with
the laser source having the largest variance.

To examine the how important this error is in the measurement of the calibration
factor, we performed the simulation described earlier, this time adding noise
to the measurement of each sister cell pair. [@Fig:model1_sim] shows the result of
this simulation. The measurement noise range was chosen to cover the extremes
of the potential variation, ranging from $10^{-3} \times \alpha$ to $10^2 \times \alpha$.


![**Measurement noise for three illumination sources.** (A) Center 10 by 10 pixel
region of a homogeneously fluorescent slide imaged under the three illumination
sources. Blue and yellow pixels correspond to low and high intensities, respectively.
(B) Empirical cumulative distributions (black points) for each illumination
source and the Gaussian approximation (red line). The pixel values from the
images above the distributions were rescaled to the raw value of the image mean
and centered at zero. This allows for direct comparison of the three distributions.](../figs/shot_noise_measurement.pdf){#fig:shot_noise_measurement}


![**Numerical error estimate from neglecting random measurement noise.**
Black points indicate estimated values of $\alpha$ from each simulated experiment.
The inset shows the range of measurement noise relevant to the three described
illumination sources.](../figs/error_est_model1.pdf){#fig:model1_sim}




### Temporal variation

![**Measurement of temporal variation for three illumination sources.** (A)
Measurement of total image intensity of a fluorescent slide over a 100 exposure
experiment. Each "pixel" represents the sum total intensity of the image at
that frame number. The rows correspond to exposures taken with an Hg lamp, LED,
and 589 nm laser illumination, respectively. (B) Normalized intensity
distributions of images summarized in (A) after normalization. Red curves are
the Gaussian approximation of the measured distribution.
](../figs/temporal_noise_measurement.pdf){#fig:temporal_noise_measurement}

# References

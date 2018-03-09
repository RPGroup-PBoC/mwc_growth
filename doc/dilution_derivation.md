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


The real test of this method is being able to faithfully reproduce results
seen before in their entirety. This means growing the same batch of cells using
the same growth conditions and methods. [@Fig:example_experiments] shows the
results of such reproduced experiments using three different types of excitation
wavelengths -- a mercury (Hg) lamp, white light emitting diode (LED), and a
589 nm laser. The top two and bottom left-hand plots in [@fig:example_experiments]
show the quantified fluctuations in intensity and the estimated calibration
factors for each set. It's notable that for both the Hg lamp and the LED system,
the fit agrees quite nicely with the means (indicated by red dots for arbitrary
bins), but disagrees with the laser illumination. This is particularly obvious
when using another method to test the accuracy of the counts. Without going
into the gory details of the model, the bottom right-hand plot in [@fig:example_experiments]
shows the predicted change in the expression of a reporter gene as the number
of repressors (which are counted using the dilution method) is modulated. For
both the LED and Hg lamp derived calibration factors, the data agrees with the
prediction. However, the repressor copy number measurements are shifted by an
 factor of approximately five to six. This is a factor that has been high repeatable
 when using a laser excitation source, even when a different fluorophore is used.

 To determine why this particular illumination method does not work, we examined
 various sources of experimental error and measured their resulting bias on
 the calibration factor estimation *in silico*.

![**Calibration factor estimation under three different illumination sources.**
Representative examples of experimental measurements using the three different
sources of excitation illumination. Top-left and top-right correspond to experiments performed
with an Hg lamp source and a LED, respectively. The bottom left corner was collected using a 589 nm laser source. The bottom right panel shows the fold-change in gene expression measurements from all three experiments shown in the other panels.](../figs/light_source_experiments.pdf){#fig:example_experiments}

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
to the measurement of each sister cell pair. [@Fig:err_est_model1] shows the result of
this simulation. The measurement noise range was chosen to cover the extremes
of the potential variation, ranging from $10^{-3} \alpha$ to $10^2 \alpha$.
These results indicate that the effect of random measurement error does not
strongly affect the resulting calculation of the calibration factor, at least
until the relative value of the error becomes large (> 1\%). Past this point,
the noise very strongly influences the measurement. Experimental measurement
of the random noise of the three illumination systems shown in [@fig:shot_noise_measurement]
are shown in the inset of [@fig:err_est_model1]. In this range, the measured calibration
factor is nearly identical to the true value. This is perhaps not surprising,
but is a satisfying confirmation that this is not the source of the repeatable
factor of five error when using a laser.


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
illumination sources.](../figs/error_est_model1.pdf){#fig:err_est_model1}

### Temporal variation

Another potential source of noise is variation in illumination intensity from
snapshot to snapshot. We often think of these light sources as static, unwavering
beams of light. In reality, there is constant and fairly predictable variations
in intensity over many time scales. As the intensity of the emitted light from
the fluorophore is dependent on the intensity of the excitation light, it's
plausible that slight changes in the illumination could skew the cellular measurements.

We can, for the purposes of this document, assume that the intensity of the
emitted light $I$ is proportional to the intensity of the excitation light $\Phi$,
$$
I = \Phi\tau,
$${#eq:temporal_i_phi_tau}
where $\tau$ is the exposure time. With this assumption, we can rephrase
[@eq:antot] to be depending on the intensity of the excitation light as
$$
I = \phi(\sigma_2)\alpha N,
$${#eq:temporal_ian}
where $\phi$ is a rescaled parameterization of $\Phi$ which is normally
distributed with a mean of one and a variance $\sigma_2^2$. If the intensity
is fixed across time, $\phi$ should be a delta function at 1, meaning all
imaged sister cells obey as is presented in the perfect case. However, as
the variation in the intensity increases, some sister cells will be exposed
to more or less excitation light, making their measured fluorescence to be

To get a handle on the degree of variation for the three illumination schemes
in question, we imaged a single field of a homogeneously fluorescent slide
for one hundred consecutive frames at intervals of 100 ms. [@Fig:temporal_noise_measurement]
shows the results of these experiments. [@Fig:temporal_noise_measurement] (A) gives a qualtitative
comparison of the measured intensity over the course of the acquisition. Each
"pixel" represents the sum total intensity of the image at that frame number.
(B) shows the distribution of the sum total intensities (black points) along
with a Gaussian approximation (red lines) which describes the distributions
reasonably well.

We then simulated the impact of temporal variation in intensity by using [@eq:temporal_ian]
to determine the single-cell intensities. In our simulations, we treated each
sister cell pair as it's own position and therefore had it's own value of $\phi$.
[@Fig:temporal_sim] summarizes the error in the estimation of the calibration
factor when this noise source is neglected. As was seen for random measurement
noise, the influence of temporal variation only becomes pronounced when the
variation becomes large. The $x$-axis of [@fig:temporal_sim] shows the standard
deviation $\sigma_2$ used to determine the intensity. The error only becomes
large when the intensity is changing by a factor of two or more for each snapshot.
The inset shows the region of parameter space where the errors measured in
[@fig:temporal_noise_measurement] are of interest. It appears that all three
methods do not introduce appreciable noise to the measurement and can thus likely
be neglected.

![**Measurement of temporal variation for three illumination sources.** (A)
Measurement of total image intensity of a fluorescent slide over a 100 exposure
experiment. Each "pixel" represents the sum total intensity of the image at
that frame number. The rows correspond to exposures taken with an Hg lamp, LED,
and 589 nm laser illumination source, respectively. (B) Normalized intensity
distributions of images summarized in (A) after normalization. Red curves are
the Gaussian approximation of the measured distribution.
](../figs/temporal_noise_measurement.pdf){#fig:temporal_noise_measurement}


![**Numerical error estimate from neglecting temporal variation.** Each black point represents the the most likely parameter value for the calibration factor when temporal variation in intensity is neglected. The inset shows the region of parameter space that is appropriate for the errors measured in [@fig:temporal_noise_measurement].](../figs/err_temporal_variation.pdf)


### Spatial variation
Much as we imagine these light sources to be stable over time, we unfortunately
often think of them as completely uniform across their width. While we are
probably all comfortable thinking of a laser beam as Gaussian in profile, it's
important to remember that there are fringes caused by interference. These
fringes can be thought of as spatially localized variations in intensity, as
we worked through in the previous section. For large objects or averaged measurements,
this is not too much of a concern. However, these fringes can be on the length
scale of a few bacterial cells, which can be very important for our system of
interest.

Incorporating this noise in to our model now becomes a little more tricky.
Rather than saying we can draw the values of our noise from a tidy distribution,
the position of these fringes will depend strongly on the collimation of the beam
and its alignment with the chip of the camera. This means that any given
sister cell pair has a probability $p$ of being in a position that rests on
one of these fringes. We can modify [@eq:antot] as
$$
I = \theta_{pos}\alpha N,
$${#eq:spatial_ian}
where $\theta_{pos}$ is on the intensity of the excitation beam at that
position and is on the range $[0, 1]$, with $0$ being completely dark and $1$
being maximally bright. Unlike the previous models of noise we've wrestled with,
there are two parameters here to consider -- the degree of intensity variation
as well as the fraction of cells affected.

The error estimates for this model can be seen in [@fig:err_spatial_variation].
Here, we've simulated the experiment varying both the fraction of the
sister cell pairs affected ($x$-axis) as well as the magnitude of the
intensity variation of the fringe (colored lines). This effect is rather striking
and, as expected, can introduce considerable error into the experiment. It
appears that  even with 10\% variation in intensity from fringe to fringe,
errors in the estimation of $\alpha$ can be off several fold, even when only
1\% of the sisters are affected. This, so far, appears to be the largest
contributor to potential measurement error.


![**Numerical error estimate from neglecting spatial intensity variation.** Each point indicates the most likely parameter value for the calibration factor when all sources of spatial variation are ignored. The two independent parameters (number of lineages affected and fractional intensity difference) were calculated in tandem.](../figs/err_spatial_variation.pdf){#fig:err_spatial_variation}

Unlike the other two sources of error, the degree of spatial variation for
our three light sources is a bit trickier to measure. To get a sense for what
the illumination field looks like, it is no longer sufficient to look at a
homogeneously fluorescent slide as the features of the illumination can be
overwhelmed from out-of-plane fluorescence. To image the illumination field,
we took a sample of highly dilute fluorescent beads and slowly but surely rastered
a single bead across the field of view. By scanning a point source across the
illumination field, we can accurately measure the degree of intensity variation
on a length scale comparable to a single cell. A maximum projection of such
an experiment can be seen in [@fig:max_projection]. This image has been flatted
to correct for large-scale non-uniformity of the laser beam profile. This means
that any small variations in intensity that are left are non-static fringes
that could be present in our experimental data sets.

![**Map of spatial variation in laser illumination.** Maximum projection of a 0.1 $\mu$m fluorescent bead scanned across the field of view using a 589nm laser excitation source. This image has been flattened using an average image of the illumination profile from a fluorescent Indium Tin Oxide coated slide. Scalebar is 10 $\mu$m.](../figs/max_projection.pdf){#fig:max_projection}

As a first pass at quantifying these variations, we examined each row of imaged
positions by calculating the total intensity of each bead and measuring the
fractional difference from it's neighbor,
$$
\delta I = {I_n - I_{n+1} \over I_n},
$${#eq:fractional_difference}
for all $n$ images in a given row or column. This approach can been seen in
[@fig:laser_spatial_variation] where we examined only the $x$-dimensional spatial
variations. Each faint blue line corresponds to a single row of beads. The
dark blue line highlights a representative trace whose  images can be seen
at the top of the plot. The variation in intensity is larger than one would have initially suspected {\text{color{red} I am still suspicious of the magnitude here -- up to a factor of 2?!}},
suggesting that the small-scale spatial fluctuations in intensity may be in
the range of plausibility for this experiment.

![**Percent variation in intensity between adjacent particles.**  Individual rows
of the bead shown in[@fig:max_projection] were isolated and the fractional difference in bead intensity was calculated between adjacent neighbors. Each row is plotted as a thin blue line. The dark line is a representative trace highlighted for clarity. The beads shown in the top panel is the row highlighted in the plot.](../figs/laser_spatial_fluctuations.pdf){#fig:laser_spatial_variation}



## The state of the union
It seems like temporal and measurement noise can be neglected, so long as they
aren't enormous. Measurement of the errors for the three different types of
illumination schemes used in this work show that we are operating far below the
these extreme limits. Spatial variation in intensity, however, can lead to
very large errors in the estimated calibration factor value. Our cursory
measurement of the spatial variation of a 589 nm laser illumination source
shows that the percent variation in intensity on the length scale of a typical
cell can result in these large overestimations of the calibration factor.

For experiments in which the absolute intensity of a given cell is of
quantitative interest, it may be worthwhile to take the hit in the signal to
noise ratio and photobleaching effects to ensure spatial uniformity.


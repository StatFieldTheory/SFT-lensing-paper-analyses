# Proposed revisions: the five items left after the external-review triage

*2026-08-29. **APPLIED** in full, minimal variants, +21/-15 lines. The six
confirmed errors were fixed earlier in `c796237`; these were the remaining
items, all wording rather than error. Kept as the record of why each was
changed and of what was deliberately not done.*

Each item gives the exact current text and the exact replacement. Items are
independent; take any subset.

---

## 1. Poisson gauge naming  (reviewer item 6) — RECOMMENDED

**The defect.** Poisson gauge is *defined* by removing the scalar part of the
shift. The metric ansatz keeps `- 2 \partial_i B d\tau dx^i` and calls itself
the Poisson gauge, and `B` is called "the Poisson-gauge scalar shift", which is
self-contradictory. The paper already knows this: cosmology.tex:661 says "In
the standard Poisson gauge the shift is set to zero by the gauge choice,
$B=0$". So only the naming at the point of definition is wrong, not the physics.

**Second, sharper defect the reviewer only half-caught.** cosmology.tex:663
offers to identify `B` with "a frame-dragging vector sector". It cannot be one:
`\partial_i B` is a pure gradient, i.e. longitudinal, whereas a frame-dragging
mode is a transverse `B_i` with `\partial^i B_i = 0`. The ansatz therefore adds
the one scalar piece Poisson gauge excludes *and* omits the genuine vector
channel entirely.

**Why "Poisson gauge" can stay downstream.** All ten later uses sit where
`B = 0` already, so they are correct as written. Only three lines need touching.

### 1a. cosmology.tex:25

> before: perturbation, `$B$ the Poisson-gauge scalar shift`, and $h_{ij}$ the

> after: perturbation, `$B$ a scalar shift, retained for generality and set to
> zero by the Poisson gauge itself`, and $h_{ij}$ the

### 1b. cosmology.tex:38

> before: `All four perturbation channels are` retained throughout this section;

> after: `All four channels ($\Psi$, $\Phi$, $B$, $h_{ij}$) are` retained
> throughout this section;

*Rationale: "perturbation channels" invites the scalar/vector/tensor reading,
under which the list is wrong, since no vector mode appears anywhere. Naming
the four fields is accurate and costs two words.*

### 1c. cosmology.tex:662-663

> before: contributions drop out; we record them here for reference `in case the
> gauge is relaxed or $B$ is identified with a frame-dragging vector sector.`

> after: contributions drop out; we record them here for reference `in case the
> gauge is relaxed. A frame-dragging mode is a transverse $B_i$ with
> $\partial^i B_i = 0$, which the scalar shift $\partial_i B$ cannot represent;
> the vector sector is not treated here.`

**Alternative, if you would rather not admit a missing sector:** delete the
`\partial_i B` term from Eq. (49), delete 1a, and drop Eqs. (eq: Phi00 B) and
(eq: Psi0 B). The numerics never use them. This is cleaner but discards work.

---

## 2. "connected" on the observable side  (reviewer item 2) — MARGINAL

**The defect, and how small it is.** The reviewer claimed five sentences
misuse "connected" on the observable side. Three of them are correct usage,
because they are scoped to *leading* order, where the free theory is a centred
Gaussian and no disconnected piece exists yet (the disc topology first appears
at Order-2, with two vertices): cosmology.tex:725, :1065, :1075 all stand.

Two remain, and both are claims about what the *formalism* resolves, which is
true of the formalism whichever object you choose to plot. So this is a
tightening, not a correction.

### 2a. intro.tex:60-61

> before: The resulting diagrammatic expansion resolves every `connected`
> $n$-point correlation function ($n$PCF) of convergence and shear into three
> ingredients:

> after: The resulting diagrammatic expansion resolves every $n$-point
> correlation function ($n$PCF) of convergence and shear into three
> ingredients:

*The second "connected" two lines later, on the driving-field cumulants, is
correct and stays.*

### 2b. conclusion.tex:9

> before: the Sachs trajectory to a path integral turns every `connected`
> $n$-point function of convergence and shear into a diagrammatic expansion,

> after: the Sachs trajectory to a path integral turns every $n$-point
> function of convergence and shear into a diagrammatic expansion,

**Judgement: I would take this.** It costs one word twice and removes the only
handle a referee has for the "you plot moments, not cumulants" objection.

---

## 3. Magnitude on logarithmic axes  (reviewer item, figures) — RECOMMENDED

**The defect, scoped honestly.** The reviewer named five figures. Only two
plot a genuinely sign-changing quantity under a bare label: Fig. 15's top row
and Fig. 17. Figs. 10 and 11 carry no axis label at all, and Fig. 14's EB panel
is on a null quantity. The convention is already disclosed six times in
captions and inside three panels, so no reader is misled — but the ordinate
really is `|xi|`, and saying so costs one clause and no re-rendering.

**One clause, applied at four places** (insights.tex:134, insights.tex:152,
conclusion.tex:99, conclusion.tex:116; appendix.tex:619 carries the same
sentence):

> before: `Markers are filled where positive, hollow where negative.`

> after: `The ordinate is the magnitude; markers are filled where positive,
> hollow where negative.`

and at insights.tex:152, the compressed variant:

> before: `Markers filled/hollow for positive/negative.`

> after: `The ordinate is the magnitude, markers filled/hollow for
> positive/negative.`

**Not proposed:** re-rendering the figures with `$|\xi_\kappa|$` axis labels.
That touches five figure products for a caption-level clarification, and the
paper's own text repeatedly discusses where the curves change sign.

---

## 4. The screen basis behind Eq. (42)  (reviewer item 3) — RECOMMENDED

**Why the reviewer's headline is wrong, and what is genuinely missing.** The
reviewer says Eq. (42)'s ordinary Legendre expansion is "valid only for
spin-0" and that the E/B results therefore rest on a wrong formula. That does
not survive: each component correlator is a bounded function of the single
angle $\gamma$, hence lies in $L^2([-1,1])$ where $P_\ell$ is complete, and the
`d^\ell_{2,\pm2}` are themselves expandable in $P_\ell$. Read as *defining*
`C_\ell^{ab} = 2\pi \int C^{(2)}_{ab}(\gamma) P_\ell(\cos\gamma) d\cos\gamma`,
Eq. (42) is exact for every $a,b$, and the substitution into Eq. (43) is exact.
Nothing downstream is corrupted.

What is genuinely absent is a statement that the components are taken in the
*pair-aligned screen basis*. Without it, a referee cannot tell whether
`C_\ell^{ab}` are spin-weighted coefficients (they are not) or Legendre
coefficients of pair-aligned components (they are), and the reviewer's
objection is exactly the confusion that gap invites.

### path_int.tex:286-289

> before: $C_\ell^{\sachi\sachj}(\lambda',\lambda'')$ `is the corresponding
> angular power spectrum between the driving-field components` $\sachi$ and
> $\sachj$ at affine parameters $\lambda'$ and $\lambda''$.

> after: $C_\ell^{\sachi\sachj}(\lambda',\lambda'')$ `collects the Legendre
> coefficients of the corresponding component correlator, with the components
> taken in the pair-aligned screen basis fixed in
> Sec.~\ref{sec: cosmological setup} ($\hat x$ along the separation). Each such
> correlator is a function of the single angle between the two directions, so
> the expansion is exact for the spin-2 components as well as the scalar one;
> the spin weight is carried by the basis, not by the expansion.` The components
> are at affine parameters $\lambda'$ and $\lambda''$.

*This is the one place where a referee could otherwise repeat the reviewer's
objection, and it is cheap to close.*

---

## 5. Readout-versus-propagation pointer  (reviewer item 1) — OPTIONAL

**Status: the reviewer's charge failed, but their order counting is right.**
With $\theta = \dot{\bar D}/\bar D$ and $X_1 = \theta - \bar\theta$, one has
$X_1 = d\ln(1-\kappa)/d\lambda$ exactly, so $-\int X_1 d\lambda = -\ln(1-\kappa)
= \kappa + \kappa^2/2 + \dots$, and the dropped pieces are $O(P^2)$, the same
order as FF and FK. Their operative claim, that the paper never says so, is
false: conclusion.tex:41-44 states it, and the legends read "Full (O0+FF+FK)",
which is self-defining. Their grep missed it because the paper writes
"reduced-shear" hyphenated and line-broken.

So no change is *required*. The optional improvement is to put the caveat in
front of the reader where the numbers are first quoted, rather than only in the
conclusion twenty pages later.

### insights.tex, after the $1.3$--$1.7\%$ sentence in Sec. V C

> add: `These are corrections to the propagation of the bundle; the mapping
> from $(\kappa,\gamma)$ to a measured galaxy shape carries its own correction
> of the same order, additive to these (Sec.~\ref{sec: conclusion}).`

**Judgement: worth it if you expect a referee to raise reduced shear.** It
pre-empts the single most likely objection at the cost of one sentence, and it
is the objection this reviewer raised first and ranked highest.

---

## Not proposed

- **Re-running to $\ell_{\max}=30720$.** Would settle the convergence question
  properly, but `driver_field_emulators/notes/finding_ell_max_resolved.md`
  records the fiducial CAMB table ($k_{\max} = 206\,h/$Mpc) exhausting above
  $\ell \sim 27000$, so it needs a longer $P(k)$ table or an explicit high-$k$
  extrapolation first. Not a text edit.
- **Reworking the observables through the full Jacobi map.** The reviewer's
  largest suggestion. It is a research task, not a revision, and the paper's
  scope statement already covers it.

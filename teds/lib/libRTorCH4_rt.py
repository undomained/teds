import sys

import numpy as np
import torch

class ISRF(object):
    """Base class for instrument spectral response function."""

    def __init__(self, *args, **kwargs):
        pass

    def set_parameters(
        self,
        obs_wls: torch.Tensor,
        obs_fwhm: torch.Tensor,
        lbl_wls: torch.Tensor,
        cache_hint: bool
    ):
        """Set parameters for the ISRF.

        If `cache_hint` is True, indicates the ISRF object may assume the
        parameters have not changed since the previous call---by the same
        callers---to `set_parameters` (if any), and may not need to reprocess
        anything.
        Advanced users should note that in cases where the ISRF is used by
        multiple callers, it is up to the user to guarantee that calls to
        `set_parameters` are mutually consistent, or that no caching is used.
        (Such cases, however, are highly nonstandard; if you do not know
        whether or not you are dealing with such a case, you are not dealing
        with one, and can ignore this notice.)

        @param obs_wls Observed wavelengths (N × N_obs_wls)
        @param obs_fwhm Full width at half maximum for observed wavelengths
        ([N ×] N_obs_wls)
        @param lbl_wls Line-by-line wavelengths (N_lbl_wls)
        @param cache_hint Flag indicating if the ISRF can use cached
        parameters.
        """
        raise NotImplemented()

    def tensor(self):
        """Output a tensor of shape (N × N_obs_wls × N_lbl_wls) corresponding
        to the convolution operation, so that, given a target tensor `t` of
        shape (N × N_lbl_wls), the following holds:

            torch.einsum(
                "Nol,Nl->No",
                self.tensor(),
                t
            ) == self.convolve(t)

        @return Tensor of shape (N × N_obs_wls × N_lbl_wls)
        """
        raise NotImplemented()

    def gradient(self):
        """Output a tensor of shape (N × N_obs_wls × N_lbl_wls) corresponding
        to the gradient of `self.tensor()` w.r.t. observed wavelength.

        @return Tensor of shape (N × N_obs_wls × N_lbl_wls)
        """
        raise NotImplemented()

    def convolve(
        self,
        target: torch.Tensor
    ):
        """Convolve a tensor with the ISRF. The target is a tensor of
        dimension (N × N_lbl_wls [× P]), and the convolution should return a
        tensor of dimension (N × N_obs_wls [× P]).

        @param target Tensor to convolve with the ISRF. (N × N_lbl_wls [× P])
        @return Tensor of shape (N × N_obs_wls [× P])
        """
        raise NotImplemented()

    def convolve_gradient(
        self,
        target: torch.Tensor
    ):
        """Convolve a tensor with the ISRF's gradient. The target is a tensor
        of dimension (N × N_lbl_wls [× P]), and the convolution should return a
        tensor of dimension (N × N_obs_wls [× P]).

        @param target Tensor to convolve with the ISRF's gradient.
        (N × N_lbl_wls [× P])
        @return Tensor of shape (N × N_obs_wls [× P])
        """
        raise NotImplemented()

class GaussianISRF(ISRF):
    """Class representing a Gaussian ISRF.

    N.B. this class uses caching."""

    def __init__(self, epsilon=1e-6):
        """Initialiser.

        Values of the Gaussian that are below `epsilon` are set to zero.
        """
        self.epsilon = epsilon
        self.N = None
        self.N_obs_wls = None
        self.N_lbl_wls = None
        self.isrf_tensor = None
        self.gradient_tensor = None
        self.device = None
        self.dtype = torch.float32

    def set_parameters(
        self,
        obs_wls: torch.Tensor,
        obs_fwhm: torch.Tensor,
        lbl_wls: torch.Tensor,
        cache_hint: bool
    ):
        if self.isrf_tensor is not None and cache_hint:
            return

        if self.device is None:
            self.device = obs_wls.device

        self.N, self.N_obs_wls = obs_wls.shape
        self.N_lbl_wls = lbl_wls.shape[0]

        if len(obs_fwhm.shape) == 1:
            obs_fwhm = torch.outer(
                torch.ones(self.N, device=self.device, dtype=self.dtype),
                obs_fwhm
            )

        # Difference between observed and lbl wavelengths
        isrf_wavediff = (
            torch.einsum(
                "No,l->Nol",
                obs_wls,
                torch.ones(self.N_lbl_wls, device=self.device, dtype=self.dtype)
            ) - torch.einsum(
                "No,l->Nol",
                torch.ones(
                    self.N, self.N_obs_wls,
                    device=self.device, dtype=self.dtype
                ),
                lbl_wls,
            )
        )

        # isrf_cst_grid: grid of standard deviation constants
        # (N, N_lbl_wls, N_obs_wls)
        isrf_cst_grid = 1/torch.einsum(
            "No,l->Nol",
            obs_fwhm**2 / 4 / np.log(2),
            torch.ones(
                self.N_lbl_wls,
                device=self.device,
                dtype=self.dtype
            )
        )

        # isrf:
        # (N × N_obs_wls × N_lbl_wls)
        isrf_num = torch.exp(
            -isrf_wavediff**2 * isrf_cst_grid
        )
        mask = torch.abs(isrf_num) < self.epsilon
        isrf_num[mask] = 0
        isrf_denom = 1/torch.sum(isrf_num, dim=2, keepdim=True)
        self.isrf_tensor = isrf_num * isrf_denom

        grd_cst = 2*isrf_cst_grid*isrf_wavediff
        self.gradient_tensor = (
            grd_cst
            - (
                torch.sum(grd_cst * isrf_num, dim=2, keepdim=True)
                * isrf_denom
            )
        ) * self.isrf_tensor

    def tensor(self):
        return self.isrf_tensor

    def gradient(self):
        return self.gradient_tensor

    def convolve(self, target):
        if len(target.shape) == 2:
            return torch.einsum(
                "Nol,Nl->No",
                self.isrf_tensor,
                target
            )
        else:
            return torch.einsum(
                "Nol,NlP->NoP",
                self.isrf_tensor,
                target
            )

    def convolve_gradient(self, target):
        if len(target.shape) == 2:
            return torch.einsum(
                "Nol,Nl->No",
                self.gradient_tensor,
                target
            )
        else:
            return torch.einsum(
                "Nol,NlP->NoP",
                self.gradient_tensor,
                target
            )

class BasicRadiativeTransfer(object):
    """Simplistic singlge-scattering radiative transfer class. Compute
    top-of-atmosphere Earth-reflected radiance given scattering geometry,
    optical depth, albedo and solar irradiance."""
    epsilon = 1e-9 # Fudge factor to prevent division by 0

    def __init__(
        self,
        sza,
        vza,
        **kwargs
    ):
        """Initialise RT for given scattering geometry.

        @param sza Solar zenith angle per pixel.
        @param vza Viewing zenith angle per pixel."""

        if kwargs is not None and len(kwargs) > 0:
            print(
                (
                    "WARNING: kwargs to BasicRadiativeTransfer "
                    "initializer will not be used."
                ),
                file=sys.stderr
            )

        self.N = sza.shape[0]
        self.sza = torch.clone(sza)
        self.vza = torch.clone(vza)

        self.mu_sol = torch.cos(self.sza*torch.pi/180)
        if torch.any(
            (torch.abs(self.mu_sol) < self.epsilon) | (self.mu_sol > 1)
        ):
            raise ValueError("Solar cosine out of range")
        self.mu_view = torch.cos(self.vza*torch.pi/180)
        if torch.any(
            (self.mu_view < -1) | (self.mu_view > 1)
            | (torch.abs(self.mu_view) < self.epsilon)
        ):
            raise ValueError("Viewing cosine out of range")

        self.mu_eff = torch.abs(1/self.mu_sol)+torch.abs(1/self.mu_view)

    def radiance_toa(self, tau, albedo, irradiance_lbl, grad=True):
        """Compute top-of-atmosphere reflected radiance using simple
        optical depth model, albedo, and solar irradiance.

        @param tau Optical depth per pixel [per atmospheric layer] per
        wavelength (N × [N_layers ×] N_lbl_wls)
        @param albedo Surface albedo per wavelength (N × N_lbl_wls)
        @param irradiance_lbl Solar irradiance (N_lbl_wls)
        @param grad If True, also compute gradients w.r.t. optical depth and
        albedo.
        @return If grad == False:
            Top-of-atmosphere Earth-reflected radiance
        (N × N_lbl_wls).
        If grad == True: (
            Top-of-atmosphere Earth-reflected radiance (N × N_lbl_wls),
            ∂R/∂(tau_total) (N × N_lbl_wls),
            ∂R/∂A (N × N_lbl_wls)
        )
        """

        dev_alb = torch.einsum(
            "N,Na,a->Na",
            self.mu_sol,
            torch.exp(
                -torch.einsum(
                    "Na,N->Na",
                    (torch.sum(tau, dim=1) if len(tau.shape) == 3 else tau),
                    self.mu_eff
                )
            ),
            irradiance_lbl
        ) / torch.pi
        rad_trans = torch.einsum(
            "Na,Na->Na",
            albedo,
            dev_alb
        )
        if grad:
            dev_tau = -torch.einsum("N,Nl->Nl", self.mu_eff, rad_trans)
            self.last_output__radiance_toa = (rad_trans, dev_tau, dev_alb)
        else:
            self.last_output__radiance_toa = rad_trans

        return self.last_output__radiance_toa

class RTForwardModel(object):
    """Basic tensorised forward model class."""

    epsilon = 1e-9 # Fudge factor to prevent division by 0

    def __init__(
        self,
        obs_wls,
        obs_fwhm,
        lbl_wls,
        tau_base,
        sun_lbl,
        sza,
        vza,
        tau_offset=None,
        device=None,
        dtype=None,
        isrf=None,
        rt_class=BasicRadiativeTransfer,
        rt_args={}
    ):
        """Initialiser. Tensor dimensions are inferred from the parameters
        passed to this function.

        An ISRF object should be given, taylored to the problem at hand. If
        `isrf` is None, a GaussianISRF will be created upon initialisation.

        By default, this class uses BasicRaditiveTransfer as the radiative
        transfer class, however, a different class may be selected using the
        `rt_class` argument. In this case, the initialiser is called as
        `rt_class(sza, vza, **rt_args)`.

        @param obs_wls Wavelengths in observation dataset [per pixel]
        ([N ×] N_obs_wls)
        @param obs_fwhm Spectral full width at half-maximum [per pixel]
        ([N ×] N_obs_wls)
        @param lbl_wls Wavelengths used in line-by-line computations
        (N_lbl_wls)
        @param tau_base Optical depth [per pixel] per species per layer per LBL
        wavelength ([N ×] N_species × N_layers × N_lbl_wls)
        @param sun_lbl Solar spectrum sampled at lbl wavelengths (N_lbl_wls)
        @param sza Solar zenith angle per pixel (N)
        @param vza Viewing zenith angle per pixel (N)
        @param tau_offset Constant optical depth term, e.g. for non-retrieved
        absorbers ([N ×] N_layers × N_lbl_wls)
        @param device PyTorch device to use. If None, use `obs_wls.device`
        @param dtype Data type to use for tensors. Should generally be float32
        or float64. If None, use `obs_wls.dtype`
        @param isrf Instrument spectral response function object
        @param rt_class Radiative transfer class to use.
        @param rt_args Dict of arguments to radiative transfer class
        initialiser."""
        self.N = sza.shape[0]
        self.device = device if device is not None else obs_wls.device
        self.dtype = dtype if dtype is not None else obs_wls.dtype

        if len(obs_wls.shape) == 1:
            self.N_obs_wls = obs_wls.shape[0]
            self.obs_wls = torch.outer(
                torch.ones(self.N, dtype=self.dtype, device=self.device),
                obs_wls
            )
        else:
            self.N_obs_wls = obs_wls.shape[1]
            self.obs_wls = torch.clone(obs_wls)

        if len(obs_fwhm.shape) == 1:
            self.obs_fwhm = torch.outer(
                torch.ones(self.N, dtype=self.dtype, device=self.device),
                obs_fwhm
            )
        else:
            self.obs_fwhm = torch.clone(obs_fwhm)


        self.lbl_wls = torch.clone(lbl_wls)
        self.N_lbl_wls = lbl_wls.shape[0]

        if len(tau_base.shape) == 3:
            self.N_species = tau_base.shape[0]
            self.N_layers = tau_base.shape[1]
        else:
            self.N_species = tau_base.shape[1]
            self.N_layers = tau_base.shape[2]

        self.sun_lbl = torch.clone(sun_lbl)
        self.inversion_exceptions = []

        self.rt_class = rt_class
        self.rt_args = rt_args

        if isrf is not None:
            self.isrf = isrf
        else:
            print(
                "Notice: no ISRF provided; will use GaussianISRF.",
                file=sys.stderr
            )
            self.isrf = GaussianISRF()

        self.isrf_cache_good = True

        self.setup_measurement(
            sza,
            vza,
            tau_base,
            tau_offset=tau_offset
        )

        # Experimental feature: pass computed albedo through the function
        #
        #   (roll_albedo_max-roll_albedo_min)*(tanh(a)+1)/2+roll_albedo_min
        #
        # to constrain its value
        self.roll_albedo = False
        self.roll_albedo_min = torch.zeros(
            self.N, self.N_lbl_wls,
            dtype=self.dtype, device=self.device
        )
        self.roll_albedo_max = torch.ones(
            self.N, self.N_lbl_wls,
            dtype=self.dtype, device=self.device
        )

    def setup_measurement(
        self,
        sza,
        vza,
        tau_base=None,
        tau_offset=None,
    ):
        """Setup measurement variables (geometry and optical depth base) while
        preserving more generic variables like wavelength definitions. This
        method may be used to switch to a different patch of the same scene
        with reduced overhead compared to re-initialising the full class.

        @param sza Solar zenith angle per pixel (N)
        @param vza Viewing zenith angle per pixel (N)
        @param tau_base Optical depth prior [per pixel] per species per layer
        per LBL wavelength ([N ×] N_species × N_layers × N_lbl_wls) or None
        @param tau_offset Optical depth constant term [per pixel] per layer per
        LBL wavelength ([N ×] N_layers × N_lbl_wls) or None
        or N)"""
        self.rt_solver = self.rt_class(sza, vza, **(self.rt_args))

        if tau_base is not None:
            if len(tau_base.shape) == 3:
                tau_base = torch.einsum(
                    "N,slz->Nslz",
                    torch.ones(self.N, device=self.device, dtype=self.dtype),
                    tau_base
                )
            self.tau_per_molec_base = torch.clone(tau_base)

            # Sum over layers
            # (N × N_species × N_layers × N_lbl_wls) ->
            #   (N × N_species × N_lbl_wls)
            self.tau_species_sum = torch.sum(self.tau_per_molec_base, dim=2)

        if tau_offset is not None:
            if len(tau_offset.shape) == 2:
                tau_offset = torch.einsum(
                    "N,lw->Nlw",
                    torch.ones(self.N, device=self.device, dtype=self.dtype),
                    tau_offset
                )
            self.tau_offset = torch.clone(tau_offset)
        else:
            self.tau_offset = torch.zeros(
                self.N, self.N_layers, self.N_lbl_wls,
                device=self.device, dtype=self.dtype
            )

        # Null out the albedo dimension and related constants, so they can be
        # recomputed the first time they are needed.
        self.N_alb = None
        self.N_params = None
        self.N_params1 = None
        self.id_params = None
        self.alb_wl_grid = None
        self.alb_wls = None

    def compute_albedo(self, x_alb):
        """Compute albedo from polynomial in wavelength.

        Assume `self.N_alb` terms. If `self.N_alb` is `None`, set `self.N_alb`
        to `x_alb.shape[1]`. This method depends on the normalised
        line-by-line wavelengths `self.alb_wls`. If this variable is not
        present, it is computed from `self.lbl_wls`.
        @param x_alb Series coefficients. (N × N_alb)
        @return Albedo expansion (N × N_lbl_wls)"""

        if not hasattr(self, "N_alb") or self.N_alb is None:
            self.N_alb = x_alb.shape[1]
        if not hasattr(self, "alb_wl_grid") or self.alb_wl_grid is None:
            if not hasattr(self, "alb_wls") or self.alb_wls is None:
                self.alb_wls = (
                    (self.lbl_wls-self.lbl_wls[0])
                    / (self.lbl_wls[-1]-self.lbl_wls[0])
                )
            # Normalised wavelengths ^ expansion coefficients
            # (N_lbl_wls), (N × N_alb) -> (N × N_lbl_wls × N_alb)
            self.alb_wl_grid = (
                torch.einsum(
                    "l,Na->Nla",
                    self.alb_wls,
                    torch.ones(
                        self.N, self.N_alb,
                        device=self.device, dtype=self.dtype
                    )
                ) ** torch.einsum(
                    "a,Nl->Nla",
                    torch.arange(
                        self.N_alb,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.ones(
                        self.N, self.N_lbl_wls,
                        device=self.device, dtype=self.dtype
                    )
                )
            )

        albedo = torch.sum(
            torch.einsum(
                "Na,l->Nla",
                x_alb,
                torch.ones(
                    self.N_lbl_wls,
                    device=self.device,
                    dtype=self.dtype
                )
            ) * self.alb_wl_grid,
            dim=2
        )
        self.last_output__compute_albedo = albedo
        return self.last_output__compute_albedo

    def compute_optical_depth(self, column_scales):
        """Compute optical depth. Multiply the precomputed unscaled optical
        depths by the given column scales per species to obtain optical depth
        per layer per species per wavelength. Sum over species and add
        tau_offset to get optical depth per layer.

        IMPORTANT NOTE: the first returned item -- optical depth per species
        per layer -- does NOT include tau_total, whereas tau_total HAS been added
        to the second returned item (optical depth per layer).

        @param column_scales Gas column scales per pixel per species
        (N × N_species)
        @return (
            Optical depth per molecule per species
                (N × N_species × N_layers × N_lbl_wls),
            Optical depth per layer (N × N_layers × N_lbl_wls)
        )"""

        # (N × N_species × N_layers × N_lbl_wls), (N × N_species) ->
        #   (N × N_species × N_layers × N_lbl_wls)
        tau_per_molec = torch.einsum(
            "Nszl,Ns->Nszl",
            self.tau_per_molec_base,
            column_scales,
        )

        # (N × N_layers × N_lbl_wls)
        tau_layers = torch.sum(tau_per_molec, dim=1) + self.tau_offset

        self.last_output__compute_optical_depth = (tau_per_molec, tau_layers)
        return self.last_output__compute_optical_depth

    def forward_lbl(self, column_scales, x_alb, grad=False):
        """Forward model, line-by-line; compute radiance and optionally
        Jacobian, but do not convolve with ISRF.

        If if grad==True, the Jacobian will have a shape of
        (N × N_lbl_wls × (N_species+N_alb)).

        @param column_scales Gas column scaling factors. (N × N_species)
        @param x_alb Albedo series coefficients. (N × N_alb)
        @param waveshift Wavelength shift. (N)
        @param grad Set to True to also compute the Jacobian.
        @return If grad=True: tuple(
            Radiance (N × N_lbl_wls),
            Jacobian (N × N_lbl_wls × (N_species+N_alb))
        ). If grad=False: Radiance (N × N_lbl_wls)
        """

        if self.N_alb is None:
            self.N_alb = x_alb.shape[1]

        if self.N_params is None:
            self.N_params = self.N_species+self.N_alb
            self.N_params1 = self.N_params+1

        tau_per_molec, tau_layers = self.compute_optical_depth(
            column_scales
        )

        orig_albedo = self.compute_albedo(x_alb)
        albedo = (
            orig_albedo if not self.roll_albedo
            else (
                (self.roll_albedo_max-self.roll_albedo_min)
                * (torch.tanh(orig_albedo)+1) / 2
                + self.roll_albedo_min
            )
        )

        if not grad:
            self.last_output__forward_lbl = self.rt_solver.radiance_toa(
                tau_layers,
                albedo,
                self.sun_lbl,
                grad=False
            )
            return self.last_output__forward_lbl

        rad_trans, dev_tau, dev_alb = self.rt_solver.radiance_toa(
            tau_layers,
            albedo,
            self.sun_lbl,
            grad=True
        )

        # Part of gradients w.r.t. albedo factors
        # (N × N_lbl_wls), (N × N_lbl_wls × N_alb)
        #   -> (N × N_lbl_wls × N_alb)
        dev_alb_out = (
            dev_alb / (
                2 * torch.cosh(orig_albedo)**2
            )
            if self.roll_albedo
            else dev_alb
        )

        # Jacobian
        J = torch.cat(
            (
                torch.einsum(
                    "Nsl,Nl->Nls",
                    self.tau_species_sum,
                    dev_tau
                ), # Gradients w.r.t. column scales
                torch.einsum(
                    "Nl,Nla->Nla",
                    dev_alb_out,
                    self.alb_wl_grid
                ) # Gradients w.r.t. albedo factors
            ),
            dim=2
        )

        self.last_output__forward_lbl = (rad_trans, J)
        return self.last_output__forward_lbl

    def forward_inst(self, column_scales, x_alb, waveshift, grad=False):
        """Forward model, simulated instrument observation: compute radiance
        and optionally Jacobian, and convolve with ISRF.

        If if grad==True, the Jacobian will have a shape of
        (N × N_obs_wls × (N_species+N_alb+1)) if waveshift is not None,
        or (N × N_obs_wls × (N_species+N_alb)) if waveshift is None.

        @param column_scales Gas column scaling factors. (N × N_species)
        @param x_alb Albedo series coefficients. (N × N_alb)
        @param waveshift Wavelength shift. (N)
        @param grad Set to True to also compute the Jacobian.
        @return If grad=True: tuple(
            Radiance (N × N_obs_wls),
            Jacobian (N × N_obs_wls × (N_species+N_alb[+1]))
        ). If grad=False: Radiance (N × N_obs_wls)
        """

        obs_wls = (
            (
                self.obs_wls.detach() +
                torch.outer(
                    waveshift,
                    torch.ones(
                        self.N_obs_wls,
                        dtype=self.dtype,
                        device=self.device
                    )
                )
            )
            if waveshift is not None else self.obs_wls
        )

        if waveshift is not None:
            self.isrf_cache_good = False

        rad_trans = None
        J_in = None
        if grad:
            rad_trans, J_in = self.forward_lbl(
                column_scales, x_alb, grad=True
            )
        else:
            rad_trans = self.forward_lbl(column_scales, x_alb, grad=False)

        self.isrf.set_parameters(
            obs_wls,
            self.obs_fwhm,
            self.lbl_wls,
            cache_hint=self.isrf_cache_good
        )
        if waveshift is None:
            self.isrf_cache_good = True

        rad_out = self.isrf.convolve(rad_trans)

        if not grad:
            self.last_output__forward_inst = rad_out
            return self.last_output__forward_inst

        N_params = self.N_params if waveshift is None else self.N_params1

        # Jacobian
        J = self.isrf.convolve(J_in)

        J_out = J if waveshift is None else torch.cat(
            (J, self.isrf.convolve_gradient(rad_trans)),
            dim=2
        )

        self.last_output__forward_inst = (rad_out, J_out)
        return self.last_output__forward_inst

class GaussNewtonRetrieval(object):
    def __init__(
        self,
        radtran,
        radiance = None,
        inverse_covariance = None,
        covariance = None,
        uncertainty = None
    ):
        self.radtran = radtran
        self.inv_cov = None
        self.radiance = None

        self.device = self.radtran.device
        self.dtype = self.radtran.dtype
        self.epsilon = self.radtran.epsilon

        self.id_params = None

        if (
            inverse_covariance is not None
            or covariance is not None
            or (radiance is not None and uncertainty is not None)
        ):
            self.setup_measurement(
                radiance,
                inverse_covariance,
                covariance,
                uncertainty
            )
        elif isinstance(radiance, torch.Tensor):
            self.radiance = torch.clone(radiance)

    def setup_measurement(
        self,
        radiance,
        inverse_covariance = None,
        covariance = None,
        uncertainty = None
    ):
        """Set up measurement parameters for retrieval.

        The inverse covariance needed to compute χ² can be provided to this
        method as-is, or computed either from the covariance itself, or from
        the radiance and an uncertainty. The parameters are checked in the
        order `inverse_covariance`, `covariance`, `uncertainty`, and the first
        that is not-None will be used.

        `inverse_covariance` and `covariance` may be tensors of 3, 2, 1, or 0
        dimensions. A two-dimensional tensor will be copied to each pixel. A
        one-dimensional tensor is promoted to a diagonal matrix and then
        copied to each pixel. A zero-dimensional tensor is multiplied by the
        identity matrix and copied to each pixel.

        If, instead of the (inverse) covariance, the uncertainty is given, it
        must be a tensor of dimension 2, 1, or 0. If it is of dimension 2, the
        first dimension is the pixel index, and the second is the uncertainty
        per wavelength for that pixel. A one-dimensional `uncertainty`
        represents the uncertainty per wavelength, which is assumed to be
        equal for each pixel. A zero-dimenisional `uncertainty` shall have the
        same uncertainty for each wavelength, for each pixel.

        @param radiance Observed radiance (N × N_obs_wls)
        @param inverse_covariance Inverse covariance matrix (see text)
        @param covariance Covariance matrix (see text)
        @param uncertainty Radiance uncertainty (see text)"""

        self.radiance = torch.clone(radiance)

        if (
            self.inv_cov is None
            and inverse_covariance is None
            and covariance is None
            and uncertainty is None
        ):
            raise ValueError(
                "Need not-None (invverse) covariance or uncertainty!"
            )

        if inverse_covariance is not None:
            if len(inverse_covariance.shape) == 3:
                self.inv_cov = torch.clone(inverse_covariance)
            elif len(inverse_covariance.shape) == 2:
                self.inv_cov = torch.einsum(
                    "N,oO->NoO",
                    torch.ones(
                        self.radtran.N,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.clone(inverse_covariance)
                )
            elif len(inverse_covariance.shape)  == 1:
                self.inv_cov = torch.einsum(
                    "N,oO->NoO",
                    torch.ones(
                        self.radtran.N,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.diag(torch.clone(inverse_covariance))
                )
            elif len(inverse_covariance.shape) == 0:
                self.inv_cov = torch.einsum(
                    "N,oO->NoO",
                    torch.ones(
                        self.radtran.N,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.clone(inverse_covariance) * torch.eye(
                        self.radtran.N_obs_wls,
                        device=self.device, dtype=self.dtype
                    )
                )
            else:
                raise ValueError(
                    "inverse_covariance must be 3, 2, 1, or 0-dimensional"
                )
        elif covariance is not None:
            cov = None
            if len(covariance.shape) == 3:
                cov = torch.clone(covariance)
            elif len(covariance.shape) == 2:
                cov = torch.einsum(
                    "N,oO->NoO",
                    torch.ones(
                        self.radtran.N,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.clone(covariance)
                )
            elif len(covariance.shape) == 1:
                cov = torch.einsum(
                    "N,oO->NoO",
                    torch.ones(
                        self.radtran.N,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.diag(torch.clone(covariance))
                )
            elif len(covariance.shape) == 0:
                cov = torch.einsum(
                    "N,oO->NoO",
                    torch.ones(
                        self.radtran.N,
                        device=self.device, dtype=self.dtype
                    ),
                    torch.clone(covariance) * torch.eye(
                        self.radtran.N_obs_wls,
                        device=self.device, dtype=self.dtype
                    )
                )
            else:
                raise ValueError(
                    "covariance must be 3, 2, 1, or 0-dimensional"
                )
            self.inv_cov = torch.zeros_like(cov)
            id_obs = torch.eye(
                self.radtran.N_obs_wls,
                dtype=self.dtype, device=self.device
            )

            for i in range(self.radtran.N):
                self.inv_cov[i,:,:] = torch.linalg.solve(
                    cov[i],
                    id_obs
                )
        elif uncertainty is not None:
            self.inv_cov = torch.einsum(
                "N,oO->NoO",
                torch.ones(
                    self.radtran.N,
                    device=self.device, dtype=self.dtype
                ),
                torch.eye(
                    self.radtran.N_obs_wls,
                    device=self.device, dtype=self.dtype
                )
            )
            if len(uncertainty.shape) == 2:
                for i in range(self.radtran.N):
                    self.inv_cov[i,:] = (
                        self.inv_cov[i,:]
                        * (uncertainty[i,:]*self.radiance[i,:]) ** -2
                    )
            elif len(uncertainty.shape) in (1, 0):
                for i in range(self.radtran.N):
                    self.inv_cov[i,:] = (
                        self.inv_cov[i,:]
                        * (uncertainty*self.radiance[i,:]) ** -2
                )
        self.inv_cov = torch.clamp(self.inv_cov, -1/self.epsilon, 1/self.epsilon)
        self.inv_cov64 = self.inv_cov.to(torch.float64)

    def inversion_step(
        self,
        column_scales,
        x_alb,
        waveshift,
        delta=False,
        ignore_pixels=None,
        damping_parameter=None,
        fwd_precomp=None,
        aux_matrices=False
    ):
        """Perform a single Gauss-Newton step.
        Inputs are the state vector components; see `forward`.
        If `damping_parameter` is specified, use a Levenberg-Marquardt-style
        inversion, whereby the Hessian J.T S J is offset by a term
        λ*diag[J.T S J].

        To save a call to `forward_inst()`, the forward model output
        corresponding to the input state vector can be supplied using the
        `fwd_precomp` argument. This should be a tuple of
        (radiance, Jacobian), of the shape
        (
            (N × N_obs_wls),
            (N × N_obs_wls × (N_species + N_alb[ + 1]))
        )
        If `fwd_precomp` is None, it will be computed using `forward_inst()`.
        N.B. the shape of the Jacobian in `fwd_precomp` should be
        (N × N_obs_wls × (N_species + N_alb + 1)) if waveshift is not None,
        or (N × N_obs_wls × (N_species + N_alb)) if waveshift is None.

        If waveshift is None, the corresponding element of the return tuple is
        also None.

        If aux_matrices is True, additionally return the state vector
        covariance Ainv = (J.T C^-1 J)^-1 and gain matrix B = Ainv (C^-1 J).T.

        @param column_scales Column scale prior (N × N_species)
        @param x_alb Albedo prior (N × N_alb)
        @param waveshift Wavelength shift prior (N)
        @param delta Set to True to return only the update to the state
        vector, rather than the updated state vector itself
        @param ignore_pixels Tensor of bools that flag pixels to be ignored
        (=not inverted) (N)
        @param damping_parameter Levenberg-Marquardt damping parameter
        (scalar or N)
        @param fwd_precomp Pre-computed forward model; see description above.
        @param gain Set to True to additionally return the gain matrix B.
        @return (
            Updated column scales/update to column scales (N × N_species)
            Updated albedo factors/update to albedo factors (N × N_alb)
            Updated wavelength shifts/update to wavelength shifts (N)
            χ² (N)
            [If aux_matrices is True:
                State vector covariance (N × N_params × N_params)
                Gain matrix (N × N_params × N_obs_wls)
            ]
        )"""
        self.inversion_exceptions = []
        y, J = (
            fwd_precomp if fwd_precomp is not None
            else self.radtran.forward_inst(
                column_scales, x_alb, waveshift, grad=True
            )
        )
        dy = self.radiance - y

        # Identity matrix for dimension of the state vector. Used for
        # torch.linalg.solve, which is preferred over torch.linalg.inv,
        # according to the documentation.
        if self.id_params is None:
            self.id_params = torch.eye(
                self.radtran.N_params,
                device=self.device,
                dtype=self.dtype
            )
            # Mask of off-diagonal components. Used for Levenberg-Marquardt.
            self.params_offdiag = ~torch.eye(
                self.radtran.N_params,
                device=self.device,
                dtype=bool
            )

            self.id_params1 = torch.eye(
                self.radtran.N_params1,
                device=self.device,
                dtype=self.dtype
            )
            self.params_offdiag1 = ~torch.eye(
                self.radtran.N_params1,
                device=self.device,
                dtype=bool
            )

        N_params = (
            self.radtran.N_params if waveshift is None
            else self.radtran.N_params1
        )
        id_params = self.id_params if waveshift is None else self.id_params1
        params_offdiag = (
            self.params_offdiag if waveshift is None else self.params_offdiag1
        )

        # Inverse covariance * Jacobian
        # (N × N_obs_wls × N_obs_wls), (N × N_obs_wls × N_params)
        #   -> (N × N_obs_wls × N_params)
        CJ = torch.einsum("NOo,Nop->NOp", self.inv_cov, J)

        # Invertible part of generalised least-squares
        # A = J.T C^-1 J
        # (N × N_obs_wls × N_params), (N × N_obs_wls × N_params)
        #   -> (N × N_params × N_params)
        A = torch.einsum("NoP,Nop->NPp", J, CJ)

        if damping_parameter is not None:
            if not isinstance(damping_parameter, torch.Tensor):
                damping_parameter = (
                    damping_parameter * torch.ones(
                        self.radtran.N, dtype=self.dtype, device=self.device
                    )
                )

            if torch.all(damping_parameter > self.epsilon):
                diag = torch.clone(A)
                diag[:,params_offdiag] = 0
                A += torch.einsum(
                    "N,NPp->NPp",
                    damping_parameter,
                    diag
                )

        Ainv = torch.zeros_like(A)

        # FIXME vectorised inversion?
        # (N × N_params × N_params) -> (N × N_params × N_params)
        for i in range(self.radtran.N):
            if ignore_pixels is None or ignore_pixels[i] == False:
                try:
                    Ainv[i,:,:] = torch.linalg.solve(A[i,:,:], id_params)
                except Exception as e:
                    self.inversion_exceptions.append((i, e))
                    Ainv[i,:,:] = torch.nan
            else:
                Ainv[i,:,:] = torch.nan

        # B = (J.T C^-1 J)^-1 (C^-1 J).T
        # (N × N_params × N_params), (N × N_obs_wls × N_params)
        #   -> (N × N_params × N_obs_wls)
        B = torch.einsum("NPp,Nop->NPo", Ainv, CJ)

        # Δx = B Δy
        # (N × N_params × N_obs_wls), (N, N_obs_wls)
        #   -> (N, N_params)
        x = torch.einsum(
            "Npo,No->Np",
            B,
            dy
        )

        # Update state vector
        if not delta:
            x[:,0:self.radtran.N_species] += column_scales
            x[
                :,
                self.radtran.N_species
                    :self.radtran.N_species+self.radtran.N_alb
            ] += (
                x_alb
            )
            if waveshift is not None:
                x[:,-1] += waveshift

        # (N × N_obs_wls), (N, N_obs_wls, N_obs_wls), (N, obs_wls) ->
        #   (N)
        chi2 = torch.einsum(
            "NO,NOo,No->N",
            dy,
            self.inv_cov,
            dy
        ) / (self.radtran.N_obs_wls - N_params)

        # Sometimes computing the chi2 goes wrong with float32, producing a
        # NaN. Although somewhat kludgy, we can just recompute it with float64
        # when that happens.
        if self.dtype != torch.float64 and torch.any(torch.isnan(chi2)):
            dy64 = dy.to(torch.float64)
            chi2 = (
                torch.einsum(
                    "NO,NOo,No->N",
                    dy64,
                    self.inv_cov64,
                    dy64
                ) / (self.radtran.N_obs_wls - N_params)
            ).to(self.dtype)

        ret_tup = (
            x[:,0:self.radtran.N_species],
            x[
                :,
                self.radtran.N_species
                    :self.radtran.N_species+self.radtran.N_alb
            ],
            x[:,-1] if waveshift is not None else None,
            chi2
        )

        if not aux_matrices:
            return ret_tup
        else:
            return ret_tup + (Ainv, B)

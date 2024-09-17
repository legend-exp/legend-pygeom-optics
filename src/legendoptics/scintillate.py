from __future__ import annotations

from typing import Literal, NamedTuple, NewType, Optional, get_args, get_type_hints

import numpy as np
import pint
from numba import njit
from pint import Quantity

u = pint.get_application_registry()


class ScintParticle(NamedTuple):
    """Configuration for the scintillation yield relative to the flat-top yield."""

    name: Literal["deuteron", "triton", "alpha", "ion", "electron"]
    yield_factor: float
    exc_ratio: Optional[float]  # noqa: UP007

    def valid_geant_particle(self) -> bool:
        return self.name.lower() in get_args(get_type_hints(ScintParticle)["name"])


class ScintConfig(NamedTuple):
    """Scintillation yield parameters, depending on the particle types."""

    flat_top: Quantity
    fano_factor: Optional[float]  # noqa: UP007
    particles: list[ScintParticle]

    def get_particle(self, name: str) -> Optional[ScintParticle]:  # noqa: UP007
        for p in self.particles:
            if p.name.upper() == name.upper():
                return p
        return None


ParticleIndex = NewType("ParticleIndex", int)
ComputedScintParams = NewType(
    "ComputedScintParams", tuple[float, float, np.ndarray, np.ndarray]
)
PARTICLE_INDEX_ELECTRON = 0
PARTICLE_INDEX_ALPHA = 1
PARTICLE_INDEX_ION = 2
PARTICLE_INDEX_DEUTERON = 3
PARTICLE_INDEX_TRITON = 4
PARTICLE_INDICES = {
    "electron": PARTICLE_INDEX_ELECTRON,
    "alpha": PARTICLE_INDEX_ALPHA,
    "ion": PARTICLE_INDEX_ION,
    "deuteron": PARTICLE_INDEX_DEUTERON,
    "triton": PARTICLE_INDEX_TRITON,
}


def particle_to_index(particle: str) -> ParticleIndex:
    """Converts the given G4 scintillation particle name to a module-internal index."""
    return PARTICLE_INDICES[particle.lower()]


def precompute_scintillation_params(
    scint_config: ScintConfig,
    time_components: tuple[Quantity, ...],
) -> ComputedScintParams:
    # validate that we have a valid configuration of scintillation parameters.
    part_with_exc_ratio = [1 for p in scint_config.particles if p.exc_ratio is not None]
    if len(part_with_exc_ratio) not in (0, len(scint_config.particles)):
        msg = "Not all particles have exc_ratio defined"
        raise ValueError(msg)
    has_2_timeconstants = len(part_with_exc_ratio) > 0
    if len(time_components) != (1 + int(has_2_timeconstants)):
        msg = "Number of time_components is not as required"
        raise ValueError(msg)

    particles = np.zeros((len(PARTICLE_INDICES), (2 + int(has_2_timeconstants))))

    electron = scint_config.get_particle("electron")
    if electron is None:
        msg = "missing electron scintillation particle component (used as fallback for all other)"
        raise ValueError(msg)
    for k, v in PARTICLE_INDICES.items():
        p = scint_config.get_particle(k)
        p = electron if p is None else p
        if has_2_timeconstants:
            particles[v] = np.array([p.yield_factor, p.exc_ratio, 1 - p.exc_ratio])
        else:
            particles[v] = np.array([p.yield_factor, 1])

    time_components = np.array([t.to(u.nanosecond).m for t in time_components])
    fano = scint_config.fano_factor if scint_config.fano_factor is not None else 1

    return scint_config.flat_top.to("1/keV").m, fano, time_components, particles


@njit
def scintillate_local(
    params: ComputedScintParams,
    particle: ParticleIndex,
    edep_keV: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generates a Poisson/Gauss-distributed number of photons according to the
    scintillation yield formula, as implemented in Geant4.

    This function only calculates the local part of scintillation.

    Parameters
    ----------
    params
        scintillation parameter tuple, as created by :meth:`precompute_scintillation_params`.
    particle
        module-internal particle index, see :meth:`particle_to_index`.
    edep_keV
        energy deposition along this step, in units pf keV.

    Returns
    -------
    array of scintillation time stamps in nanoseconds, relative to the point in
    time of energy deposition.
    """
    flat_top, fano, time_components, particles = params

    # get the particle scintillation info.
    if particle < 0 or particle > particles.shape[0]:
        msg = "unknown particle index"
        raise IndexError(msg)
    part = particles[particle]
    yield_factor = part[0]
    yields = part[1:]

    mean_num_phot = flat_top * yield_factor * edep_keV

    # derive the actual number of photons generated in this step.
    if mean_num_phot > 10:
        sigma = np.sqrt(fano * mean_num_phot)
        num_photons = int(rng.normal(mean_num_phot, sigma) + 0.5)
    else:
        num_photons = rng.poisson(mean_num_phot)

    if num_photons <= 0:
        return np.empty((0,))

    # derive number of photons for all time components.
    yields = (num_photons * yields).astype(np.int64)
    yields[-1] = num_photons - np.sum(yields[0:-1])  # to keep the sum constant.

    # now, calculate the timestamps of each generated photon.
    times = np.log(rng.uniform(size=num_photons))
    start = 0
    for num_phot, scint_t in zip(yields, time_components):
        times[start : start + num_phot] *= -scint_t
        start += num_phot

    return times


@njit
def scintillate(
    params: ComputedScintParams,
    x0_m: np.ndarray,
    x1_m: np.ndarray,
    v0_mpns: float,
    v1_mpns: float,
    t0_ns: float,
    particle: ParticleIndex,
    particle_charge: int,
    edep_keV: float,
    rng: np.random.Generator,
):
    """Generates a Poisson/Gauss-distributed number of photons according to the
    scintillation yield formula, as implemented in Geant4, along the line segment
    between x0 and x1.

    Parameters
    ----------
    params
        scintillation parameter tuple, as created by :meth:`precompute_scintillation_params`.
    x0_m
        three-vector of pre step point position (in units of meter).
    x1_m
        three-vector of post step point position (in units of meter).
    v0_mpns
        velocity of particle before step (in units of meter per nanosecond).
    v1_mpns
        velocity of particle after step (in units of meter per nanosecond).
    t0_ns
        global time offset of step start in scintillator, in nanoseconds.
    particle
        module-internal particle index, see :meth:`particle_to_index`.
    particle_charge
        charge of the particle, in units of the elementary charge.
    edep_keV
        energy deposition along this step, in units pf keV.

    Returns
    -------
    array of four-vectors containing time stamps (in nanoseconds) and global scintillation
    positions (in meter).
    """
    delta_t_scint = scintillate_local(params, particle, edep_keV, rng)
    # emission position for each single photon.
    if particle_charge != 0:
        λ = rng.uniform(size=delta_t_scint.shape[0])
    else:
        λ = np.ones(shape=delta_t_scint.shape[0])

    x = np.empty((delta_t_scint.shape[0], 4))
    # spatial components along the line segment between x0 and x1.
    x[:, 1] = x0_m[0] + λ * (x1_m[0] - x0_m[0])
    x[:, 2] = x0_m[1] + λ * (x1_m[1] - x0_m[1])
    x[:, 3] = x0_m[2] + λ * (x1_m[2] - x0_m[2])

    # time component along the velocity decrease between x0 and x1.
    x[:, 0] = np.linalg.norm(x1_m - x0_m) / (v0_mpns + λ * (v1_mpns - v0_mpns) / 2)
    # add the global time offset and the emission time offset.
    x[:, 0] += t0_ns + delta_t_scint

    return x

"""Gibbs distributions and exact heat-bath transition matrices."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .ising import IsingModel

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class MarkovChain:
    beta: float
    probabilities: FloatArray
    transition: FloatArray
    discriminant: FloatArray
    eigenvalues: FloatArray
    eigenvectors: FloatArray
    spectral_gap: float
    walk_phase_gap: float
    row_residual: float
    stationarity_residual: float
    detailed_balance_residual: float


def gibbs_distribution(energies: FloatArray, beta: float) -> FloatArray:
    log_weights = -float(beta) * np.asarray(energies, dtype=np.float64)
    log_weights -= np.max(log_weights)
    weights = np.exp(log_weights)
    return weights / weights.sum()


def heat_bath_transition(model: IsingModel, beta: float) -> FloatArray:
    """Exact random-scan single-site heat-bath transition matrix."""

    if not np.isfinite(beta):
        raise ValueError("beta must be finite")
    n_states = model.n_states
    transition = np.zeros((n_states, n_states), dtype=np.float64)
    energies = model.energies()
    rows = np.arange(n_states)
    for site in range(model.n_sites):
        minus = model.flipped_indices(site, -1)
        plus = model.flipped_indices(site, 1)
        e_minus = energies[minus]
        e_plus = energies[plus]
        shift = np.minimum(e_minus, e_plus)
        w_minus = np.exp(-beta * (e_minus - shift))
        w_plus = np.exp(-beta * (e_plus - shift))
        norm = w_minus + w_plus
        np.add.at(transition, (rows, minus), w_minus / (model.n_sites * norm))
        np.add.at(transition, (rows, plus), w_plus / (model.n_sites * norm))
    return transition


def lazify(transition: FloatArray, laziness: float) -> FloatArray:
    if not np.isfinite(laziness) or not 0.0 <= laziness <= 1.0:
        raise ValueError("laziness must lie in [0, 1]")
    if transition.ndim != 2 or transition.shape[0] != transition.shape[1]:
        raise ValueError("transition matrix must be square")
    n_states = transition.shape[0]
    return laziness * np.eye(n_states) + (1.0 - laziness) * transition


def symmetrized_discriminant(
    transition: FloatArray, probabilities: FloatArray
) -> FloatArray:
    if not np.all(np.isfinite(probabilities)) or np.any(probabilities <= 0.0):
        raise ValueError(
            "Gibbs probabilities must be strictly positive and representable"
        )
    roots = np.sqrt(probabilities)
    discriminant = roots[:, None] * transition / roots[None, :]
    return (discriminant + discriminant.T) / 2.0


def build_chain(
    model: IsingModel,
    beta: float,
    laziness: float = 0.5,
) -> MarkovChain:
    probabilities = gibbs_distribution(model.energies(), beta)
    transition = lazify(heat_bath_transition(model, beta), laziness)
    discriminant = symmetrized_discriminant(transition, probabilities)
    eigenvalues, eigenvectors = np.linalg.eigh(discriminant)
    eigenvalues = np.clip(eigenvalues, -1.0, 1.0)
    lambda_two = float(eigenvalues[-2])
    spectral_gap = 1.0 - lambda_two
    walk_phase_gap = float(np.arccos(lambda_two))
    row_residual = float(np.max(np.abs(transition.sum(axis=1) - 1.0)))
    propagated = np.sum(probabilities[:, None] * transition, axis=0)
    stationarity_residual = float(np.max(np.abs(propagated - probabilities)))
    detailed_balance_residual = float(
        np.max(
            np.abs(
                probabilities[:, None] * transition
                - probabilities[None, :] * transition.T
            )
        )
    )
    return MarkovChain(
        beta=float(beta),
        probabilities=probabilities,
        transition=transition,
        discriminant=discriminant,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        spectral_gap=spectral_gap,
        walk_phase_gap=walk_phase_gap,
        row_residual=row_residual,
        stationarity_residual=stationarity_residual,
        detailed_balance_residual=detailed_balance_residual,
    )


def squared_qsample_overlap(left: FloatArray, right: FloatArray) -> float:
    amplitude = float(np.dot(np.sqrt(left), np.sqrt(right)))
    return amplitude * amplitude

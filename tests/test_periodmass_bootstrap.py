import pytest
import numpy as np

import tess_young
from tess_young.tau_sq.periodmass import PeriodMassDistribution, PeriodMassBootstrap
from tess_young.tau_sq.spinmodel import SpinModel

pytestmark = pytest.mark.parametrize('max_q,include_blends,include_lit,mass_limits', 
                         [(0,True,True,None), 
                          (1,False,True,None),
                          (0,False,True,[0.05,0.65]), 
                          (1,True,True,[0.05,0.65]), 
                          (0,False,True,[0.65,1.35]), 
                          (1,True,True,[0.65,1.35])])

def test_bs_pmd(max_q,include_blends,include_lit,mass_limits):
    """
    Ensure that the masses produced by the bootstrap resampling
    do NOT equal the original masses
    """
    pmd_obs = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)

    pmd = PeriodMassBootstrap(pmd_obs,mass_limits)

    obs_mass = np.sort(pmd_obs.mass_raw[pmd_obs.qmask])
    mass = np.sort(pmd.mass_raw[pmd.qmask])

    assert np.all(obs_mass!=pytest.approx(mass))

def test_select_obs(max_q,include_blends,include_lit,mass_limits):
    """
    Ensure that the masses produced by the bootstrap resampling
    do NOT equal the original masses
    """
    pmd_obs = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)

    pmd = PeriodMassBootstrap(pmd_obs,mass_limits)

    sm = SpinModel("WideHat8Myr_Mattea2022",80,"log","kde")

    pmd_obs.select_obs(sm)
    pmd.select_obs(sm)

    obs_mass = np.sort(pmd_obs.mass)
    mass = np.sort(pmd.mass)

    assert np.all(obs_mass!=pytest.approx(mass))
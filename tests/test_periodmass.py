import pytest

import tess_young
from tess_young.tau_sq.periodmass import PeriodMassDistribution, PeriodMassModel

pytestmark = pytest.mark.parametrize('max_q,include_blends,include_lit,mass_limits', 
                         [(0,True,True,None), 
                          (0,True,False,None), 
                          (0,False,True,None), 
                          (1,True,True,None), 
                          (1,True,False,None), 
                          (1,False,True,None),
                          (0,True,True,[0.05,0.65]), 
                          (0,True,False,[0.05,0.65]), 
                          (0,False,True,[0.05,0.65]), 
                          (1,True,True,[0.05,0.65]), 
                          (1,True,False,[0.05,0.65]), 
                          (1,False,True,[0.05,0.65]),
                          (0,True,True,[0.65,1.35]), 
                          (0,True,False,[0.65,1.35]), 
                          (0,False,True,[0.65,1.35]), 
                          (1,True,True,[0.65,1.35]), 
                          (1,True,False,[0.65,1.35]), 
                          (1,False,True,[0.65,1.35]),])

def test_create_pmd(max_q,include_blends,include_lit,mass_limits):
    pmd = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)

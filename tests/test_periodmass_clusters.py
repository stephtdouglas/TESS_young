import pytest

import tess_young
from tess_young.tau_sq.periodmass import PeriodMassDistribution, PeriodMassModel

def test_IC_2391_pmd():
    pmd = PeriodMassDistribution(0,True,False,None,"IC_2391")
    assert len(pmd.prot)>0

def test_IC_2602_pmd():
    pmd = PeriodMassDistribution(0,True,False,None,"IC_2602")
    assert len(pmd.prot)>0

def test_NGC_2547_pmd():
    pmd = PeriodMassDistribution(0,True,False,None,"NGC_2547")
    assert len(pmd.prot)>0

def test_NGC_2451A_pmd():
    pmd = PeriodMassDistribution(0,True,False,None,"NGC_2451A")
    assert len(pmd.prot)>0

def test_Collinder_135_pmd():
    pmd = PeriodMassDistribution(0,True,False,None,"Collinder_135")
    assert len(pmd.prot)>0


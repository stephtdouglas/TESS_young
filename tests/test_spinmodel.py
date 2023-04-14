import pytest

import tess_young
from tess_young.tau_sq.spinmodel import SpinModel

pytestmark = pytest.mark.parametrize('model,init_type', 
                         [("UpSco_Mattea2015","cluster"), 
                           ("UpSco_Mattea2022","cluster"), 
                           ("UpSco_ZeroTorque","cluster"), 
                           ("WideHat8Myr_Mattea2015","tophat"),
                           ("WideHat8Myr_Mattea2022","tophat"), 
                           ("WideHat8Myr_ZeroTorque","tophat")])
def test_spinmodel_creation_linear(model,init_type):
    sm = SpinModel(model,80,"linear",init_type)

def test_spinmodel_creation_log(model,init_type):
    sm = SpinModel(model,80,"log",init_type)

def test_period_scale_linear(model,init_type):
    sm = SpinModel(model,80,"linear",init_type)
    assert sm.period_scale=="linear"

def test_are_linear_bins_even(model,init_type):
    sm = SpinModel(model,80,"linear",init_type)
    l1 = sm.period_bins[-1] - sm.period_bins[-2]
    l2 = sm.period_bins[1] - sm.period_bins[0]
    assert l1==pytest.approx(l2)

def test_period_scale_log(model,init_type):
    sm = SpinModel(model,80,"log",init_type)
    assert sm.period_scale=="log"

def test_are_log_bins_uneven(model,init_type):
    sm = SpinModel(model,80,"log",init_type)
    l1 = sm.period_bins[-1] - sm.period_bins[-2]
    l2 = sm.period_bins[1] - sm.period_bins[0]
    assert l1!=l2
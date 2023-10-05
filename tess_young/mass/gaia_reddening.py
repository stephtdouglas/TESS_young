"""
Convert reddening from E(B-V) to Gaia colors

Stephanie T. Douglas, including snippets from Jason L. Curtis
"""

import pathlib, os
import numpy as np

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent

# AG might be DR2, not sure if I updated it. can check -JLC
def AG_Law(colors, Av):
  # 0.9761 −0.1704 0.0086 0.0011 −0.0438 0.0013 0.0099
  coeff = np.array([0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099])
  return coeff[0] + coeff[1]*colors + coeff[2]*colors**2 + coeff[3]*colors**3 + coeff[4]*Av + coeff[5]*Av**2 + coeff[6]*Av*colors

# -JLC
# there are two lines for BP and RP because the first was for DR2 and the second is for DR3. 
# Didn't code up a keyword to opt in, didn't comment it out.
def EBR_Law(colors, Av):
  # coeff_BP = np.array([1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043])
  coeff_BP = np.array([1.15363197483424,-0.0814012991657388,-0.036013023976704,0.0192143585568966,-0.022397548243016,0.000840562680547171,-1.31018008013549e-05,0.00660124080271006,-0.000882247501989453,-0.000111215755291684])

  # coeff_RP = np.array([0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006])
  coeff_RP = np.array([0.66320787941067,-0.0179847164933981,0.000493769449961458,-0.00267994405695751,-0.00651422146709376,3.30179903473159e-05,1.57894227641527e-06,-7.9800898337247e-05,0.000255679812110045,1.10476584967393e-05])
  ABP = coeff_BP[0] + coeff_BP[1]*colors + coeff_BP[2]*colors**2 + coeff_BP[3]*colors**3 + coeff_BP[4]*Av + coeff_BP[5]*Av**2 + coeff_BP[6]*Av*colors
  ARP = coeff_RP[0] + coeff_RP[1]*colors + coeff_RP[2]*colors**2 + coeff_RP[3]*colors**3 + coeff_RP[4]*Av + coeff_RP[5]*Av**2 + coeff_RP[6]*Av*colors
  return (ABP-ARP)

def calc_BP_RP0(bp_rp,eb_v,r_v=3.1):
    """
    Convert from a B-V extinction to Gaia colors

    Inputs:
    -------
        bp_rp: star's BP-RP color
        eb_v: E(B-V) value
        r_v: assumed R_V value (default=3.1)
    """

    # so here's how to do it.
    # my_star_bp_rp = 0.9
    # my_star_A_V = 0.1
    # (1) ebr_coeff = EBR_Law(my_star_bp_rp - 0.42*my_star_A_V, my_star_A_V)
    # (2) my_star_bp_rp_0 = my_star_bp_rp - ebr_coeff * my_star_A_V
    # so the answer to my example i think is 0.8578614766263766

    a_v = r_v * eb_v

    approx_bp_rp0 = bp_rp - 0.42*a_v

    ebr_coeff = EBR_Law(approx_bp_rp0, a_v)

    bp_rp0 = bp_rp - ebr_coeff * a_v

    # don't extrapolate beyond -0.06 < (GBP − GRP )0 < 2.5 mag
    bad = (bp_rp0<-0.06) | (bp_rp0>2.5) | np.isnan(eb_v)
    bp_rp0[bad] = np.nan

    return bp_rp0

if __name__=="__main__":

    a_v = 0.1
    eb_v = a_v/3.1

    print(calc_BP_RP0(bp_rp=0.9,eb_v=eb_v))
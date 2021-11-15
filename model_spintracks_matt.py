import os, glob

import astropy.io.ascii as at
import numpy as np

def make_spintrack(mass,mass_tol,per,model):
    """
    mass should be array-like, run multiple at once
    ditto with per
    """
    mod_dir = os.path.join(os.path.expanduser("~/Dropbox/Models/"), model)
    mod_list = np.sort(glob.glob(os.path.join(mod_dir,"*Myr.txt")))
    nmod = len(mod_list)
    if isinstance(mass,float):
        float_mass=True
        spintrack = np.zeros(nmod)
    else:
        float_mass=False
        spintrack = np.zeros((nmod,len(mass)))
    ages = np.zeros(nmod)

    # Use the first file to find the mass index we want
    # This isn't actually what we want. Because this doesn't control the period
    # So need to find all the stars within my mass range, and then find the
    # period closest to what I want
    if float_mass is True:
        mfile = mod_list[0]
        modl = at.read(mfile,names=["mass","per"])
        mass_idx = np.where((modl["mass"]<=(mass+mass_tol)) &
                            (modl["mass"]>=(mass-mass_tol)))[0]
        # print(modl[mass_idx])

        per_midx = np.argmin(np.abs(modl["per"][mass_idx]-per))
        per_idx = mass_idx[per_midx]

        minit = modl["mass"][per_idx]
        pinit = modl["per"][per_idx]

        print(modl[per_idx])

        for i,mfile in enumerate(mod_list):
            ages[i] = float(mfile.split("/")[-1].split("_")[2].split(".")[0][:5])
            modl = at.read(mfile,names=["mass","per"])
            spintrack[i] = modl["per"][per_idx]

        at.write({"per":spintrack,"age(Myr)":ages},f"models/{model}_M{minit:.2f}_P{pinit:.2f}.csv",
                 delimiter=",",overwrite=True)

    # return spintrack



if __name__=="__main__":

    for model in ["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque"]:
        for pinit in [0.826,2.391,4.027]:
            print(model,pinit)
            make_spintrack(1.0,0.05,pinit,model)

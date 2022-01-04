import os, sys, time

import astropy.io.ascii as at
import numpy as np
import astropy.units as u
from astropy.table import join, Column
from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch

if __name__=="__main__":
    # res19 = at.read("master_table_compZorro.csv")
    gem_dir = os.path.expanduser("~/proposals/Gemini/2021A/")
    obs19 = at.read(os.path.join(gem_dir,"Douglas_Speckle_Observations_withTIC.csv"))
    det19 = at.read(os.path.join(gem_dir,"Douglas_Binaries_2019AB.txt"),delimiter="\s")
    res19 = join(obs19,det19,join_type="left",keys=["ID"])
    # res19.dtype

    rotfile = os.path.expanduser("jose_master_table.csv")
    rot19 = at.read(rotfile)
    print(rot19.dtype)
    rot19.rename_column("tic_id","TIC")
    res_rot = join(res19,rot19,keys=["TIC"],join_type="right")
    print(res_rot.dtype)
    res_rot.add_columns([Column(name="ra_tic",data=np.ones(len(res_rot))*-99),
                         Column(name="dec_tic",data=np.ones(len(res_rot))*-99)])

    for i,row in enumerate(res_rot):
        print(i)
        result = Vizier.query_region("TIC {0}".format(row["TIC"]),
                                     radius=2*u.arcsec, catalog='IV/38/tic')
    #     print(result[0])
        if len(result[0])==1:
            if row["TIC"]!=result[0]["TIC"]:
                print("UH OH")
                print(result[0])
            else:
                row["ra_tic"] = result[0]["RAJ2000"]
                row["dec_tic"] = result[0]["DEJ2000"]
        else:
            print("UH OH")
            print(result[0])
            # match on TIC id
            loc = row["TIC"]==result[0]["TIC"]
            print(loc,result[0]["TIC"][loc])
            row["ra_tic"] = result[0]["RAJ2000"][loc]
            row["dec_tic"] = result[0]["DEJ2000"][loc]

        time.sleep(2)

    table2 = XMatch.query(cat1=res_rot,
                         cat2='vizier:I/350/gaiaedr3',
                         max_distance=1 * u.arcsec,
                         colRA1='ra_tic',colDec1='dec_tic')

    tic_id_uniq, ct_uniq = np.unique(table2["TIC"],
                                     return_counts=True)
    print(tic_id_uniq[ct_uniq>1])

    for tid in tic_id_uniq[ct_uniq>1]:
        print(table2["angDist","TIC","phot_g_mean_mag","ruwe"][table2["TIC"]==tid])

    for tid in tic_id_uniq[ct_uniq>1]:
        to_del = np.where(table2["TIC"]==tid)[0]
        print(to_del)
        table2.remove_row(to_del[1])
    rr_edr3 = table2



    plt.figure(figsize=(8,7))
    ax = plt.subplot(111)
    # ax.plot(res_rot["j-k"],res_rot["period"],'k.')
    def_bad = rr_edr3["period"]>27.5
    sc = ax.scatter(rr_edr3["bp_rp"][~def_bad],
                    rr_edr3["period"][~def_bad],
                    c=rr_edr3["ruwe"][~def_bad],
                    vmin=0.6,vmax=2,label="Preliminary TESS period",
                    alpha=0.4)

    obs = (rr_edr3["ID"].mask==False) & (def_bad==False)
    sc = ax.scatter(rr_edr3["bp_rp"][obs],rr_edr3["period"][obs],
                    c=rr_edr3["ruwe"][obs],vmin=0.6,vmax=2,
                    label="Speckle observation")
    plt.colorbar(sc,label="RUWE - EDR3")

    det = (rr_edr3["notes"].mask==False) & (def_bad==False)
    ax.plot(rr_edr3["bp_rp"][det],rr_edr3["period"][det],
            'ko',mfc="none",ms=11,label="Speckle detection")

    high_ruwe = (rr_edr3["ruwe"]>=1.2) & (def_bad==False)
    ax.plot(rr_edr3["bp_rp"][high_ruwe],rr_edr3["period"][high_ruwe],
            'ks',mfc="none",ms=12,label=r"RUWE > 1.4")


    # ax.set_yscale("log")
    # ax.set_ylim(0.1,25)
    # ax.set_yticklabels([0,0.1,1,10])
    # ax.set_xlabel("J-Ks")
    # ax.set_ylabel("Period (d)")
    # plt.xlim(0,1)
    # ax.axhline(28,zorder=-10)


    ax.set_ylim(0.1,50)
    ax.set_yscale('log')
    ax.set_ylabel('Period (d)',fontsize=18)
    ax.tick_params(labelsize=16,width=2)
    ax.tick_params(which="minor",width=2)
    ax.set_yticklabels(["","0.1","1","10"])
    ax.set_xlabel("(BP - RP)",fontsize=18)
    plt.xlim(0.5,3.75)

    ques_limit = plt.Rectangle((0.5,12),3.25,38,fill=True,
                          color="C1",alpha=0.25)
    ax.add_patch(ques_limit)

    fontdict = {"weight":"bold", "size":14}

    ax.text(0.55,22,"Likely spurious",color="C1",fontdict=fontdict)


    plt.legend(loc=4)
    plt.title("IC 2391",fontsize=16)

    plt.show()
    # plt.savefig("IC2391_periodcolor_EDR3RUWE.png",bbox_inches="tight",dpi=600)

#python3
import sys
from altaipony.lcio import from_mast

if __name__ == "__main__":
    EPIC = sys.argv[1]
    campaign = int(sys.argv[2])
    i = sys.argv[3]
    print("Fetch EPIC {}, C {}".format(EPIC, campaign, i))
    path = '/home/eilin/FlaresInClusters/data/injrec/'
    flc = from_mast(EPIC, mission="K2", mode="TPF", c=campaign)
    print("Run injection-recovery iteration {}:\n -------------------------------------".format(i))
    flc, fake_flc = flc.sample_flare_recovery(inject_before_detrending=True, mode="k2sc",
                                              iterations=100, fakefreq=1e-3, ampl=[1e-4, 0.5],
                                              dur=[.001/6., 0.2/6.], save=True,
                                              path="{}_{}.csv".format(path, i, flc.targetid))

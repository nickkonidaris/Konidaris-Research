
import os
import numpy as np
import parse

try:
    f = open("2021-09-17.txt")
    lines = f.readlines()
    f.close()
except:
    print("BAD")



section = -1
alldat = []
for line in lines:
    new = False

    r = parse.parse("Ghost reflection off surface {} then {}. {}", line)
    if r is not None:
        section = 0
        res = [int(r[0]), int(r[1])]
        continue
    
    if section == -1: continue # Burn a few lines

    section += 1
    #print(section, line.rstrip())
    if section == 1: continue
    if section == 2:
        sp = line.split()
        spn = map(float,sp)
        for sp in spn: res.append(sp)
    if (section >= 3) and (section <= 7):
        sp = line.split(":")
        it = float(sp[1])
        res.append(it)
    
    if section == 8:
        alldat.append(res)

alldat = np.array(alldat)


s1, s2, _, marg, fno, rms, mrh, crh, dgp, dgf, efl = alldat.T
s1 = int(s1)
s2 = int(s2)






    


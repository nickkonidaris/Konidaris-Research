import numpy as np

from matplotlib.patches import Rectangle

import matplotlib.pyplot as plt

R = 150.
# Delta lambda / lambda = d ln lambda = 1 / R --> d log lambda = ln 10 / R
dloglam = 1 / R / np.log(10)
# boundaries of bins:
wavearr = np.logspace(np.log10(0.60), np.log10(2.40), np.round(np.log10(2.40 / 0.60) / dloglam).astype(int))
wavecen = (wavearr[1:] + wavearr[:-1]) / 2

def process(filename):
   w, f = np.loadtxt(filename, unpack=True)
   w *= 1e6 # to microns
   binned = np.zeros(wavecen.size)
   for i in range(wavecen.size):
      binned[i] = np.mean(f[(w > wavearr[i]) & (w < wavearr[i+1])])
   return binned

def plotwindows(ax): # b/w MOSFIRE filters
   #ax.axvspan(0.8, 0.95, alpha=0.3, color='gray')
   ax.axvspan(1.124, 1.153, alpha=0.3, color='gray')
   ax.axvspan(1.352, 1.466, alpha=0.3, color='gray')
   ax.axvspan(1.708, 1.921, alpha=0.3, color='gray')

   r = Rectangle([0.6, 8100.], 1.8, 500, alpha=.5)


# H mag = 9.7 Vega
# Transit is about 2 hour
scale = 10**(-0.4*(9.7-9.0))
poi = 1/np.sqrt(251000. * 2 * 3600. * scale) * 1e6
print("Poisson limit (ppm) per transit per channel: ",poi)

sc = 0.8
fig = plt.figure(figsize=(8*sc,4*sc))
ax = plt.gca()
#ax.set_aspect(.00025)

import glob
colors = ['b', 'r']
sys = 50. # ppm floor
np.random.seed(8675309)
for i, f in enumerate(['model.dat']):
   fl = process(f)
   fl *= 1e6
   # Now bin up the I band and remove points between bands
   inI = (wavecen < 0.8)
   simobs = np.concatenate([[np.mean(fl[inI])], fl[~inI]])
   w = np.concatenate([[0.7], wavecen[~inI]])
   noise = np.repeat(np.sqrt(sys**2 + poi**2), w.size)
   noise[0] = sys # broadband, poisson noise negligible


   offset = 0
   plt.plot(wavecen, fl + offset, color=colors[i], alpha=0.7)
   noisereal = noise * np.random.normal(size=simobs.size)
   if i != 0: continue
   for bands in [[0.89, 1.124], [1.153, 1.352], [1.466, 1.708], [1.921, 2.4]]:
      inband = (w > bands[0]) & (w < bands[1])
      plt.errorbar(w[inband], (simobs + noisereal)[inband] + offset, yerr=noise[inband], marker='o', color='0.2', linestyle='none',
          ms=3, label=('MIRMOS (50 ppm)' if bands[0]==0.89 else None))
   plt.errorbar([w[0]], [(simobs + noisereal)[0] + offset], yerr=[noise[0]], xerr=[0.05], marker='o', color='0.1', ms=3)

   sys = 60.
   offset = 700.
   noise = np.repeat(np.sqrt(sys**2 + poi**2), w.size)
   noise[0] = sys # broadband, poisson noise negligible
   plt.plot(wavecen, fl + offset, color=colors[i], alpha=0.7)
   noisereal = noise * np.random.normal(size=simobs.size)
   if i != 0: continue
   col = '0.0'
   #for bands in [[0.89, 1.124], [1.153, 1.352], [1.466, 1.708], [1.921, 2.4]]:
   oob = ((w>1.124) & (w<1.153)) | ((w>1.352) & (w<1.466)) | ((w>1.708) &
                                                            (w<1.921))
   noise[oob] *= 5
   for bands in [[0.89, 1.80]]:
      inband = (w > bands[0]) & (w < bands[1])
      plt.errorbar(w[inband], (simobs + noisereal)[inband] + offset, yerr=noise[inband], marker='o', color=col, linestyle='none',
          ms=2, label=('Henrietta (50 ppm)' if bands[0]==0.89 else None))
   #plt.errorbar([w[0]], [(simobs + noisereal)[0] + offset], yerr=[noise[0]], xerr=[0.1], marker='o', color=col, ms=3)

   # Plot HST data simulated...
   lam, dlam, radius, error = np.loadtxt('wakeford.dat', unpack=True)
   model = np.array([np.mean(fl[(wavecen > lam[i]-dlam[i]/2) & (wavecen < lam[i]+dlam[i]/2)]) for i in range(lam.size)])
   noise = ((radius + error)**2 - radius**2) # in fraction of flux
   noise /= 1e-6 # in ppm
   print(noise)
   model += np.random.normal(size=model.size) * noise
   plt.errorbar(lam, model - 600, yerr=noise, marker='s', color='r', linestyle='none', label='HST')
   plt.plot(wavecen, fl - 600, 'b')

   # Plot data for JWST. Note Poisson noise will be the same as MIRMOS scaled by sqrt(throughput / MIRMOS throughput).
   # While MIRMOS throughput is similar ~30% with wavelength (actually can be higher but assumed 30% in the above)
   # JWST SOSS is variable with wavelength.
   # http://jwst.astro.umontreal.ca/?page_id=51
   tput = np.loadtxt('sosstput.dat')
   jwst_poi = poi / np.sqrt(np.interp(wavecen, tput[:,0]/1000., np.clip(tput[:,1], 0.01, 1)) / 0.3)
   print(np.min(jwst_poi), " JWST min Poisson noise")
   jwst_sys = 20. # ppm, Greene et al 2016
   noise = np.sqrt(jwst_poi**2 + jwst_sys**2)
   noisereal = noise * np.random.normal(size=noise.size)
   plt.plot(wavecen, fl - 1200, 'b')
   plt.errorbar(wavecen, fl - 1200 + noisereal, yerr=noise, marker='o', color='g', linestyle='none', ms=3, label='JWST NIRISS')

plotwindows(plt.gca())
plt.xlabel(r'Wavelength [$\mu$m]')
plt.ylabel('Transit depth [ppm]')
#plt.legend()
for wave, name in zip([(0.95+1.124)/2, (1.153+1.352)/2, (1.466+1.708)/2, 2.2], ['Y', 'J', 'H', 'K']):
   plt.gca().text(wave, 7600, name, ha='center', fontsize=(14 if not name.startswith('A') else 8))
#plt.ylim([])
from matplotlib.ticker import MultipleLocator
plt.gca().yaxis.set_minor_locator(MultipleLocator(250))
plt.xlim([0.6,2.4])
plt.legend(loc='lower left', labelspacing=0.4, fontsize=9, framealpha=0)
plt.ylim([3700,7500 + 600])

def plotfeature(name, w1, w2, yoff=0):
   y0, y1 = 4300+yoff-200, 4400+yoff-200
   plt.plot([w1, w1, w2, w2], [y1, y0, y0, y1], 'k')
   plt.text((w1+w2)/2, y0-150, name, ha='center', va='center')
plotfeature(r'${\rm H}_2{\rm O}$', 1.05, 1.21)
plotfeature(r'${\rm H}_2{\rm O}$', 1.3, 1.55)
plotfeature(r'${\rm H}_2{\rm O}$', 1.7, 2.1)
plotfeature(r'CO, ${\rm CH}_4$, ${\rm CH}_3$', 2.15, 2.35)
plotfeature(r'${\rm NH}_3$', 1.2, 1.3, yoff=350)
plotfeature(r'${\rm NH}_3$', 1.85, 2.1, yoff=350)
plotfeature(r'${\rm NH}_3$, HCN, ${\rm CH}_4$', 1.5, 1.7, yoff=350)



#r = Rectangle((0.6,8300), 1.8, 200, alpha=0.2, color='black', clip_on=False)
#plt.gca().add_patch(r)
plt.tight_layout()
plt.savefig('sim.pdf')
plt.show()

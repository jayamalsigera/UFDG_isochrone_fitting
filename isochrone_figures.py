#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 09:08:52 2021
@author: veelochana
"""

from astropy.io import fits
from astropy.table import Table
from pylab import *
import operator
from matplotlib.pyplot import *
from astropy.utils.data import download_file

# 1-Bootes
event_filename_1 = download_file('https://archive.stsci.edu/hlsps/fhufd/hlsp_fhufd_hst_acs_bootesi_multi_v2.0_cat.fits',
                                 cache=True)
hdu_list_1 = fits.open(event_filename_1, memmap=True)
evt_data_1 = Table(hdu_list_1[1].data)
x_1 = evt_data_1['M606'] - evt_data_1['M814']
y_1 = evt_data_1['M814']
L_1 = sorted(zip(x_1, y_1), key=operator.itemgetter(0))
xp_1, yp_1 = zip(*L_1)
data_1 = xp_1, yp_1

fname1_1 = "bootes_he33.iso"
V1_1, I1_1 = np.genfromtxt(fname1_1, dtype=float, comments="#", usecols=(10, 15), unpack=True)

dmod1_1 = 22.1
E1_1 = -0.74
dmod2_1 = 22.1
E2_1 = -0.74

# extinction
Av1_1 = 3.1 * E1_1
Ai1_1 = 1.8 * E1_1

Vfin1_1 = V1_1 + dmod1_1 + Av1_1
Ifin1_1 = I1_1 + dmod1_1 + Ai1_1

X1_1 = Vfin1_1 - Ifin1_1
Y1_1 = Vfin1_1

fname2_1 = "bootes.iso"
V2_1, I2_1 = np.genfromtxt(fname2_1, dtype=float, comments="#", usecols=(10, 15), unpack=True)

# extinction
Av2_1 = 3.1 * E2_1
Ai2_1 = 1.8 * E2_1

Vfin2_1 = V2_1 + dmod2_1 + Av2_1
Ifin2_1 = I2_1 + dmod2_1 + Ai2_1

X2_1 = Vfin2_1 - Ifin2_1
Y2_1 = Vfin2_1

# 2-Canes Venatici II
event_filename_2 = download_file(
    'https://archive.stsci.edu/hlsps/fhufd/hlsp_fhufd_hst_acs_canesvenaticiii_multi_v2.0_cat.fits', cache=True)
hdu_list_2 = fits.open(event_filename_2, memmap=True)
evt_data_2 = Table(hdu_list_2[1].data)
x_2 = evt_data_2['M606'] - evt_data_2['M814']
y_2 = evt_data_2['M814']
L_2 = sorted(zip(x_2, y_2), key=operator.itemgetter(0))
xp_2, yp_2 = zip(*L_2)
data_2 = xp_2, yp_2

fname1_2 = "canes_venatici_he33.iso"
V1_2, I1_2 = np.genfromtxt(fname1_2, dtype=float, comments="#", usecols=(10, 15), unpack=True)

dmod1_2 = 23.9
E1_2 = -0.76
dmod2_2 = 23.9
E2_2 = -0.76

# extinction
Av1_2 = 3.1 * E1_2
Ai1_2 = 1.8 * E1_2

Vfin1_2 = V1_2 + dmod1_2 + Av1_2
Ifin1_2 = I1_2 + dmod1_2 + Ai1_2

X1_2 = Vfin1_2 - Ifin1_2
Y1_2 = Vfin1_2

fname2_2 = "canes_venatici.iso"
V2_2, I2_2 = np.genfromtxt(fname2_2, dtype=float, comments="#", usecols=(10, 15), unpack=True)

# extinction
Av2_2 = 3.1 * E2_2
Ai2_2 = 1.8 * E2_2

Vfin2_2 = V2_2 + dmod2_2 + Av2_2
Ifin2_2 = I2_2 + dmod2_2 + Ai2_2

X2_2 = Vfin2_2 - Ifin2_2
Y2_2 = Vfin2_2

# 3-Coma Berenices
event_filename_3 = download_file(
    'https://archive.stsci.edu/hlsps/fhufd/hlsp_fhufd_hst_acs_comaberenices_multi_v2.0_cat.fits', cache=True)
hdu_list_3 = fits.open(event_filename_3, memmap=True)
evt_data_3 = Table(hdu_list_3[1].data)
x_3 = evt_data_3['M606'] - evt_data_3['M814']
y_3 = evt_data_3['M814']
L_3 = sorted(zip(x_3, y_3), key=operator.itemgetter(0))
xp_3, yp_3 = zip(*L_3)
data_3 = xp_3, yp_3

fname1_3 = "coma_berenices_he33.iso"
V1_3, I1_3 = np.genfromtxt(fname1_3, dtype=float, comments="#", usecols=(10, 15), unpack=True)

dmod1_3 = 20.99
E1_3 = -0.744
dmod2_3 = 20.99
E2_3 = -0.744

# extinction
Av1_3 = 3.1 * E1_3
Ai1_3 = 1.8 * E1_3

Vfin1_3 = V1_3 + dmod1_3 + Av1_3
Ifin1_3 = I1_3 + dmod1_3 + Ai1_3

X1_3 = Vfin1_3 - Ifin1_3
Y1_3 = Vfin1_3

fname2_3 = "coma_berenices.iso"
V2_3, I2_3 = np.genfromtxt(fname2_3, dtype=float, comments="#", usecols=(10, 15), unpack=True)

# extinction
Av2_3 = 3.1 * E2_3
Ai2_3 = 1.8 * E2_3

Vfin2_3 = V2_3 + dmod2_3 + Av2_3
Ifin2_3 = I2_3 + dmod2_3 + Ai2_3

X2_3 = Vfin2_3 - Ifin2_3
Y2_3 = Vfin2_3

# 4-Hercules
event_filename_4 = download_file(
    'https://archive.stsci.edu/hlsps/fhufd/hlsp_fhufd_hst_acs_hercules_multi_v2.0_cat.fits', cache=True)
hdu_list_4 = fits.open(event_filename_4, memmap=True)
evt_data_4 = Table(hdu_list_4[1].data)
x_4 = evt_data_4['M606'] - evt_data_4['M814']
y_4 = evt_data_4['M814']
L_4 = sorted(zip(x_4, y_4), key=operator.itemgetter(0))
xp_4, yp_4 = zip(*L_4)
data_4 = xp_4, yp_4

fname1_4 = "hercules_he33.iso"
V1_4, I1_4 = np.genfromtxt(fname1_4, dtype=float, comments="#", usecols=(10, 15), unpack=True)

dmod1_4 = 23.45
E1_4 = -0.71
dmod2_4 = 23.45
E2_4 = -0.71

# extinction
Av1_4 = 3.1 * E1_4
Ai1_4 = 1.8 * E1_4

Vfin1_4 = V1_4 + dmod1_4 + Av1_4
Ifin1_4 = I1_4 + dmod1_4 + Ai1_4

X1_4 = Vfin1_4 - Ifin1_4
Y1_4 = Vfin1_4

fname2_4 = "hercules.iso"
V2_4, I2_4 = np.genfromtxt(fname2_4, dtype=float, comments="#", usecols=(10, 15), unpack=True)

# extinction
Av2_4 = 3.1 * E2_4
Ai2_4 = 1.8 * E2_4

Vfin2_4 = V2_4 + dmod2_4 + Av2_4
Ifin2_4 = I2_4 + dmod2_4 + Ai2_4

X2_4 = Vfin2_4 - Ifin2_4
Y2_4 = Vfin2_4

# 5-Leo IV
event_filename_5 = download_file('https://archive.stsci.edu/hlsps/fhufd/hlsp_fhufd_hst_acs_leoiv_multi_v2.0_cat.fits',
                                 cache=True)
hdu_list_5 = fits.open(event_filename_5, memmap=True)
evt_data_5 = Table(hdu_list_5[1].data)
x_5 = evt_data_5['M606'] - evt_data_5['M814']
y_5 = evt_data_5['M814']
L_5 = sorted(zip(x_5, y_5), key=operator.itemgetter(0))
xp_5, yp_5 = zip(*L_5)
data = xp_5, yp_5

fname1_5 = "leo_he33.iso"
V1_5, I1v = np.genfromtxt(fname1_5, dtype=float, comments="#", usecols=(10, 15), unpack=True)

dmod1_5 = 23.85
E1_5 = -0.735
dmod2_5 = 23.85
E2_5 = -0.735

# extinction
Av1_5 = 3.1 * E1_5
Ai1_5 = 1.8 * E1_5

Vfin1_5 = V1_5 + dmod1_5 + Av1_5
Ifin1_5 = I1v + dmod1_5 + Ai1_5

X1_5 = Vfin1_5 - Ifin1_5
Y1_5 = Vfin1_5

fname2_5 = "leo.iso"
V2_5, I2_5 = np.genfromtxt(fname2_5, dtype=float, comments="#", usecols=(10, 15), unpack=True)

# extinction
Av2_5 = 3.1 * E2_5
Ai2_5 = 1.8 * E2_5

Vfin2_5 = V2_5 + dmod2_5 + Av2_5
Ifin2_5 = I2_5 + dmod2_5 + Ai2_5

X2_5 = Vfin2_5 - Ifin2_5
Y2_5 = Vfin2_5

# 6-Ursa Major I
event_filename_6 = download_file(
    'https://archive.stsci.edu/hlsps/fhufd/hlsp_fhufd_hst_acs_ursamajori_multi_v2.0_cat.fits', cache=True)
hdu_list_6 = fits.open(event_filename_6, memmap=True)
evt_data_6 = Table(hdu_list_6[1].data)
x_6 = evt_data_6['M606'] - evt_data_6['M814']
y_6 = evt_data_6['M814']
L_6 = sorted(zip(x_6, y_6), key=operator.itemgetter(0))
xp_6, yp_6 = zip(*L_6)
data_6 = xp_6, yp_6

fname1_6 = "ursa_major_he33.iso"
V1_6, I1_6 = np.genfromtxt(fname1_6, dtype=float, comments="#", usecols=(10, 15), unpack=True)

dmod1_6 = 23
E1_6 = -0.76
dmod2_6 = 23
E2_6 = -0.76

# extinction
Av1_6 = 3.1 * E1_6
Ai1_6 = 1.8 * E1_6

Vfin1_6 = V1_6 + dmod1_6 + Av1_6
Ifin1_6 = I1_6 + dmod1_6 + Ai1_6

X1_6 = Vfin1_6 - Ifin1_6
Y1_6 = Vfin1_6

fname2_6 = "ursa_major.iso"
V2_6, I2_6 = np.genfromtxt(fname2_6, dtype=float, comments="#", usecols=(10, 15), unpack=True)

# extinction
Av2_6 = 3.1 * E2_6
Ai2_6 = 1.8 * E2_6

Vfin2_6 = V2_6 + dmod2_6 + Av2_6
Ifin2_6 = I2_6 + dmod2_6 + Ai2_6

X2_6 = Vfin2_6 - Ifin2_6
Y2_6 = Vfin2_6

# ..................................fitting...........................................

fig, ax = subplots()

# 1-Bootes
plt.subplot(2, 3, 1)
plt.scatter(xp_1, yp_1, s=3, marker='D', facecolors='black', edgecolors='none', label="Data")
plt.plot(X1_1, Y1_1, color='r', label="Y = 0.33", linewidth=1.5)
plt.plot(X2_1, Y2_1, color='b', label="Y = 0.245 + 1.5 * Z", linewidth=1.5)
plt.ylabel('M814')
plt.xlabel('M606 - M814')
plt.title('Bootes I')
plt.ylim([18.0, 27.0])
plt.xlim([-0.8, 0.0])
plt.gca().invert_yaxis()

# 2-Canes Venatici II
plt.subplot(2, 3, 2)
plt.scatter(xp_2, yp_2, s=3, marker='D', facecolors='black', edgecolors='none', label="Data")
plt.plot(X1_2, Y1_2, color='r', label="Y = 0.33", linewidth=1.5)
plt.plot(X2_2, Y2_2, color='b', label="Y = 0.245 + 1.5 * Z", linewidth=1.5)
plt.ylabel('M814')
plt.xlabel('M606 - M814')
plt.title('Canes Venatici II')
plt.ylim([18.0, 27.0])
plt.xlim([-0.8, 0.0])
plt.gca().invert_yaxis()

# 3-Coma Berenices
plt.subplot(2, 3, 3)
plt.scatter(xp_3, yp_3, s=3, marker='D', facecolors='black', edgecolors='none', label="Data")
plt.plot(X1_3, Y1_3, color='r', label="Y = 0.33", linewidth=1.5)
plt.plot(X2_3, Y2_3, color='b', label="Y = 0.245 + 1.5 * Z", linewidth=1.5)
plt.ylabel('M814')
plt.xlabel('M606 - M814')
plt.title('Coma Berenices')
plt.ylim([18.0, 27.0])
plt.xlim([-0.8, 0.0])
plt.gca().invert_yaxis()

# 4-Hercules
plt.subplot(2, 3, 4)
plt.scatter(xp_4, yp_4, s=3, marker='D', facecolors='black', edgecolors='none', label="Data")
plt.plot(X1_4, Y1_4, color='r', label="Y = 0.33", linewidth=1.5)
plt.plot(X2_4, Y2_4, color='b', label="Y = 0.245 + 1.5 * Z", linewidth=1.5)
plt.ylabel('M814')
plt.xlabel('M606 - M814')
plt.title('Hercules')
plt.ylim([18.0, 27.0])
plt.xlim([-0.8, 0.0])
plt.gca().invert_yaxis()

# 5-Leo IV
plt.subplot(2, 3, 5)
plt.scatter(xp_5, yp_5, s=3, marker='D', facecolors='black', edgecolors='none', label="Data")
plt.plot(X1_5, Y1_5, color='r', label="Y = 0.33", linewidth=1.5)
plt.plot(X2_5, Y2_5, color='b', label="Y = 0.245 + 1.5 * Z", linewidth=1.5)
plt.ylabel('M814')
plt.xlabel('M606 - M814')
plt.title('Leo IV')
plt.ylim([18.0, 27.0])
plt.xlim([-0.8, 0.0])
plt.gca().invert_yaxis()

# 6-Ursa Major I
plt.subplot(2, 3, 6)
plt.scatter(xp_6, yp_6, s=3, marker='D', facecolors='black', edgecolors='none', label="Data")
plt.plot(X1_6, Y1_6, color='r', label="Y = 0.33", linewidth=1.5)
plt.plot(X2_6, Y2_6, color='b', label="Y = 0.245 + 1.5 * Z", linewidth=1.5)
plt.ylabel('M814')
plt.xlabel('M606 - M814')
plt.title('Ursa Major I')
plt.ylim([18.0, 27.0])
plt.xlim([-0.8, 0.0])
plt.gca().invert_yaxis()

lines, labels = fig.axes[-1].get_legend_handles_labels()

fig.legend(lines, labels, loc='lower center')
plt.show()

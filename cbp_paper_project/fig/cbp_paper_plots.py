# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # CBP paper plots
#
# This notebook produces plots for the CBP paper.

# +
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import savgol_filter
import sys
import scipy
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '../../../stardice-master/analysis/cbp_paper/')

from cbp_dataset import CBPDataSet, SolarCellRun #, calculate_npulses, get_info_from_filename, from_detected_to_calibrated_wavelengths

from cbp_dr4 import *
data_release = "dr4/"

# %matplotlib notebook

#plt_params["figure.dpi"] = 100
plt.rcParams.update(plt_params)

# -

aa_plot(twocols=True)

# ## Solar cell QE

# +
aa_plot(twocols=False, height=3)

ref_data_root = "/data/STARDICE/cbp/solarcell/refCalData/"
solarcell_QE_file = 'SC_QE_NewCell_20220401_LongQECurve-2.txt'
SC_QE_data_file = os.path.join(ref_data_root, solarcell_QE_file)
SC_QE = np.loadtxt(SC_QE_data_file, skiprows=1, delimiter=",").T

solarcell_QE_file_reduced = os.path.join(ref_data_root, 'SC_QE_NewCell_20220401_LongQECurve-2_reduced_filtered.txt')
SC_QE_reduced = np.loadtxt(solarcell_QE_file_reduced, skiprows=1, delimiter=",").T

gs_kw = dict(width_ratios=[1], height_ratios=[3,1])

fig, ax = plt.subplots(2,1, sharex="all", gridspec_kw=gs_kw)
ax[0].errorbar(SC_QE[0], SC_QE[1], yerr=SC_QE[2], linestyle="", color="k",label="raw data")
p = ax[0].plot(SC_QE_reduced[0], SC_QE_reduced[1], label="reduced data")
p = ax[0].fill_between(SC_QE_reduced[0], SC_QE_reduced[1]+SC_QE_reduced[2], SC_QE_reduced[1]-SC_QE_reduced[2], color=p[0].get_color(), alpha=0.5)
#ax[0].grid()
ax[0].set_ylim(0.57, 1)
ax[1].set_xlabel("Calibrated $\lambda$ [nm]")
ax[0].set_ylabel("Solar cell quantum efficiency")
ax[1].set_ylabel("Uncertainties")
#ax[1].errorbar(SC_QE[0], SC_QE[1]-np.interp(SC_QE[0], SC_QE_reduced[0], SC_QE_reduced[1]), yerr=SC_QE[2], linestyle="", color="k",label="raw data")
p = ax[1].plot(SC_QE_reduced[0], SC_QE_reduced[1]-SC_QE_reduced[1], label="reduced data")
p = ax[1].fill_between(SC_QE_reduced[0], SC_QE_reduced[2], -SC_QE_reduced[2], color=p[0].get_color(), alpha=0.5)
ax[1].set_ylim(-8e-2, 8e-2)
ax[1].set_yscale('symlog', linthresh=1e-3)
#ax[1].grid()
ax[0].yaxis.set_label_coords(-0.11,0.55)
fig.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
fig.savefig("solarcell_qe.pdf")
plt.show()

# -

# ## Typical charge sequence

datapath='/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_04_solar_cell_5mm/'
fits_filename = os.path.join(datapath, "NOIMG_0095971.fits")  # 469nm
d = CBPDataSet(fits_filename=fits_filename)
print(d)

# +
aa_plot(twocols=False, height=1.2*len(d.subsystems))

fig, axes = plt.subplots(len(d.subsystems), 1)

d.plot_data_set(axes=axes)
for ax in axes:
    ax.grid(False)
fig.tight_layout()
fig.savefig(f"sc_dataset_{d.wavelength:.0f}.pdf")
# -

fits_filename = os.path.join(datapath, "NOIMG_0095977.fits")  # 966nm
# fits_filename = os.path.join(datapath, "NOIMG_0095974.fits")  # 551nm
d = CBPDataSet(fits_filename=fits_filename)
print(d)

# +
aa_plot(twocols=False, height=1.2*len(d.subsystems))

fig, axes = plt.subplots(len(d.subsystems), 1)

d.plot_data_set(axes=axes)
for ax in axes:
    ax.grid(False)
fig.tight_layout()
fig.savefig(f"sc_dataset_{d.wavelength:.0f}.pdf")

# +
aa_plot(twocols=False, height=2.5)

time = d.timing.data["timing"][d.timing.data["pinstate"]==4] * d.timing.factor - d.pd.keithleytimer
d.pd.data["time"][:len(time)] = time
pd_rows = d.pd.extraction(time_breaks=d.get_time_breaks(), pd_with_digital_analyzer=True, 
                          plot=True, verbose=True)
axes = plt.gcf().axes
for ax in axes:
    ax.grid(False)
plt.gcf().savefig(f"pd_reduc_{d.wavelength:.0f}.pdf")
plt.show()

# +
aa_plot(twocols=False, height=2)

timing_rows = d.timing.extraction()
sc_rows = d.sc.extraction(time_breaks=d.get_time_breaks() - timing_rows["acq_start"][0], 
                          plot=True,
                          fix_time_breaks=True, min_len=10, verbose=True)
axes = plt.gcf().axes
for ax in axes:
    ax.grid(False)
plt.gcf().savefig(f"sc_reduc_{d.wavelength:.0f}.pdf")
plt.show()
# -

fig =  plt.gcf()
ax = fig.axes
ax[0].set_xlim(0.62, 1.06)
ax[0].set_ylim(0.49e-7, 0.53e-7)
ax[1].set_ylim(-2.1e-12, 2.2e-12)
plt.savefig(f"sc_reduc_{d.wavelength:.0f}_zoom.pdf")
#plt.show()

# ## Spectrograph reduction

datapath='/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_04_solar_cell_5mm/'
fits_filename = os.path.join(datapath, "NOIMG_0096008.fits")
d = CBPDataSet(fits_filename=fits_filename)
print(d)

# +
aa_plot(twocols=False, height=2)

d.spectro.get_laser_flux(plot=True)
ax = plt.gca()
ax.set_xlim(510,670)
ax.grid(False)
plt.gcf().savefig(f"spectro_reduc_{d.wavelength:.0f}.pdf")
plt.show()
# -

bkgd = d.spectro.get_stacked_spectrum()[500:]
bkgd_err = d.spectro.get_stacked_spectrum_err()[500:]
print(np.mean(bkgd), np.std(bkgd))
res = (bkgd - np.mean(bkgd)) / bkgd_err
print(np.mean(res), np.std(res))

# ## Wavelength calibration stability

data_path = "/data/STARDICE/cbp/cbp_bench2/golden_sample/"
catalogs = os.listdir(os.path.join(data_path, data_release))
#catalogs = ["/data/STARDICE/cbp/cbp_bench2/golden_sample/dr2/2022_03_04_solar_cell_5mm.npy", "/data/STARDICE/cbp/cbp_bench2/golden_sample/dr2/2022_03_05_stardice_5mm.npy", 
#    "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_05_stardice_5mm/catalog_method_C3_spectro_test_lineunweightedfit.npy", 
#            "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_04_solar_cell_5mm/catalog_method_C3_spectro_test_lineunweightedfit.npy"]
catalogs.sort()
cats = {}
for c in catalogs:
    if '.npy' not in c or "spectro" in c or "2022" not in c:
        continue
    date = '/'.join(c.split('_')[:3])
    print(c)
    cat = np.load(os.path.join(data_path, data_release, c), allow_pickle=True)
    #cat = np.load(c, allow_pickle=True)
    #cats[os.path.basename(c)] = cat
    cats[c] = cat


# +
aa_plot(twocols=False, height=2.2)

fig  = plt.figure()
counter = 0
for c in cats.keys():
    print(c)
    date = '/'.join(os.path.basename(c).split('_')[:3])
    #cat = cats[os.path.basename(c)]
    cat = cats[c]
    #print(cat.dtype.names)
    ind = (cat["spectro_laser_wl_err"] < 1) & (cat["spectro_laser_wl"] > 0)
    label = date
    if "stardice" in c:
        label += " stardice"
    else:
        label += " solar cell"
    plt.errorbar(cat["set_wl"][ind], cat["spectro_laser_wl_cal"][ind]-cat["set_wl"][ind], 
                 yerr=np.sqrt(cat["spectro_laser_wl_cal_stat"][ind]**2+cat["spectro_laser_wl_cal_stat"][ind]**2), 
                 marker="+", linestyle="none", label=label)
    counter += len(cat["spectro_laser_wl"])
print("Number of data points:" , counter)
#plt.grid()
#plt.title(f"Data release: {data_release}")
plt.xlabel("$\lambda_L$ [nm]")
plt.ylabel("$\lambda_c - \lambda_L$ [nm]")
plt.legend(ncol=2)
plt.xlim(350, 1100)
plt.ylim(-2.2, 1.5)
#fig.tight_layout()
fig.savefig("wavelength_stability.pdf")
fig.savefig("wavelength_stability.png")
plt.show()
# -

fig  = plt.figure()
for c in cats.keys():
    print(c)
    date = '/'.join(os.path.basename(c).split('_')[:3])
    cat = cats[c]
    ind = (cat["spectro_laser_wl_err"] < 1) & (cat["spectro_laser_wl"] > 0)
    stats_wl = {}
    stats_wl_err = {}
    for wl in np.unique(cat["set_wl"][ind]):
        key = wl  # f"{wl:.0f}"
        ind_wl = ind & (cat["set_wl"] == wl)
        # if key in stats.keys():
        #     stats[key] += list(cat["spectro_laser_wl"][ind_wl]-cat["set_wl"][ind_wl])
        #     stats_err[key] += list(cat["spectro_laser_wl_err"][ind_wl]+0.01)
        # else:
        #     stats[key] = list(cat["spectro_laser_wl"][ind_wl]-cat["set_wl"][ind_wl])
        #     stats_err[key] = list(cat["spectro_laser_wl_err"][ind_wl]+0.01)
        if key in stats_wl.keys():
            stats_wl[key] += list(cat["spectro_laser_wl_cal"][ind_wl]-cat["set_wl"][ind_wl])
            stats_wl_err[key] += list(cat["spectro_laser_wl_cal_syst"][ind_wl])
        else:
            stats_wl[key] = list(cat["spectro_laser_wl_cal"][ind_wl]-cat["set_wl"][ind_wl])
            stats_wl_err[key] = list(cat["spectro_laser_wl_cal_syst"][ind_wl])
    rms = []
    label = date
    if "stardice" in c:
        label += " stardice"
    else:
        label += " solar cell"
    for key in stats_wl.keys():
        rms.append(np.std((stats_wl[key]-np.mean(stats_wl[key]))/np.array(stats_wl_err[key])))
    #rms = savgol_filter(rms, 51, 1)
    plt.plot(list(stats_wl.keys()), rms, marker="+", linestyle="none", label=label)
counter += len(cat["spectro_laser_wl"])
print("Number of data points:" , counter)
plt.grid()
plt.title(f"Data release: {data_release}")
plt.xlabel("$\lambda_L$ [nm]")
plt.ylabel("RMS of $(\lambda_c - \lambda_L)/\sigma_\lambda$ [units of $\sigma_\lambda$]")
plt.legend(ncol=2)
plt.xlim(350, 1100)
plt.ylim(0.0, 3)
#plt.yscale("log")
#fig.savefig("wavelength_stability_pull.pdf")
plt.show()

stats_wl = {}
stats_wl_err = {}
for c in cats.keys():
    print(c)
    date = '/'.join(os.path.basename(c).split('_')[:3])
    cat = cats[c]
    ind = (cat["spectro_laser_wl_err"] < 1) & (cat["spectro_laser_wl"] > 0)
    for wl in np.unique(cat["set_wl"][ind]):
        key = wl  # f"{wl:.0f}"
        ind_wl = ind & (cat["set_wl"] == wl)
        if key in stats_wl.keys():
            stats_wl[key] += list(cat["spectro_laser_wl"][ind_wl]-cat["set_wl"][ind_wl])
            stats_wl_err[key] += list(cat["spectro_laser_wl_err"][ind_wl])
        else:
            stats_wl[key] = list(cat["spectro_laser_wl"][ind_wl]-cat["set_wl"][ind_wl])
            stats_wl_err[key] = list(cat["spectro_laser_wl_err"][ind_wl])


fig  = plt.figure()
counter = 0
for key in stats_wl.keys():
    #plt.errorbar(key, np.std((np.array(stats[key])-np.mean(stats[key]))/np.array(stats_err[key])), marker="+", linestyle="none", color="r")
    plt.errorbar(key, np.std(np.array(stats_wl[key])), marker="+", linestyle="none", color="r")
    counter += len(cat["spectro_laser_wl"])
print("Number of data points:" , counter)
plt.grid()
plt.title(f"Data release: {data_release}")
plt.xlabel("$\lambda_L$ [nm]")
plt.ylabel("$\lambda_d - \lambda_L$ RMS [nm]")
#plt.legend(ncol=2)
plt.xlim(350, 1100)
#plt.ylim(-10, 10)
plt.yscale("log")
#fig.savefig("wavelength_stability.pdf")
plt.show()

fig  = plt.figure()
for c in cats.keys():
    print(c)
    date = '/'.join(os.path.basename(c).split('_')[:3])
    #cat = cats[os.path.basename(c)]
    cat = cats[c]
    ind = (cat["spectro_laser_wl_err"] < 1) & (cat["spectro_laser_wl"] > 0)
    plt.plot(cat["set_wl"][ind], cat["spectro_laser_wl_err"][ind], "+", label=c)
plt.grid()
plt.legend()
plt.yscale("log")
plt.show()

# ### Vizualize the error budget

cats.keys()

# +
titles = ["2022/03/01 Stardice 75$\,\mathrm{\mu}$m run", "2022/03/04 Solar cell $5\,$mm run", "2022/03/05 Stardice $5\,$mm run"]
datapaths = ['2022_03_01_stardice_transmission_radius.npy',
             '2022_03_04_solar_cell_5mm.npy',
             '2022_03_05_stardice_5mm.npy'
             ]
aa_plot(twocols=False, height=3.5)

fig, ax = plt.subplots(len(datapaths), 1)
for k, path in enumerate(datapaths):
    print(c)
    date = '/'.join(os.path.basename(c).split('_')[:3])
    #cat = cats[os.path.basename(c)]
    run = cats[c]
    # run = SolarCellRun(directory_path=path)
    #wls_calib, wls_calib_stat, wls_calib_syst = from_detected_to_calibrated_wavelengths(run.data["spectro_laser_wl"], 
    #                                                run.data["spectro_laser_wl_err"],
    #                                                spectro_calib_filename=spectro_calib_filename, 
    #                                                spectro_calib_cov_filename=spectro_calib_cov_filename)
    wls_calib = run["spectro_laser_wl_cal"]
    wls_calib_stat = run["spectro_laser_wl_cal_stat"]
    wls_calib_syst = run["spectro_laser_wl_cal_syst"]
    ax[k].plot(wls_calib, wls_calib_stat, marker=".", linestyle="none", color="gray", label=r"$\sigma_\lambda^{\mathrm{stat}}$")
    ax[k].plot(wls_calib, wls_calib_syst, marker=".", linestyle="none", color="k", label=r"$\sigma_\lambda^{\mathrm{cal}}$")
    ax[k].plot(wls_calib, np.sqrt(wls_calib_stat**2+wls_calib_syst**2), marker=".", linestyle="none", color="r", label=r"$\sigma_\lambda$")
    ax[k].set_xlim(349, 1101)
    ax[k].set_yscale("log")
    ax[k].set_ylim(0.01,0.2)
    #ax[k].grid()
    ax[k].set_title(titles[k])
    if k==0:
        ax[k].legend(ncol=3)
    if k==len(datapaths) - 1:
        ax[k].set_xlabel("$\lambda_c$ calibrated wavelength [nm]")
    ax[k].set_ylabel(r"$\sigma_{\lambda}$ [nm]")
fig.tight_layout()
fig.savefig("spectrograph_error_budget.pdf")
fig.savefig("spectrograph_error_budget.png")
plt.show()

# +
datapath='2022_03_04_solar_cell_5mm.npy'
#datapath='/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_05_stardice_5mm/'
run = cats[datapath]
#run.load_from_file(os.path.join(datapath, f"catalog_method_C3_spectro_test_lineunweightedfit_newstaterr.npy"))
#
indices = (run["spectro_laser_wl"] > 560) & (run["spectro_laser_wl"] < 640) & (run["spectro_532_wl_err"]<1)  & (run["spectro_532_wl"] > 530.1) & (run["spectro_532_wl"] <532)

wls = run["spectro_532_wl"][indices]
wls_calib, wls_calib_stat, wls_calib_syst = from_detected_to_calibrated_wavelengths(wls, 
                                                                                   run["spectro_532_wl_err"][indices])
#wls_calib, wls_calib_err = from_detected_to_calibrated_wavelengths(wls, run.data["spectro_532_wl_err"][indices])
wls_calib_mean = np.mean(wls_calib)
wls_calib_std = np.std(wls_calib)
wls_calib_mad = stats.median_abs_deviation(wls_calib)
print(wls_calib_mean, wls_calib_std)


fig = plt.figure()
plt.errorbar(run["expnum"][indices], wls_calib, np.sqrt(wls_calib_stat**2+0*wls_calib_syst**2), marker=".", linestyle="none", color="r")
plt.axhline(wls_calib_mean, color="k", label="mean")
plt.axhline(wls_calib_mean+wls_calib_std, color="k", linestyle="--", label = "RMS")
plt.axhline(wls_calib_mean-wls_calib_std, color="k", linestyle="--")
#plt.axhline(wls_calib_mean+wls_calib_mad)
#plt.axhline(wls_calib_mean-wls_calib_mad)
plt.xlabel("Expnum")
plt.ylabel("$\lambda_c$ calibrated [nm]")
plt.grid()
plt.legend()
fig.tight_layout()
plt.show()
#plt.ylim(250,1200)

# +
res = (wls_calib-wls_calib_mean)/wls_calib_stat
x = stats.norm.rvs(size=100, scale=1, random_state=123456)
mad_gauss = stats.median_abs_deviation(x)
mad_res = stats.median_abs_deviation(res, scale=mad_gauss)

fig = plt.figure()
plt.hist(res, bins=50, label=f"mean={np.mean(res):.4f}$\pm${np.sqrt(np.sum(wls_calib_stat**2))/len(indices):.4f} (stat) $\pm$ {np.mean(wls_calib_syst):.4f} (syst)\nRMS={np.std(res):.3f}\nMAD={mad_res:.3f}", color='r', density=True, alpha=0.7)
x_gauss = np.linspace(-5, 5, 100)
y_gauss = np.exp(-0.5*x_gauss**2) / (np.sqrt(2*np.pi))
plt.plot(x_gauss, y_gauss, label="Normal distribution")
plt.grid()
plt.xlabel("Residuals normalized to stat uncertainties")
plt.legend()
# +
data_path = "/data/STARDICE/cbp/cbp_bench2/golden_sample/"
#data_release = ""
dirs = os.listdir(os.path.join(data_path, data_release))
dirs.sort()

aa_plot(twocols=False, height=3)

fig, ax  = plt.subplots(2, 1) #, figsize=(10,8))
counter = 0
all_wls = []
# for direc in dirs:
#     cat_filenames = []
#     if not os.path.isdir(os.path.join(data_path, direc)):
#         continue
#     for f in os.listdir(os.path.join(data_path, direc)):
#         if ".npy" not in f:
#             continue
#         if 'C3_spectro' in f:
#             cat_filenames.append(os.path.join(data_path, direc, f))
#     # try for subdirectories
#     for direc2 in os.listdir(os.path.join(data_path, direc)):
#         if not os.path.isdir(os.path.join(data_path, direc, direc2)):
#             continue            
#         for f in os.listdir(os.path.join(data_path, direc, direc2)):
#             if ".npy" not in f:
#                 continue
#             if 'C3_spectro' in f:
#                 cat_filenames.append(os.path.join(data_path, direc, direc2, f))
#     for cat_filename in cat_filenames:

for c in cats.keys():
    if "stardice" in c: continue
    print(c)
    date = '/'.join(os.path.basename(c).split('_')[:3])
    #cat = cats[os.path.basename(c)]
    run = cats[c]
    if len(run) < 5000:  # very small run
        continue
    #cat = np.load(os.path.join(data_path, data_release, c))
    #print(cat.dtype.names)
    #run = SolarCellRun(directory_path=os.path.join(data_path, direc))
    #run.load_from_file(cat_filename)
    try:
        indices = (run["spectro_laser_wl"] > 580) & (run["spectro_laser_wl"] < 640) & (run["spectro_532_wl_err"]<1)  & (run["spectro_532_wl"] > 530.1) & (run["spectro_532_wl"] <532)
    except ValueError:
        continue
    # XXXX TO REMOVE IN DEVLOPEMENT
    flux = run["spectro_532_flux_per_pulse"] #* run["set_npulses"]  * run["set_nbursts"]
    flux_err = run["spectro_532_flux_per_pulse_err"] #* run["set_npulses"]  * run["set_nbursts"]
    new_wl = np.copy(run["spectro_532_wl"])
    new_wl_err = np.copy(run["spectro_532_wl_err"])**2
    #for key in pvals_filtered_func.keys():
    #    new_wl -= pvals_filtered_func[key](run["set_wl"]) * flux **(deg-key)
    #    new_wl_err += (flux**(deg-key) * pvals_filtered_err_func[key](run["set_wl"]))**2 +  ((deg-k)*flux_err*flux**(deg-key-1) * pvals_filtered_func[key](run["set_wl"]))**2
    new_wl_err = np.sqrt(new_wl_err)
    # XXXXXXXXXXXXXXXX
    wls_calib, wls_calib_stat, wls_calib_syst = from_detected_to_calibrated_wavelengths(new_wl[indices], 
                                                        new_wl_err[indices], spectro_wl_psf_error=spectro_wl_psf_error,
                                                        spectro_calib_filename=spectro_calib_filename, 
                                                        spectro_calib_cov_filename=spectro_calib_cov_filename)


    all_wls +=  list(wls_calib)

    # XXXX HERE WARNING : check below formula XXXXX
    wls_calib_err = np.sqrt(wls_calib_stat**2+wls_calib_syst**2)
    p = ax[0].errorbar(run["expnum"][indices], wls_calib, wls_calib_err, marker=".", linestyle="none", label=date, alpha=0.3, zorder=-42)

    res = (wls_calib-np.mean(wls_calib))/wls_calib_stat
    print("mean", np.mean(wls_calib),"std", np.std(wls_calib))
    x = stats.norm.rvs(size=100, scale=1, random_state=123456)
    mad_gauss = stats.median_abs_deviation(x)
    mad_res = stats.median_abs_deviation(res, scale=mad_gauss)

    #ax[1].hist(res, bins=50, label=f"mean={np.mean(res):.3f}$\pm${np.sqrt(np.sum(wls_calib_stat**2))/len(indices):.3f} (stat) $\pm$ {np.mean(wls_calib_syst):.3f} (syst), RMS={np.std(res):.3f}", color=p[0].get_color(), density=True, alpha=0.2)
    ax[1].hist(res, bins=50, label=f"RMS={np.std(res):.3f}", color=p[0].get_color(), density=True, alpha=0.2)

    counter += len(indices)
print("Number of data points:" , counter)
wls_calib_mean = np.mean(all_wls)
wls_calib_std = np.std(all_wls)
ax[0].axhline(wls_calib_mean, color="k", label=f"mean={wls_calib_mean:.2f}$\,$nm")
ax[0].axhline(wls_calib_mean+wls_calib_std, color="k", linestyle="--", label = f"RMS={wls_calib_std:.2f}$\,$nm")
ax[0].axhline(wls_calib_mean-wls_calib_std, color="k", linestyle="--")
#ax[0].grid()
ax[0].set_xlabel("Exposure index")
ax[0].set_ylabel("532$\,$nm pump line\nwavelength [nm]")
ax[0].legend(ncol=2)
x_gauss = np.linspace(-5, 5, 100)
y_gauss = np.exp(-0.5*x_gauss**2) / (np.sqrt(2*np.pi))
ax[1].plot(x_gauss, y_gauss, label="Normal distribution")
#ax[1].grid()
ax[1].set_xlabel("Residuals normalized to $\sigma_\lambda^{\mathrm{stat}}$ [$\sigma$ units]")
ax[1].legend(ncol=1)

ax[0].set_ylim(531.9, 532.4)
#plt.title(f"Data release: {data_release}")
fig.tight_layout()
fig.savefig("wavelength_error_model_consistency.pdf")
fig.savefig("wavelength_error_model_consistency.png")
plt.show()
# +
data_path = "/data/STARDICE/cbp/cbp_bench2/golden_sample/"
#data_release = ""
dirs = os.listdir(os.path.join(data_path, data_release))
dirs.sort()

fig, ax  = plt.subplots(1, 2, figsize=(15,5))
counter = 0
all_alphas = []
for c in cats.keys():
        if "stardice" in c: continue
        print(c)
        date = '/'.join(os.path.basename(c).split('_')[:3])
        #cat = cats[os.path.basename(c)]
        run = cats[c]
        indices = (run["spectro_laser_wl"] > 540) & (run["spectro_laser_wl"] < 644) #& (run.data["spectro_532_wl_err"]<1)  & (run.data["spectro_532_wl"] > 530.1) & (run.data["spectro_532_wl"] <532)

        #wls_calib, wls_calib_stat, wls_calib_syst = from_detected_to_calibrated_wavelengths(run.data["spectro_532_wl"][indices], 
        #                                        run.data["spectro_532_wl_err"][indices],
        #                                        spectro_calib_filename="../../../analysis/cbp_paper/stardice_analysis/spectro_calibration_golden_sample_combined_poly_values.npy", 
        #                                        spectro_calib_cov_filename="../../../analysis/cbp_paper/stardice_analysis/spectro_calibration_golden_sample_combined_poly_cov.npy")

        
        alpha = run["spectro_laser_flux_per_pulse"][indices] / run["spectro_laser_max_fit"][indices]
        ind_538 = indices & (run["set_wl"] == 554)
        #alpha /= np.mean(alpha)
        #print(direc, run.data["expnum"][ind_538])
        if np.any(alpha > 0.15):
            toto = np.argmax(alpha)
            print(direc, run["expnum"][indices][toto], run["spectro_532_flux_per_pulse"][indices][toto], run["spectro_laser_flux_per_pulse"][indices][toto])
        alpha_stat = np.abs(alpha) * np.sqrt((run["spectro_532_flux_per_pulse_err"][indices]/run["spectro_532_flux_per_pulse"][indices])**2 + (run["spectro_laser_flux_per_pulse_err"][indices]/run["spectro_laser_flux_per_pulse"][indices])**2 ) 
        all_alphas +=  list(alpha)

        # XXXX HERE WARNING : check below formula XXXXX
        #wls_calib_err = np.sqrt(wls_calib_stat**2+wls_calib_syst**2)
        p = ax[0].errorbar(run["set_wl"][indices], alpha, 
                           alpha_stat, marker=".", linestyle="none", label=date, alpha=0.3, zorder=-42)

        res = (alpha-np.mean(alpha))/alpha_stat
        x = stats.norm.rvs(size=100, scale=1, random_state=123456)
        mad_gauss = stats.median_abs_deviation(x)
        mad_res = stats.median_abs_deviation(res, scale=mad_gauss)

        ax[1].hist(res, bins=50, label=f"mean={np.mean(res):.3f}$\pm${np.sqrt(np.sum(alpha_stat**2))/len(indices):.3f} (stat), RMS={np.std(res):.3f}", color=p[0].get_color(), density=True, alpha=0.2)

        counter += len(indices)
print("Number of data points:" , counter)
wls_calib_mean = np.mean(all_alphas)
wls_calib_std = np.std(all_alphas)
ax[0].axhline(wls_calib_mean, color="k", label=f"mean={wls_calib_mean:.2f}")
ax[0].axhline(wls_calib_mean+wls_calib_std, color="k", linestyle="--", label = f"RMS={wls_calib_std:.2f}")
ax[0].axhline(wls_calib_mean-wls_calib_std, color="k", linestyle="--")
ax[0].grid()
ax[0].set_xlabel("Exposure index")
ax[0].set_ylabel("Measured wavelength from 532nm line\n after wavelength calibration [nm]")
ax[0].legend(ncol=3)
x_gauss = np.linspace(-5, 5, 100)
y_gauss = np.exp(-0.5*x_gauss**2) / (np.sqrt(2*np.pi))
ax[1].plot(x_gauss, y_gauss, label="Normal distribution")
ax[1].grid()
ax[1].set_xlabel("Residuals normalized to stat uncertainties [$\sigma$ units]")
#ax[1].legend(ncol=3)

#ax[0].set_ylim(529, 534)
plt.title(f"Data release: {data_release}")
#fig.tight_layout()
#fig.savefig("error_model_consistency.pdf")
plt.show()
# +
data_dir = "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_04_solar_cell_5mm"  
run5mm= SolarCellRun(directory_path=data_dir)
run5mm.load_from_file(os.path.join(data_dir, f"catalog_method_C3_spectro.npy"))

run_cal = apply_calibration_to_solarcell_catalog(run5mm.data, sc_532=False, sc_dark=False, pd_532=False, 
                                                  sc_1064=False, pd_1064=False)


# -

wls = np.arange(350, 1101, 1)
sigma_clip = 5
cbp_lambdas, cbp_transmission, cbp_transmission_err = compute_cbp_response(run_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="spectro_laser_wl_cal", plot=True)
# cbp_lambdas_qsw298, cbp_transmission_qsw298, cbp_transmission_qsw298_err = compute_cbp_response(run_cal, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=True)


# +
run_cal = np.lib.recfunctions.append_fields(run_cal, data=[run_cal["spectro_laser_wl_cal"]+run_cal["spectro_laser_wl_cal_syst"]], names=["wl_+"], usemask=False)
run_cal = np.lib.recfunctions.append_fields(run_cal, data=[run_cal["spectro_laser_wl_cal"]-run_cal["spectro_laser_wl_cal_syst"]], names=["wl_-"], usemask=False)


cbp_lambdas_plus, cbp_transmission_plus, cbp_transmission_plus_err = compute_cbp_response(run_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="wl_+", plot=False)
cbp_lambdas_minus, cbp_transmission_minus, cbp_transmission_minus_err = compute_cbp_response(run_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="wl_-", plot=False)
# -

run_cal["wl_+"]

run_cal["wl_-"]

fig = plt.figure()
#plt.errorbar(cbp_lambdas, cbp_transmission, yerr=cbp_transmission_err, marker="+", linestyle="none")
plt.errorbar(cbp_lambdas_plus, (cbp_transmission_plus-cbp_transmission)/cbp_transmission, yerr=0*cbp_transmission_plus_err, marker="+", linestyle="none")
plt.errorbar(cbp_lambdas_minus, (cbp_transmission_minus-cbp_transmission)/cbp_transmission, yerr=0*cbp_transmission_minus_err, marker="+", linestyle="none")
plt.xlim(350, 1100)
plt.xlabel(f'$\lambda_c$ [nm]')
plt.ylabel(r"$Q_{\mathrm{solar}}^{\mathrm{dark}}/Q_{\mathrm{phot}}^{\mathrm{dark}}$")
fig.tight_layout()
#plt.savefig("sc_dark_qswmax.pdf")
plt.show()



# ## Solar cell dark current

datapath="/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_02_22_2022_02_22_stardice_transmission_75um/"
#datapath="/data/STARDICE/cbp/cbp_bench2/solarcell/2022_02_24_solarcell_goldensample/"
#catalog_C3 = np.load(os.path.join(datapath, 'catalog_method_C3.npy'))
# %ls {datapath}

# +
filenames = os.listdir(datapath)
filenames = [f for f in filenames if ".fits" in f]
filenames.sort()

#fits_filename = os.path.join(datapath, "IMG_0025000.fits")

datasets = []
for f in filenames[:20]:
    d = CBPDataSet(fits_filename=os.path.join(datapath, f))
    datasets.append(d)
    #print(d)
    # d.plot_data_set()

# +
aa_plot(twocols=False, height=2)

fig = plt.figure()
mean_PSD = []
for d in datasets:
    x = np.copy(d.sc.data["time"])[:550]
    y = np.copy(d.sc.data["CHAR"])[:550]
    pval = np.polyfit(x, y, deg=1)
    #yp = (y[1:] - y[:-1]) / 2e-3
    y -= np.polyval(pval, x)

    #y = np.copy(d.sc.data[mode])
    N = len(y)
    inter = 2e-3
    
    yf = scipy.fft.fft(y)
    yf = scipy.fft.fftshift(yf)
    yPSD = np.abs(yf) ** 2

    xf = scipy.fft.fftfreq(N, inter)
    xf = scipy.fft.fftshift(xf)

    plt.plot(xf, yPSD, 'gray')
    mean_PSD.append(yPSD)
plt.plot(xf, np.mean(mean_PSD, axis=0), 'b-')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.7, 250)
plt.gca().axvspan(1/0.2, 250, alpha=0.5, color='lightgray')
plt.ylim(1e-21, 1e-14)
plt.xlabel(r"$f$ [Hz]")
plt.ylabel("Solar cell dark power spectrum [C$^2$/Hz]")
fig.tight_layout()
plt.savefig("sc_dark_ps.pdf")
plt.show()
# -

# ## Number of pulses

# +
aa_plot(twocols=False, height=1.5)

data_dir = "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_06_solar_cell_5mm"  
run_sc= SolarCellRun(directory_path=data_dir)
run_sc.load_from_file(os.path.join(data_dir, f"catalog_method_C3_spectro.npy"))
data_dir = "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_05_stardice_5mm"  
run_sd= SolarCellRun(directory_path=data_dir)
run_sd.load_from_file(os.path.join(data_dir, f"catalog_method_C3_spectro.npy"))

fig = plt.figure()
plt.plot(run_sd.data["set_wl"], run_sd.data["set_npulses"], 'r+', label="StarDICE run")
plt.plot(run_sc.data["set_wl"], run_sc.data["set_npulses"], 'b+', label="Solar cell run")
plt.xlabel(r"$\lambda_L$ [nm]")
plt.ylabel(r"Number of laser pulses")
plt.ylim(0.8, 220)
plt.yscale('log')
plt.legend()
plt.savefig("npulses.pdf")
plt.show()
# -


# ## Sensor QEs

# +
pd_qe_data_file = "/data/STARDICE/cbp/solarcell/refCalData/SM05PD1B_QE.csv"
spectro_qe_data_file = "/data/STARDICE/cbp/solarcell/refCalData/ocean65000_qe.npy"
pd_qe = np.loadtxt(pd_qe_data_file, skiprows=1, delimiter=",").T
pd_qe[1] *= (const.h * const.c / (pd_qe[0] * 1e-9 * u.meter * const.e.value)).value
pd_qe_f = interp1d(pd_qe[0], pd_qe[1], bounds_error=False, fill_value=np.min(pd_qe[1]))
spectro_qe = np.load(spectro_qe_data_file)
sp_qe_f = interp1d(spectro_qe['wavelengths'], spectro_qe['spectro_qe'], bounds_error=False)

wl_tot = np.arange(350, 1101, 1)

aa_plot(twocols=False, height=2)

fig = plt.figure() #figsize=(15,7))
plt.plot(wl_tot, pd_qe_f(wl_tot), '+', label='Spectrograph', color=instrument_colors['spectro'])
plt.plot(wl_tot, sp_qe_f(wl_tot), '+', label='Photodiode', color=instrument_colors['photodiode'])
plt.ylabel('Quantum efficiency') #, fontsize=20)
plt.xlabel('Wavelength [nm]') #, fontsize=20)
#plt.xticks(fontsize=18)
#plt.yticks(fontsize=18)
plt.legend() #fontsize=18)
plt.tight_layout()
fig.savefig("qe_phototiode_spectro.pdf")
plt.show()

# -

# ## CBP response

data_dir = "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_04_solar_cell_5mm"  
run5mm= SolarCellRun(directory_path=data_dir)
run5mm.load_from_file(os.path.join(data_dir, f"catalog_method_C3_spectro.npy"))

catalog = apply_calibration_to_solarcell_catalog(run5mm.data, sc_532=False, sc_dark=False, pd_532=False, 
                                                  sc_1064=False, pd_1064=False)

# +
wls = np.arange(350, 1101, 1)
sigma_clip = 5

good_indices = preliminary_cuts(catalog)
good_indices = good_indices & (catalog["qsw"] == "MAX") & (catalog["SC_distance"] == "short")
wl_bins, tr_cbp, tr_cbp_err = compute_binned_values(catalog, wls, key_x="set_wl", key_y="ratio_cbp_cal",
                                                    key_y_stat="ratio_cbp_cal_stat", key_y_syst="ratio_cbp_cal_syst",
                                                    indices=good_indices, sigma_clip=sigma_clip)

# +
from saunerie import bspline

aa_plot(twocols=False, height=3)

fig, ax = plt.subplots(3, 1, sharex='all')
ax[0].errorbar(catalog["set_wl"][good_indices], catalog["ratio_cbp_cal"][good_indices], 
             yerr=np.sqrt(catalog["ratio_cbp_cal_stat"][good_indices]**2+catalog["ratio_cbp_cal_syst"][good_indices]**2), color="k", alpha=0.1, marker='+',
             linestyle="none", label="per burst ratios", zorder=-42)
ax[0].plot(wl_bins, tr_cbp, color="b", alpha=1, marker='+', linestyle="none", label=rf"${sigma_clip}\sigma$ clip average")
ax[0].fill_between(wl_bins, tr_cbp - tr_cbp_err, tr_cbp + tr_cbp_err, color="b", alpha=0.3)
ax[0].set_xlabel(f'$\lambda_L$ [nm]')
ax[0].set_ylabel(r"$Q_{\mathrm{solar}}/Q_{\mathrm{phot}}$")
ax[0].set_ylim(30,170)
ax[0].set_xlim(350, 1100)
#plt.grid()
#plt.title("CBP response")
ax[0].legend()

wl = run5mm.data["set_wl"][good_indices]
opentr = run5mm.data["sc_charge_per_burst"][good_indices] / run5mm.data["pd_charge_per_burst"][good_indices] 
opentr_err = opentr * np.sqrt( (run5mm.data["sc_charge_per_burst_err"][good_indices]/run5mm.data["sc_charge_per_burst"][good_indices])**2 + (run5mm.data["pd_charge_per_burst_err"][good_indices]/run5mm.data["pd_charge_per_burst"][good_indices])**2)
B = bspline.CardinalBSpline(n=300, x_range=(wl.min(), wl.max()+0.1))
A = B.eval(wl)
P = B.linear_fit(wl, opentr)
residuals = opentr - (A*P)

ax[1].errorbar(wl, 100*opentr_err/opentr, marker='.', linestyle='', color='k', alpha=0.3)
ax[1].errorbar(wl_bins, 100*tr_cbp_err/tr_cbp, marker='.', linestyle='', color='b', alpha=1)
ax[1].set_xlabel(f'$\lambda_L$ [nm]')
ax[1].set_ylabel("Relative statistical\nuncertainties [%]")
ax[1].set_yscale("log")
ax[1].set_ylim(0.7e-2,5)
#ax[1].grid()

ax[2].errorbar(wl, residuals/opentr_err, yerr=opentr_err/opentr_err, marker='.', linestyle='', color='k', alpha=0.3)
#ax[2].axhline(1, color='r', linestyle='--')
#ax[2].axhline(-1, color='r', linestyle='--')
ax[2].set_xlabel(f'$\lambda_L$ [nm]')
ax[2].set_ylabel("$\sigma$-normalized\nresiduals")
ax[2].set_ylim(-15, 15)
#ax[2].grid()

fig.tight_layout()
fig.subplots_adjust(hspace=0)
plt.savefig("cbp_charge_ratio.pdf")
plt.savefig("cbp_charge_ratio.png")
plt.show()

# -

fig = plt.figure()
ind_blue = catalog[good_indices]["set_wl"] < 669
ind_red = catalog[good_indices]["set_wl"] > 669
plt.hist(residuals[ind_blue]/opentr_err[ind_blue], bins=50, color="b", alpha=0.5)
plt.hist(residuals[ind_red]/opentr_err[ind_red], bins=50, color="r", alpha=0.5)
plt.show()

# ## 532 and 1064nm contamination

cat = np.load("/data/STARDICE/cbp/cbp_bench2/golden_sample/dr3/2022_03_06_solar_cell_5mm.npy")

# +
fig = plt.figure()

indices = cat["qsw"] == "MAX"
plt.errorbar(cat["set_wl"][indices], cat["alpha_532"][indices], yerr=cat["alpha_532_err"][indices], 
             color="b", linestyle="none", marker='+', label="$\alpha(\lambda)$ at QSW=MAX")
plt.errorbar(cat["set_wl"][indices], cat["beta_1064"][indices], yerr=cat["beta_1064_err"][indices], 
             color="r", linestyle="none", marker='+', label="$\beta(\lambda)$ at QSW=MAX")
indices = cat["qsw"] == "298"
plt.errorbar(cat["set_wl"][indices], cat["alpha_532"][indices], yerr=cat["alpha_532_err"][indices], 
             color="c", linestyle="none", marker='+', label="$\alpha(\lambda)$ at QSW=298")
plt.errorbar(cat["set_wl"][indices], cat["beta_1064"][indices], yerr=cat["beta_1064_err"][indices], 
             color="lightcoral", linestyle="none", marker='+', label="$\beta(\lambda)$ at QSW=298")
plt.xlim(350, 1100)
plt.legend()
plt.show()
# -
# ## 532nm correction

# +
full_cat_path = f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/2022_03_03_stardice_repeatability_tot.npy"
pinhole='75um'
wl_range = np.arange(350,1101,1)
pinhole='75um'
sigma_clip = 5
filt = 'g'


catalog=np.load(full_cat_path, allow_pickle=True)
q_pd = catalog["pd_charge_total"]
q_sd = catalog["sd_charge_total"]
# compute ratio and uncertainties
ratio = q_sd / q_pd
ratio_stat = ratio * np.sqrt((catalog["sd_charge_total_err"]/q_sd)**2 + (catalog["pd_charge_total_err"]/q_pd)**2)
# append new columns to catalog
catalog = np.lib.recfunctions.append_fields(catalog, data=[ratio, ratio_stat],
                                            names=['ratio_sd', 'ratio_sd_err'], usemask=False)


ind_filter = catalog['filter'] == filt.encode('utf-8')
sd_wl = catalog[ind_filter]['set_wl']
good_indices = preliminary_cuts(catalog)
good_indices = good_indices & (catalog["cbppinhole"] == pinhole) 
sd_wl_empty, SD_tr_empty, SD_tr_empty_err = compute_binned_values(catalog, wl_range, key_x="set_wl", key_y="ratio_sd_cal",
                                                    key_y_stat="ratio_sd_cal_stat", key_y_syst="ratio_sd_cal_syst",
                                                    indices=good_indices & (catalog["filter"] == b'EMPTY'), sigma_clip=sigma_clip)
sd_wl, SD_tr, SD_tr_err = compute_binned_values(catalog, wl_range, key_x="set_wl", key_y="ratio_sd_cal",
                                                    key_y_stat="ratio_sd_cal_stat", key_y_syst="ratio_sd_cal_syst",
                                                    indices=good_indices & (catalog["filter"] == filt.encode("utf-8")), sigma_clip=sigma_clip)
sd_wl_532, SD_tr_532, SD_tr_532_err = compute_binned_values(catalog, wl_range, key_x="set_wl", key_y="ratio_sd",
                                                    key_y_stat="ratio_sd_err", key_y_syst=None,
                                                    indices=good_indices & (catalog["filter"] == filt.encode("utf-8")), sigma_clip=sigma_clip)

SD_tr_norm = 100 * SD_tr/np.interp(sd_wl,sd_wl_empty,SD_tr_empty)
SD_tr_532_norm = 100 * SD_tr_532/np.interp(sd_wl,sd_wl_empty,SD_tr_empty)

SD_tr_norm_err = np.abs(SD_tr_norm) * np.sqrt((SD_tr_err/SD_tr)**2 + (np.interp(sd_wl,sd_wl_empty,SD_tr_empty_err)/np.interp(sd_wl,sd_wl_empty,SD_tr_empty))**2)
SD_tr_532_norm_err = np.abs(SD_tr_532_norm) * np.sqrt((SD_tr_532_err/SD_tr_532)**2 + (np.interp(sd_wl,sd_wl_empty,SD_tr_empty_err)/np.interp(sd_wl,sd_wl_empty,SD_tr_empty))**2)


# +
aa_plot(twocols=False, height=2)

fig, ax = plt.subplots(1,1) #,figsize=(7,4))
ax.errorbar(sd_wl_532, SD_tr_532_norm, SD_tr_532_norm_err, marker = '', linestyle='--', label=f"no $532\,$nm correction", color=filter_colors[filt.encode("utf-8")])
ax.errorbar(sd_wl, SD_tr_norm, SD_tr_norm_err, marker = '', linestyle='-', label=f"with $532\,$nm correction", color=filter_colors[filt.encode("utf-8")])
#ax.grid()
ax.legend(ncol=1)
ax.set_xlabel('$\lambda_L$ [nm]') #, fontsize = 13)
ax.set_ylabel(f'{filt} filter transmission [%]') #, fontsize = 13)

# inset axes....
axins = ax.inset_axes([0.45, 0.3, 0.47, 0.2])
axins.errorbar(sd_wl_532, SD_tr_532_norm, SD_tr_532_norm_err, marker = '', linestyle='--', label=f"{filt} with no $532\,$nm correction", color=filter_colors[filt.encode("utf-8")])
axins.errorbar(sd_wl, SD_tr_norm, SD_tr_norm_err, marker = '', linestyle='-', label=f"{filt} with $532\,$nm correction", color=filter_colors[filt.encode("utf-8")])
# subregion of the original image
x1, x2, y1, y2 = 535, 655, -3e-1, 3
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
#axins.set_xticklabels([])
#axins.set_yticklabels([])
ax.indicate_inset_zoom(axins, edgecolor="black")

plt.tight_layout()
plt.savefig("g_filter_532.pdf")
plt.show()
# -

# ## Solar cell dark

data_dir = "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_07_solar_cell_5mm_with_cap"  
run5mm= SolarCellRun(directory_path=data_dir)
run5mm.load_from_file(os.path.join(data_dir, f"catalog_method_C3_spectro.npy"))
run_cal = apply_calibration_to_solarcell_catalog(run5mm.data, sc_532=False, sc_dark=False, pd_532=False, 
                                                 pd_1064=False, sc_1064=False)

# +
aa_plot(twocols=False, height=2)

wls = np.arange(350, 1101, 1)
sigma_clip = 5
cbp_lambdas_qswmax, cbp_transmission_qswmax, cbp_transmission_qswmax_err = compute_cbp_response(run_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=True)
# cbp_lambdas_qsw298, cbp_transmission_qsw298, cbp_transmission_qsw298_err = compute_cbp_response(run_cal, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=True)

fig = plt.gcf()
ax = plt.gca()
ax.grid(False)
plt.xlim(350, 1100)
ax.set_xlabel(f'$\lambda_L$ [nm]')
ax.set_ylabel(r"$Q_{\mathrm{solar}}^{\mathrm{dark}}/Q_{\mathrm{phot}}^{\mathrm{dark}}$")
fig.tight_layout()
plt.savefig("sc_dark_qswmax.pdf")
plt.savefig("sc_dark_qswmax.png")
plt.show()
# -

# ## Solar cell QSW ratios

data_dir = "/data/STARDICE/cbp/cbp_bench2/golden_sample/2022_03_04_solar_cell_5mm"  
run5mm= SolarCellRun(directory_path=data_dir)
run5mm.load_from_file(os.path.join(data_dir, f"catalog_method_C3_spectro.npy"))
run_raw = apply_calibration_to_solarcell_catalog(run5mm.data, sc_532=False, sc_dark=False, pd_532=False, 
                                                 pd_1064=False, sc_1064=False)
#run_cal = apply_calibration_to_solarcell_catalog(run5mm.data, sc_532=True, sc_dark=True, pd_532=True, 
#                                                 pd_1064=False, sc_1064=False)
run_cal = np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/2022_03_04_solar_cell_5mm.npy")

wls = np.arange(350, 1101, 1)
# wls = np.arange(1000, 1101, 1)
sigma_clip = 5
cbp_lambdas_qsw298_raw, cbp_transmission_qsw298_raw, cbp_transmission_qsw298_err_raw = compute_cbp_response(run_raw, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qswmax_raw, cbp_transmission_qswmax_raw, cbp_transmission_qswmax_err_raw = compute_cbp_response(run_raw, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qsw298, cbp_transmission_qsw298, cbp_transmission_qsw298_err = compute_cbp_response(run_cal, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qswmax, cbp_transmission_qswmax, cbp_transmission_qswmax_err = compute_cbp_response(run_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)

# +
aa_plot(twocols=False, height=3)

fig, ax = plt.subplots(2, 1, sharex="all", height_ratios=[3, 1])

alpha=0.5
lw=3


cbp_298 = np.interp(cbp_lambdas_qswmax_raw, cbp_lambdas_qsw298_raw, cbp_transmission_qsw298_raw)
cbp_298_err = np.interp(cbp_lambdas_qswmax_raw, cbp_lambdas_qsw298_raw, cbp_transmission_qsw298_err_raw)
ratio = cbp_298 / cbp_transmission_qswmax_raw
ratio_err = np.abs(ratio) * np.sqrt((cbp_298_err/cbp_298)**2 + (cbp_transmission_qswmax_err_raw/cbp_transmission_qswmax_raw)**2)
ax[0].errorbar(cbp_lambdas_qswmax_raw, ratio, yerr=ratio_err, color="b", marker="+", linestyle="none", label=r"before any correction", alpha=1)

ind = cbp_lambdas_qswmax_raw < 400
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw < 532) & (cbp_lambdas_qswmax_raw >= 400)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw >= 532) & (cbp_lambdas_qswmax_raw < 669)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw >= 669) & (cbp_lambdas_qswmax_raw < 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw >= 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)

label = r"after $r_{\mathrm{CBP}}^{\mathrm{dark}}$ subtraction and"+"\nlaser contamination correction"
cbp_298 = np.interp(cbp_lambdas_qswmax, cbp_lambdas_qsw298, cbp_transmission_qsw298)
cbp_298_err = np.interp(cbp_lambdas_qswmax, cbp_lambdas_qsw298, cbp_transmission_qsw298_err)
ratio = cbp_298 / cbp_transmission_qswmax
ratio_err = np.abs(ratio) * np.sqrt((cbp_298_err/cbp_298)**2 + (cbp_transmission_qswmax_err/cbp_transmission_qswmax)**2)
ax[0].errorbar(cbp_lambdas_qswmax, ratio, yerr=ratio_err, color="red", marker="+", linestyle="none", label=label, alpha=1)

ind = cbp_lambdas_qswmax_raw < 400
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw < 532) & (cbp_lambdas_qswmax_raw >= 400)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax >= 532) & (cbp_lambdas_qswmax < 669)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax >= 669) & (cbp_lambdas_qswmax < 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax >= 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)


ax[0].axhline(1+1e-3, color='k', linestyle='--')
ax[0].axhline(1, color='k', linestyle='-')
ax[0].axhline(1-1e-3, color='k', linestyle='--')
#ax[0].set_ylabel("CBP charge ratios\n (QSW=298)/(QSW=MAX)")
ax[0].set_ylabel(r"$r_{\mathrm{CBP}}^{\mathrm{QSW}=298} / r_{\mathrm{CBP}}^{\mathrm{QSW}=MAX}$")

ax[0].set_xlabel("$\lambda_L$ [nm]")
ax[0].set_title("Solar cell run 2022/03/04")
ax[0].legend()
ax[0].set_xlim(350, 1100)
ax[0].set_ylim(0.98, 1.02)

ax[1].axhline(1+1e-3, color='k', linestyle='--')
ax[1].axhline(1, color='k', linestyle='-')
ax[1].axhline(1-1e-3, color='k', linestyle='--')
ax[1].set_ylabel("Binned")
ax[1].set_xlabel("$\lambda_L$ [nm]")
#ax[1].legend()
ax[1].set_xlim(350, 1100)
ax[1].set_ylim(0.993, 1.007)



fig.tight_layout()
fig.subplots_adjust(hspace=0)
plt.savefig("sc_qsw_ratios.pdf")
plt.show()
# -

# ## Repeatability

run1_cal = np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/2022_03_04_solar_cell_5mm.npy")
run2_cal = np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/2022_03_06_solar_cell_5mm.npy")
run3_cal = np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/2022_03_08_solar_cell_distance.npy")

wls = np.arange(350, 1101, 1)
sigma_clip = 5
cbp_lambdas_qsw298, cbp_transmission_qsw298, cbp_transmission_qsw298_err = compute_cbp_response(run1_cal, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qswmax, cbp_transmission_qswmax, cbp_transmission_qswmax_err = compute_cbp_response(run1_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qsw298_2, cbp_transmission_qsw298_2, cbp_transmission_qsw298_err_2 = compute_cbp_response(run2_cal, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qswmax_2, cbp_transmission_qswmax_2, cbp_transmission_qswmax_err_2 = compute_cbp_response(run2_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qsw298_3, cbp_transmission_qsw298_3, cbp_transmission_qsw298_err_3 = compute_cbp_response(run3_cal, wls, qsw="298", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)
cbp_lambdas_qswmax_3, cbp_transmission_qswmax_3, cbp_transmission_qswmax_err_3 = compute_cbp_response(run3_cal, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False)

cbp_response_mean = np.zeros_like(cbp_lambdas_qswmax).astype(float)
cbp_response_mean_err = np.zeros_like(cbp_lambdas_qswmax).astype(float)
for k, lbda in enumerate(cbp_lambdas_qswmax):
    n = 1
    response = cbp_transmission_qswmax[k]
    response_err = cbp_transmission_qswmax_err[k]**2
    if lbda in cbp_lambdas_qswmax_2:
        response += cbp_transmission_qswmax_2[cbp_lambdas_qswmax_2==lbda][0]
        response_err += cbp_transmission_qswmax_err_2[cbp_lambdas_qswmax_2==lbda][0]**2
        n += 1
    if lbda in cbp_lambdas_qswmax_3:
        response += cbp_transmission_qswmax_3[cbp_lambdas_qswmax_3==lbda][0]
        response_err += cbp_transmission_qswmax_err_3[cbp_lambdas_qswmax_3==lbda][0]**2
        n += 1
    # print(k, lbda, np.sqrt(response_err)/response, n)
    cbp_response_mean[k] = response / n
    cbp_response_mean_err[k] = np.sqrt(response_err) / n

# +
aa_plot(twocols=False, height=2.5)
fig, ax = plt.subplots(2, 1, sharex="all", height_ratios=[3, 1])

alpha=0.5
lw=3

label = r"run 1"
ratio = cbp_transmission_qswmax/cbp_response_mean
ratio_err = cbp_transmission_qswmax/cbp_response_mean  * np.sqrt((cbp_transmission_qswmax_err/cbp_transmission_qswmax)**2 + (cbp_response_mean_err/cbp_response_mean)**2)
ax[0].errorbar(cbp_lambdas_qswmax_raw, ratio, yerr=ratio_err, color="b", marker="+", linestyle="none", label=label, alpha=1)

ind = cbp_lambdas_qswmax_raw < 400
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw < 532) & (cbp_lambdas_qswmax_raw >= 400)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw >= 532) & (cbp_lambdas_qswmax_raw < 669)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw >= 669) & (cbp_lambdas_qswmax_raw < 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw >= 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "b-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="b", alpha=alpha)

label = r"run 2"
ratio = cbp_transmission_qswmax_2/cbp_response_mean
ratio_err = cbp_transmission_qswmax_2/cbp_response_mean  * np.sqrt((cbp_transmission_qswmax_err_2/cbp_transmission_qswmax_2)**2 + (cbp_response_mean_err/cbp_response_mean)**2)
ax[0].errorbar(cbp_lambdas_qswmax, ratio, yerr=ratio_err, color="red", marker="+", linestyle="none", label=label, alpha=1)

ind = cbp_lambdas_qswmax_raw < 400
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax_raw < 532) & (cbp_lambdas_qswmax_raw >= 400)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax >= 532) & (cbp_lambdas_qswmax < 669)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax >= 669) & (cbp_lambdas_qswmax < 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)
ind = (cbp_lambdas_qswmax >= 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax[ind], ratio[ind], "r-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="r", alpha=alpha)

label = r"run 3"
ratio = cbp_transmission_qswmax_3/np.interp(cbp_lambdas_qswmax_3, cbp_lambdas_qswmax, cbp_response_mean)
ratio_err = np.abs(ratio) * np.sqrt((cbp_transmission_qswmax_err_3/cbp_transmission_qswmax_3)**2 + (np.interp(cbp_lambdas_qswmax_3, cbp_lambdas_qswmax, cbp_response_mean_err)/np.interp(cbp_lambdas_qswmax_3, cbp_lambdas_qswmax, cbp_response_mean))**2)
ax[0].errorbar(cbp_lambdas_qswmax_3, ratio, yerr=ratio_err, color="green", marker="+", linestyle="none", label=label, alpha=1)

ind = cbp_lambdas_qswmax_3 < 400
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax_3[ind], ratio[ind], "g-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax_3[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="g", alpha=alpha)
ind = (cbp_lambdas_qswmax_3 < 532) & (cbp_lambdas_qswmax_3 >= 400)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax_3[ind], ratio[ind], "g-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax_3[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="g", alpha=alpha)
ind = (cbp_lambdas_qswmax_3 >= 532) & (cbp_lambdas_qswmax_3 < 669)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax_3[ind], ratio[ind], "g-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax_3[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="g", alpha=alpha)
ind = (cbp_lambdas_qswmax_3 >= 669) & (cbp_lambdas_qswmax_3 < 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax_3[ind], ratio[ind], "g-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax_3[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="g", alpha=alpha)
ind = (cbp_lambdas_qswmax_3 >= 1064)
ratio_err[ind] = np.std(ratio[ind]) / np.sum(ind)
ratio[ind] = np.mean(ratio[ind])
ax[1].plot(cbp_lambdas_qswmax_3[ind], ratio[ind], "g-", lw=lw, alpha=1)
ax[1].fill_between(cbp_lambdas_qswmax_3[ind], ratio[ind]+ratio_err[ind], ratio[ind]-ratio_err[ind],  color="g", alpha=alpha)


ax[0].axhline(1+1e-3, color='k', linestyle='--')
ax[0].axhline(1, color='k', linestyle='-')
ax[0].axhline(1-1e-3, color='k', linestyle='--')
ax[0].set_ylabel("$r_{\mathrm{CBP}}^{\mathrm{run}\ i}\;/\;\overline{r_{\mathrm{CBP}}}$")
ax[0].set_xlabel("$\lambda_L$ [nm]")
#ax[0].set_title("Solar cell run 2022/03/04")
ax[0].legend()
ax[0].set_xlim(350, 1100)
ax[0].set_ylim(0.988, 1.012)

ax[1].axhline(1+1e-3, color='k', linestyle='--')
ax[1].axhline(1, color='k', linestyle='-')
ax[1].axhline(1-1e-3, color='k', linestyle='--')
ax[1].set_ylabel("Binned")
ax[1].set_xlabel("$\lambda_L$ [nm]")
#ax[1].legend()
ax[1].set_xlim(350, 1100)
ax[1].set_ylim(0.9975, 1.0025)



fig.tight_layout()
fig.subplots_adjust(hspace=0)
plt.savefig("sc_runi_ratios.pdf")
plt.savefig("sc_runi_ratios.png")
plt.show()
# -

# ## Distance

run_dist = np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/2022_03_08_solar_cell_distance.npy")

wls = np.arange(350, 1102, 1)
indices = run_dist["set_wl"]!=620
sigma_clip = 5
cbp_lambdas_short, cbp_transmission_short, cbp_transmission_short_err = compute_cbp_response(run_dist, wls, qsw="MAX", SC_distance="short", sigma_clip=sigma_clip, wl_key="set_wl", plot=False, indices=indices)
cbp_lambdas_long, cbp_transmission_long, cbp_transmission_long_err = compute_cbp_response(run_dist, wls, qsw="MAX", SC_distance="long", sigma_clip=sigma_clip, wl_key="set_wl", plot=False, indices=indices)

lbdas = np.intersect1d(cbp_lambdas_long, cbp_lambdas_short)
lbdas_long = [lbda for lbda in cbp_lambdas_long if lbda in lbdas]
lbdas_short = [lbda for lbda in cbp_lambdas_short if lbda in lbdas]

# +
cbp_tr_short = np.interp(cbp_lambdas_long, cbp_lambdas_short, cbp_transmission_short)
cbp_tr_short_err = np.interp(cbp_lambdas_long, cbp_lambdas_short, cbp_transmission_short_err)

cbp_tr_short = np.interp(lbdas_short, cbp_lambdas_short, cbp_transmission_short)
cbp_tr_short_err = np.interp(lbdas_short, cbp_lambdas_short, cbp_transmission_short_err)
cbp_tr_long = np.interp(lbdas_long, cbp_lambdas_long, cbp_transmission_long)
cbp_tr_long_err = np.interp(lbdas_long, cbp_lambdas_long, cbp_transmission_long_err)

ratio = cbp_tr_long/cbp_tr_short
ratio_err = np.abs(ratio) * np.sqrt((cbp_tr_long_err/cbp_tr_long)**2 + (cbp_tr_short_err/cbp_tr_short)**2)

pval = np.polyfit(lbdas, ratio, deg=1, w=1/ratio_err)

aa_plot(twocols=False, height=2)

fig = plt.figure()
#plt.errorbar(cbp_lambdas_long, ratio_uncorr, yerr=ratio_uncorr_err, marker="+", linestyle="", label="long/short")
p = plt.errorbar(lbdas, ratio, yerr=ratio_err, marker="+", linestyle="", color="k", alpha=1, label="long/short binned data")
plt.plot(cbp_lambdas_long, np.polyval(pval, cbp_lambdas_long), color="b", lw=3, 
                 label=f'fit: ${pval[0]:.8f}\\times\lambda+{pval[1]:.4f}$', alpha=1)
#plt.grid()
plt.legend()
plt.xlabel("$\lambda_L$ [nm]")
plt.ylabel(r"$r_{\mathrm{CBP}}^{\mathrm{long\ distance}} / r_{\mathrm{CBP}}^{\mathrm{short\ distance}}$")
fig.tight_layout()
fig.savefig("sc_distance.pdf")
plt.show()
# -



# ## CBP response

full_cat_cbp = np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/cbp_response_full_cat.npy", allow_pickle=True)

ind = (full_cat_cbp["qsw"]=="MAX") & (full_cat_cbp["SC_distance"]=="short")
wls = np.arange(350, 1102, 1)
good_indices = preliminary_cuts(full_cat_cbp) & ind
wl_key="set_wl"
wl_bins_stat, tr_cbp, tr_cbp_stat = compute_binned_values(full_cat_cbp, wls, key_x=wl_key, key_y="ratio_cbp_cal",
                                                    key_y_stat="ratio_cbp_cal_stat", key_y_syst=None,
                                                    indices=good_indices, sigma_clip=sigma_clip)
wl_bins_syst, tr_cbp, tr_cbp_syst = compute_binned_values(full_cat_cbp, wls, key_x=wl_key, key_y="ratio_cbp_cal",
                                                    key_y_stat=None, key_y_syst="ratio_cbp_cal_syst",
                                                    indices=good_indices, sigma_clip=sigma_clip)
wl_bins, tr_cbp, tr_cbp_syst_alpha_beta = compute_binned_values(full_cat_cbp, wls, key_x=wl_key, key_y="ratio_cbp_cal",
                                                    key_y_stat=None, key_y_syst="ratio_cbp_cal_syst_alpha_beta",
                                                    indices=good_indices, sigma_clip=sigma_clip)

# +
aa_plot(twocols=False, height=2.5)

fig = plt.figure()
plt.plot(wl_bins_stat, tr_cbp_stat/tr_cbp, '.', label="statistics")
plt.plot(full_cat_cbp["set_wl"][ind], full_cat_cbp["cbp_response_cal_syst_SCeff"][ind]/full_cat_cbp["cbp_response_cal"][ind], '.', label="$\eta_{SC}$")
plt.plot(full_cat_cbp["set_wl"][ind], full_cat_cbp["ratio_cbp_cal_syst_distance"][ind]/full_cat_cbp["ratio_cbp_cal"][ind], '.', label="scattered light")
plt.plot(full_cat_cbp["set_wl"][ind], full_cat_cbp["ratio_cbp_cal_syst_linearity"][ind]/full_cat_cbp["ratio_cbp_cal"][ind], '.', label="linearity")
plt.plot(full_cat_cbp["set_wl"][ind], full_cat_cbp["ratio_cbp_cal_syst_repeatability"][ind]/full_cat_cbp["ratio_cbp_cal"][ind], '.', label="repeatability")
plt.plot(wl_bins, tr_cbp_syst_alpha_beta/tr_cbp, '.', label=r"$532\,$nm and $1064\,$nm"+"\ncontamination")
plt.plot(full_cat_cbp["set_wl"][ind], full_cat_cbp["ratio_cbp_cal_syst"][ind]/full_cat_cbp["ratio_cbp_cal"][ind], 'k.', label="total systematics")
plt.yscale("log")
plt.xlim(350, 1100)
plt.xlabel("$\lambda_L$ [nm]")
plt.ylabel("Relative uncertainties on the CBP response")
plt.legend(ncol=2, loc="lower right")
fig.tight_layout()
plt.savefig("cbp_error_budget.pdf")
plt.savefig("cbp_error_budget.png")
plt.show()

# +
#wl_range = np.arange(350,1102,1)
#sigma_clip = 5
#cbp_wl, cbp_tr, cbp_tr_err = compute_cbp_response(full_cat_cbp, wl_range, qsw='MAX', SC_distance="short", 
#                                                  sigma_clip=sigma_clip, wl_key="spectro_laser_wl_cal")

#tr = cbp_tr/SC_QE_f(cbp_wl)/const.e.value
#tr_err = tr * np.sqrt(0*(SC_QE_f_err(cbp_wl)/SC_QE_f(cbp_wl))**2 + (cbp_tr_err/cbp_tr)**2)
#cbp_tr = tr.copy()
#cbp_tr_err = tr_err.copy()

#cbp_response = np.rec.fromarrays([cbp_wl, cbp_tr, cbp_tr_err], dtype=[('cbp_wl', '<f8'), ('cbp_tr', '<f8'), ('cbp_tr_err', '<f8')])
cbp_response = np.array(np.load(f"/data/STARDICE/cbp/cbp_bench2/golden_sample/{data_release}/cbp_response.npy"))
cbp_wl = cbp_response["cbp_wl"]
cbp_tr = cbp_response["cbp_tr"]
cbp_tr_err = cbp_response["cbp_tr_err"]

# +
aa_plot(twocols=False, height=3)

fig, ax = plt.subplots(2,1, sharex=True)
ax[0].errorbar(cbp_wl, cbp_tr, cbp_tr_err, marker = '+', color='crimson')
ax[0].set_ylabel('CBP response [$\gamma$/C]')
#ax[0].set_xlabel('Wavelengths [nm]')
props = dict(boxstyle='square', facecolor='white', alpha=0.8)
#ax[0].text(0.05, 0.95, r"$R_{\mathrm{CBP}} = \frac{Q_{\mathrm{solar}}}{Q_{\mathrm{phot}} \times \epsilon_{\mathrm{solar}} \times e}$", 
           #transform=ax[0].transAxes, fontsize=18, verticalalignment='top', bbox=props, color='crimson')
#ax[0].grid()
ax[0].set_ylim(0, 1.5e21)

# ax[1].plot(full_cat_cbp["spectro_laser_wl_cal"][ind], 100*full_cat_cbp["cbp_response_cal_syst"][ind]/full_cat_cbp["cbp_response_cal"][ind], marker='+', linestyle='none', color='g', label="total error")
ax[1].plot(cbp_wl, 100*cbp_tr_err/cbp_tr, marker='none', linestyle='-', color='crimson', label="total error")
ax[1].plot(wl_bins_stat, 100*tr_cbp_syst/tr_cbp, marker='none', linestyle='-', color='orange', label="total systematics")
ax[1].plot(wl_bins_syst, 100*tr_cbp_stat/tr_cbp, marker='none', linestyle='-', color='grey', label="total statistics")
ax[1].set_yscale("log")
ax[1].set_ylim(2e-3, 0.1*100)
ax[1].set_xlim(350, 1050)
ax[1].set_ylabel(r'Relative uncertainties [$\%$]')
ax[1].set_xlabel('$\lambda_c$ [nm]')
ax[1].legend(loc="upper right")
plt.tight_layout()
#fig.savefig("cbp_response.pdf")
plt.savefig("cbp_response.pdf")
plt.show()
# -







# # BACKUP
#
#
# ## Flux-wavelength dependence

# +
wls = {400: "b", 550: "g", 750: "r", 900: "gray"}

fluxes = {}
fluxes_532 = {532: None}

fig  = plt.figure()
counter = 0
for c in cats.keys():
    print(c)
    cat = cats[os.path.basename(c)]
    for wl in wls.keys():
        flux = cat["spectro_laser_flux_per_pulse"] * cat["set_npulses"] * cat["set_nbursts"]
        ind = (cat["set_wl"] == wl) & ~np.isnan(flux) & ~np.isnan(cat["spectro_laser_wl_cal"])
        if wl not in fluxes.keys():
            fluxes[wl] = [flux[ind], cat["spectro_laser_wl"][ind]-cat["set_wl"][ind], cat["spectro_laser_wl_err"][ind]]
        else:
            fluxes[wl][0] = np.concatenate([fluxes[wl][0], flux[ind]])
            fluxes[wl][1] = np.concatenate([fluxes[wl][1], cat["spectro_laser_wl"][ind]-cat["set_wl"][ind]])
            fluxes[wl][2] = np.concatenate([fluxes[wl][2], cat["spectro_laser_wl_err"][ind]])
        ind = (cat["set_wl"] > 560) &  (cat["set_wl"] < 640) &  (cat["spectro_532_wl"] > 525) &  (cat["spectro_532_wl"] < 545) & ~np.isnan(flux)
    if fluxes_532[532] is None:
        fluxes_532[532] = [flux[ind], cat["spectro_532_wl"][ind]-530.8, cat["spectro_532_wl_err"][ind]]
    else:
        fluxes_532[532][0] = np.concatenate([fluxes_532[532][0], flux[ind]])
        fluxes_532[532][1] = np.concatenate([fluxes_532[532][1], cat["spectro_532_wl"][ind]-530.8])
        fluxes_532[532][2] = np.concatenate([fluxes_532[532][2], cat["spectro_532_wl_err"][ind]])
        #plt.plot(flux[ind], cat["spectro_laser_wl_cal"][ind]-cat["set_wl"][ind], "+", color=wls[wl])
        #pval = np.polyfit(flux[ind], cat["spectro_laser_wl_cal"][ind]-cat["set_wl"][ind], deg=1, w=1/cat["spectro_laser_wl_cal_stat"][ind])
        #plt.plot(flux, np.polyval(pval, flux), "-", color=wls[wl]) 
        #ind = (cat["set_wl"] > 560) &  (cat["set_wl"] < 640) &  (cat["spectro_532_wl"] > 525)
        #plt.plot(flux[ind], cat["spectro_532_wl"][ind]-530.8, "+", color="orange")
        
slopes = []
wl=532
deg = 1
pval, pcov = np.polyfit(fluxes_532[wl][0], fluxes_532[wl][1], deg=deg, cov=True)  #, w=1/np.array(fluxes_532[wl][2]), 
p = plt.errorbar(fluxes_532[wl][0], fluxes_532[wl][1], yerr=0*fluxes_532[wl][2], linestyle="none", marker="+", label=f"532: slope={pval[0]:.2g}", color="orange")
x = np.linspace(0, 1.2*np.max(fluxes_532[wl][0]), 50)
plt.plot(x, np.polyval(pval, x), "-", color=p[0].get_color()) #, color=wls[wl]) 
slopes.append(pval[0])
for wl in wls.keys():
    #print(wl, fluxes[wl][0], fluxes[wl][1])
    pval, pcov = np.polyfit(fluxes[wl][0], fluxes[wl][1], deg=deg, cov=True)  #,w=1/np.array(fluxes[wl][2]),
    p = plt.errorbar(fluxes[wl][0], fluxes[wl][1], yerr=fluxes[wl][2], linestyle="none", marker="+", label=f"{wl}: slope={pval[0]:.2g}", color=wls[wl])
    x = np.linspace(0, 1.2*np.max(fluxes[wl][0]), 50)
    plt.plot(x, np.polyval(pval, x), "-", color=p[0].get_color()) #, color=wls[wl]) 
    slopes.append(pval[0])

plt.grid()
plt.xlabel("Stacked spectrum peak maximum $F_p$ [ADU]")
plt.ylabel("$\lambda_g - \lambda_L$ [nm]")
plt.legend(ncol=2)
#plt.ylim(-2.2, 1.5)
plt.title(f"Data release: {data_release}")
#fig.savefig("wavelength_stability.pdf")
plt.show()


# +
wls = np.arange(370,1064,1)

counter = 0
for c in cats.keys():
    #if '.npy' not in c:
    #    continue
    #date = '/'.join(c.split('_')[:3])
    print(c)
    #cat = np.load(os.path.join(data_path, data_release, c), allow_pickle=True)
    cat = cats[os.path.basename(c)]
    for wl in wls:  #.keys():
        flux = cat["spectro_laser_flux_per_pulse"] * cat["set_npulses"]  * cat["set_nbursts"]
        ind = (cat["set_wl"] == wl) & ~np.isnan(flux) & ~np.isnan(cat["spectro_laser_wl"])
        if wl not in fluxes.keys():
            fluxes[wl] = [flux[ind], cat["spectro_laser_wl"][ind]-cat["set_wl"][ind], cat["spectro_laser_wl_err"][ind]]
        else:
            fluxes[wl][0] = np.concatenate([fluxes[wl][0], flux[ind]])
            fluxes[wl][1] = np.concatenate([fluxes[wl][1], cat["spectro_laser_wl"][ind]-cat["set_wl"][ind]])
            fluxes[wl][2] = np.concatenate([fluxes[wl][2], cat["spectro_laser_wl_err"][ind]])
        #plt.plot(flux[ind], cat["spectro_laser_wl_cal"][ind]-cat["set_wl"][ind], "+", color=wls[wl])
        #pval = np.polyfit(flux[ind], cat["spectro_laser_wl_cal"][ind]-cat["set_wl"][ind], deg=1, w=1/cat["spectro_laser_wl_cal_stat"][ind])
        #plt.plot(flux, np.polyval(pval, flux), "-", color=wls[wl]) 
        #ind = (cat["set_wl"] > 560) &  (cat["set_wl"] < 640) &  (cat["spectro_532_wl"] > 525)
        #plt.plot(flux[ind], cat["spectro_532_wl"][ind]-530.8, "+", color="orange")


# +
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy import stats
        
pvals = []
pvals_err = []
deg = 1
for wl in wls:  #.keys():
    #print(wl, fluxes[wl][0], fluxes[wl][1])
    pval, pcov = np.polyfit(fluxes[wl][0], fluxes[wl][1], deg=deg, cov=True)  #, w=1/np.array(fluxes[wl][2]), 
    pvals.append(pval)
    pvals_err.append([np.sqrt(pcov[k,k]) for k in range(deg+1)])
pvals = np.array(pvals)
pvals_err = np.array(pvals_err)

pvals_filtered_func = {}
pvals_filtered_err_func = {}

for k in range(deg):
    fig  = plt.figure()
    params = pvals[:, k]
    params_err = pvals_err[:, k]
    # filt = ~np.isnan(params) &  (params_err < 0.01 * np.max(params_err))
    filt = ~np.isnan(params) &  (params_err < 3*stats.median_abs_deviation(params_err))
    print("coefficient:", k, "filtered data:", 100*np.sum(filt) / len(params), "%")
    params_filtered = savgol_filter(params[filt], 21, 3)
    params_filtered_err = np.std(params_filtered - params[filt]) * np.ones_like(params_filtered)
    #params_filtered = params[filt]
    print(np.mean(params_filtered_err))

    params_filtered_func = interp1d(wls[filt], params_filtered, bounds_error=False, fill_value=0)
    pvals_filtered_func[k] = params_filtered_func
    pvals_filtered_err_func[k] = interp1d(wls[filt], params_filtered_err, bounds_error=False, fill_value=np.mean(params_filtered_err))

    plt.errorbar(wls, params, yerr=params_err, linestyle="none", marker="+", color="k")
    plt.plot(wls[filt], params_filtered, "r-")
    plt.fill_between(wls[filt], params_filtered-params_filtered_err, params_filtered+params_filtered_err, color="r", alpha=0.5)
    plt.xlabel("Set wavelength $\lambda_L$ [nm]")
    plt.ylabel("Wavelength shift slope with peak amplitude [nm/ADU]")
    plt.grid()
    plt.show()
# -

wls = np.arange(350, 1100)
np.save("../../../analysis/cbp_paper/spectro_calibration_flux_dependance_slopes.npy", np.array([wls, pvals_filtered_func[0](wls), pvals_filtered_err_func[0](wls)]).T)

# +
a = np.load("../../../analysis/cbp_paper/spectro_calibration_flux_dependance_slopes.npy")

fig = plt.figure()
plt.plot(a.T[0], a.T[1])
plt.show()
# -

fig  = plt.figure()
counter = 0
for c in catalogs:
    cat = cats[os.path.basename(c)]
    date = '/'.join(c.split('_')[:3])

    flux = cat["spectro_laser_flux_per_pulse"] * cat["set_npulses"]  * cat["set_nbursts"]
    flux_err = cat["spectro_laser_flux_per_pulse_err"] * cat["set_npulses"]  * cat["set_nbursts"]
    new_wl = np.copy(cat["spectro_laser_wl"])
    new_wl_err = np.copy(cat["spectro_laser_wl_err"])**2
    for key in pvals_filtered_func.keys():
        new_wl -= pvals_filtered_func[key](cat["set_wl"]) * flux **(deg-key)
        # TODO: add full cov matrix for error propagation if needed
        new_wl_err += (flux**(deg-key) * pvals_filtered_err_func[key](cat["set_wl"]))**2 +  ((deg-k)*flux_err*flux**(deg-key-1) * pvals_filtered_func[key](cat["set_wl"]))**2
    new_wl_err = np.sqrt(new_wl_err)
    if "set_wl" in cat.dtype.names:
        plt.errorbar(cat["set_wl"], new_wl-cat["set_wl"], yerr=new_wl_err, marker="+", linestyle="none", label=date)
    print(c, "median wl err", np.nanmedian(new_wl_err))
    counter += len(cat["spectro_laser_wl"])
print("Number of data points:" , counter)
plt.grid()
plt.xlabel("Set wavelength $\lambda_L$ [nm]")
plt.ylabel("$\lambda_{det}-\lambda_L$ [nm]")
plt.legend(ncol=2)
plt.ylim(-2.2, 1.5)
plt.xlim(350, 1100)
plt.title(f"Data release: {data_release}")
#fig.savefig("wavelength_stability.pdf")
plt.show()






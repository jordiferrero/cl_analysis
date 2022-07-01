###########
import os, glob
import lumispy as lum
import hyperspy.api as hs

# PARAMS TO MODIFY
# root = os.path.abspath(r'G:\My Drive\PhD\projects\cl_pl_correlation')
# session = os.path.relpath(r'20220212_CL_SED_nFTIR_correlation_day1')

root = os.path.abspath(r'F:\HYPCards_PROD')
session = r'20220630-JORDI'
folder = os.path.join(root, session)

extension = '*\*Card.sur'

# Use the bkg subtracted .sur file by Attolight instead of the RAW data
use_bkg_removed_sur = False

# Bkg remove using lumispy
bkg_removal = True
bkg_name_wildcard = 'BG*s.txt'

spike_removal = True

# Panchromatic ROIs
plot_filtered_pan_im = True
import hyperspy.api as hs
roi1 = hs.roi.SpanROI(left=485., right=520.)      # sets a digitalbandfilter for PbI2 / n=1
roi2 = hs.roi.SpanROI(left=520., right=560.)      # sets a digitalbandfilter for broader PbI2 peak
roi3 = hs.roi.SpanROI(left=750., right=850.)      # sets a digitalbandfilter for 3d perovskite blue shift
roi4 = hs.roi.SpanROI(left=750., right=850.)      # sets a digitalbandfilter for 3d perovskite
rois = [roi1, roi2, roi3, roi4]

cmaps = ['Blues', 'Greens', 'Oranges', 'Reds',]
roi_phase_name = ['PbI2 / n=1', 'broader PbI2', '3D pvk blue-shift', '3D perovskite']

#######################

# LOAD HYPMAP FILES

def get_paths_list(ses_path, bkg_included):
    # Load HYPMaps
    if not bkg_included:
        paths = [p for p in glob.glob(ses_path, recursive=True) if 'HYPCard.sur' in p]
    else:
        paths = [p for p in glob.glob(ses_path, recursive=True) if 'HYPCard' in p and 'HYPCard.sur' not in p]
    return paths


def check_if_file_already_been_processed(base_folder, overwrite=False):
    # Returns True if already processed, else False
    if overwrite:
        return False
    else:
        if os.path.exists(os.path.join(base_folder, 'HYPCard_processed.hspy')):
            return True
        else:
            return False


session_path = os.path.join(folder, extension)
paths_hypmap = get_paths_list(session_path, use_bkg_removed_sur)

if len(paths_hypmap) == 0:
    extension = '*.sur'
    session_path = os.path.join(folder, extension)
    paths_hypmap = get_paths_list(session_path, use_bkg_removed_sur)
    
paths_hypmap.sort()
print([('i={}'.format(i), os.path.split(os.path.dirname(f))[-1]) for i,f in enumerate(paths_hypmap)])

# Metadata pandas list
import pandas as pd
df_names = ['md_ccd', 'md_scan', 'md_sem', 'md_site_image', 'md_spectrometer']
dfs = {key: pd.DataFrame() for key in df_names}

for f in paths_hypmap[:]:
    
    base_folder = os.path.dirname(f)
    bool_file_already_checked = check_if_file_already_been_processed(base_folder, overwrite=False)
    # If file already processed, then skip iteration:
    if bool_file_already_checked:
        continue
    
    cl = hs.load(f, signal_type='CL_SEM')
    name = os.path.split(os.path.dirname(f))[-1].split('_')[0]
    cl.metadata.General.title = name

    # PRE-PROCESSING
    # Background subtraction
    if bkg_removal:
        def bkg_subtraction_single_spectrum(spectrum, background_intensities):
            return spectrum - background_intensities

        dirname = os.path.dirname(f)

        bkg_path = glob.glob(os.path.join(dirname, bkg_name_wildcard))
        if len(bkg_path) != 1:
            import warnings
            warnings.warn("Warning...........Either too many '*bkg-*.txt' files OR none exist in the folder {}".format(os.path.split(dirname)[-1]))
        else:
            import numpy as np

            bkg = np.loadtxt(bkg_path[0])
            bkg_intensity = bkg[1]

            cl.map(bkg_subtraction_single_spectrum, background_intensities=bkg_intensity,
                   show_progressbar=True, parallel=True, inplace=True)

    # Grating shift correction
    calibration_factor = 131072 # From Gunnar

    grating = int(cl.original_metadata.Object_0_Channel_0.Parsed.SPECTROMETER.Grating__Groove_Density)
    magnification = cl.original_metadata.Object_0_Channel_0.Parsed.SEM.Real_Magnification

    if grating == 150:
        correction_factor_grating = 2.73E-04 # 150 gr/mm grating
    elif grating == 600:
        correction_factor_grating = 6.693659836087227e-05 # 600 gr/mm grating
    else:
        raise ImportError('Grating correction not available')

    grating_calibrations = {
        'cal_factor_x_axis' : calibration_factor,
        'corr_factor_grating' : correction_factor_grating,
        'sem_magnification' : magnification,
    }

    cl.correct_grating_shift(**grating_calibrations)


    # Crop edges
    cl_processed = cl.crop_edges(crop_px=5)
    
    # Remove spikes
    if spike_removal:
        cl_processed = cl_processed.remove_spikes()

    # Get eV datafiles
    cl_processed_ev = cl_processed.to_eV(inplace=False)

    # PLOTING
    import matplotlib.pyplot as plt
    
    cl_processed.T.mean().plot(cmap='Greys') # Panchromatic
    plt.savefig(os.path.join(base_folder, 'im_panchromatic_{}.tiff'.format(name)))
    plt.close()

    # Spectral figure
    title = cl_processed.metadata.General.title
    x = cl_processed.axes_manager.signal_axes[0].axis
    y = cl_processed.mean().data
    x_ev = cl_processed_ev.axes_manager.signal_axes[0].axis
    y_ev = cl_processed_ev.mean().data
    f, axs = plt.subplots(ncols=2, figsize=(10,4))
    axs[0].plot(x, y, c='C0')
    axs[1].plot(x_ev, y_ev, c='C1')
    axs[0].set_xlabel("Wavelength (nm)")
    axs[1].set_xlabel("Energy (eV)")
    axs[0].set_ylabel("Mean CL intensity (cnts)")
    f.suptitle(f"Mean spectrum for {title}")
    plt.savefig(os.path.join(base_folder, 'im_spectrum_{}.tiff'.format(name)))
    plt.close()
    
    # Panchromatic filtered
    if plot_filtered_pan_im:
        im = cl_processed.T

        n_rois = len(rois)
        from matplotlib_scalebar.scalebar import ScaleBar
        fig = plt.figure(figsize=(4*n_rois, 5),)
        gridsize = (1, n_rois)
        plt.matplotlib.gridspec.GridSpec(gridsize[0], gridsize[1])

        for i in range(n_rois):
            ax = plt.subplot2grid(gridsize, (0,i))
            roi = rois[i]
            if len(cmaps) == len(rois):
                cmap_i = cmaps[i]
            else:
                cmap_i = 'viridis'
            
            if len(roi_phase_name) == len(rois):
                phase_name_i = roi_phase_name[i]
            else:
                phase_name_i = ''
                
            im_roi = ax.imshow(roi(im).mean().data, cmap=cmap_i)
            roi_width = roi.right - roi.left
            roi_centre = roi.left + 0.5*roi_width
            title = '{} \n{:.0f} $\pm$ {:.0f} nm'.format(phase_name_i, roi_centre, roi_width/2)
            ax.set_title(title, color='k')
            plt.axis('off')
            plt.colorbar(im_roi, shrink=0.75, pad=0.02)

        units = roi(im).mean().axes_manager[0].units
        scale = roi(im).mean().axes_manager[0].scale
        scalebar = ScaleBar(scale, units, location='lower right') # 1 pixel = 0.2 meter
        plt.gca().add_artist(scalebar)
        plt.suptitle(roi(im).mean().metadata.General.title + ' bandpass filter')
        plt.axis('off')

        plt.tight_layout(rect=(0,0.03,1,1))
        plt.savefig(os.path.join(base_folder, 'im_bandpass_filter_{}.tiff'.format(name)))
        plt.close()

    # Save processed file
    cl_processed.save(os.path.join(base_folder, 'HYPCard_processed.hspy'), overwrite=True)

    # Save metadata
    md = cl.original_metadata.Object_0_Channel_0.Parsed

    md_all = {'md_ccd': md.CCD.as_dictionary(),
              'md_scan': md.SCAN.as_dictionary(),
              'md_sem': md.SEM.as_dictionary(),
              'md_site_image': md.SITE_IMAGE.as_dictionary(),
              'md_spectrometer': md.SPECTROMETER.as_dictionary(),
              }
    for key, value in md_all.items():
        num = name.split('P')[-1]
        if num.isnumeric():
            num = int(num)
        else:
            num = name

        df = pd.DataFrame(value, index=[name])
        df['map_number'] = num

        # Add the calibration of the axis
        ax = cl.axes_manager.navigation_axes
        nx, ny = ax[0].scale, ax[1].scale
        units = ax[0].units

        df['pixel_size_nx'] = nx
        df['pixel_size_ny'] = ny
        df['pixel_size_unit'] = units

        # Append to the dfs
        dfs[key] = dfs[key].append(df)

    print('{} completed'.format(name))

# Save pandas metadata files
for key, df in dfs.items():
    try:
        path_pd = os.path.join(folder, 'metadata_dataframes', (key + '.csv'))
        df.to_csv(path_pd)
    except FileNotFoundError:
        path_pd = os.path.join(folder, (key + '.csv'))
        df.to_csv(path_pd)

print("Metadata csv's saved successfully!")

""" RUN THE NOTEBOOK FIRST TO GET AN IDEA OF:

- If any of the cl maps need to be cropped beforehand
- To know the ideal INITIAL fitting parameters

"""
# PARAMS
import os
root = os.path.abspath(r'G:\My Drive\PhD\projects\external_measurements')
session = os.path.relpath(r'20220118-CsAgBiBr3_2D_3D_films')
extension = '*\*processed*.hspy'

i_interest = []

rebin_nav = 1
rebin_sig = 2

# Initial conditions



###############
#CODE

import lumispy as lum
import hyperspy.api as hs
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='darkgrid')


# ## Loading HYPCard files

# GO TO AUTOMATISED TO GET GENERAL PROCESSED DATA
import os, glob
folder = os.path.join(root, session)
session_path = os.path.join(root, session, extension)
# For HYPMaps
paths_hypmap = [p for p in glob.glob(session_path, recursive=True)]

# Filter out i not of interest
if len(i_interest) != 0:
    paths_hypmap = [p for i, p in enumerate(paths_hypmap) if i in i_interest]

paths_hypmap.sort()
print([('i={}'.format(i), os.path.split(os.path.dirname(f))[-1]) for i,f in enumerate(paths_hypmap)])

cls = hs.load(paths_hypmap, signal_type='CL_SEM')
print(cls)

#Start for loop to iterate over all the files in the folder
for cl in cls[:]:

    # FITTING PART
    def fit_gaussians(cl_object):

        # Create model (with 2 Gaussian)
        m = cl_object.create_model()

        g_pvk = hs.model.components1D.Expression(
            expression="height * exp(-(x - x0) ** 2 * 4 * log(2)/ fwhm ** 2)",
            name="Perovskite",
            position="x0",
            height=1,
            fwhm=1,
            x0=1,
            module="numpy")

        g_pbi2 = hs.model.components1D.Expression(
            expression="height * exp(-(x - x0) ** 2 * 4 * log(2)/ fwhm ** 2)",
            name="broad_red_peak",
            position="x0",
            height=1,
            fwhm=1,
            x0=1,
            module="numpy")

        # g_decay = hs.model.components1D.Expression(
        #     expression="height * exp(-(x - x0) ** 2 * 4 * log(2)/ fwhm ** 2)",
        #     name="Decay peak",
        #     position="x0",
        #     height=1,
        #     fwhm=1,
        #     x0=1,
        #     module="numpy")

        bkg = hs.model.components1D.Polynomial(order=1)

        # m.extend([g_pvk, g_pbi2, g_decay, bkg])
        m.extend([g_pvk, g_pbi2, bkg])

        # EDIT HERE
        ######################
        # COPY HERE
        # Set limits to the gausians
        g_pvk.x0.value = 2.60
        g_pvk.x0.bmin, g_pvk.x0.bmax = 2.4, 2.8
        g_pvk.height.value, g_pvk.height.bmin, g_pvk.height.bmax = 10, 0, 500
        g_pvk.fwhm.value, g_pvk.fwhm.bmin, g_pvk.fwhm.bmax = 0.5, 0.005, 2

        g_pbi2.x0.value = 1.93
        g_pbi2.x0.bmin, g_pbi2.x0.bmax = 1.7, 2.1
        g_pbi2.height.value, g_pbi2.height.bmin, g_pbi2.height.bmax = 50, 0, 1000
        g_pbi2.fwhm.value, g_pbi2.fwhm.bmin, g_pbi2.fwhm.bmax = 0.5, 0.005, 2

#         bkg.a.value, bkg.a.bmin, bkg.a.bmax = -1, -100, 100
#         bkg.b.value, bkg.b.bmin, bkg.b.bmax = 5, 0, 100

        # g_decay.x0.value = 650
        # g_decay.x0.bmin, g_decay.x0.bmax = 510, 750
        # g_decay.height.value, g_decay.height.bmin, g_decay.height.bmax = 10, -10, 100
        # g_decay.fwhm.value, g_decay.fwhm.bmin, g_decay.fwhm.bmax = 20, 0, 100
        # END COPY HERE
        ######################

        #Fit for all the positions
        m.assign_current_values_to_all()
        m.fit(bounded=True)
        m.multifit(bounded=True, show_progressbar=True, iterpath='serpentine')

        # Plot
        import matplotlib.pyplot as plt
        m.plot(plot_components=True)
        plt.savefig(os.path.join(folder, 'fit_imgs', '{}_spectrum.tiff'.format(cl_object.metadata.General.title)))
        plt.close()

        # return m, g_pvk, g_pbi2, g_decay, bkg
        return m, g_pvk, g_pbi2, bkg


    # Make new folder
    try:
        os.mkdir(os.path.join(folder, 'fit_imgs'))
    except:
        pass

    # Rebin
    scale = (rebin_nav, rebin_nav, rebin_sig)
    cl = cl.rebin(scale=scale)
    
    # Convert to eV
    cl = cl.to_eV(inplace=False)

    # Run the fit gaussian function
    # m, g_pvk, g_pbi2, g_decay, bkg = fit_gaussians(cl)
    m, g_pvk, g_pbi2, bkg = fit_gaussians(cl)

    path_m = os.path.join(folder, 'fit_imgs', '{}_binned_{}{}_{}_fitted.hspy'.format(cl.metadata.General.title, rebin_nav, rebin_nav, rebin_sig))
    m.save(path_m, "gaus_fit", overwrite=True)
    print('{} map fitted.'.format(cl.metadata.General.title))

print('Fitting finished and saved.')



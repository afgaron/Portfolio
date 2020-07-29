'''plots all of the desired plots for grismProcessing'''

import numpy as np
import matplotlib
matplotlib.use('AGG') #non-interactive backend for pngs
from matplotlib import pyplot as plt
from matplotlib.colors import SymLogNorm

def plotCompressedGrism(gimg, pRange, gColumnTotalFine, wpa, minimization, plotDir):
    '''for use after finding the offsets of the compressed grism'''
    f, (dataPlot, residPlot) = plt.subplots(2)
    dataPlot.plot(gimg, alpha=0.1)
    dataPlot.plot(pRange, gColumnTotalFine, label='mean')
    outputModel = np.zeros(len(pRange))
    for profNum, prof in enumerate(wpa):
        dataPlot.plot(pRange, prof, label='model %i' % profNum)
        outputModel += prof
    
    dataPlot.legend(loc='best')
    dataPlot.set_ylabel('intensity')
    dataPlot.set_title('Compressed grism: %s\nweights = %s' % ('success' if minimization['success'] else 'failure', minimization['x'][0::2]))
    residPlot.plot(pRange, gColumnTotalFine - outputModel, label='residual')
    residPlot.legend(loc='best')
    residPlot.set_xlabel('y (pixels)')
    residPlot.set_ylabel('residual')
    plt.savefig('%s/compressedGrism.png' % plotDir)
    plt.close(f)
    return

def plotIndividualWavelengths(pRange, gColumnFine, weightedProfArrays, minimization, x, plotDir):
    '''for use after finding the amplitudes at specific wavelengths'''
    f, (dataPlot, residPlot) = plt.subplots(2)
    dataPlot.plot(pRange, gColumnFine, label='observed')
    model = np.zeros(len(pRange))
    for profNum, prof in enumerate(weightedProfArrays):
        dataPlot.plot(pRange, prof, label='model %i' % profNum)
        model += prof
    
    dataPlot.legend(loc='best')
    dataPlot.set_ylabel('intensity')
    dataPlot.set_title('At x = %i: %s\nweights = %s' % (x, 'success' if minimization['success'] else 'failure', minimization['x']))
    residPlot.plot(pRange, gColumnFine - model, label='residual')
    residPlot.legend(loc='best')
    residPlot.set_xlabel('y (pixels)')
    residPlot.set_ylabel('residual')
    plt.savefig('%s/grismAt%i.png' % (plotDir, x))
    plt.close(f)
    return

def plotFinalResults(entry, img, totalSuccess, profile, c_profiles, pRange, gimgMasked, gColumnTotalFine, contamGimg, wpa, subtractGimg, stageDir):
    '''for use after fully subtracting the contaminants'''
    gdx = gimgMasked.shape[1]
    gdy = gimgMasked.shape[0]
    xmin, xmax, ymin, ymax = max(entry['X_IMAGE']-gdx/2., 0), min(entry['X_IMAGE']+gdx/2., img.shape[1]-1), \
                             max(entry['Y_IMAGE']-gdy/2., 0), min(entry['Y_IMAGE']+gdy/2., img.shape[0]-1)
    
    f, ((directImagePlot, directImageProfilePlot), (grismPlot, grismProfilePlot), (contamPlot, contamProfilePlot), \
        (subtractedPlot, subtractedProfilePlot)) = plt.subplots(4, 2, sharey=True)

    #plot a grism-sized piece of the direct image
    crop = img[ymin:ymax, xmin:xmax]
    directImagePlot.imshow(crop, norm=SymLogNorm(0.1))
    if totalSuccess:
        directImagePlot.set_title('object %i' % entry['NUMBER'])
    else:
        directImagePlot.set_title('object %i (minimization unsuccessful)' % entry['NUMBER'])
    
    directImagePlot.set_ylabel('direct')
    directImagePlot.axis([0, gdx, 0, gdy])
    directImageProfile = profile
    for profNum in c_profiles:
        directImageProfile += c_profiles[profNum]
    
    directImageProfilePlot.plot(directImageProfile, pRange)
    directImageProfilePlot.set_title('profiles')
    
    #plot the measured grism and its intensity profile
    vmin = np.ma.min(gimgMasked)
    vmax = np.ma.max(gimgMasked)
    grismPlot.imshow(gimgMasked, norm=SymLogNorm(0.1, vmin=vmin, vmax=vmax))
    grismPlot.set_ylabel('grism')
    grismProfilePlot.plot(gColumnTotalFine, pRange)
    
    #plot the contamination model and its intensity profile
    contamPlot.imshow(contamGimg, norm=SymLogNorm(0.01))
    contamPlot.set_ylabel('model')
    contamProfilePlot.plot(sum(wpa), pRange)
    
    #plot the subtracted grism and its intensity profile
    subtractedPlot.imshow(subtractGimg, norm=SymLogNorm(0.1, vmin=vmin, vmax=vmax))
    subtractedPlot.set_ylabel('subtracted')
    subtractedProfile = np.ma.average(subtractGimg, axis=-1)
    subtractedProfilePlot.plot(subtractedProfile, xrange(gdy))
    
    plt.savefig('%s/stages%i.png' % (stageDir, entry['NUMBER']))
    plt.close(f)
    return

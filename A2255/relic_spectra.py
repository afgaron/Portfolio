import numpy as np
from pprint import pprint
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rc, path
from matplotlib.ticker import AutoMinorLocator
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import wcs
from scipy import special, integrate, optimize
from scipy.interpolate import interp1d
from astropy.modeling.models import Linear1D, Gaussian1D
from astropy.modeling.fitting import LevMarLSQFitter
import operator, itertools
from collections import OrderedDict
plt.ion()
rc('text', usetex=True)

# dictionary of filepaths
files = {'P16':'final_maps_P/A2255_P_regridL_16.im.fits',
		 'P30':'final_maps_P/A2255_P_regridL_30.im.fits',
		 'P30x16':'final_maps_P/A2255_P_regridL_30x16.im.fits',
		 'P180x16':'in_progress_fits/A2255_P_regridL_180x16.im.fits',
		 'L16':'final_maps_P/A2255_L_16.image.tt0.fits',
		 'L30':'final_maps_P/A2255_L_30.image.tt0.fits',
		 'L30x16':'final_maps_P/A2255_L_30x16.image.tt0.fits',
		 'L180x16':'in_progress_fits/A2255_L_180x16.image.tt0.fits',
		 'PL16':'final_maps_P/A2255_PL_16_3sigma.alpha.fits',
		 'LL16':'final_maps_P/A2255_L_16.alpha.fits',
		 'PL30':'final_maps_P/A2255_PL_30_3sigma.alpha.fits',
		 'LL30':'final_maps_P/A2255_L_30.alpha.fits',
		 'LL30x16':'final_maps_P/A2255_L_30x16.alpha.fits',
		 'LL180x16':'in_progress_fits/A2255_L_180x16.alpha.fits',
		 'L_uvtaper':'final_maps_P/A2255_L_uvtaper.image.tt0.fits',
		 'L_hi-res':'FINAL_MAPS_L/A2255C-allms-poldata.pbcor.image.tt0.4.fits'}

# global variables
pm = np.array([1., -1.])
T, F = True, False
nu1, nu3, nuL, nuP = 1.26e9, 1.78e9, 1.52e9, 3.68e8 # in Hz
scale = 1.563 # kpc/arcsec
rms_dict = {'L': {'16':3.80432e-5, '30':1.09882e-4, '30x16':6.32815e-5,                  '180x16':2.61971e-4},
			'P': {'16':1.40804e-4, '30':1.77871e-4, '30x16':1.51664e-4,      '180x16':5.06084e-4}} # ra=17:14:30, dec=+64:07:25, radius=175"

# manually define crop size and axes labels
axes_dict = {'full': {'x1':1, 'y1':1, 'x2':401, 'y2':521,
					  'ticks_x':6, 'ticks_y':6,
					  'ra_ticks':['17:15:30', '17:14:30', '17:13:30', '17:12:30', '17:11:30'],
					  'ra_labels':['$17^h15^m30^s$', '$14^m30^s$', '$13^m30^s$', '$12^m30^s$', '$11^m30^s$'],
					  'dec_ticks':['+63:48:00', '+63:54:00', '+64:00:00', '+64:06:00', '+64:12:00', '+64:18:00', '+64:24:00'],
					  'dec_labels':["$+63^\circ48'$", "$54'$", "$+64^\circ00'$", "$06'$", "$12'$", "$18'$", "$24'$"]},
			 'central': {'x1':178, 'y1':201, 'x2':310, 'y2':286,
					  'ticks_x':5, 'ticks_y':5,
					  'ra_ticks':['17:13:30', '17:13:10', '17:12:50', '17:12:30', '17:12:10'],
					  'ra_labels':['$17^h13^m30^s$', '$13^m10^s$', '$12^m50^s$', '$12^m30^s$', '$12^m10^s$'],
					  'dec_ticks':['+64:01:30', '+64:02:45', '+64:04:00', '+64:05:15', '+64:06:30'],
					  'dec_labels':["$+64^\circ01'30''$", "$2'45''$", "$4'00''$", "$5'15''$", "$6'30''$"]},
			 'embryo': {'x1':22, 'y1':175, 'x2':81, 'y2':257,
					  'ticks_x':5, 'ticks_y':4,
					  'ra_ticks':['17:15:20', '17:15:10', '17:15:00', '17:14:50'],
					  'ra_labels':['$17^h15^m20^s$', '$15^m10^s$', '$15^m00^s$', '$14^m50^s$'],
					  'dec_ticks':['+64:00:00', '+64:01:00', '+64:02:00', '+64:03:00', '+64:04:00', '+64:05:00'],
					  'dec_labels':["$+64^\circ00'$", "$1'$", "$2'$", "$3'$", "$4'$", "$5'$"]},
			 'beaver': {'x1':171, 'y1':12, 'x2':236, 'y2':105,
					  'ticks_x':4, 'ticks_y':6,
					  'ra_ticks':['17:13:40', '17:13:20', '17:13:00'],
					  'ra_labels':['$17^h13^m40^s$', '$13^m20^s$', '$13^m00^s$'],
					  'dec_ticks':['+63:47:00', '+63:48:30', '+63:50:00', '+63:51:30', '+63:53:00'],
					  'dec_labels':["$+63^\circ47'00''$", "$48'30''$", "$50'00''$", "$51'30''$", "$53'00''$"]},
			 'trg': {'x1':250, 'y1':204, 'x2':292, 'y2':255,
					  'ticks_x':5, 'ticks_y':6,
					  'ra_ticks':['17:12:40', '17:12:30', '17:12:20'],
					  'ra_labels':['$17^h12^m40^s$', '$12^m30^s$', '$12^m20^s$'],
					  'dec_ticks':['+64:02:00', '+64:03:00', '+64:04:00', '+64:05:00'],
					  'dec_labels':["$+64^\circ02'$", "$3'$", "$4'$", "$5'$"]},
			 'trg_hi-res': {'x1':1117, 'y1':893, 'x2':1315, 'y2':1133,
					  'ticks_x':5, 'ticks_y':6,
					  'ra_ticks':['17:12:40', '17:12:30', '17:12:20'],
					  'ra_labels':['$17^h12^m40^s$', '$12^m30^s$', '$12^m20^s$'],
					  'dec_ticks':['+64:02:00', '+64:03:00', '+64:04:00', '+64:05:00'],
					  'dec_labels':["$+64^\circ02'$", "$3'$", "$4'$", "$5'$"]},
			 'ne_relic': {'x1':167, 'y1':342, 'x2':319, 'y2':443,
					  'ticks_x':5, 'ticks_y':6,
					  'ra_ticks':['17:13:30', '17:13:15', '17:13:00', '17:12:45', '17:12:30', '17:12:15', '17:12:00'],
					  'ra_labels':['$17^h13^m30^s$', '$13^m15^s$', '$13^m00^s$', '$12^m45^s$', '$12^m30^s$', '$12^m15^s$', '$12^m00^s$'],
					  #'ra_ticks':['17:13:15','17:12:45', '17:12:15'],
					  #'ra_labels':['$17^\circ13^m15^s$', '$12^m45^s$', '$12^m15^s$'],
					  #'dec_ticks':['+64:11:30', '+64:13:00', '+64:14:30', '+64:16:00', '+64:17:30', '+64:19:00'],
					  #'dec_labels':["$+64^\circ11'30''$", "$13'00''$", "$14'30''$", "$16'00''$", "$17'30''$", "$19'00''$"]
					  'dec_ticks':['+64:12:00', '+64:13:00', '+64:14:00', '+64:15:00', '+64:16:00', '+64:17:00', '+64:18:00', '+64:19:00'],
					  'dec_labels':["$+64^\circ12'$", "$13'$", "$14'$", "$15'$", "$16'$", "$17'$", "$18'$", "$19'$"]},
			 'goldfish': {'x1':185, 'y1':230, 'x2':236, 'y2':290,
					  'ticks_x':5, 'ticks_y':6,
					  'ra_ticks':['17:13:30', '17:13:20', '17:13:10', '17:13:00'],
					  'ra_labels':['$17^h13^m30^s$', '$13^m20^s$', '$13^m10^s$', '$13^m00^s$'],
					  'dec_ticks':['+64:03:30', '+64:04:30', '+64:05:30', '+64:06:30'],
					  'dec_labels':["$+64^\circ03'30''$", "$4'30''$", "$5'30''$", "$6'30''$"]}}

####################################
def fitTRG_linear():
	''' fitting linear fit to bright TRG for error determination'''
	ix = data['flux']>2
	x = 1000/rms*data['flux'][ix]
	y = data['alpha_mean'][ix]

	# best fit parameters (any poly fit)
	coeff = np.polyfit(x, y, 1)
	x0 = np.linspace(min(x), max(x), 100)
	y0 = sum([c*x0**ix for ix,c in enumerate(np.flip(coeff))])

	# error calc (lin fit)
	S = np.sqrt( sum((y-coeff[0]*x-coeff[1])**2) / (len(x)-2) )
	err_m = S * np.sqrt( len(x) / (len(x)*sum(x**2)-sum(x)**2) ) / np.abs(coeff[0])
	err_b = S * np.sqrt( sum(x**2) / (len(x)*sum(x**2)-sum(x)**2) ) / np.abs(coeff[1])

	# plot best fit line
	plt.figure()
	plt.errorbar(1000/rms*data['flux'], data['alpha_mean'], fmt='.', label='TRG data', color='C0')
	plt.plot(x0, y0, label='best fit', color='C1')
	plt.fill_between(x0, (1-err_m)*coeff[0]*x0+(1-err_b)*coeff[1], (1+err_m)*coeff[0]*x0+(1+err_b)*coeff[1], alpha=.5, color='C1')
	plt.xlabel('Flux [std dev]')
	plt.ylabel('alpha')
	plt.legend(loc='lower right')
	plt.tight_layout()

####################################
def fitTRG_power():
	''' fitting alpha-flux power law for TRG'''
	ix = np.hstack([np.arange(3,8), np.arange(9, 16)])
	x = np.log10(1000/rms*data['flux'][ix])
	y = np.log10(-1*data['alpha_mean'][ix])

	# best fit parameters (any poly fit)
	coeff = np.polyfit(x, y, 1)
	x0 = np.linspace(min(x), max(x), 100)
	y0 = sum([c*x0**ix for ix,c in enumerate(np.flip(coeff))])

	# error calc (lin fit)
	S = np.sqrt( sum((y-coeff[0]*x-coeff[1])**2) / (len(x)-2) )
	err_m = S * np.sqrt( len(x) / (len(x)*sum(x**2)-sum(x)**2) ) / np.abs(coeff[0])
	err_b = S * np.sqrt( sum(x**2) / (len(x)*sum(x**2)-sum(x)**2) ) / np.abs(coeff[1])

	# plot best fit line
	plt.figure()
	plt.errorbar(x, y, fmt='.', label='TRG data', color='C0')
	plt.plot(x0, y0, label='best fit', color='C1')
	plt.fill_between(x0, (1-err_m)*coeff[0]*x0+(1-err_b)*coeff[1], (1+err_m)*coeff[0]*x0+(1+err_b)*coeff[1], alpha=.5, color='C1')
	plt.xlabel('log10(Flux [std dev])')
	plt.ylabel('log10(-alpha)')
	plt.legend(loc='lower left')
	plt.tight_layout()

	# convert best fit back to linear space
	x_prime = 1000/rms*data['flux']
	y_prime = data['alpha_mean']
	x1 = np.linspace(min(x_prime), max(x_prime), 100)

	# plot best fit power law
	plt.figure()
	plt.errorbar(x_prime, y_prime, fmt='.', label='TRG data', color='C0')
	plt.plot(x1, -10**coeff[1] * x1**coeff[0], label='best fit', color='C1')
	plt.fill_between(x1, -10**((1-err_b)*coeff[1]) * x1**((1-err_m)*coeff[0]), -10**((1+err_b)*coeff[1]) * x1**((1+err_m)*coeff[0]), alpha=.5, color='C1')
	plt.xlabel('Flux (std dev)')
	plt.ylabel('Mean alpha')
	plt.legend(loc='lower right')
	plt.tight_layout()
	
####################################
def read_regions(regfile, fitsimage, n=1):
	'''given a region file and a fits image, determine the appropriate shape of the region and read the corresponding values from the fits image
	returns the center pixel of circles; every nth pixel of boxes or polygons; and every pixel along a line'''

	# load fits image
	with fits.open(fitsimage) as hdul:
		data = hdul[0].data[0][0]
		w = wcs.WCS(hdul[0])
	
	# load shape from region file
	with open(regfile, 'r') as f:
		lines = f.readlines()

	shape, circles, polygons, boxes, profiles = '', [], [], [], []
	for line in lines:
		if line[:6]=='circle':
			shape = 'circle'
			_line = line[7:].split(',')
			circles.append((_line[:2], _line[2].split('"')[0]))
		
		elif line[:7]=='polygon':
			shape = 'polygon'
			polygons.append(line.split('(')[1].split(')')[0].split(','))
		
		elif line[:3]=='box':
			shape = 'box'
			boxes.append(line.split('(')[1].split(')')[0].split(','))
		
		elif line[:4]=='line':
			shape = 'profile'
			profiles.append(line.split('(')[1].split(')')[0].split(','))

	if shape=='':
		raise ValueError('Region file must contain polygons, boxes, circles, or lines')
	else:
		print 'Input shape is:', shape

	# convert shape to pixel values
	coords_pix = []
	if shape=='circle':
		c_wcs = []
		for circle in circles:
			c_wcs.append(SkyCoord(*circle[0], unit=(u.hourangle if circle[0][0][2]==':' else u.deg, u.deg)))
		
		for c in c_wcs:
			c_pix = w.wcs_world2pix(c.ra, c.dec, 0, 0, 1)
			coords_pix.append([c_pix[0], c_pix[1]])
	
	elif shape=='polygon':
		for polygon in polygons:
			_coords = []
			for ix in range(0, len(polygon), 2):
				c_wcs = SkyCoord(polygon[ix], polygon[ix+1], unit=(u.hourangle if polygon[0][2]==':' else u.deg, u.deg))
				c_pix = np.array(w.wcs_world2pix(c_wcs.ra, c_wcs.dec, 0, 0, 1))
				_coords.append([c_pix[0], c_pix[1]])
			_coords.append(_coords[0])
			coords_pix.append(_coords)
		
	elif shape=='box':
		for box in boxes:
			c_wcs = SkyCoord(box[0], box[1], unit=(u.hourangle if box[0][2]==':' else u.deg, u.deg))
			c_pix = np.array(w.wcs_world2pix(c_wcs.ra, c_wcs.dec, 0, 0, 1)[:2])
			xdelt = pm* float(box[2][:-1])/(2*w.to_header()['CDELT1']*3600)
			ydelt = pm* float(box[3][:-1])/(2*w.to_header()['CDELT2']*3600)
			_coords = np.array([c_pix + np.array([xdelt[0], ydelt[0]]),
								c_pix + np.array([xdelt[0], ydelt[1]]),
								c_pix + np.array([xdelt[1], ydelt[1]]),
								c_pix + np.array([xdelt[1], ydelt[0]])])
			
			theta = float(box[4])*np.pi/180.
			rot_matrix = np.array([[np.cos(theta),  np.sin(theta)],
								   [-np.sin(theta), np.cos(theta)]])
			trans_matrix = np.array([-c_pix]*4)
			
			_coords = np.matmul(_coords+trans_matrix, rot_matrix) - trans_matrix
			coords_pix.append(np.vstack([_coords, _coords[0]]))
	
	elif shape=='profile':
		for profile in profiles:
			c0 = SkyCoord(profile[0], profile[1], unit=(u.hourangle if profile[0][2]==':' else u.deg, u.deg))
			c1 = SkyCoord(profile[2], profile[3], unit=(u.hourangle if profile[0][2]==':' else u.deg, u.deg))
			c0_pix = w.wcs_world2pix(c0.ra, c0.dec, 0, 0, 1)
			c1_pix = w.wcs_world2pix(c1.ra, c1.dec, 0, 0, 1)
			x0 = int( np.floor(c0_pix[0]) if c0_pix[0]<c1_pix[0] else np.ceil(c0_pix[0]) )
			x1 = int( np.floor(c1_pix[0]) if c0_pix[0]>c1_pix[0] else np.ceil(c1_pix[0]) )
			y0 = int( np.floor(c0_pix[1]) if c0_pix[1]<c1_pix[1] else np.ceil(c0_pix[1]) )
			y1 = int( np.floor(c1_pix[1]) if c0_pix[1]>c1_pix[1] else np.ceil(c1_pix[1]) )
			l = int(round(np.sqrt((x1-x0)**2+(y1-y0)**2)))
			xx = np.linspace(x0, x1, l)
			yy = interp1d([x0, x1], [y0, y1])(xx)
			xx, yy = np.floor(xx), np.floor(yy)
			xx[-1], yy[-1] = xx[-1]-1, yy[-1]-1
			coords_pix.append(np.vstack([xx, yy]).T)
	
	# read the values for every nth pixel in the regions, or the center of each circle
	if shape == 'circle':

		shape_obj = path.Path(coords_pix)
		bbox = shape_obj.get_extents()

		vals, samples = [], []
		for i, j in shape_obj.vertices:
			vals.append(data[int(np.round(j)), int(np.round(i))])
			samples.append([int(np.round(i)), int(np.round(j))])
		
		region_list = [{'vals':np.array(vals), 'samples':np.array(samples), 'ix':0}]
	
	else: # box, polygon, or line
	
		region_list = []
		for ix, region in enumerate(coords_pix):
			
			shape_obj = path.Path(region)
			bbox = shape_obj.get_extents()
			
			vals, samples = [], []
			if shape=='box' or shape=='polygon':
				for i in np.arange(np.floor(bbox.xmin), np.ceil(bbox.xmax)):
					for j in np.arange(np.floor(bbox.ymin), np.ceil(bbox.ymax)):
						if not (i+j)%n and shape_obj.contains_point([i, j]):
							vals.append(data[int(j), int(i)])
							samples.append([i, j])
			
			else: # line
				for i, j in shape_obj.vertices:
					vals.append(data[int(j), int(i)])
					samples.append([i, j])
			
			region_list.append({'vals':np.array(vals), 'samples':np.array(samples), 'ix':ix})
		
	return region_list

####################################
def color_I(regfile, freq, spec, beam=16, n=1, color_regions=True):
	'''given a region file, plot alpha vs flux
	freq is "P" or "L"
	spec is "PL" or "LL"
	beam is 16 or 30
	color_regions=False plots each point a different color'''
	
	assert freq in ['P', 'L'], 'freq must be P or L'
	assert spec in ['PL', 'LL'], 'spec must be PL or LL'
	assert beam in [16, 30], 'beam must be 16 or 30'

	# input files, all with same beam and gridding
	fluxfile = files[freq+str(beam)]
	alphafile = files[spec+str(beam)]
	with fits.open(fluxfile) as hdul:
		data = hdul[0].data[0][0]

	# read region files
	fluxes = read_regions(regfile, fluxfile, n)
	alphas = read_regions(regfile, alphafile, n)

	# Strip infs and plot
	fig, (ax0, ax1) = plt.subplots(1,2)
	
	for flux, alpha in zip(fluxes, alphas):
		color = 'C%i'%flux['ix'] if color_regions else None
		vals_f = flux['vals']
		vals_a = alpha['vals']
		locs = flux['samples']
		mask = np.isfinite(vals_a)
		ax0.scatter(vals_f[mask], vals_a[mask], s=2, c=color)
		ax1.scatter(*locs.T, s=2, c=color)
	
	ax0.set_xlabel('Flux (Jy/beam)')
	ax0.set_ylabel('alpha')
	ax1.axis('scaled')
	
	xlim = ax1.get_xlim()
	ax1.set_xlim(xlim -pm* 0.05*(xlim[1]-xlim[0]))
	ylim = ax1.get_ylim()
	ax1.set_ylim(ylim -pm* 0.05*(ylim[1]-ylim[0]))
	
	ax1.imshow(np.log10(data))
	fig.suptitle('%s: %s vs %s' % (regfile.split('_')[0], freq, spec))
	fig.tight_layout(rect=[0.03, 0.03, 1, 0.97])

####################################
def get_alpha(S1, S2, nu1, nu3):
	return (np.log(S2)-np.log(S1)) / (np.log(nu3)-np.log(nu1))

def color_color(regfile, beam=16, flux_lim=-np.inf, n=1, color_regions=True):
	'''given a region file containing a bunch of circles, get the flux, alpha, and error for each circle and plot alpha-alpha
	flux_lim is the minimum log10(flux) to include
	color_regions=False plots each point a different color'''
	
	assert beam in [16, 30], 'beam must be 16 or 30'
	
	# input files, all with 16" beam and same gridding
	fluxfile_P = files['P'+str(beam)]
	alphafile_PL = files['PL'+str(beam)]
	fluxfile_L = files['L'+str(beam)]
	alphafile_LL = files['LL'+str(beam)]
	with fits.open(fluxfile_L) as hdul:
		data = hdul[0].data[0][0]
	
	# read region files
	fluxes_P = read_regions(regfile, fluxfile_P, n)
	alphas_PL = read_regions(regfile, alphafile_PL, n)
	fluxes_L = read_regions(regfile, fluxfile_L, n)
	alphas_LL = read_regions(regfile, alphafile_LL, n)
	
	fig, (ax0, ax1) = plt.subplots(1,2)
	for flux_P, alpha_PL, flux_L, alpha_LL in zip(fluxes_P, alphas_PL, fluxes_L, alphas_LL):
		
		color = 'C%i'%alpha_PL['ix'] if color_regions and len(fluxes_P)<=10 else None
		
		# calculate error bars
		rms_L, rms_P = rms_dict['L'][str(beam)], rms_dict['P'][str(beam)]
		flux_1 = flux_L['vals'] * (nu1/nuL)**alpha_LL['vals']
		flux_2 = flux_L['vals'] * (nu3/nuL)**alpha_LL['vals']
		sigma_PL = np.sqrt( (rms_L/flux_L['vals'])**2 + (rms_P/flux_P['vals'])**2 ) / (np.log(nuL)-np.log(nuP))
		sigma_LL = np.sqrt( (np.sqrt(2)*rms_L/flux_2)**2 + (np.sqrt(2)*rms_L/flux_1)**2 ) / (np.log(nu3)-np.log(nu1))
		
		# plot alpha vs alpha
		log_flux = np.log10(flux_L['vals'])
		mask = np.isfinite(log_flux) & (log_flux>flux_lim) & np.isfinite(alpha_PL['vals']) & np.isfinite(alpha_LL['vals'])
		for x,y,xerr,yerr in np.array(zip(alpha_PL['vals'], alpha_LL['vals'], sigma_PL, sigma_LL))[mask]:
			ax0.errorbar(x, y, xerr=xerr, yerr=yerr, marker='.', color=color, zorder=6, alpha=0.8)

		# plot spatial data
		for sample in flux_L['samples'][mask]:
			ax1.scatter(*sample, color=color, s=2)

	ax0.set_xlabel('alpha_PL')
	ax0.set_ylabel('alpha_LL')
	ax0.axis('scaled')
	ax0.plot([-5, 5], [-5, 5], label='power law', c='black', lw=1, zorder=10)
	ax0.legend(loc='best')
	ax1.axis('scaled')
	
	xlim = ax1.get_xlim()
	ax1.set_xlim(xlim -pm* 0.05*(xlim[1]-xlim[0]))
	ylim = ax1.get_ylim()
	ax1.set_ylim(ylim -pm* 0.05*(ylim[1]-ylim[0]))
	
	ax1.imshow(np.log10(data), origin='lower', cmap='gray')
	fig.suptitle('%s color-color' % regfile.split('_')[0])
	fig.tight_layout(rect=[0.03, 0.03, 1, 0.97])

	# generate model color-color line
	model = 'JP'
	
	if model=='JP':
		x0, y0 = JP_plot(alpha_inj=-0.5)
		ax0.plot(x0, y0, label='JP', zorder=1)
		ax0.legend(loc='best')
	
	return fig, ax0, ax1

def JP(nu, nu_break, alpha_inj):
		return nu**alpha_inj * np.exp(-nu/nu_break)

def JP_plot(alpha_inj=-0.5):
	nu_break = np.arange(1,15)*1e9
	S1 = JP(nu1, nu_break, alpha_inj)
	S2 = JP(nu3, nu_break, alpha_inj)
	SL = JP(nuL, nu_break, alpha_inj)
	SP = JP(nuP, nu_break, alpha_inj)
	alpha_LL = get_alpha(S1, S2, nu1, nu3)
	alpha_PL = get_alpha(SP, SL, nuP, nuL)
	
	# fit line to model data
	p = np.polyfit(alpha_PL, alpha_LL, 1)
	x0 = np.linspace(-5, alpha_inj, 100)
	y0 = sum([c*x0**ix for ix,c in enumerate(np.flip(p))])
	
	return x0, y0

####################################
def KGJP():
	'''trying to get KGJP model to work'''

	def N(t_on, t_off, alpha_inj, nu, x, B, z):
	
		# time scales, Myr -> sec
		t_on_s = t_on * np.pi * 1e13
		t_off_s = t_off * np.pi * 1e13
	
		# power loss mechanism, B fields in Gauss
		B_G = B*1e-6
		B_CMB = 3.25e-6*(1+z)**2
		B_IC = np.sqrt(2./3) * B_CMB
		b = 2.37e-3 * B_G**2 * (2./3 + (B_IC/B_G)**2)
	
		# particle energy distribution, N(E) \propto E**-gamma
		gamma = 2*alpha_inj+1
		E = 4e-10 * np.sqrt(nu/x) / np.sqrt(B_G)
	
		if E < 1./(b*t_on_s):
			return E**(-gamma-1) / ( b*(gamma-1)*((1-b*E*(t_off_s-t_on_s))**(gamma-1) - (1-b*E*t_off_s)**(gamma-1)) )
		elif E > 1./(b*(t_off_s-t_on_s)):
			return 0.
		else:
			return E**(-gamma-1) / ( b*(gamma-1)*(1-b*E*(t_off_s-t_on_s))**(gamma-1) )

	def F(x):
		return x * integrate.quad(lambda z: special.kv(5./3, z), x, np.inf)[0]

	def S_KGJP(nu, t_on, t_off, alpha_inj, B, z):
		'''t_on is active phase duration (Myr)
		t_off is time elapsed since source shut-down (Myr)
		alpha_inj is injection spectral index
		B is the local magnetic field strength (uG)
		z is the redshift of the source'''
		return np.sqrt(nu) * integrate.quad(lambda x: F(x) * x**-1.5 * N(t_on, t_off, alpha_inj, nu, x, B, z), 0, np.inf)[0]

	# derived from Shulevski et al. (2015)
	model = 'KGJP'
	if model=='KGJP':
		dt = 33
		t_on = np.arange(10, 121, dt)
		t_off = np.arange(1, 2, dt)
		alpha_inj, B, z = -0.7, 7.4, 0.0806 # default B (in uG) from Keshet & Loeb 2010
		tons, toffs, alpha_LL, alpha_PL = [], [], [], []
	
		# for plotting spectra
		markers = np.array(['o', 's', '*', '+', 'x', 'v', '^', '>', '<', 'p', 'D'])
		labels = []
		plt.scatter(t_on, t_on, c=t_on)
		cb0 = plt.colorbar()
		plt.clf()
	
		print 't_on (Myr)\tt_off (Myr)'
		for t0 in t_on:
			for t1, m in zip(t_off, markers):
				print '%.0f\t\t%.0f' % (t0, t1)
				S1 = S_KGJP(nu1, t0, t1, alpha_inj, B, z)
				S2 = S_KGJP(nu3, t0, t1, alpha_inj, B, z)
				SL = S_KGJP(nuL, t0, t1, alpha_inj, B, z)
				SP = S_KGJP(nuP, t0, t1, alpha_inj, B, z)
			
				if m not in labels and not np.isnan(x) and not np.isnan(y):
					labels.append(m)
					label = True
				else:
					label = False
				plt.plot([nuP, nu1, nuL, nu3], [SP, S1, SL, S2], color=cb0.to_rgba(t0), marker=m, label=(t1 if label else None))
			
				alpha_LL.append(get_alpha(S1, S2, nu1, nu3))
				alpha_PL.append(get_alpha(SP, SL, nuP, nuL))
				tons.append(t0)
				toffs.append(t1)
	
		plt.xlabel('Frequency (Hz)')
		plt.ylabel('Intensity (arbitrary units)')
		plt.legend(loc='best', title='t_off (Myr)')
		sm = plt.cm.ScalarMappable(cmap=cb0.cmap, norm=mpl.colors.BoundaryNorm(np.hstack([t_on, t_on[-1]+dt]), cb0.cmap.N))
		sm.set_array([])
		cb1 = plt.colorbar(sm)
		cb1.set_label('t_on (Myr)')
		cb1.set_ticks(t_on + dt/2.)
		cb1.set_ticklabels(t_on)
		plt.tight_layout()

	# plot on color-color diagram
	alpha_LL = np.array(alpha_LL)
	alpha_PL = np.array(alpha_PL)
	tons = np.array(tons)
	toffs = np.array(toffs)
	plt.scatter(alpha_PL, alpha_LL, c=tons)
	cb0 = plt.colorbar()
	plt.clf()
	markers = np.array(['o', 's', '*', '+', 'x', 'v', '^', '>', '<', 'p', 'D'])

	labels = []
	for t0 in np.unique(tons):
		for t1, m in zip(np.unique(toffs), shapes):
			x = alpha_PL[np.logical_and(tons==t0, toffs==t1)]
			y = alpha_LL[np.logical_and(tons==t0, toffs==t1)]
			if m not in labels and not np.isnan(x) and not np.isnan(y):
				labels.append(m)
				label = True
			else:
				label = False
			plt.scatter(x, y, color=cb0.to_rgba(t0), marker=m, label=(t1 if label else None))

	# finesse the color bar
	sm = plt.cm.ScalarMappable(cmap=cb0.cmap, norm=mpl.colors.BoundaryNorm(np.hstack([t_on, t_on[-1]+dt]), cb0.cmap.N))
	sm.set_array([])
	cb1 = plt.colorbar(sm)
	cb1.set_label('t_on (Myr)')
	cb1.set_ticks(t_on + dt/2.)
	cb1.set_ticklabels(t_on)

	plt.xlabel('alpha_PL')
	plt.ylabel('alpha_LL')
	plt.axis('scaled')
	plt.tight_layout()
	plt.legend(loc='best', title='t_off (Myr)')

####################################
def error_approx(freq, beam, flux_0, alpha, plot=False):
	'''estimating the error for in-band alpha map via Monte Carlo
	fitsimage is corresponding flux map'''
	
	freq = freq.strip().upper()
	assert freq in ['L', 'P'], 'freq must be L or P'
	assert beam in [16, '16', 30, '30', '30x16', '180x16'], 'beam must be 16, 30, or 30x16'
	
	fitsimage = files[freq+str(beam)]
	rms = rms_dict[freq][str(beam)]
	
	if freq=='L':
		nu_0, nu_min, nu_max = nuL, 1.01e9, 2.03e9
	elif freq=='P':
		nu_0, nu_min, nu_max = nuP, 2.88e8, 4.48e8
	
	nu = np.linspace(nu_min, nu_max, 20, retstep=True, endpoint=False)
	nu = nu[0] + nu[1]/2.
	lognu = np.log10(nu)
	
	flux = flux_0*(nu/nu_0)**alpha
	logflux = np.log10(flux)

	alpha_obs = []
	n = []
	try:
		for i in range(1000):
			err = np.random.normal(0, rms*np.sqrt(len(flux)), len(flux))
			logfluxerr = np.log10(flux+err)
			mask = np.isfinite(logfluxerr)
			#n.append(sum(mask))
			if sum(mask):
				linfit = np.polyfit(lognu[mask], logfluxerr[mask], 1)
				alpha_obs.append(linfit[0])
				#plt.plot(lognu, linfit[1]+linfit[0]*lognu, lw=1)
	except TypeError:
		print flux, alpha, mask
		raise

	alpha_obs = np.array(alpha_obs)[np.isfinite(alpha_obs)]
	alpha_std = np.std(alpha_obs)
		
	if plot:
		fig, ax = plt.subplots(1)
		ax.hist(np.array(alpha_obs), bins=20)
		ax.axvline(alpha, color='k')
	
	return alpha, alpha_std

####################################
def solve_profile(y, x=None, n=None, slope_0=0., intercept_0=0., amplitude_1=1., mean_1=0., stddev_1=1., amplitude_2=1., mean_2=0., stddev_2=1., amplitude_3=1., mean_3=0., stddev_3=1., plot=True, scale=1, legend=None, linestyle=None):
	
	# Initialize models
	one_gauss = Linear1D + Gaussian1D
	two_gauss = Linear1D + Gaussian1D + Gaussian1D
	three_gauss = Linear1D + Gaussian1D + Gaussian1D + Gaussian1D
	fit = LevMarLSQFitter()
	
	# Fit model to data
	if x is None:
		x = range(len(y))
	x_ = np.linspace(x[0], x[-1], 1000)
	y = 1e3*y # convert units to mJy
	
	if n is None or n==1: # single Gaussian
		m_init1 = one_gauss(slope_0=slope_0, intercept_0=intercept_0, amplitude_1=amplitude_1, mean_1=mean_1, stddev_1=stddev_1)
		m1 = fit(m_init1, x, y)
		if plot:
			plt.plot(x_, scale*m1(x_), label=('%.1e'%sum((y-m1(x))**2) if legend is None else legend), linestyle=linestyle)
	
	if n is None or n==2: # double Gaussian
		m_init2 = two_gauss(slope_0=slope_0, intercept_0=intercept_0, amplitude_1=amplitude_1, mean_1=mean_1, stddev_1=stddev_1, amplitude_2=amplitude_2, mean_2=mean_2, stddev_2=stddev_2)
		m2 = fit(m_init2, x, y)
		if plot:
			plt.plot(x_, scale*m2(x_), label=('%.1e'%sum((y-m2(x))**2) if legend is None else legend), linestyle=linestyle)
	
	if n==3: # triple Gaussian
		m_init3 = three_gauss(slope_0=slope_0, intercept_0=intercept_0, amplitude_1=amplitude_1, mean_1=mean_1, stddev_1=stddev_1, amplitude_2=amplitude_2, mean_2=mean_2, stddev_2=stddev_2, amplitude_3=amplitude_3, mean_3=mean_3, stddev_3=stddev_3)
		m3 = fit(m_init3, x, y)
		if plot:
			plt.plot(x_, scale*m3(x_), label=('%.1e'%sum((y-m3(x))**2) if legend is None else legend), linestyle=linestyle)
	
	if plot:
		# Plot the data and the best fit
		plt.xlabel('Pixels along profile')
		plt.ylabel('Flux (mJy/beam)')
		plt.scatter(x, scale*y, s=5)
		plt.legend(loc='best', title=('Sum sqr err' if legend is None else None))
	
	# Return the prefered fit
	if n==1:
		return m1
	elif n==2:
		return m2
	elif n==3:
		return m3

def bg_subtract(region, n=None, slope_0=0., intercept_0=0., amplitude_1=1., mean_1=0., stddev_1=1., amplitude_2=1., mean_2=0., stddev_2=1., amplitude_3=1., mean_3=0., stddev_3=1., outfile=None, outprint=False, scale=1, legend=None, linestyle=None):

	# Inspect the profiles and fit Gaussian(s)
	if n is None:
		profile = solve_profile(region['vals'], scale=scale, legend=legend, linestyle=linestyle) # inspect plot
	else:
		profile = solve_profile(region['vals'], n=n, slope_0=slope_0, intercept_0=intercept_0, amplitude_1=amplitude_1, mean_1=mean_1, stddev_1=stddev_1, amplitude_2=amplitude_2, mean_2=mean_2, stddev_2=stddev_2, amplitude_3=amplitude_3, mean_3=mean_3, stddev_3=stddev_3, plot=operator.not_(outprint), scale=scale, legend=legend, linestyle=linestyle) # refined fit
	
	# Once fit, save coords and values of peaks to file
	if n is not None:
		if len(profile.parameters) >= 5:
			x1 = profile.mean_1.value
			y1 = profile[1](x1)
			if np.round(x1)>=0 and np.round(x1)<len(region['samples']):
				loc1 = region['samples'][int(np.round(x1))]
				if outprint:
					with open(outfile, 'a') as f:
						print >> f, loc1, y1/1e3
		
		if len(profile.parameters) >= 8:
			x2 = profile.mean_2.value
			y2 = profile[2](x2)
			if np.round(x2)>=0 and np.round(x2)<len(region['samples']):
				loc2 = region['samples'][int(np.round(x2))]
				if outprint:
					with open(outfile, 'a') as f:
						print >> f, loc2, y2/1e3
		
		if len(profile.parameters) >= 11:
			x3 = profile.mean_3.value
			y3 = profile[3](x3)
			if np.round(x3)>=0 and np.round(x3)<len(region['samples']):
				loc3 = region['samples'][int(np.round(x3))]
				if outprint:
					with open(outfile, 'a') as f:
						print >> f, loc3, y3/1e3
	
	return profile

def bg_subtract_manual():
	
	# NE relic
	regfile = 'NE_relic_30x16_profiles.reg'
	regions_L = rs.read_regions(regfile, files['L30x16'])
	outfile_L = 'NE_relic_L30x16_peaks.txt'
	regions_P = rs.read_regions(regfile, files['P30x16'])
	outfile_P = 'NE_relic_P30x16_peaks.txt'
	region = regions_L[0] # etc
	
	# TRG filaments
	regfile = 'TRG_filament_profiles.reg'
	regions_L = rs.read_regions(regfile, files['L16'])
	outfile_L = 'TRG_filament_L16_peaks.txt'
	regions_P = rs.read_regions(regfile, files['P16'])
	outfile_P = 'TRG_filament_P16_peaks.txt'
	region = regions_L[0] # etc
	
	# TRG tail
	regfile = 'TRG_16_ridge_profiles.reg'
	regions_L = rs.read_regions(regfile, files['L16'])
	outfile_L = 'TRG_ridge_L16_peaks.txt'
	regions_P = rs.read_regions(regfile, files['P16'])
	outfile_P = 'TRG_ridge_P16_peaks.txt'
	region = regions_L[0] # etc
	
	rs.bg_subtract(region, outfile=outfile_L, legend=None, outprint=F, n=None, mean_1=5, mean_2=10)
	
	profiles_L = ['TRG_filament_L16_peaks.txt', 'TRG_ridge_L16_peaks.txt']
	profiles_P = ['TRG_filament_P16_peaks.txt', 'TRG_ridge_P16_peaks.txt']
	beam = 16
	
	profiles_L = 'NE_relic_L30x16_peaks.txt'
	profiles_P = 'NE_relic_P30x16_peaks2.txt'
	beam = '30x16'

def bg_subtract_color_color(profiles_L, profiles_P, beam=16, label=None, hires=False):
	
	assert beam in [16, '16', 30, '30', '30x16', '180x16'], 'beam must be 16, 30, or 30x16'
	assert type(profiles_L) in [str, np.string_, list, np.ndarray] and type(profiles_P) in [str, np.string_, list, np.ndarray], 'profiles must be string, list of strings, or array of strings'
	
	# read map files
	im_file = files['L'+str(beam)]	
	alpha_file = files['LL'+str(beam)]

	with fits.open(im_file) as hdul:
		im_data = hdul[0].data[0][0]
		w = wcs.WCS(hdul[0])
	
	with fits.open(alpha_file) as hdul:
		alpha_data = hdul[0].data[0][0]
	
	if hires:
		with fits.open(files['L_hi-res']) as hdul:
			im_data = hdul[0].data[0][0]
			del hdul[0].header['ORIGIN']
			w = wcs.WCS(hdul[0])

	# read profile files
	if type(profiles_L) in [str, np.string_]:
		profiles_L = [profiles_L]

	if type(profiles_P) in [str, np.string_]:
		profiles_P = [profiles_P]
	
	# read flux values and coordinates along each profile
	profile_dict = {}
	for frequency, profiles in zip(['P', 'L'], [profiles_P, profiles_L]):
		for ix, (freq, profile) in enumerate(itertools.product(frequency, profiles)):
			
			with open(profile, 'r') as f:
				lines = f.readlines()
			
			xx, yy, vals = [], [], []
			for line in lines:
				xx.append(float(line.split(' ')[0][1:]))
				yy.append(float(line.split(' ')[1][:-1]))
				vals.append(float(line.split(' ')[2].strip()))
			
			coords = np.array(zip(xx, yy))
			profile_dict['%s%i'%(freq,ix)] = {'vals':np.array(vals), 'coords':coords, 'filename':profile}
	
	# identify matching coordinates
	for ix in range(len(profile_dict)/2):
		ix_L, ix_P = [], []
		for i, val1 in enumerate(profile_dict['L%i'%(ix)]['coords']):
			for j, val2 in enumerate(profile_dict['P%i'%(ix)]['coords']):
				if val1[0]==val2[0] and val1[1]==val2[1]:
					ix_L.append(i)
					ix_P.append(j)
		profile_dict['L%i'%(ix)]['ix'] = np.array(ix_L)
		profile_dict['P%i'%(ix)]['ix'] = np.array(ix_P)
	
	# apply manual masking and color-coding
	region_dict = {'TRG_filament_P16_peaks.txt':np.array([2, 2, 3, -1, 2, 1, -1, -1]), 
				   'TRG_ridge_P16_peaks.txt':np.zeros(7).astype(int), 
				   'TRG_filament_L16_peaks.txt':np.array([2, 2, 3, -1, 2, 1, -1, -1]), 
				   'TRG_ridge_L16_peaks.txt':np.zeros(7).astype(int),
				   'NE_relic_L30x16_peaks.txt':np.array(3*[0]+[1, -1]+13*[2]+3*[-1]),
				   'NE_relic_P30x16_peaks2.txt':np.array(3*[0]+[1, -1]+13*[2]+3*[-1]),
				   'NE_relic_L180x16_peaks.txt':np.array([0, 1, 2]),
				   'NE_relic_P180x16_peaks2.txt':np.array([0, 1, 2]),
				   'Goldfish_L16_peaks.txt':np.array([0]*11+[1]*4),
				   'Goldfish_P16_peaks.txt':np.array([1]*4+[0]*11)}
	for region, profile in itertools.product(region_dict, profile_dict):
		if profile_dict[profile]['filename'] == region:
			profile_dict[profile]['ix'] = profile_dict[profile]['ix'][region_dict[region]>=0]
			profile_dict[profile]['regions'] = region_dict[region][region_dict[region]>=0]
			for key in ['coords', 'vals']:
				profile_dict[profile][key] = profile_dict[profile][key][profile_dict[profile]['ix']]
	
	# initialize figure
	source = 'trg' if profile_dict['L0']['filename'][:3]=='TRG' else (
	         'ne_relic' if profile_dict['L0']['filename'][:2]=='NE' else 
			 'goldfish' )
	if hires:
		_source = source
		source += '_hi-res'
	fig, (ax0, ax1) = plt.subplots(1,2)

	for i in range(len(profile_dict)/2):
		for (xx, yy), c in zip(profile_dict['L%i'%i]['coords'], profile_dict['L%i'%i]['regions']):
			if not hires:
				ax1.scatter(xx-axes_dict[source]['x1']+1, yy-axes_dict[source]['y1']+1, s=5, color='C%i'%c)
			else:
				ratio = 1. * axes_dict[source]['x1'] / axes_dict[_source]['x1']
				ax1.scatter(ratio*xx-axes_dict[source]['x1']+1, ratio*yy-axes_dict[source]['y1']+1, s=5, color='C%i'%c)
	
	# calculate error bars 
	rms_L, rms_P = rms_dict['L'][str(beam)], rms_dict['P'][str(beam)]
	for i in range(len(profile_dict)/2):
		flux_P, flux_L = profile_dict['P%i'%i]['vals'], profile_dict['L%i'%i]['vals']
		alpha_PL = get_alpha(flux_P, flux_L, nuP, nuL)
		alpha_LL = []
		for (x, y) in profile_dict['L%i'%i]['coords']:
			alpha_LL.append(alpha_data[int(y),int(x)])
		
		alpha_LL = np.array(alpha_LL)
		flux_1 = flux_L * (nu1/nuL)**alpha_LL
		flux_2 = flux_L * (nu3/nuL)**alpha_LL
		sigma_PL = np.sqrt( (rms_L/flux_L)**2 + (rms_P/flux_P)**2 ) / (np.log(nuL)-np.log(nuP))
		#sigma_LL = np.sqrt( (np.sqrt(2)*rms_L/flux_2)**2 + (np.sqrt(2)*rms_L/flux_1)**2 ) / (np.log(nu3)-np.log(nu1))
		sigma_LL = []
		for flux, alpha in zip(flux_L, alpha_LL):
			sigma_LL.append(error_approx('L', beam, flux, alpha)[1])
		
		sigma_LL = np.array(sigma_LL)
		
		# plot alpha vs alpha
		log_flux = np.log10(flux_L)
		mask = np.isfinite(log_flux) & np.isfinite(alpha_PL) & np.isfinite(alpha_LL)
		for x,y,xerr,yerr,color in np.array(zip(alpha_PL, alpha_LL, sigma_PL, sigma_LL, profile_dict['L%i'%i]['regions']))[mask]:
			ax0.errorbar(x, y, xerr=xerr, yerr=yerr, marker=None, color='C%i'%int(color), zorder=6, alpha=0.5)
			ax0.scatter(x, y, marker='o', color='C%i'%int(color), zorder=6, s=6)
	
	# prettify plot
	ax0.plot([-5, 5], [-5, 5], label='Power law', c='black', lw=1, zorder=10)
	alpha_inj = -0.5
	x0, y0 = JP_plot(alpha_inj)
	ax0.plot(x0, y0, label='JP, $\\alpha_{inj}=%.1f$'%alpha_inj, zorder=1, linestyle='--')
	
	origin_wcs = w.wcs_pix2world(axes_dict[source]['x1'], axes_dict[source]['y1'], 0, 0, 1)[0:2]
	im_crop = im_data[axes_dict[source]['y1']-1:axes_dict[source]['y2'], axes_dict[source]['x1']-1:axes_dict[source]['x2']]
	im = ax1.imshow(np.log10(im_crop), origin='lower', cmap='gray')	
	#cb = fig.colorbar(im)
	#cb.set_label('log10(Lband flux [mJy])')
	
	# build x-axis
	ra_degs = ra_to_deg(axes_dict[source]['ra_ticks'])
	x_ticks = w.wcs_world2pix(ra_degs, origin_wcs[1], 0, 0, 1)[0] - axes_dict[source]['x1']
	ax1.set_xticks(x_ticks)
	ax1.set_xticklabels(axes_dict[source]['ra_labels'])
	ax1.xaxis.set_minor_locator(AutoMinorLocator(axes_dict[source]['ticks_x']))

	# build y-axis
	dec_degs = dec_to_deg(axes_dict[source]['dec_ticks'])
	y_ticks = w.wcs_world2pix(origin_wcs[0], dec_degs, 0, 0, 1)[1] - axes_dict[source]['y1']
	ax1.set_yticks(y_ticks)
	ax1.set_yticklabels(axes_dict[source]['dec_labels'])
	ax1.yaxis.set_minor_locator(AutoMinorLocator(axes_dict[source]['ticks_y']))
	
	ax0.axis('scaled')
	ax1.axis('scaled')
	ax0.set_xlabel('$\\alpha_{PL}$')
	ax0.set_ylabel('$\\alpha_{LL}$')
	ax0.legend(loc='best')
	#fig.suptitle('%s color-color diagram' % (source))
	fig.tight_layout() #(rect=[0.03, 0.03, 1, 0.97])

####################################

def tdiff(t0, t1):
	'''find difference in seconds in two times written as HH:MM:SS'''
	
	h0, m0, s0 = np.array(t0.split(':')).astype(int)
	h1, m1, s1 = np.array(t1.split(':')).astype(int)
	t0_sec = s0+m0*60+h0*3600
	t1_sec = s1+m1*60+h1*3600
	return np.abs(t1_sec-t0_sec)

def get_obs_times(filename):
	'''read a listr file from AIPS and calculate the total time spent observing each source (in seconds)'''
	
	with open(filename) as f:
		lines = f.readlines()
	
	sources = {}
	for line in lines:
		if line=='Source summary\n':
			break
		
		if unicode(line[:4].strip()).isnumeric():
			source = line[5:22].strip()
			t0 = line[43:51]
			t1 = line[58:66]
			if source not in sources:
				sources[source] = tdiff(t0, t1)
			else:
				sources[source] += tdiff(t0, t1)
	
	return sources
	
def ra_to_hhmmss(ra):
	'''converts RA value from degrees to HH:MM:SS'''
	
	try:
		ra = float(ra)
		h = int(ra/15.)
		m = int((ra-h*15)*4)
		s = int(round((ra-h*15-m/4.)*240))
		if m==60:
			h += 1
			m = 0
		if s==60:
			m += 1
			s = 0
		return '%02i:%02i:%02i' % (h, m, s)
	
	except:
		if type(ra) in [list, np.ndarray]:
			ras = []
			for ra0 in ra:
				ras.append(ra_to_hhmmss(ra0))
			return np.array(ras)
		else:
			raise TypeError, 'ra must be number, list, or array'

def ra_to_deg(ra):
	'''converts RA value from HH:MM:SS to degrees'''
	
	try:
		h, m, s = ra.split(':')
		return int(h)*15+int(m)/4.+int(s)/240.
	
	except:
		if type(ra) in [list, np.ndarray]:
			ras = []
			for ra0 in ra:
				ras.append(ra_to_deg(ra0))
			return np.array(ras)
		else:
			raise TypeError, 'ra must be "HH:MM:SS", list, or array'

def dec_to_ddmmss(dec):
	'''converts Dec value from degrees to DD:MM:SS'''
	
	try:
		sign = np.abs(dec)/dec
		dec = np.abs(dec)
		d = int(dec)
		m = int((dec - d)*60)
		s = int(round((dec - d - m/60.)*3600))
		if m==60:
			h += 1
			m = 0
		if s==60:
			m += 1
			s = 0
		return '%+02i:%02i:%02i' % (sign*d, m, s)
	
	except:
		if type(dec) in [list, np.ndarray]:
			decs = []
			for dec0 in dec:
				decs.append(dec_to_hhmmss(dec0))
			return np.array(decs)
		else:
			raise TypeError, 'dec must be number, list, or array'

def dec_to_deg(dec):
	'''converts Dec value from DD:MM:SS to degrees'''
	
	try:
		sign = -1. if dec[0]=='-' else 1.
		d, m, s = dec.split(':')
		return sign * (np.abs(int(d))+int(m)/60.+int(s)/3600.)
	
	except:
		if type(dec) in [list, np.ndarray]:
			decs = []
			for dec0 in dec:
				decs.append(dec_to_deg(dec0))
			return np.array(decs)
		else:
			raise TypeError, 'dec must be "DD:MM:SS", list, or array'

def print_map(source, freq='P', beam=16):
	'''Create a frame for printing maps of sources within A2255
	source options are TRG, Goldfish, Beaver, Embryo, NE_relic, full (for the entire field), and central (for the Goldfish+TRG)
	freq and beam set the example map plotted'''
	
	freq = freq.strip().upper()
	assert freq in ['L', 'P', 'LL', 'PL'], 'freq must be L, P, LL, or PL'
	assert beam in [16, '16', 30, '30', '30x16', '180x16'], 'beam must be 16, 30, or 30x16'
	source = source.lower()
	
	# load fits image
	fitsimage = files[freq+str(beam)]
	with fits.open(fitsimage) as hdul:
		data = hdul[0].data[0][0]
		w = wcs.WCS(hdul[0])
	
	# prepare subplot
	crop = data[axes_dict[source]['y1']-1:axes_dict[source]['y2'], axes_dict[source]['x1']-1:axes_dict[source]['x2']]
	origin_wcs = w.wcs_pix2world(axes_dict[source]['x1'], axes_dict[source]['y1'], 0, 0, 1)[0:2]
	fig, ax = plt.subplots(1)
	im = ax.imshow(np.log10(crop), origin='lower', cmap='gray')
	ax.set_xlabel('Right Ascension (J2000)')
	ax.set_ylabel('Declination (J2000)')
	
	# build x-axis
	ra_degs = ra_to_deg(axes_dict[source]['ra_ticks'])
	x_ticks = w.wcs_world2pix(ra_degs, origin_wcs[1], 0, 0, 1)[0] - axes_dict[source]['x1']
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(axes_dict[source]['ra_labels'])
	ax.xaxis.set_minor_locator(AutoMinorLocator(axes_dict[source]['ticks_x']))

	# build y-axis
	dec_degs = dec_to_deg(axes_dict[source]['dec_ticks'])
	y_ticks = w.wcs_world2pix(origin_wcs[0], dec_degs, 0, 0, 1)[1] - axes_dict[source]['y1']
	ax.set_yticks(y_ticks)
	ax.set_yticklabels(axes_dict[source]['dec_labels'])
	ax.yaxis.set_minor_locator(AutoMinorLocator(axes_dict[source]['ticks_y']))

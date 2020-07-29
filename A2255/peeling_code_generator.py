import glob, os, sys
if os.getcwd().split('/')[-1] != 'peeling':
	sys.exit('Not in peeling directory; exiting')

for submit in glob.glob('submit_jobs*'):
	os.remove(submit)

for obs_n in [1,2]:
	times = ['09:34:00~10:46:00', '10:46:00~11:58:00', '11:58:00~13:10:00', '14:12:00~15:10:00'] if obs_n==1 else ['10:30:00~11:50:00', '11:50:00~13:10:00', '14:20:00~15:10:00', '15:10:00~16:00:00']
	for ix, timerange in enumerate(times):
		for stokes in ['XX', 'YY']:
			
			# Add jobs to bash script
			
			for label in ['A', 'BC', 'D', 'E', 'F', 'G']:
				with open('submit_jobs_%s.sh' % label, 'a') as f:
					print >> f, 'submit -f ~/CASA/peeling/batch_obs%i_%s_%i_%s.req' % (obs_n, stokes, ix, label)
			
			# Create job files
			
			for label in ['A', 'BC', 'D', 'E', 'F', 'G']:
				with open('batch_obs%i_%s_%i_%s.req' % (obs_n, stokes, ix, label), 'w') as f:
					print >> f, '''# These are required
VERS="%i_%s_%i"
PMEM="27gb"   # physmem used by process. Wont kill job exceeded
WORK_DIR="/lustre/aoc/observers/nm-10124/CASA"
CASAPATH="/home/casa/packages/RHEL7/release/current"
COMMAND="xvfb-run -d /opt/local/bin/mpicasa ${CASAPATH}/bin/casa --nogui -c /lustre/aoc/observers/nm-10124/CASA/peeling/params_obs${VERS}_%s.py"

# These are optional
NUM_NODES="1"     # default is 1
NUM_CORES="7"     # default is 1
QUEUE="batch"     # default is the batch queue
STDOUT="peeling/batch_obs${VERS}.out"   # file relative to WORK_DIR.  default is no output
STDERR="peeling/batch_obs${VERS}_err.out"   # file relative to WORK_DIR.  default is no output
MAIL_OPTIONS="abe"   # default is "n" therefore no email
JOB_NAME="A2255_P${VERS}_%s"   # default is <username>_qsub.<procID>''' % (obs_n, stokes, ix, label, label)
			
			# Create CASA scripts
			
			# Step 1: rotate, clean, and selfcal
			with open('params_obs%i_%s_%i_A.py' % (obs_n, stokes, ix), 'w') as f:
				print >> f, '''# Step 1: rotate, clean, and selfcal

import glob, os, shutil

# uv-peeling steps (for each polarization and time step)
obs_n = %i
ix, timerange = %i, '%s'
stokes = '%s'
imagename = 'output/A2255_obs{0}_peeled_{1}_{2}'.format(obs_n, stokes, ix)
vis0='A2255_obs{0}_{1}_{2}.ms'.format(obs_n, stokes, ix)

# Reset the corrected and model data columns in the MS
clearcal(vis=vis0)

# Flag additional data from visual inspection
flagdata(vis=vis0, spw='3:102~105')

flagdata(vis=vis0,
	mode='clip',
	antenna='4',
	clipminmax=[0,10],
	channelavg=True,
	chanbin=4)

if obs_n==2 and ix>0:
	flagdata(vis=vis0, spw='8')

# Delete existing maps
old_maps = glob.glob(imagename+'_A.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Rotate phase center to bright source and generate model for self-cal
tclean(vis=vis0,
    datacolumn='corrected',
    imagename=imagename+'_A',
    timerange=timerange,
    stokes=stokes,
    imsize=4800,
    cell='4.0arcsec',
    phasecenter='J2000 260.000deg +64.077deg', # peak of bright source
    interpolation='linear',
    gridder='widefield',
    nterms=2,
    wprojplanes=-1,
    aterm=False,
    psterm=True,
    wbawp=False,
    conjbeams=False,
    computepastep=360.0,
    rotatepastep=10.0,
    pblimit=0.01,
    deconvolver='mtmfs',
    scales=[0,3,9],
    weighting='briggs',
    robust=0,
    cyclefactor=4,
    niter=20000, # just generating an initial model
    mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_4.mask',
    savemodel='modelcolumn',
    calcres=True,
    calcpsf=True,
    parallel=True)

# Perform gain calibration
caltable = 'peeling/pcal_obs{}_{}_{}_A.tab'.format(obs_n, stokes, ix)
rmtables(caltable)
gaincal(vis=vis0,
    caltable=caltable,
    spw='1~4,7~10', # some spws are fully flagged already
    timerange=timerange,
    gaintype='T',
    refant='11',
    calmode='p',
    solint='1min',
    minsnr=3.0,
    minblperant=6)''' % (obs_n, ix, timerange, stokes)
			
			# Step 2: model bright source and subtract it
			# Step 3: rotate back, clean, and selfcal (or crosscal)
			with open('params_obs%i_%s_%i_BC.py' % (obs_n, stokes, ix), 'w') as f:
				print >> f, '''# Step 2: model bright source and subtract it
# Step 3: rotate back, clean, and selfcal (or crosscal)

import glob, os, shutil

# uv-peeling steps (for each polarization and time step)
obs_n = %i
ix, timerange = %i, '%s'
stokes = '%s'
imagename = 'output/A2255_obs{0}_peeled_{1}_{2}'.format(obs_n, stokes, ix)
vis0='A2255_obs{0}_{1}_{2}.ms'.format(obs_n, stokes, ix)
vis1='A2255_obs{0}_{1}_{2}_peeled.ms'.format(obs_n, stokes, ix)

# Apply calibration
caltable = 'peeling/pcal_obs{}_{}_{}_A.tab'.format(obs_n, stokes, ix)
applycal(vis=vis0,
    gaintable=[caltable],
    gainfield='',
    calwt=False,
    flagbackup=True,
    interp='linear')

# Delete existing maps
old_maps = glob.glob(imagename+'_B.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Clean just the bright source
tclean(vis=vis0,
    datacolumn='corrected',
    imagename=imagename+'_B',
    timerange=timerange,
    stokes=stokes,
    imsize=4800,
    cell='4.0arcsec',
    phasecenter='J2000 260.000deg +64.077deg', # peak of bright source
    interpolation='linear',
    gridder='widefield',
    nterms=2,
    wprojplanes=-1,
    aterm=False,
    psterm=True,
    wbawp=False,
    conjbeams=False,
    computepastep=360.0,
    rotatepastep=10.0,
    pblimit=0.01,
    deconvolver='mtmfs',
    scales=[0,3,9],
    weighting='briggs',
    robust=0,
    cyclefactor=4,
    niter=5000, # it won't converge on its own since only one source is being cleaned 
    mask='circle[[259.998deg,+64.077deg],30arcsec]', # clean box around the bright source
    savemodel='modelcolumn',
    calcres=True,
    calcpsf=True,
    parallel=True)

# Subtract the model of the bright source from the uv data
uvsub(vis=vis0)

# Initialize new ms
if os.path.exists(vis1):
	shutil.rmtree(vis1)

if os.path.exists(vis1+'.flagversions'):
	shutil.rmtree(vis1+'.flagversions')

# Split the peeled data into the new ms
split(vis=vis0,
	outputvis=vis1,
	datacolumn='corrected',
	keepflags=False)

# Delete existing maps
old_maps = glob.glob(imagename+'_C.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Return to pointing center and clean
tclean(vis=vis1,
    datacolumn='corrected',
    imagename=imagename+'_C',
    timerange=timerange,
    stokes=stokes,
    imsize=4800,
    cell='4.0arcsec',
    phasecenter='J2000 17h12m42.0 +64d09m44.0', # original pointing center
    interpolation='linear',
    gridder='widefield',
    nterms=2,
    wprojplanes=-1,
    aterm=False,
    psterm=True,
    wbawp=False,
    conjbeams=False,
    computepastep=360.0,
    rotatepastep=10.0,
    pblimit=0.01,
    deconvolver='mtmfs',
    scales=[0,3,9,27],
    weighting='briggs',
    robust=0,
    cyclefactor=4,
    niter=50000, # deep clean now
    mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_4.mask',
    savemodel='modelcolumn',
    calcres=True,
    calcpsf=True,
    parallel=True)

# Copy day 1 model over for crosscal if necessary
do_crosscal = obs_n == 2
if do_crosscal:
	old_maps = glob.glob(imagename+'_temp.*')
	for f in old_maps:
		if os.path.isdir(f):
			shutil.rmtree(f)
		else:
			os.remove(f)
	
	tclean(vis=vis1,
		datacolumn='corrected',
		imagename=imagename+'_temp',
		timerange=timerange,
		stokes=stokes,
		imsize=4800,
		cell='4.0arcsec',
		phasecenter='J2000 17h12m42.0 +64d09m44.0', # original pointing center
		interpolation='linear',
		gridder='widefield',
		nterms=2,
		wprojplanes=-1,
		aterm=False,
		psterm=True,
		wbawp=False,
		conjbeams=False,
		computepastep=360.0,
		rotatepastep=10.0,
		pblimit=0.01,
		deconvolver='mtmfs',
		scales=[0,3,9,27],
		weighting='briggs',
		robust=0,
		cyclefactor=4,
		niter=0, # no actual cleaning
		mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_4.mask',
		startmodel=['output/A2255_obs1_peeled_boxed_2.model.tt0', 'output/A2255_obs1_peeled_boxed_2.model.tt1'], # good day 1 model
		savemodel='modelcolumn',
		calcres=True,
		calcpsf=True,
		parallel=True)
	
	old_maps = glob.glob(imagename+'_temp.*')
	for f in old_maps:
		if os.path.isdir(f):
			shutil.rmtree(f)
		else:
			os.remove(f)

# Perform self-cal on the peeled data
caltable = 'peeling/pcal_obs{}_{}_{}_C.tab'.format(obs_n, stokes, ix)
rmtables(caltable)
gaincal(vis=vis1,
    caltable=caltable,
    spw='1~4,7~10', # some spws are fully flagged already
    timerange=timerange,
    gaintype='T',
    refant='11',
    calmode='p',
    solint='1min',
    minsnr=3.0,
    minblperant=6)''' % (obs_n, ix, timerange, stokes)
			
			# Step 4: clean again, test phase selfcal (should be ~0) and amp selfcal (should be ~1)
			with open('params_obs%i_%s_%i_D.py' % (obs_n, stokes, ix), 'w') as f:
				print >> f, '''# Step 4: clean again, test phase selfcal (should be ~0) and amp selfcal (should be ~1)

import glob, os, shutil

# uv-peeling steps (for each polarization and time step)
obs_n = %i
ix, timerange = %i, '%s'
stokes = '%s'
imagename = 'output/A2255_obs{0}_peeled_{1}_{2}'.format(obs_n, stokes, ix)
vis1='A2255_obs{0}_{1}_{2}_peeled.ms'.format(obs_n, stokes, ix)
vis2='A2255_obs{0}_{1}_{2}_selfcal.ms'.format(obs_n, stokes, ix)

# Apply new calibration
caltable = 'peeling/pcal_obs{}_{}_{}_C.tab'.format(obs_n, stokes, ix)
applycal(vis=vis1,
    gaintable=[caltable],
    gainfield='',
    calwt=False,
    flagbackup=True,
    interp='linear')

# Initialize new ms for second round of selfcal
if os.path.exists(vis2):
	shutil.rmtree(vis2)

#if os.path.exists(vis2+'.flagversions'):
	shutil.rmtree(vis2+'.flagversions')

# Split the selfcal'ed data into the new ms
split(vis=vis1,
	outputvis=vis2,
	datacolumn='corrected',
	keepflags=False)

# Delete existing maps
old_maps = glob.glob(imagename+'_D.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Clean again to check
tclean(vis=vis2,
    datacolumn='corrected',
    imagename=imagename+'_D',
    timerange=timerange,
    stokes=stokes,
    imsize=4800,
    cell='4.0arcsec',
    phasecenter='J2000 17h12m42.0 +64d09m44.0', # original pointing center
    interpolation='linear',
    gridder='widefield',
    nterms=2,
    wprojplanes=-1,
    aterm=False,
    psterm=True,
    wbawp=False,
    conjbeams=False,
    computepastep=360.0,
    rotatepastep=10.0,
    pblimit=0.01,
    deconvolver='mtmfs',
    scales=[0,3,9,27],
    weighting='briggs',
    robust=0,
    cyclefactor=4,
    niter=50000, # deep clean
    mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_4.mask',
    savemodel='modelcolumn',
    calcres=True,
    calcpsf=True,
    parallel=True)

# Perform self-cal again; solutions should all be around zero
caltable = 'peeling/pcal_obs{}_{}_{}_D.tab'.format(obs_n, stokes, ix)
rmtables(caltable)
gaincal(vis=vis2,
    caltable=caltable,
    spw='1~4,7~10', # some spws are fully flagged already
    timerange=timerange,
    gaintype='T',
    refant='11',
    calmode='p',
    solint='1min',
    minsnr=3.0,
    minblperant=6)

# Also perform amplitude self-cal; if phase self-cal has converged, see what the effects of amp cal are
caltable = 'peeling/acal_obs{}_{}_{}_D.tab'.format(obs_n, stokes, ix)
rmtables(caltable)
gaincal(vis=vis2,
    caltable=caltable,
    spw='1~4,7~10', # some spws are fully flagged already
    timerange=timerange,
    gaintype='T',
    refant='11',
    calmode='ap',
    solint='1min',
    minsnr=3.0,
    minblperant=6,
    solnorm=True)''' % (obs_n, ix, timerange, stokes)
			
			# Step 4b: if selfcal hadn't converged, apply new solutions; this step shouldn't be necessary
			with open('params_obs%i_%s_%i_E.py' % (obs_n, stokes, ix), 'w') as f:
				print >> f, '''# Step 4b: if selfcal hadn't converged, apply new solutions; this step shouldn't be necessary

import glob, os, shutil

# uv-peeling steps (for each polarization and time step)
obs_n = %i
ix, timerange = %i, '%s'
stokes = '%s'
imagename = 'output/A2255_obs{0}_peeled_{1}_{2}'.format(obs_n, stokes, ix)
vis2='A2255_obs{0}_{1}_{2}_selfcal.ms'.format(obs_n, stokes, ix)
vis3='A2255_obs{0}_{1}_{2}_ampcal.ms'.format(obs_n, stokes, ix)

# Initialize new ms for experimental round of ampcal
if os.path.exists(vis3):
	shutil.rmtree(vis3)

if os.path.exists(vis3+'.flagversions'):
	shutil.rmtree(vis3+'.flagversions')

# Split the selfcal'ed data into the new ms
split(vis=vis2,
	outputvis=vis3,
	datacolumn='data',
	keepflags=False)

# Apply new calibration
caltable = 'peeling/acal_obs{}_{}_{}_D.tab'.format(obs_n, stokes, ix)
applycal(vis=vis3,
    gaintable=[caltable],
    gainfield='',
    calwt=False,
    flagbackup=True,
    interp='linear')

# Delete existing maps
old_maps = glob.glob(imagename+'_E.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Clean with ampcal
tclean(vis=vis3,
    datacolumn='corrected',
    imagename=imagename+'_E',
    timerange=timerange,
    stokes=stokes,
    imsize=4800,
    cell='4.0arcsec',
    phasecenter='J2000 17h12m42.0 +64d09m44.0', # original pointing center
    interpolation='linear',
    gridder='widefield',
    nterms=2,
    wprojplanes=-1,
    aterm=False,
    psterm=True,
    wbawp=False,
    conjbeams=False,
    computepastep=360.0,
    rotatepastep=10.0,
    pblimit=0.01,
    deconvolver='mtmfs',
    scales=[0,3,9,27],
    weighting='briggs',
    robust=0,
    cyclefactor=4,
    niter=50000, # deep clean
    mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_4.mask',
    savemodel='modelcolumn',
    calcres=True,
    calcpsf=True,
    parallel=True)''' % (obs_n, ix, timerange, stokes)
			
			# Step 5: combine all ms's and tclean together, then manually update clean mask (see ../tclean_params.py)
			
			# Step 6: copy combined model into binned ms's and amp selfcal
			with open('params_obs%i_%s_%i_F.py' % (obs_n, stokes, ix), 'w') as f:
				print >> f, '''# Step 6: copy combined model into binned ms's and selfcal

import glob, os, shutil

# uv-peeling steps (for each polarization and time step)
obs_n = %i
ix, timerange = %i, '%s'
stokes = '%s'
imagename = 'output/A2255_obs{0}_peeled_{1}_{2}'.format(obs_n, stokes, ix)
vis2='A2255_obs{0}_{1}_{2}_selfcal.ms'.format(obs_n, stokes, ix)
vis3='A2255_obs{0}_{1}_{2}_combined_crosscal.ms'.format(obs_n, stokes, ix)

# Initialize new ms for experimental round of ampcal
if os.path.exists(vis3):
	shutil.rmtree(vis3)

if os.path.exists(vis3+'.flagversions'):
	shutil.rmtree(vis3+'.flagversions')

# Split the selfcal'ed data into the new ms
split(vis=vis2,
	outputvis=vis3,
	datacolumn='data',
	keepflags=False)

# Copy combined model over for crosscal
old_maps = glob.glob(imagename+'_temp.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

tclean(vis=vis3,
	datacolumn='corrected',
	imagename=imagename+'_temp',
	timerange=timerange,
	stokes=stokes,
	imsize=4800,
	cell='4.0arcsec',
	phasecenter='J2000 17h12m42.0 +64d09m44.0', # original pointing center
	interpolation='linear',
	gridder='widefield',
	nterms=2,
	wprojplanes=-1,
	aterm=False,
	psterm=True,
	wbawp=False,
	conjbeams=False,
	computepastep=360.0,
	rotatepastep=10.0,
	pblimit=0.01,
	deconvolver='mtmfs',
	scales=[0,3,9,27],
	weighting='briggs',
	robust=0,
	cyclefactor=4,
	niter=0, # no actual cleaning
	mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_5.mask', # deeper combined mask
	startmodel=['output/A2255_combined_selfcal.model.tt0', 'output/A2255_combined_selfcal.model.tt1'], # good combined model
	savemodel='modelcolumn',
	calcres=True,
	calcpsf=True,
	parallel=True)

old_maps = glob.glob(imagename+'_temp.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Perform phase cross-cal
caltable = 'peeling/pcal_obs{}_{}_{}_F.tab'.format(obs_n, stokes, ix)
rmtables(caltable)
gaincal(vis=vis3,
    caltable=caltable,
    spw='1~4,7~10', # some spws are fully flagged already
    timerange=timerange,
    gaintype='T',
    refant='11',
    calmode='p',
    solint='1min',
    minsnr=3.0,
    minblperant=6)

# Perform amplitude cross-cal
caltable = 'peeling/acal_obs{}_{}_{}_F.tab'.format(obs_n, stokes, ix)
rmtables(caltable)
gaincal(vis=vis3,
    caltable=caltable,
    spw='1~4,7~10', # some spws are fully flagged already
    timerange=timerange,
    gaintype='T',
    refant='11',
    calmode='ap',
    solint='1min',
    minsnr=3.0,
    minblperant=6,
    solnorm=True)''' % (obs_n, ix, timerange, stokes)

			# Step 7: apply calibrations if needed and tclean again
			with open('params_obs%i_%s_%i_G.py' % (obs_n, stokes, ix), 'w') as f:
				print >> f, '''# Step 7: apply calibrations if needed

import glob, os, shutil

# uv-peeling steps (for each polarization and time step)
obs_n = %i
ix, timerange = %i, '%s'
stokes = '%s'
imagename = 'output/A2255_obs{0}_peeled_{1}_{2}'.format(obs_n, stokes, ix)
vis3='A2255_obs{0}_{1}_{2}_combined_crosscal.ms'.format(obs_n, stokes, ix)

# Apply new calibration
clearcal(vis=vis3)
caltable = 'peeling/acal_obs{}_{}_{}_F.tab'.format(obs_n, stokes, ix)
applycal(vis=vis3,
    gaintable=[caltable],
    gainfield='',
    calwt=False,
    flagbackup=True,
    interp='linear')

# Delete existing maps
old_maps = glob.glob(imagename+'_G.*')
for f in old_maps:
	if os.path.isdir(f):
		shutil.rmtree(f)
	else:
		os.remove(f)

# Clean again to check
tclean(vis=vis3,
    datacolumn='corrected',
    imagename=imagename+'_G',
    timerange=timerange,
    stokes=stokes,
    imsize=4800,
    cell='4.0arcsec',
    phasecenter='J2000 17h12m42.0 +64d09m44.0', # original pointing center
    interpolation='linear',
    gridder='widefield',
    nterms=2,
    wprojplanes=-1,
    aterm=False,
    psterm=True,
    wbawp=False,
    conjbeams=False,
    computepastep=360.0,
    rotatepastep=10.0,
    pblimit=0.01,
    deconvolver='mtmfs',
    scales=[0,3,9,27],
    weighting='briggs',
    robust=0,
    cyclefactor=4,
    niter=50000, # deep clean
	mask='/lustre/aoc/observers/nm-10124/CASA/output/clean_boxes_5.mask', # deeper combined mask
    savemodel='modelcolumn',
    calcres=True,
    calcpsf=True,
    parallel=True)''' % (obs_n, ix, timerange, stokes)

# Make files executable

for submit in glob.glob('submit_jobs*.sh'):
	os.chmod(submit, 0o755)

for params in glob.glob('params*.py'):
	os.chmod(params, 0o755)

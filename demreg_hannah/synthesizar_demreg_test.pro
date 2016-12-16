pro synthesizar_demreg_test
  ;~~~~~~~~
  ;Read in all channel data from FITS files into the array aia_cube for a given timestamp
  time_stamp='000900'
  aia_files = findfile('/data/datadrive2/ar_viz/synthesizar_v01demo/SDO_AIA/*/map_t'+time_stamp+'.fits')
  aia_files = shift(aia_files,1)
  read_sdo,aia_files,ind,data
  
  ;A Few tests to make sure we're reading things in correctly
  ;aia_lct, rr, gg, bb, wavelnth=94, /load
  ;plot_image,alog(aia_cube[*,*,0]),min=alog(.1),max=alog(50)
 
  ;~~~~~~~~~
  ;Set some basic parameters about the shape of the data
  nx=n_elements(data[*,0,0])
  ny=n_elements(data[0,*,0])
  nf=n_elements(data[0,0,*])

  xc=ind[0].xcen
  yc=ind[0].ycen
  dx=ind[0].cdelt1
  dy=ind[0].cdelt2
 
  ; ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Filter out the "bad" pixels
  ; Only need to set the actual bad ones to 0 (not all filters at that location)
  ; as the DEM code will only perform the calc if all filters are non-zero at that pixel location

  ; If the value is above sat_lvl get rid of it
  ; assuming bad from saturation
  sat_lvl=1.5e4
  id=where(data ge sat_lvl,nid)
  if (nid gt 1) then data[id]=0.0

  ; Anything less than 0 just set to 0
  id=where(data le 0,nid)
  if (nid gt 1) then data[id]=0.0

  ; Set the edge few pixels to 0 to minimise problems - from prep'd data that has been shifted (?)
  edg0=10
  data[0:edg0-1,*,*]=0.0
  data[*,0:edg0-1,*]=0.0
  data[nx-edg0:nx-1,*,*]=0.0
  data[*,ny-edg0:ny-1,*]=0.0

  ; Work out the errors
  ; The synoptic data is in DN/px (ours is in DN/px/s, is this ok?)
  edata=fltarr(nx,ny,nf)
  ; Proper way but a bit slow?
  ;  for i=0,nf-1 do edata[*,*,i]=aia_bp_estimate_error(reform(data[*,*,i]),replicate(wavenum[i],nx,ny),n_sample=16)
  ; So instead quick approx based on aia_bp_estimate_error.pro
  npix=4096.^2/(nx*ny)
  edata=fltarr(nx,ny,nf)
  gains=[18.3,17.6,17.7,18.3,18.3,17.6]
  dn2ph=gains*[94,131,171,193,211,335]/3397.
  rdnse=1.15*sqrt(npix)/npix
  drknse=0.17
  qntnse=0.288819*sqrt(npix)/npix
  ; error in DN/px
  serr_per=0.0
  min_snr=3.0
  for i=0, nf-1 do begin
    etemp=sqrt(rdnse^2.+drknse^2.+qntnse^2.+(dn2ph[i]*abs(data[*,*,i]))/(npix*dn2ph[i]^2))
    esys=serr_per*data[*,*,i]/100.
    edata[*,*,i]=sqrt(etemp^2. + esys^2.)
  endfor

  ; Get rid of data with too large an uncertaintity
  id=where(data/edata le min_snr,nid)
  if (nid gt 1) then data[id]=0.0
  
  ;~~~~~~~~~~~~~~~~~~
  ; Response functions and temperature binning
  logtemps=5.7+findgen(17)*0.1
  temps=10d^logtemps

  ; This is is the temperature bin mid-points
  mlogt=get_edges(logtemps,/mean)
  nt=n_elements(mlogt)
  ; Restore AIA response functions
  restore,file='/home/wtb2/Documents/demreg/idl/aia_resp.dat'
  
  ; Only want the coronal ones without 304A
  idc=[0,1,2,3,4,6]

  tr_logt=tresp.logte
  ; Don't need the response outside of the T range we want for the DEM
  gdt=where(tr_logt ge min(logtemps) and tr_logt le max(logtemps),ngd)
  tr_logt=tr_logt[gdt]
  TRmatrix=tresp.all[*,idc]
  TRmatrix=TRmatrix[gdt,*]
  
  ;~~~~~~~~~~~~~~~~~~~~~
  ;DEM calculation
  ; Do DEM calculation without the specified intial weighting
  dn2dem_pos_nb, data,edata,TRmatrix,tr_logt,temps,dem,edem,elogt,chisq,dn_reg,/timed

  ; Convert the DEM output into maps
  dem_maps=replicate(make_map(fltarr(nx,ny),xc=xc,yc=yc,dx=dx,dy=dy,date=date),nt)
  
  ; Show EM rather than DEM (units of cm^-5)
  dt=fltarr(nt)
  for i=0,nt-1 do dt[i]=10d^(logtemps[i+1])-10d^(logtemps[i])
  em=dblarr(nx,ny,nt)
  for xx=0,nx-1 do begin
    for yy=0,ny-1 do begin
      em[xx,yy,*]=dem[xx,yy,*]*dt
    endfor
  endfor
  dem_maps.data=em
  idd='EM '
  
  ; Plotting
  temp_ids=strarr(nt)
  for i=0,nt-1 do temp_ids[i]=string(logtemps[i],format='(f3.1)')+' - '+string(logtemps[i+1],format='(f3.1)')
  dem_maps.id=idd+'logT: '+temp_ids

  ; plot them - assuming using default 16 temp bins and no interpolating
  WINDOW, 0, XSIZE=1000, YSIZE=1000
  !p.multi=[0,4,4]
  !p.Color = '000000'xL
  !p.Background = 'FFFFFF'xL
  loadct,5,/silent

  drang=[1d25,1d30]
  for i=0,nt-1 do plot_map,dem_maps[i],title=dem_maps[i].id,/log,dmin=drang[0],dmax=drang[1],chars=2.0,Background=cgColor('white')
  ;xyouts,10,10,date,chars=1.2,/device
  
  stop
 
end
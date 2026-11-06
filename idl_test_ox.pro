; Test of slit_func_ox (libslitdec_ox_idl.so) from IDL on the swath saved by mkslitf.pro.
; Calls the REDUCE library (slit_func_2d) as mkslitf.pro does, then slit_func_ox.
;   IDL> .run idl_test_ox
; or from the shell:  idl -e "idl_test_ox"   (with this folder on !PATH or as cwd)

function sfox_reduce, lib, d, sp, sfsm, msk, sfbin
  ; slit_func_2d of the REDUCE library with mkslitf.pro's argument list
  nc = d.nc & nysf = d.nysf & osamp = d.osamp
  sp = dblarr(nc) & sp_old = sp & unc = dblarr(nc)
  sfbin = dblarr(nc, nysf)
  ny = (nysf + 1L) * osamp + 1L
  sfsm = dblarr(ny)
  msk = d.msk
  xi = bytarr(16L, 4L, ny, nc)
  zeta = bytarr(16L, 3L * (osamp + 1L), nysf, nc)
  m_zeta = lonarr(nysf, nc)
  l_Aij = dblarr(ny * (4 * osamp + 1))
  l_bj = dblarr(ny)
  p_Aij = dblarr(((4 * d.idelta_x + 1) > 3) * nc)
  p_bj = dblarr(nc)
  return, CALL_EXTERNAL(lib, 'slit_func_2d', long(nc), long(nysf), long(osamp), d.sf $
                       , d.ycen_fract, long(d.ycen_int), long(d.y_lower_lim), d.PSF_curve $
                       , long(d.idelta_x), double(d.lambda_sp), double(d.lambda_sf), msk $
                       , sfbin, sp, sp_old, sfsm, unc, xi, zeta, m_zeta, l_Aij, l_bj $
                       , p_Aij, p_bj, long(0))
end

function sfox_call, lib, d, p, sp, sfsm, msk, sfbin, sp_fine, info
  ; slit_func_ox with all 24 arguments; p holds the parameters
  nc = d.nc & nysf = d.nysf & osamp = d.osamp
  ny = (nysf + 1L) * osamp + 1L
  msk = d.msk
  sfbin = dblarr(nc, nysf) & sp = dblarr(nc) & sp_fine = dblarr(nc * long(p.osample_x))
  sfsm = dblarr(ny) & unc = dblarr(nc) & info = dblarr(11)
  return, CALL_EXTERNAL(lib, 'slit_func_ox', long(nc), long(nysf), long(osamp) $
                       , long(p.osample_x), double(d.sf), double(d.ycen_fract), long(d.ycen_int) $
                       , long(d.y_lower_lim), double(d.PSF_curve) $
                       , double(p.mu1_x), double(p.mu2_x), double(p.mu3_x) $
                       , double(p.mu1_y), double(p.mu2_y), long(p.maxiter), double(p.kappa) $
                       , long(0), msk, sfbin, sp, sp_fine, sfsm, unc, info)
end

pro sfox_report, name, d, sp, sfbin, msk, ref=ref
  dx = d.idelta_x
  i = lindgen(d.nc - 2 * dx) + dx
  res = (sfbin - d.sf)[i, *]
  good = where(msk[i, *] gt 0)
  s = string(name, median(sp[i]), stddev(res[good]), total(msk eq 0), $
             format='(A-44, " median sP=", F10.1, "  rms(model-im)=", F7.2, "  masked=", I4)')
  if n_elements(ref) gt 0 then s += string(stddev(sp[i] / ref[i] - 1), format='("  rms(sP/ref-1)=", E9.2)')
  print, s
end

pro idl_test_ox, savfile, old_lib=old_lib, new_lib=new_lib
  here = file_dirname((routine_info('idl_test_ox', /source)).path)
  if n_elements(savfile) eq 0 then savfile = here + '/../test.sav'
  if n_elements(old_lib) eq 0 then old_lib = expand_path('~/IDL/REDUCE') + '/reduce.so.darwin.arm64.64'
  if n_elements(new_lib) eq 0 then new_lib = here + '/libslitdec_ox_idl.so'
  restore, savfile
  d = {nc:nc, nysf:nysf, osamp:osamp, sf:sf, ycen_fract:ycen_fract, ycen_int:ycen_int, $
       y_lower_lim:y_lower_lim, PSF_curve:PSF_curve, idelta_x:idelta_x, lambda_sp:lambda_sp, $
       lambda_sf:lambda_sf, msk:msk}
  print, savfile, nc, nysf, osamp, y_lower_lim, $
         format='(A, ": nc=", I0, " nysf=", I0, " osamp=", I0, " y_lower_lim=", I0)'
  print, 'old: ', old_lib
  print, 'new: ', new_lib

  ; 1. REDUCE library as mkslitf.pro calls it
  r = sfox_reduce(old_lib, d, sp0, sfsm0, msk0, sfbin0)
  sfox_report, 'reduce.so (lambda_sf='+strtrim(lambda_sf, 2)+')', d, sp0, sfbin0, msk0

  ; 2. slitdec_ox, recommended parameters
  p = {osample_x:8L, mu1_x:0d, mu2_x:0d, mu3_x:1.64d-4, mu1_y:1d-2, mu2_y:0d, maxiter:100L, kappa:10d}
  r = sfox_call(new_lib, d, p, sp8, sfsm8, msk8, sfbin8, fine, info)
  print, 'return code', r
  sfox_report, 'slitdec_ox ox=8 mu3_x=1.64e-4 mu1_y=1e-2', d, sp8, sfbin8, msk8, ref=sp0
  print, info[3], info[4], info[7], info[8], info[9], info[10], $
         format='("    iterations ", I0, "  delta_x ", I0, "  damping x ", G10.3, "  res x ", F5.3, " px  damping y ", G10.3, "  res y ", F5.3, " px")'
  i = lindgen(nc - 2 * idelta_x) + idelta_x
  fine = reform(fine, 8, nc)
  print, max(abs(total(fine[*, i], 1) / 8 - sp8[i]) / abs(sp8[i])), $
         format='("    max |mean(sP_fine) - sP| / sP = ", E9.2)'

  ; 3. slit-function penalty: spectrum and slit function
  foreach mu1_y, [1d-3, 1d-1] do begin
    q = p & q.mu1_y = mu1_y
    r = sfox_call(new_lib, d, q, spy, sfsmy, msky, sfbiny, finey, infoy)
    sfox_report, string('slitdec_ox mu1_y=', mu1_y, format='(A, G0.3)'), d, spy, sfbiny, msky, ref=sp8
    print, infoy[10], max(sfsmy) / max(sfsm8), format='("    res y ", F5.3, " px  max(sL) / max(sL, mu1_y=1e-2) ", F6.4)'
  endforeach

  ; 4. historic model: osample_x=1, no spectrum penalty, mu1_y equivalent to lambda_sf.
  ;    Differs from reduce.so only by the tilt reference point (y-subpixel centre).
  q = p & q.osample_x = 1L & q.mu3_x = 0d & q.mu1_y = 1.927d-3   ; lambda_sf=2 on test.sav
  r = sfox_call(new_lib, d, q, sp1, sfsm1, msk1, sfbin1, fine1, info1)
  sfox_report, 'slitdec_ox ox=1 mu_x=0 mu1_y=1.927e-3', d, sp1, sfbin1, msk1, ref=sp0

  ; 5. wrong calls are refused with a message
  q = p & q.mu1_y = 0d
  print, 'mu1_y = mu2_y = 0: return code', sfox_call(new_lib, d, q, a, b, c, e, f, g)
  print, 'too few arguments: return code', CALL_EXTERNAL(new_lib, 'slit_func_ox', long(nc))
end

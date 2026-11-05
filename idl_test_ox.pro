; Test of libslitdec_ox_idl.so from IDL on the swath saved by mkslitf.pro.
; Calls slit_func_2d exactly as mkslitf.pro does (25 arguments) with the REDUCE
; library and with slitdec_ox, then with the optional arguments.
;   IDL> .run idl_test_ox
; or from the shell:  idl -e "idl_test_ox"   (with this folder on !PATH or as cwd)

function sfox_call, lib, d, sp, sfsm, msk, sfbin, unc, extra=extra
  ; One slit_func_2d call with mkslitf.pro's argument list; d holds the restored swath.
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
  if n_elements(extra) eq 0 then begin
    iret = CALL_EXTERNAL(lib, 'slit_func_2d', long(nc), long(nysf), long(osamp), d.sf $
                        , d.ycen_fract, long(d.ycen_int), long(d.y_lower_lim), d.PSF_curve $
                        , long(d.idelta_x), double(d.lambda_sp), double(d.lambda_sf), msk $
                        , sfbin, sp, sp_old, sfsm, unc, xi, zeta, m_zeta, l_Aij, l_bj $
                        , p_Aij, p_bj, long(0))
  endif else begin
    e = extra
    sp_fine = e.sp_fine & info = e.info   ; variables, so the outputs come back
    iret = CALL_EXTERNAL(lib, 'slit_func_2d', long(nc), long(nysf), long(osamp), d.sf $
                        , d.ycen_fract, long(d.ycen_int), long(d.y_lower_lim), d.PSF_curve $
                        , long(d.idelta_x), double(d.lambda_sp), double(d.lambda_sf), msk $
                        , sfbin, sp, sp_old, sfsm, unc, xi, zeta, m_zeta, l_Aij, l_bj $
                        , p_Aij, p_bj, long(0) $
                        , long(e.maxiter), long(e.osample_x), double(e.mu1_x), double(e.mu2_x) $
                        , double(e.mu3_x), double(e.mu1_y), double(e.mu2_y), sp_fine, info)
    e.sp_fine = sp_fine & e.info = info
    extra = e
  endelse
  return, iret
end

pro sfox_report, name, d, sp, sfbin, msk, ref=ref
  dx = d.idelta_x
  i = lindgen(d.nc - 2 * dx) + dx
  res = (sfbin - d.sf)[i, *]
  good = where(msk[i, *] gt 0)
  s = string(name, median(sp[i]), stddev(res[good]), total(msk eq 0), $
             format='(A-36, " median sP=", F10.1, "  rms(model-im)=", F7.2, "  masked=", I4)')
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
  r = sfox_call(old_lib, d, sp0, sfsm0, msk0, sfbin0, unc0)
  sfox_report, 'reduce.so (as mkslitf)', d, sp0, sfbin0, msk0

  ; 2. slitdec_ox with the 25 mkslitf arguments: defaults (osample_x 8, mu3_x 1.64e-4)
  setenv, 'SLITDEC_OX_OSAMPLE_X=' & setenv, 'SLITDEC_OX_MU3_X=' & setenv, 'SLITDEC_OX_MU1_Y='
  r = sfox_call(new_lib, d, sp8, sfsm8, msk8, sfbin8, unc8)
  print, 'return code', r
  sfox_report, 'slitdec_ox defaults (as mkslitf)', d, sp8, sfbin8, msk8, ref=sp0

  ; 3. historic model through env vars: must match reduce.so up to the
  ;    half-subpixel tilt correction (t_offset), see NOTES.md
  setenv, 'SLITDEC_OX_OSAMPLE_X=1' & setenv, 'SLITDEC_OX_MU3_X=0'
  r = sfox_call(new_lib, d, sp1, sfsm1, msk1, sfbin1, unc1)
  sfox_report, 'slitdec_ox env osample_x=1 mu=0', d, sp1, sfbin1, msk1, ref=sp0
  setenv, 'SLITDEC_OX_T_OFFSET=0'
  r = sfox_call(new_lib, d, sp1, sfsm1, msk1, sfbin1, unc1)
  sfox_report, '  same, t_offset=0 (slitdec.c tilt)', d, sp1, sfbin1, msk1, ref=sp0
  print, max(abs(sp1 - sp0)), max(abs(sfsm1 - sfsm0)), max(abs(sfbin1 - sfbin0)), total(msk1 ne msk0), $
         format='("    max |dsP| ", G10.3, "  |dsL| ", G10.3, "  |dmodel| ", G10.3, "  mask diffs ", I0)'
  setenv, 'SLITDEC_OX_OSAMPLE_X=' & setenv, 'SLITDEC_OX_MU3_X=' & setenv, 'SLITDEC_OX_T_OFFSET='

  ; 4. optional arguments, info and sP_fine; recommended slit-function penalty
  for k = 0, 1 do begin
    mu1_y = k eq 0 ? 0d : 1d-2
    e = {maxiter:100L, osample_x:8L, mu1_x:0d, mu2_x:0d, mu3_x:1.64d-4, mu1_y:mu1_y, mu2_y:0d, $
         sp_fine:dblarr(nc * 8L), info:dblarr(11)}
    r = sfox_call(new_lib, d, spx, sfsmx, mskx, sfbinx, uncx, extra=e)
    sfox_report, k eq 0 ? 'slitdec_ox args, legacy lambda_sf' : 'slitdec_ox args, mu1_y=1e-2', $
                 d, spx, sfbinx, mskx, ref=sp8
    print, e.info[3], e.info[7], e.info[8], e.info[9], e.info[10], $
           format='("    iterations ", I0, "  damping x ", G10.3, "  res x ", F5.3, " px  damping y ", G10.3, "  res y ", F5.3, " px")'
    fine = reform(e.sp_fine, 8, nc)
    i = lindgen(nc - 2 * idelta_x) + idelta_x
    print, max(abs(total(fine[*, i], 1) / 8 - spx[i]) / abs(spx[i])), $
           format='("    max |mean(sP_fine) - sP| / sP = ", E9.2)'
  endfor
  print, 'slit function: legacy vs mu1_y=1e-2 at the profile maximum: ', max(sfsm8), max(sfsmx)
end

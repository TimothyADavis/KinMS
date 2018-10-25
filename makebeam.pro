function makebeam, beam_pix,rot=rot,npixel=npixel
	if not keyword_set(rot) then rot=0.0
  	if not keyword_set(npixel) then begin
		npixel=max(beam_pix)*[5,5]
  		if not odd(npixel[0]) then npixel+=1
	endif
	cents=FLOOR((npixel-1)/2)
	beam_pix_std=beam_pix/2.355d
	dirfac=(tan(!dtor*rot) gt 0)*2 - 1 ;; get the sign of the rotatation
	x=transpose(findgen(npixel[1],npixel[0]) mod npixel[1]) - cents[0]
	y=(findgen(npixel[0],npixel[1]) mod npixel[0]) - cents[1]
	a=(cos(!dtor*(rot))^2)/(2.0*(beam_pix_std[0]^2)) + (sin(!dtor*(rot))^2)/(2.0*(beam_pix_std[1]^2))
	b=((1*dirfac)*(sin(2.0*!dtor*(rot))^2)/(4.0*(beam_pix_std[0]^2))) + ((-1*dirfac)*(sin(2.0*!dtor*(rot))^2)/(4.0*(beam_pix_std[1]^2)))
	c=(sin(!dtor*(rot))^2)/(2.0*(beam_pix_std[0]^2)) + (cos(!dtor*(rot))^2)/(2.0*(beam_pix_std[1]^2))
	psf=exp(-1*(a*(x^2) - 2.0*b*(x*y) + c*(y^2)))
	error = check_math(MASK=32)
    return,psf
end
function makebeam, npixel,st_dev,rot=rot
	if not keyword_set(rot) then rot=0.0
	cents=npixel/2.
    st_dev=st_dev/2.355	
	if tan(!dtor*rot) eq 0 then dirfac=1 else dirfac=tan(!dtor*rot)/abs(tan(!dtor*rot))
	x=transpose(findgen(npixel[1],npixel[0]) mod npixel[1]) - cents[0]
	y=(findgen(npixel[0],npixel[1]) mod npixel[0]) - cents[1]   
    a=(cos(!dtor*(rot))^2)/(2.0*(st_dev[0]^2)) + (sin(!dtor*(rot))^2)/(2.0*(st_dev[1]^2))
    b=((1*dirfac)*(sin(2.0*!dtor*(rot))^2)/(4.0*(st_dev[0]^2))) + ((-1*dirfac)*(sin(2.0*!dtor*(rot))^2)/(4.0*(st_dev[1]^2)))
    c=(sin(!dtor*(rot))^2)/(2.0*(st_dev[0]^2)) + (cos(!dtor*(rot))^2)/(2.0*(st_dev[1]^2))
    psf=exp(-1*(a*(x^2) - 2.0*b*(x*y) + c*(y^2)))         
	error = check_math(MASK=32)
    return,psf
end
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin

	resolution = "40x40x40"
	file_path = "~/Desktop/wrfout_d01_0001-01-01_00_00_00"

	; We would like to normalize the space spans.
	real_x_span = 3900.0
	real_y_span = 3900.0
	real_z_span = 1.0

	f_in = addfile(file_path, "r")

	dim_sizes = dimsizes(f_in->Times)
	num_time = dim_sizes(0)
	delete(dim_sizes)

	dim_sizes = dimsizes(f_in->U)
	num_z_full = dim_sizes(1)
	num_y_full = dim_sizes(2)
	num_x_half = dim_sizes(3)-1 ; exclude virtual grids
	delete(dim_sizes)

	dim_sizes = dimsizes(f_in->V)
	num_y_half = dim_sizes(2)-1 ; exclude virtual grids
	num_x_full = dim_sizes(3)
	delete(dim_sizes)

	dim_sizes = dimsizes(f_in->WW)
	num_z_half = dim_sizes(1)
	delete(dim_sizes)

	dx = 1.0/num_x_full
	dy = 1.0/num_y_full

	time = WRF_Times2Udunits_c(f_in->Times, 0)

	x_full = new(num_x_full, double)
	do i = 0, num_x_full-1
		x_full(i) = dx*0.5+i*dx
	end do
	x_half = new(num_x_half, double)
	do i = 0, num_x_half-1
		x_half(i) = i*dx
	end do

	y_full = new(num_y_full, double)
	do j = 0, num_y_full-1
		y_full(j) = dy*0.5+j*dy
	end do
	y_half = new(num_y_half, double)
	do j = 0, num_y_half-1
		y_half(j) = j*dy
	end do

	z_full = new(num_z_full, double)
	z_full = (/f_in->ZNU(0,::-1)/)
	z_half = new(num_z_half, double)
	z_half = (/f_in->ZNW(0,::-1)/)

	global_att = True
	global_att@time_step_size_in_seconds = f_in@DT
	global_att@domain_type = "Cartesian 3d"
	global_att@domain_x_span = "0-1"
	global_att@domain_real_x_span = "0-"+real_x_span
	global_att@domain_x_boundary_condition = "periodic-periodic"
	global_att@domain_y_span = "0-1"
	global_att@domain_real_y_span = "0-"+real_y_span
	global_att@domain_y_boundary_condition = "periodic-periodic"
	global_att@domain_z_span = "0-1"
	global_att@domain_z_boundary_condition = "rigid-rigid"
	global_att@num_tracer = 2

	do l = 0, num_time-1
		file_path = "wrf_les."+resolution+"."+f_in->Times(l,:)+".nc"
		system("echo 'writing "+file_path+" ...'")
		system("rm -rf "+file_path)
		f_out = addfile(file_path, "c")

		setfileoption(f_out, "DefineMode", True)

		fileattdef(f_out, global_att)

		dim_names = (/"x",       "x_bnds",  "y",       "y_bnds",  "z",       "z_bnds",  "time"  /)
		dim_sizes = (/num_x_full,num_x_half,num_y_full,num_y_half,num_z_full,num_z_half,num_time/)
		dim_unlim = (/False,     False,     False,     False,     False,     False,     True    /)
		filedimdef(f_out, dim_names, dim_sizes, dim_unlim)

		filevardef(f_out, "x",      "double", (/"x"     /))
		filevardef(f_out, "x_bnds", "double", (/"x_bnds"/))
		filevardef(f_out, "y",      "double", (/"y"     /))
		filevardef(f_out, "y_bnds", "double", (/"y_bnds"/))
		filevardef(f_out, "z",      "double", (/"z"     /))
		filevardef(f_out, "z_bnds", "double", (/"z_bnds"/))
		filevardef(f_out, "time",   "double", (/"time"  /))
		filevardef(f_out, "u",      "double", (/"time","z",     "y_bnds","x"     /))
		filevardef(f_out, "v",      "double", (/"time","z",     "y",     "x_bnds"/))
		filevardef(f_out, "w",      "double", (/"time","z_bnds","y",     "x"     /))
		filevardef(f_out, "q0",     "double", (/"time","z",		"y",	 "x"	 /))
		filevardef(f_out, "q1",     "double", (/"time","z",		"y",	 "x"	 /))

		f_out->time = (/time(l)/)
        f_out->time@units = "hours since 0001-01-01 00:00:00"
		f_out->x = (/x_full/)
		f_out->x@long_name = "x"
		f_out->x@units = "1"
		delete(f_out->x@_FillValue)
		f_out->x_bnds = (/x_half/)
		f_out->x_bnds@long_name = "x bounds"
		f_out->x_bnds@units = "1"
		delete(f_out->x_bnds@_FillValue)
		f_out->y = (/y_full/)
		f_out->y@long_name = "y"
		f_out->y@units = "1"
		delete(f_out->y@_FillValue)
		f_out->y_bnds = (/y_half/)
		f_out->y_bnds@long_name = "y bounds"
		f_out->y_bnds@units = "1"
		delete(f_out->y_bnds@_FillValue)
		f_out->z = (/z_full/)
		f_out->z@long_name = "z"
		f_out->z@units = "1"
		delete(f_out->z@_FillValue)
		f_out->z_bnds = (/z_half/)
		f_out->z_bnds@long_name = "z bounds"
		f_out->z_bnds@units = "1"
		delete(f_out->z_bnds@_FillValue)

		f_out->u(0,:,:,:) = (/f_in->U(l,0:num_z_full-1:-1,0:num_y_full-1,0:num_x_half-1)/real_x_span/)
		f_out->u@long_name = "velocity x component"
		f_out->u@units = "s-1"
		f_out->v(0,:,:,:) = (/f_in->V(l,0:num_z_full-1:-1,0:num_y_half-1,0:num_x_full-1)/real_y_span/)
		f_out->v@long_name = "velocity y component"
		f_out->v@units = "s-1"
		mu = f_in->MU(l,0:num_y_full-1,0:num_x_full-1)+f_in->MUB(l,0:num_y_full-1,0:num_x_full-1)
		do k = 0, num_z_half-1
			f_out->w(0,k,:,:) = (/f_in->WW(l,num_z_half-1-k,0:num_y_full-1,0:num_x_full-1)/mu/)
		end do
		f_out->w@long_name = "velocity z component"
		f_out->w@units = "s-1"
		do k = 0, num_z_full-1
			f_out->q0(0,k,:,:) = 1.0
		end do
		f_out->q0@long_name = "background tracer"
		f_out->q0@units = "kg"
		do k = 0, num_z_full-1
			f_out->q1(0,k,:,:) = (/f_in->tr17_1(l,num_z_full-1-k,0:num_y_full-1,0:num_x_full-1)*mu/)
		end do
		f_out->q1@long_name = "test tracer"
		f_out->q1@units = "kg"
	end do

end

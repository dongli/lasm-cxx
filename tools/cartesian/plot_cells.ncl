load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	; f = addfile("lasm.00001.nc", "r")
	; f = addfile("fake_data.00000.nc", "r")
	f = addfile("diff.nc", "r")

	wks = gsn_open_wks("pdf", "cells")

	gsn_define_colormap(wks, "rainbow")

	res = True
	res@cnLinesOn = False
	res@cnFillOn = True
	; res@cnFillMode = "RasterFill"
	res@gsnSpreadColors = True
	; res@cnLevelSelectionMode = "ManualLevels"
	; res@cnMinLevelValF = 0.0
	; res@cnMaxLevelValF = 1.0
	; res@cnLevelSpacingF = 0.05
	; res@cnMinLevelValF = -0.1
	; res@cnMaxLevelValF = 0.1
	; res@cnLevelSpacingF = 0.01

	plot = gsn_csm_contour(wks, f->q0(0,0,:,:), res)

end

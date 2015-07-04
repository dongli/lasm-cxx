load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function plot_cell(wks, test_case, var, time, data)
begin

	if (isStrSubset(test_case, "wrf_les")) then
        gsn_define_colormap(wks, "WhBlGrYeRe")
    end if

    ; Set plot resources.
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    ; res@gsnLeftStringOrthogonalPosF = 0
    ; res@gsnRightStringOrthogonalPosF = 0
    ; res@gsnLeftString = "Test tracer: "+var
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnConstFEnableFill = True
    res@cnConstFLabelOn = False
    res@cnLevelSelectionMode = "ManualLevels"
    if (isStrSubset(test_case, "wrf_les")) then
        res@trYReverse = True
        if (var .eq. "w") then
            res@cnMinLevelValF = -0.001
            res@cnMaxLevelValF = 0.001
            res@cnLevelSpacingF = 0.0001
        end if
        if (var .eq. "q0") then
            res@cnMinLevelValF = -0.1
            res@cnMaxLevelValF = 1.1
            res@cnLevelSpacingF = 0.1
        end if
        if (var .eq. "q1") then
            res@cnMinLevelValF = 0.0
            res@cnMaxLevelValF = 21000.0
            res@cnLevelSpacingF = 1000.0
        end if
    end if

    plot = gsn_csm_contour(wks, data(0,:,:), res)

    return plot

end

begin
    ; wrf_les
    root = "."
    test_case = "wrf_les"
    var = "q1"
    resolution = "40x40x40"
    time_stamp = "%5.5d"
    start_time = 7200
    end_time = 7200
    time_step = 100

    file_prefix = "lasm."+test_case+"."+resolution

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    num_time = dimsizes(fs)

    if (var .ne. "") then
        wks = gsn_open_wks("pdf", file_prefix+"."+var+".cells")
    else
        wks = gsn_open_wks("pdf", file_prefix+".cells")
    end if

    do l = 0, num_time-1
        system("echo file: "+fs(l))
        f = addfile(fs(l), "r")
        if (isStrSubset(test_case, "wrf_les")) then
            q = (/f->$var$(:,:,20,:)/)
        end if

        copy_VarCoords(f->$var$, q)
        plot = plot_cell(wks, test_case, var, f->time, q)
        draw(plot)
        frame(wks)
    end do

end

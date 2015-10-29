;; ==================================================================
;; Manipulate errors from Glamdring
;; ==================================================================
;; warrenj 20151028 Routine to manipulate and organise the results
;; from the glamdring Monte Carlo.

pro man_errors
galaxies = ['ngc3557', 'ic1459', 'ic1531', 'ic4296', 'ngc0612', 'ngc1399', 'ngc3100', 'ngc7075', 'pks0718-34', 'eso443-g024']
galaxy = galaxies[1]

glamdring_file = "/Data/vimosindi/analysis/" + galaxy + $
    "/glamdring_result.txt"
glamdring_file2 = "/Data/vimosindi/analysis/" + galaxy + $
    "/glamdring_result2.txt"

RDFLOAT, glamdring_file, bins, rep, vel, sig, COLUMNS=[2,3,4,5]
RDFLOAT, glamdring_file2, bins2, rep2, h3s, h4s, COLUMNS=[2,3,4,5]

;v = make_array(max(bins))
;v_s = make_array(max(bins))
;s = make_array(max(bins))
;s_s = make_array(max(bins))
;h3 = make_array(max(bins2))
;h3_s = make_array(max(bins2))
;h4 = make_array(max(bins2))
;h4_s = make_array(max(bins2))








return
end


pro save
;repeat over all bins:
for bin=0, max(bins)-1 do begin
i_bin = where(bins eq bin)
v[bin] = mean(vel[i_bin])
v_s[bin] = stddev(vel[i_bin])
s[bin] = mean(sig[i_bin])
s_s[bin] = stddev(sig[i_bin])
endfor
;repeat over all bins:
for bin=0, max(bins2)-1 do begin
i_bin = where(bins2 eq bin2)
h3[bin] = mean(h3s[i_bin])
h3_s[bin] = stddev(h3s[i_bin])
h4[bin] = mean(h4s[i_bin])
h4_s[bin] = stddev(h4s[i_bin])
endfor


v_file = "/Data/vimosindi/analysis/" + galaxy + "/results/glamdring_v.dat"
s_file = "/Data/vimosindi/analysis/" + galaxy + "/results/glamdring_s.dat"
h3_file = "/Data/vimosindi/analysis/" + galaxy + "/results/glamdring_h3.dat"
h4_file = "/Data/vimosindi/analysis/" + galaxy + "/results/glamdring_h4.dat"


forprint, v, v_s, textout=v_file
forprint, s, s_s, textout=s_file
forprint, h3, h3_s, textout=h3_file
forprint, h4, h4_s, textout=h4_file



return
end

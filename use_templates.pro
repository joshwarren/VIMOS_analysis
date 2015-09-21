;; ==================================================================
;; Using a limited number of templates
;; ==================================================================
;; warrenj 20150529 Rewriting the routine to which specifies which
;; templates from the full library are to be used.

pro use_templates, galaxy, templatesToUse

;templatesToUse = [29, 30, 44, 64 ,71, 75, 103, 104, 177, 183,243, $
;246, 251, 252, 276, 286, 409, 452, 455, 459, 472, 498, 508, 544, $
;548, 587, 684, 697, 720, 722, 732, 741, 764, 781, 789, 824, 838, $
;873, 884, 890, 927, 931, 932, 940, 941, 949, 965] 

template_weighting = '/Data/vimosindi/analysis/' + galaxy + $
	'/templates.txt' 

READCOL, template_weighting, templatesToUse, FORMAT = 'I', /SILENT



return 
end

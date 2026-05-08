#
# convert magnetic form factor table
# @author Tobias Weber <tweber@ill.fr>
# @date 8-may-2026
# @license see 'LICENSE' file
#

using Printf

# get table from here: https://github.com/SunnySuite/Sunny.jl/blob/main/src/FormFactor.jl
include("FormFactor.jl");


f = stdout
outfile = "magffacts.xml"
try
	global f
	f = open(outfile, "w")
catch err
	@printf("Could not open file \"%s\" for writing.\n", outfile)
	exit()
end


# header
@printf(f, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
@printf(f, "<magnetic_form_factors>\n")


# iterate form factors
for ffact in magnetic_ion_data
	name = ffact[1]
	cfg = ffact[2]
	term = ffact[3]
	coeffs_0 = ffact[4]
	coeffs_2 = ffact[5]

	@printf(f, "\t<ion name=\"%s\" electrons=\"%s\" terms=\"%s\">\n", ffact[1], ffact[2], ffact[3])

	# coefficients
	@printf(f, "\t\t<coefficients_0> ")
	for coeff in coeffs_0
		@printf(f, "%f ", coeff)
	end
	@printf(f, "</coefficients_0>\n")

	# coefficients
	@printf(f, "\t\t<coefficients_2> ")
	for coeff in coeffs_2
		@printf(f, "%f ", coeff)
	end
	@printf(f, "</coefficients_2>\n")

	@printf(f, "\t</ion>\n")
end


@printf(f, "</magnetic_form_factors>\n")

close(f)
@printf("Wrote magnetic form factors to \"%s\".\n", outfile)

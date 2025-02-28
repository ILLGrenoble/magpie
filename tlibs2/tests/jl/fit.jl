#
# tlibs2 julia module test
# @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
# @date 2017 -- 2018
# @license GPLv3, see <<LICENSE>> file
# @desc Forked on 7-Nov-2018 from my privately and TUM-PhD-developed <<tlibs>> project (https://github.com/t-weber/tlibs).
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

include("tl2.jl")


# load data file
#a = tl2.loadinstr("/Users/tw/Measurements/mnsi/mira-18/data/13511_00035079.dat")
a = tl2.loadinstr("/home/tw/Measurements/mnsi/mira_18/data/13511_00035079.dat")
cols = a[1]
data = a[2]
#println(cols)
#println(typeof(data))

idxx = findfirst(isequal("h"), cols)
idxCtr = findfirst(isequal("ctr1"), cols)

xs = data[:,idxx]
cts = data[:,idxCtr]
cts_err = tl2.poisson_err.(cts)
#println(typeof(xs))
#println(typeof(cts))
#println(typeof(cts_err))


# print data
println(lpad(cols[idxx], 12), " ", lpad(cols[idxCtr], 12))
for (x, ct, err) in zip(xs, cts, cts_err)
	println(lpad(x, 12), " ", lpad(ct, 12), " ", lpad(err, 12))
end
println()


# fit gaussian
fitresult = tl2.fit(tl2.gauss_model_amp, xs, cts, cts_err, fixed = [ "offs" ],
	values = Dict("x0" => 1., "sigma" => 0.003, "amp" => 20000., "offs" => 5.),
	errors = Dict("x0" => 0.05, "sigma" => 0.002, "amp" => 2000., "offs" => 2.5))
println(fitresult)
println()

data_fine = tl2.eval_func(fitresult, minimum(xs), maximum(xs))


# write results
f = open("tst_pts.dat", "w")
for (x,y,yerr) in zip(xs, cts, cts_err)
	println(f, rpad(x, 20), " ", rpad(y, 20), " ", rpad(yerr, 20))
end
close(f)

f = open("tst_fit.dat", "w")
for (x,y) in zip(data_fine[1], data_fine[2])
	println(f, rpad(x, 20), " ", rpad(y, 20))
end
close(f)


# plot results
run(`gnuplot -p -e "plot \"tst_pts.dat\" w yerrorbars pt 7 lc rgb \"#ff0000\",
	\"tst_fit.dat\" with lines lw 2 lc rgb \"#0000ff\""`)

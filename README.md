IsingModels.jl
==============

Draw inexact samples from the Ising model using a Gibbs sampler.

	using IsingModels

	d = Ising(20, isferromagnetic = true)

	for t in 1.0:1.0:5.0
	    d.temperature = t
	    show(rand(d))
	end

	d.isferromagnetic = false

	for t in 1.0:1.0:5.0
	    d.temperature = t
	    show(rand(d))
	end

	writecsv("draw.csv",
	         rand(Ising(200,
	                    isferromagnetic = true,
	                    temperature = 3.0)))

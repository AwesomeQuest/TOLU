### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 9bbadf02-b08c-11ef-36df-b3ec229764a4
begin
	using CairoMakie
	using LinearAlgebra, IterTools
	using Typstry, MakieTeX

	Iter = Iterators
end

# ╔═╡ ceb342a1-0d9a-45b9-83a4-5b55bf65e471
function yf(y0, h, N, f, step)
	iterated((x->(x[1] + step(f, x[1] ,h, x[2]*h),x[2]+1)), (y0,0)) |> x->Iter.take(x,N) .|> first
end

# ╔═╡ f89cb49c-727b-41d1-8db9-55b48661407c
function eulerstep(f,y,h,t)
	h*f(t,y)
end

# ╔═╡ 096fb1ac-0abb-4a50-85a6-d6e1f62fce2c
α, β, γ, δ = 0.5, 0.01, 0.005, 0.2

# ╔═╡ 5fb77950-b7b1-4ef6-a478-5fe94405d884
y(t,yp) = [
		α*yp[1] - β*yp[1]*yp[2]
		γ*yp[1]*yp[2] - δ*yp[2]
	]

# ╔═╡ 782d5143-7657-46b1-bcd9-19de24aeb4cd
let
	N = 35
	xs = range(0.,100, N)
	ys = range(0.,100, N)
	ps = [Point2(x,y) for x in xs for y in ys]
	ns = map(x->y(0,x) |> Vec2, ps) 
	ls = norm.(ns)
	fig,ax,l = arrows(
		ps, ns, 
		normalize = true,
		arrowcolor = ls .|> log10, 
		lengthscale=2,
		axis=(
			xlabel=typst"$y_1$",
			ylabel=L"y_2",
		),
	)
	Colorbar(fig[1,2], l, tickformat = vals->[L"10^{%$(val)}" for val in vals])
	fig
end

# ╔═╡ 1b07f4be-afbe-46af-a648-908d336e1b27
let
	y0 = [40,9]
	y1 = [30,10]
	t0 = 0; tend=15; N = 100; h = tend/N;
	output = yf(y0, h, N, y, eulerstep)
	fig,ax = lines(t0+h:h:tend, output .|> x->x[1], label=L"y_1 \text{ með }")
	lines!(t0+h:h:tend, output .|> x->x[2], label=L"y_2 \text{ með }")
	output2 = yf(y1, h, N, y, eulerstep)
	lines!(t0+h:h:tend, output2 .|> x->x[1], label=L"y_1 \text{ með }")
	lines!(t0+h:h:tend, output2 .|> x->x[2], label=L"y_1 \text{ með }")
	axislegend(ax)
	current_figure()
end

# ╔═╡ Cell order:
# ╠═9bbadf02-b08c-11ef-36df-b3ec229764a4
# ╠═ceb342a1-0d9a-45b9-83a4-5b55bf65e471
# ╠═f89cb49c-727b-41d1-8db9-55b48661407c
# ╠═096fb1ac-0abb-4a50-85a6-d6e1f62fce2c
# ╠═5fb77950-b7b1-4ef6-a478-5fe94405d884
# ╠═782d5143-7657-46b1-bcd9-19de24aeb4cd
# ╠═1b07f4be-afbe-46af-a648-908d336e1b27

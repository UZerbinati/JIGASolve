module JuliaIGA

using Plots #Plotting library
using LinearAlgebra

include("CAD.jl");

greet() = print("Hello World, I'm JuliaIGA release 0.0.1")
# Variational Formulation

struct VarFormulation
	a::String
	trial::String
	test::String
	f::String
	n::Int32
end

# EXAMPLES

function BezierEx()

	C = BezierCurve2D(2,[[0 0],[0.5 1],[1 0]])
	BezierPlot(C,0.01);
end
function BSplineEx()
	C1 = BSplineCurve(2,[0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,1.0,2.0,3.0,3.0,3.0]);
	BSplinePlot(C1,0.1);

	C2 = BSplineCurve2D(3,[[-4.0 -4.0],[-2.0 4.0],[2.0 -4.0],[4.0 4.0]],[0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]);
	BSplinePlot(C2,0.01);

end
function NURBSEx()
	S = NURBSCurve2D(2,[[1.0 0.0],[1.0 1.0],[0.0 1.0],[-1.0 1.0],[-1.0 0.0],[-1.0 -1.0],[0.0 -1.0],[1.0 -1.0],[1.0 0.0]],[1.0,0.5*sqrt(2),1.0,0.5*sqrt(2),1.0,0.5*sqrt(2),1.0,0.5*sqrt(2),1.0],[0.0,0.0,0.0,0.5*π,0.5*π,π,π,1.5*π,1.5*π,2*π,2*π,2*π]);
	NURBSPlot(S,0.01);
end
function SurfaceEx()
	S = BezierSurface(2,2,[[1.0 2.0 1.0] [3.0 4.0 2.0];[2.0 1.0 1.0] [4.0 1.0 2.0]])
	println(DeCasteljau(S,[0.3, 0.1]))
	BezierPlot(S,0.01)
end
end # module

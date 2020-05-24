module JuliaIGA

using Plots #Plotting library
using LinearAlgebra

include("CAD.jl");
include("NM/BEF.jl");
greet() = print("Hello World, I'm JuliaIGA release 0.0.1")
# Variational Formulation
BlaBla = true;

mutable struct VariationalFormulation
	a::String
	trial::String
	test::String
	f::String
	Ω
	n::Int32

end
function Assemble(VarForm)
	#We switch among the different type of numerical methods.
	Ω = VarForm.Ω;
	println(typeof(Ω));
	if isa(Ω,NURBSurface)	
		BS, C = BezierExtraction(Ω);
		Ndof = (BS.p*BS.q)
		A = zeros(Ndof,Ndof);
		@static if BlaBla 
			println("Assembling Stifness Matrix Using Bezier FE");
			println("Ndof: ",Ndof,", dim(Ω): ",VarForm.n);
		end
		return A;
	end
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

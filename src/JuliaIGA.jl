module JuliaIGA

using Plots #Plotting library
using LinearAlgebra

include("Geo/CAD.jl");
include("Geo/Tri.jl");
include("NM/BEF.jl");
greet() = print("Hello World, I'm JuliaIGA release 0.0.1")
# Variational Formulation
BlaBla = true;

mutable struct VariationalFormulation
	a::String
	f::String
	trial::String
	test::String
	Ω
	n::Int32

end
function Assemble(VarForm)
	#We switch among the different type of numerical methods.
	Ω = VarForm.Ω;
	println(typeof(Ω));
	bformflag = 0; #Flag for the type of problem
	if isa(Ω,NURBSurface)	
		BS, C, Wb, W = BezierExtraction(Ω);
		Ndof = ((Ω.p+1)*(Ω.q+1))^2;
		A = zeros(Ndof,Ndof);
		@static if BlaBla 
			println("Assembling Stifness Matrix Using Bezier FE ... ");
			println("Ndof: ",Ndof,", dim(Ω): ",VarForm.n);
		end
		if VarForm.a == "∇"*VarForm.test*"∇"*VarForm.trial || VarForm.a == "∇"*VarForm.trial*"∇"*VarForm.test
			@static if BlaBla
				println("Bilinear form of the Poisson problem.");
				bformflag = 1;
			end
		end
		if bformflag == 1
			for i in 1:Ndof
				for j in 1:Ndof 
					e = i%9;
					if e == 0
						e = 9;
					end
					k1 = i%9;
					k2 = j%9;
					if k1 == 0
						k1 = 9;
					end
					if k2 == 0
						k2 = 9;
					end
					#we use a Gaussian quadrature of second order.
					# (-1/sq(3),1/sq(3)) ------------ (1/sq(3),1/sq(3))
					# 	|                                 |
					# 	|                                 |
					# 	|                                 |
					# 	|                                 |
					# (-1/sq(3),-1/sq(3)) ----------- (1/sq(3),-1/sq(3))
					#
					d = size(Ω.ω)[1]
					R1, ∂R1,J1 = JuliaIGA.BEFShape2D([-1/sqrt(3),-1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
					R2, ∂R2,J2 = JuliaIGA.BEFShape2D([-1/sqrt(3),1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
					R3, ∂R3,J3 = JuliaIGA.BEFShape2D([1/sqrt(3),1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
					R4, ∂R4,J4 = JuliaIGA.BEFShape2D([1/sqrt(3),-1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
					#∂φ(i,1)∂φ(j,1)+∂φ(i,2)∂φ(j,2) 
					if cld(i,((Ω.p+1)*(Ω.q+1))) == cld(j,((Ω.p+1)*(Ω.q+1)))
						A[i,j] = J1*(∂R1[k1,1]*∂R1[k2,1] + ∂R1[k1,2]*∂R1[k2,2])+J2*(∂R2[k1,1]*∂R2[k2,1] + ∂R2[k1,2]*∂R2[k2,2])+J3*(∂R3[k1,1]*∂R3[k2,1] + ∂R3[k1,2]*∂R3[k2,2])+J4*(∂R4[k1,1]*∂R4[k2,1] + ∂R4[k1,2]*∂R4[k2,2]);
					end

				end
			end
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

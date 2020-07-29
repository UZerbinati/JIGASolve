module JuliaIGA

using Plots #Plotting library
using LinearAlgebra

include("Geo/CAD.jl");
include("Geo/Tri.jl");
include("NM/BEF.jl");
greet() = print("Hello World, I'm JuliaIGA release 0.0.1")
# Variational Formulation
BlaBla = true;

function fcnFromString(s)
    f = eval(Meta.parse("x -> " * s))
    return x -> Base.invokelatest(f, x)
end
function IntBern(Ω,BS,C,e,k1,k2)

	#we use a Gaussian quadrature of second order.
	# (-1/sq(3),1/sq(3)) ------------ (1/sq(3),1/sq(3))
	# 	|                                 |
	# 	|                                 |
	# 	|                                 |
	# 	|                                 |
	# (-1/sq(3),-1/sq(3)) ----------- (1/sq(3),-1/sq(3))
	d = size(Ω.ω)[1];
	R1, ∂R1,J1 = JuliaIGA.BEFShape2D([-1/sqrt(3),-1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
	R2, ∂R2,J2 = JuliaIGA.BEFShape2D([-1/sqrt(3),1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
	R3, ∂R3,J3 = JuliaIGA.BEFShape2D([1/sqrt(3),1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
	R4, ∂R4,J4 = JuliaIGA.BEFShape2D([1/sqrt(3),-1/sqrt(3)],e,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
	return J1*(∂R1[k1,1]*∂R1[k2,1] + ∂R1[k1,2]*∂R1[k2,2])+J2*(∂R2[k1,1]*∂R2[k2,1] + ∂R2[k1,2]*∂R2[k2,2])+J3*(∂R3[k1,1]*∂R3[k2,1] + ∂R3[k1,2]*∂R3[k2,2])+J4*(∂R4[k1,1]*∂R4[k2,1] + ∂R4[k1,2]*∂R4[k2,2]);
end
mutable struct VariationalFormulation
	a::String
	f::String
	trial::String
	test::String
	Ω
	n::Int32

end
function AssembleBLF(VarForm)
	#We switch among the different type of numerical methods.
	Ω = VarForm.Ω;
	println(typeof(Ω));
	bformflag = 0; #Flag for the type of problem
	if isa(Ω,NURBSurface)	
		BS, C, Wb, W = BezierExtraction(Ω);
		Ndof = (Ω.p+1)*(Ω.q+1);
		A = zeros((BS.p-1)*(BS.q-1),(BS.p-1)*(BS.q-1));
		@static if BlaBla 
			println("Assembling Stifness Matrix Using Bezier FE ... ");
			println("Ndof: ",Ndof^2,", dim(Ω): ",VarForm.n);
		end
		if VarForm.a == "∇"*VarForm.test*"∇"*VarForm.trial || VarForm.a == "∇"*VarForm.trial*"∇"*VarForm.test
			@static if BlaBla
				println("Bilinear form of the Poisson problem.");
				bformflag = 1;
			end
		end
		if bformflag == 1 
			IENM =  IEN(Ω.p,Ω.q);
			for k in 1:Ndof
				CI = inv(C[k]);
				Ke = zeros(Ndof,Ndof);
				for i in 1:Ndof
					for j in 1:Ndof
						#we use a Gaussian quadrature of second order.
						# (-1/sq(3),1/sq(3)) ------------ (1/sq(3),1/sq(3))
						# 	|                                 |
						# 	|                                 |
						# 	|                                 |
						# 	|                                 |
						# (-1/sq(3),-1/sq(3)) ----------- (1/sq(3),-1/sq(3))
						#
						d = size(Ω.ω)[1]
						R1, ∂R1,J1 = JuliaIGA.BEFShape2D([-1/sqrt(3),-1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
						R2, ∂R2,J2 = JuliaIGA.BEFShape2D([-1/sqrt(3),1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
						R3, ∂R3,J3 = JuliaIGA.BEFShape2D([1/sqrt(3),1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
						R4, ∂R4,J4 = JuliaIGA.BEFShape2D([1/sqrt(3),-1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
						Ke[i,j] = J1*(∂R1[i,1]*∂R1[j,1] + ∂R1[i,2]*∂R1[j,2])+J2*(∂R2[i,1]*∂R2[j,1] + ∂R2[i,2]*∂R2[j,2])+J3*(∂R3[i,1]*∂R3[j,1] + ∂R3[i,2]*∂R3[j,2])+J4*(∂R4[i,1]*∂R4[j,1] + ∂R4[i,2]*∂R4[j,2]);
					end
				end
				CIKe = CI*Ke;
				for i in 1:Ndof
					for j in 1:Ndof
						A[Int(IENM[i,k]),Int(IENM[j,k])] = A[Int(IENM[i,k]),Int(IENM[j,k])]+CIKe[i,j];
					end
				end
			end
		end
		return A;
	end
end
function AssembleLF(VarForm)
	@static if BlaBla
		println("The linear functional is, ",VarForm.f);
	end
	Ω = VarForm.Ω;
	f = fcnFromString(replace(replace(VarForm.f,"x"=>"x[1]"),"y"=>"x[2]"))
	if isa(Ω,NURBSurface)	
		BS, C, Wb, W = BezierExtraction(Ω);
		Ndof = (Ω.p+1)*(Ω.q+1);
		F = zeros((BS.p-1)*(BS.q-1),1);
		@static if BlaBla 
			println("Assembling Linear form... ");
			println("Ndof: ",Ndof,", dim(Ω): ",VarForm.n);
		end
		IENM =  IEN(Ω.p,Ω.q);
		for k in 1:Ndof
			Fe = zeros(Ndof,1);
			for i in 1:Ndof
			
				d = size(Ω.ω)[1]
				R1, ∂R1,J1 = JuliaIGA.BEFShape2D([-1/sqrt(3),-1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
				R2, ∂R2,J2 = JuliaIGA.BEFShape2D([-1/sqrt(3),1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
				R3, ∂R3,J3 = JuliaIGA.BEFShape2D([1/sqrt(3),1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
				R4, ∂R4,J4 = JuliaIGA.BEFShape2D([1/sqrt(3),-1/sqrt(3)],k,BS,[2,2],C,reshape(Ω.ω',1,d*d),Wb,(Ω.p+1)*(Ω.q+1))
				Fe[i] = J1*(R1[i]*f([-1/sqrt(3),-1/sqrt(3)]))+J2*(R2[i]*f([-1/sqrt(3),1/sqrt(3)]))+J3*(R3[i]*f([1/sqrt(3),1/sqrt(3)]))+J4*(R4[i]*f([1/sqrt(3),-1/sqrt(3)]));

			end
			for i in 1:Ndof
				F[Int(IENM[i,k])]=F[Int(IENM[i,k])]+Fe[i];
			end
		end

	end
	return F;
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

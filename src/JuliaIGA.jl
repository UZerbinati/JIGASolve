module JuliaIGA

using Plots #Plotting library
using LinearAlgebra

include("CAD.jl");

greet() = print("Hello World, I'm JuliaIGA release 0.0.1")
# BEZIER EXTRACTION
function IEN(p,q)
	M = zeros((p+1)^2,(q+1)^2);
	for m in 0:p
		for k in 0:p
			σ=(k+m)*(p+q+1)+1;
			for i in 0:p
				for j in 0:p
					M[k*(p+1)+1+j,m*(q+1)+1+i] = σ+i+j;
				end
			end
		end
	end
	return M;
end
function LocalExtractionOperators(ξ,p)
	C = [];
	#Inizializations
	m = length(ξ) #Size of Knot Veector
	a = p+1;
	b = a+1;
	nb = 1;
	C = push!(C,Matrix{Float64}(I, a, a));
	K = ξ[1:a];
	while b < m
		#Adding next extraction operation
		C = push!(C,Matrix{Float64}(I, p+1, p+1));
		i = b;
		#Counting multiplicity of knot
		K = push!(K,ξ[b]); #Filling new vector knot 
		while b < m && ξ[b+1]==ξ[b]
			K = push!(K,ξ[b+1])#Filling new vector knot
			b=b+1
		end
		K = push!(K,ξ[b]); #Filling new vector knot
		multiplicity = b-i+1;
		if multiplicity < p
			#Computing alpha for knot inserction
			N = ξ[b]-ξ[a];
			A = zeros(1,p+1);
			for j in p:-1:multiplicity+1
				A[j-multiplicity] = N / (ξ[a+j]-ξ[a]);	
			end
			r = p-multiplicity;
			# Updating matrix coefficenets
			for j in 1:r
				save=r-j+1;
				s=multiplicity+j;
				for k in p+1:-1:s+1
					α = A[k-s];
					C[nb][:,k] = α*C[nb][:,k]+(1.0-α)*C[nb][:,k-1];
				end
				if b < m
					#Updating overlaping coefficents
					C[nb+1][save:j+save,save] = C[nb][p-j+1:p+1,p+1];
				end
			end
			nb = nb+1;
			if b<m
				#Update Index a,b
				a=b;
				b=b+1;
			end
		end
	end
	return C[1:end-1],K[1:end-1];
end
function BezierExtraction(S)
	n = length(S.K)-S.p-1;
	m = length(S.H)-S.q-1;
	B = reshape(S.B,n,3,m); #Reshape S.B to be accesed
	#Creating local Bezier extraction operator
	Cξ,K = LocalExtractionOperators(S.K,S.p); 
	Cη,H = LocalExtractionOperators(S.H,S.q);
	Q = [];
	Qx = [];
	Qy = [];
	Qz = [];
	T = [];
	#Renumbering the control points (IMPROVABLE)
	P = [];
	for i in 1:n
		for j in 1:m
			P = push!(P,vcat(cat(B[i,:,j],dims=1),S.ω[i,j]));
		end
	end
	if typeof(S) == NURBSurface
		for i in 1:S.p+1
			for j in 1:S.q+1
				e = 3*(i-1)+j;
				A = IEN(S.p,S.q); #Building the IEN elements
				#Using the IEN matrix to build P
				Pe = zeros((S.p+1)*(S.q+1),4);
				for a in 1:(S.p+1)*(S.q+1)
					k = A[a,e];
					Pe[a,1] = P[Int(k)][1]
					Pe[a,2] = P[Int(k)][2]
					Pe[a,3] = P[Int(k)][3]
					Pe[a,4] = P[Int(k)][4]
				end
				C = kron(Cη[i],Cξ[j]);
				"""
				println("------|Pe|------");
				println(Pe);
				println("------|W|------");	
				"""
				Wb = diagm(transpose(C)*Pe[:,4]);
				W = diagm(Pe[:,4])
				"""
				println(W);
				println("------|Wb|------");	
				println(Wb);
				println("------|Qe|------");
				"""
				Qe = inv(Wb)*transpose(C)*W*Pe;
				"""
				println(Qe);
				println("------|-|------");
				"""

				Q = push!(Q,Qe);
				for l in 1:(S.p+1)*(S.q+1)
					Qx=push!(Qx,Qe[l,1]);
					Qy=push!(Qy,Qe[l,2]);
					Qz=push!(Qz,Qe[l,3]);
					T = push!(T,"("*string(e)*","*string(l)*")");
				end
			end
		end
	end
	"""
	println("******|Q|******");
	println(Q)
	"""
	A =  IEN(S.p+1,S.q+1);
	#Creating Bezier Control Points
	BCP = zeros(S.p+S.q+3,3,S.p+S.q+3);
	#BCP[i,:,j] = Q[Int(k)][1:3];
	for e = 1:(S.p+1)*(S.q+1)
		for l = 1:(S.p+1)*(S.q+1)
			i = div(e,(S.p+1))*S.p+div(l,(S.p+1))+1;
			if (e%(S.p+1))==0
				i=i-S.p;
			end
			if (l%(S.p+1))==0
				i=i-1;
				j = ((e-1)%(S.p+1))*S.p+(S.p+1);
			else
				j = (l%(S.p+1))+((e-1)%(S.p+1))*S.p
			end
			"""
			println("(",i,",",j,") (",e,",",l,")");
			println(Q[e][l,1:3]);
			"""
			BCP[i,:,j] = Q[e][l,1:3];
		end
	end
	BCM = reshape(BCP,S.p+S.q+3,3*(S.p+S.q+3));#Converting in a format for Bezier Surface definition
	
	#Creating the Bezier Surface	
	BS = BezierSurface((S.p+S.q+2),(S.p+S.q+2),[[1 0] [0 0]; [0 0] [0 1]]);
	BS.B = BCM;
	return BS;
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

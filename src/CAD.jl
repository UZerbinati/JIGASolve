#-----------------
#| BEZIER CURVES |
#-----------------
struct BezierCurve3D
	p::Int32
	V::Array{Array{Real,2},1}
end
struct BezierCurve2D
	p::Int32
	V::Array{Array{Real,2},1}
end
struct BezierCurve
	p::Int32
	V::Array{Real,1}
end
function DeCasteljau(C,t)
	if typeof(C) == BezierCurve
		n = length(C.V);
		P = zeros(n,n);
		for i in 1:n
			P[1,i] = C.V[i]
		end
		for k in 2:n
			for j in 1:n-k+1
				P[k,j] = (1-t)*P[k-1,j]+t*P[k-1,j+1];
			end
		end
		return P[n,1]
	elseif typeof(C) == BezierCurve2D
		X = [];
		Y = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
		end
		Cx = BezierCurve(C.p,X);
		Cy = BezierCurve(C.p,Y);
		return [DeCasteljau(Cx,t) DeCasteljau(Cy,t)]
	elseif typeof(C) == BezierCurve3D
		X = [];
		Y = [];
		Z = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
			push!(Z,C.V[k][3])
		end
		Cx = BezierCurve(C.p,X);
		Cy = BezierCurve(C.p,Y);
		Cz = BezierCurve(C.p,Z);
		return [DeCasteljau(Cx,t) DeCasteljau(Cy,t) DeCasteljau(Cz,t)]
	elseif typeof(C) == BezierSurface
		#C.B C.p righe C.q Colonne
		B = reshape(C.B,C.p+1,3,C.q+1);
		#println(B);
		Q = [];
		for j in 1:C.p+1	
			V = [];
			for i in 1:C.q+1
				#println(cat(B[j,:,i],dims=2))
				push!(V,cat(B[j,:,i],dims=2))
			end
			#println(V)
			X = BezierCurve3D(C.q,V);
			push!(Q,DeCasteljau(X,t[2]));
		end
		#println(Q)
		Y = BezierCurve3D(C.p,Q);
		return DeCasteljau(Y,t[1]);
	else
		error("DeCaseteljau algorithm only work with Bezier curves");
	end
end
function BezierPlot(C,h,Poly=true)
	if typeof(C) == BezierCurve2D
		T = 0:h:1;
		BeX = [DeCasteljau(C,t)[1] for t in T]
		BeY = [DeCasteljau(C,t)[2] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		plot(BeX,BeY,label="Curve")
		if Poly
			scatter!(Vx,Vy,label="Control Points")
		end
	elseif typeof(C) == BezierCurve3D
		T = 0:h:1;
		BeX = [DeCasteljau(C,t)[1] for t in T]
		BeY = [DeCasteljau(C,t)[2] for t in T]
		BeZ = [DeCasteljau(C,t)[3] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		Vz = [C.V[i][3] for i in 1:length(C.V)]
		plot(BeX,BeY,BeZ,label="Curve")
		if Poly
			scatter!(Vx,Vy,Vz,label="Control Points")
		end
	elseif typeof(C) == BezierSurface
		x = 0:h:1
		n = length(x);
		X =  zeros(n,n)
		Y =  zeros(n,n)
		Z =  zeros(n,n)
		for i in 1:n
			for j in 1:n
				X[i,j]=DeCasteljau(C,[x[i] x[j]])[1]
				Y[i,j]=DeCasteljau(C,[x[i] x[j]])[2]
				Z[i,j]=DeCasteljau(C,[x[i] x[j]])[3]
			end
		end
		return(X,Y,Z)
	else
		error("De Casteljau algorithm only work with Bezier curves");
	end
end

#-------------------
#| BSPLINES CURVES |
#-------------------

struct BSplineCurve4D
	p::Int32
	V::Array{Array{Real,2},1}
	K::Array{Real,1}
end
struct BSplineCurve3D
	p::Int32
	V::Array{Array{Real,2},1}
	K::Array{Real,1}
end
struct BSplineCurve2D
	p::Int32
	V::Array{Array{Real,2},1}
	K::Array{Real,1}
end
struct BSplineCurve
	p::Int32
	V::Array{Real,1}
	K::Array{Real,1}
end
function BSplinePlot(C,h,Poly=true)
	if typeof(C) == BSplineCurve
		T = 0:h:maximum(C.K);
		BeY = [deBoor(C,t) for t in T]
		I = [];
		for i in C.K
			if !(i in I)
				push!(I,i);
			end
		end
		plot(T,BeY,label="Curve")
	elseif typeof(C) == BSplineCurve2D
		T = 0:h:maximum(C.K);
		BeX = [deBoor(C,t)[1] for t in T]
		BeY = [deBoor(C,t)[2] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		plot(BeX,BeY,label="Curve")
		if Poly
			scatter!(Vx,Vy,label="Control Points")
		end
	elseif typeof(C) == BSplineCurve3D
		T = 0:h:maximum(C.K);
		BeX = [deBoor(C,t)[1] for t in T]
		BeY = [deBoor(C,t)[2] for t in T]
		BeZ = [deBoor(C,t)[3] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		Vz = [C.V[i][3] for i in 1:length(C.V)]
		plot(BeX,BeY,BeZ,label="Curve")
		if Poly
			scatter!(Vx,Vy,Vz,label="Control Points")
		end
	end
end
function deBoor(C,t)
	##print(t)

	#print(" ",i);
	if typeof(C) == BSplineCurve
		i = 0;
		for l in 1:length(C.K)-1
			if C.K[l]<=t<C.K[l+1]
				i = l-1;
			end
		end
		if t == maximum(C.K)
			return C.V[end]
		end
		n = length(C.V);
		m = length(C.K);
		if  m != n+C.p+1
			error("Number of control points and length of knot vector don't match the degree of the spline");
		end
		D = [C.V[j+i-C.p] for j in 1:C.p+1]
		for r in 1:C.p
			for j in C.p:-1:r
				α = (t - C.K[j+i-C.p+1])/(C.K[j+1+i-r+1]-C.K[j+i-C.p+1]);
				D[j+1] = (1.0-α)*D[j]+α*D[j+1];
			end
		end
		return D[C.p+1];
	elseif typeof(C) == BSplineCurve2D
		X = [];
		Y = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
		end
		Cx = BSplineCurve(C.p,X,C.K);
		Cy = BSplineCurve(C.p,Y,C.K);
		return [deBoor(Cx,t) deBoor(Cy,t)]
	elseif typeof(C) == BSplineCurve3D
		X = [];
		Y = [];
		Z = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
			push!(Z,C.V[k][3])
		end
		Cx = BSplineCurve(C.p,X,C.K);
		Cy = BSplineCurve(C.p,Y,C.K);
		Cz = BSplineCurve(C.p,Z,C.K);
		return [deBoor(Cx,t) deBoor(Cy,t) deBoor(Cz,t)]
	elseif typeof(C) == BSplineCurve4D
		X = [];
		Y = [];
		Z = [];
		W = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
			push!(Z,C.V[k][3])
			push!(W,C.V[k][4])
		end
		Cx = BSplineCurve(C.p,X,C.K);
		Cy = BSplineCurve(C.p,Y,C.K);
		Cz = BSplineCurve(C.p,Z,C.K);
		Cw = BSplineCurve(C.p,W,C.K);
		return [deBoor(Cx,t) deBoor(Cy,t) deBoor(Cz,t) deBoor(Cw,t)]
	elseif typeof(C) == BSplineManifold
		n = length(C.K)-C.p-1;
		m = length(C.H)-C.q-1;
		B = reshape(C.B,n,4,m);
		Q = [];
		for j in 1:n	
			V = [];
			for i in 1:m
				push!(V,cat(B[j,:,i],dims=2))
			end
			X = BSplineCurve4D(C.q,V,C.K);
			push!(Q,deBoor(X,t[2]));
		end
		Y = BSplineCurve4D(C.p,Q,C.H);
		return deBoor(Y,t[1]);
	end
end

#---------------
#| NURBS CURVE |
#---------------

struct NURBSCurve3D
	p::Int32
	V::Array{Array{Real,2},1}
	ω::Array{Real,1}
	K::Array{Real,1}
end
struct NURBSCurve2D
	p::Int32
	V::Array{Array{Real,2},1}
	ω::Array{Real,1}
	K::Array{Real,1}
end
struct NURBSCurve
	p::Int32
	V::Array{Real,1}
	ω::Array{Real,1}
	K::Array{Real,1}
end
function HomoNURBS(C,t)
	#Mapping coordinates from R^d to R^(d+1)
	if typeof(C) == NURBSCurve
		V = [];
		for i in 1:length(C.V)
			push!(V,[C.ω[i]*C.V[i] C.ω[i]]);
		end
		BS = BSplineCurve2D(C.p,V,C.K);
		x = deBoor(BS,t);
		if (x[end]==0)
			return x[1:end-1]
		else
			return (1/x[end])*x[1:end-1]
		end
	elseif typeof(C) == NURBSCurve2D
		X = [];
		Y = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
		end
		Cx = NURBSCurve(C.p,X,C.ω,C.K);
		Cy = NURBSCurve(C.p,Y,C.ω,C.K);
		return [HomoNURBS(Cx,t) HomoNURBS(Cy,t)]
	elseif typeof(C) == NURBSCurve3D
		X = [];
		Y = [];
		Z = [];
		for k in 1:length(C.V)
			push!(X,C.V[k][1])
			push!(Y,C.V[k][2])
			push!(Z,C.V[k][3])
		end
		Cx = NURBSCurve(C.p,X,C.ω,C.K);
		Cy = NURBSCurve(C.p,Y,C.ω,C.K);
		Cz = NURBSCurve(C.p,Z,C.ω,C.K);
		return [HomoNURBS(Cx,t) HomoNURBS(Cy,t) HomoNURBS(Cz,t)]
	elseif typeof(C) == NURBSurface
		n = length(C.K)-C.p-1;
		m = length(C.H)-C.q-1;
		B = reshape(C.B,n,3,m);
		Q = zeros(n,4,m);
		for i in 1:n
			for j in 1:m
				P = zeros(4,1);
				P[1:3] = cat(B[i,:,j],dims=2)
				P[4] = 1.0;
				P = C.ω[i,j]*P;
				Q[i,:,j]=P;
			end
		end	
		BS = BSplineManifold(C.p,C.q,C.K,C.H,[[1 0] [0 0]; [0 0] [0 1]]);
		BS.B = reshape(Q,n,4*m);
		x = deBoor(BS,t);
		if (x[end]==0)
			return x[1:end-1]
		else
			return (1/x[end])*x[1:end-1]
		end

	end
end
function NURBSEval(S,t)
	if typeof(S) == NURBSurface
		n = length(S.K)-S.p-1;
		m = length(S.H)-S.q-1;
		Q = zeros(n,m);
		B = reshape(S.B,n,3,m);
		for i in 1:n
			for j in 1:m
				V1 = []
				for k in 1:length(S.K)-S.p-1
					if k ==i
						push!(V1,1.0)
					else
						push!(V1,0.0)
					end
				end
				C1=BSplineCurve(S.p,V1,S.K);
				V2 = []
				for k in 1:length(S.H)-S.q-1
					if k ==j
						push!(V2,1.0)
					else
						push!(V2,0.0)
					end
				end
				C2=BSplineCurve(S.q,V2,S.H);
				Q[i,j] = deBoor(C1,t[1])*deBoor(C2,t[2]);
			end
		end
		Σ = zeros(3,1);
		for i = 1:n
			Z=zeros(3,1);
			for j =1:m
				D = 0;
				for α = 1:n
					T=0;
					for β = 1:m
						T = T+S.ω[α,β]*Q[α,β];
					end
					D = D+T;
				end
				Z = Z+(B[i,:,j]*Q[i,j]*S.ω[i,j])/D;
			end
			Σ = Σ+Z;
		end
		return Σ;
	end
end
function NURBSPlot(C,h;Poly=true,opt="Projection")
	if typeof(C) == NURBSCurve
		T = 0:h:maximum(C.K);
		BeY = [HomoNURBS(C,t) for t in T]
		I = [];
		for i in C.K
			if !(i in I)
				push!(I,i);
			end
		end
		plot(T,BeY,label="Curve")
	elseif typeof(C) == NURBSCurve2D
		T = 0:h:maximum(C.K);
		BeX = [HomoNURBS(C,t)[1] for t in T]
		BeY = [HomoNURBS(C,t)[2] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		plot(BeX,BeY,label="Curve")
		if Poly
			scatter!(Vx,Vy,label="Control Points")
		end
	elseif typeof(C) == NURBSCurve3D
		T = 0:h:maximum(C.K);
		BeX = [HomoNURBS(C,t)[1] for t in T]
		BeY = [HomoNURBS(C,t)[2] for t in T]
		BeZ = [HomoNURBS(C,t)[3] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		Vz = [C.V[i][3] for i in 1:length(C.V)]
		plot(BeX,BeY,BeZ,label="Curve")
		if Poly
			scatter!(Vx,Vy,Vz,label="Control Points")
		end
	elseif typeof(C) == NURBSurface
		x = 0:h:maximum(C.K);
		y = 0:h:maximum(C.H);
		n = length(x);
		X =  zeros(n,n)
		Y =  zeros(n,n)
		Z =  zeros(n,n)
		for i in 1:n
			for j in 1:n
				if opt=="Projection"
					X[i,j]=HomoNURBS(C,[x[i] y[j]])[1]
					Y[i,j]=HomoNURBS(C,[x[i] y[j]])[2]
					Z[i,j]=HomoNURBS(C,[x[i] y[j]])[3]
				elseif opt=="Direct"
					X[i,j]=NURBSEval(C,[x[i] y[j]])[1]
					Y[i,j]=NURBSEval(C,[x[i] y[j]])[2]
					Z[i,j]=NURBSEval(C,[x[i] y[j]])[3]
				end
			end
		end
		return(X,Y,Z)
	end
end

#-------------
# | SURFACES |
#-------------

mutable struct BezierSurface
	p::Int64
	q::Int64
	B::Array{Real,2}
end
mutable struct BSplineManifold
	p::Int64
	q::Int64
	K::Array{Real,1}
	H::Array{Real,1}
	B::Array{Real,2}
end
struct NURBSurface
	p::Int64 #Grado 1
	q::Int64 #Grado 2
	K::Array{Real,1} #Nodi 1
	H::Array{Real,1} #Nodi 2
	B::Array{Real,2} #Punti di Controllo
	ω::Array{Real,2} #Pesi
end
#---------------------
#| BEZIER EXTRACTION |
#---------------------

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
				Wb = diagm(transpose(C)*Pe[:,4]);
				W = diagm(Pe[:,4])
				Qe = inv(Wb)*transpose(C)*W*Pe;
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
	A =  IEN(S.p+1,S.q+1);
	#Creating Bezier Control Points
	BCP = zeros(S.p+S.q+3,3,S.p+S.q+3);
	#Mapping Bezier local Bezier elements to global mesh
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
			BCP[i,:,j] = Q[e][l,1:3];
		end
	end
	BCM = reshape(BCP,S.p+S.q+3,3*(S.p+S.q+3));#Converting in a format aceptable for Bezier Surface definition
	
	#Creating the Bezier Surface	
	BS = BezierSurface((S.p+S.q+2),(S.p+S.q+2),[[1 0] [0 0]; [0 0] [0 1]]);
	BS.B = BCM;
	return BS;
end
#------------------
#| CONTROL POINTS |
#------------------

function ControlPlot(C)
	B = reshape(C.B,C.p+1,3,C.q+1);
	X = [];
	Y = [];
	Z = [];
	if typeof(C) == BezierSurface
		for i in 1:C.p+1
			for j in 1:C.q+1
				P = cat(B[i,:,j],dims=2);
				X = push!(X,P[1]);
				Y = push!(Y,P[2]);
				Z = push!(Z,P[3]);
			end
		end
	end
	scatter!(X,Y,Z);

end


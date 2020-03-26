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
		T = 0+h:h:maximum(C.K)-h;
		BeY = [deBoor(C,t) for t in T]
		I = [];
		for i in C.K
			if !(i in I)
				push!(I,i);
			end
		end
		plot(T,BeY,label="Curve")
	elseif typeof(C) == BSplineCurve2D
		T = 0+h:h:maximum(C.K)-h;
		BeX = [deBoor(C,t)[1] for t in T]
		BeY = [deBoor(C,t)[2] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		plot(BeX,BeY,label="Curve")
		if Poly
			scatter!(Vx,Vy,label="Control Points")
		end
	elseif typeof(C) == BSplineCurve3D
		T = 0+h:h:maximum(C.K)-h;
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
	elseif typeof(C) == BSplineSurface
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
	end
end
function deBoor(C,t)
	##print(t)
	i = 0;
	for l in 1:length(C.K)-1
		if C.K[l]<=t<C.K[l+1]
			i = l-1;
		end
	end
	#print(" ",i);
	if typeof(C) == BSplineCurve
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
			#print(" ",D);
		end
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
	end
	#println(" ",D[C.p+1]);
	return D[C.p+1];
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
		#print(V)
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
	end
end
function NURBSPlot(C,h,Poly=true)
	if typeof(C) == NURBSCurve
		T = 0+h:h:maximum(C.K)-h;
		BeY = [HomoNURBS(C,t) for t in T]
		I = [];
		for i in C.K
			if !(i in I)
				push!(I,i);
			end
		end
		plot(T,BeY,label="Curve")
	elseif typeof(C) == NURBSCurve2D
		T = 0+h:h:maximum(C.K)-h;
		BeX = [HomoNURBS(C,t)[1] for t in T]
		BeY = [HomoNURBS(C,t)[2] for t in T]
		Vx = [C.V[i][1] for i in 1:length(C.V)]
		Vy = [C.V[i][2] for i in 1:length(C.V)]
		plot(BeX,BeY,label="Curve")
		if Poly
			scatter!(Vx,Vy,label="Control Points")
		end
	elseif typeof(C) == NURBSCurve3D
		T = 0+h:h:maximum(C.K)-h;
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
	end
end

#-------------
# | SURFACES |
#-------------

struct BezierSurface
	p::Int64
	q::Int64
	B::Array{Real,2}
end
struct BSplineSurface
	p::Int64
	q::Int64
	B::Array{Real,2}
end


